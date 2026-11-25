using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases.EntrapmentGeneration;

namespace Readers.SpectralLibrary;

/// <summary>
/// Gives a decoy-free <c>.msl</c> library its decoys, each derived from its target's own entry, so an experimentally derived
/// library (a DDA search's, for example) can be searched with a target-decoy FDR. No spectrum is predicted: a predicted decoy
/// beside an empirical target would differ from it in kind, and a classifier could learn that instead of what is wrong.
/// </summary>
/// <remarks>
/// <para>Each target or entrapment entry gets one decoy:</para>
/// <list type="bullet">
///   <item>its sequence permuted by <c>EntrapmentPeptideGenerator.Create</c>, cleavage sites held in place, never equal to a
///   library sequence up to I/L, and the same for every charge of a peptide;</item>
///   <item>its modifications moved with their residues, terminal ones kept on their termini, so the precursor m/z is
///   unchanged;</item>
///   <item>its fragments kept (type, number, charge, neutral loss, intensity), each at the decoy's m/z;</item>
///   <item>its retention time, charge, entrapment flag and other properties kept; each accession named <c>DECOY_</c> plus the
///   target's.</item>
/// </list>
/// <para>Add decoys after any retention-time update: <c>MslRetentionTimeUpdater</c> moves targets, not decoys.</para>
/// </remarks>
public static class MslLibraryAugmenter
{
	/// <summary>
	/// The decoys of <paramref name="targets"/>, one per entry that has a possible partner.
	/// </summary>
	/// <param name="targets">The library's entries. None may be a decoy.</param>
	/// <param name="motifs">The library's digestion motifs, whose residues stay in place (<c>DigestionAgent.DigestionMotifs</c>).</param>
	/// <param name="result">What was made, and what was not.</param>
	/// <param name="knownModifications">Modifications the full sequences name; null reads mzLib's catalog.</param>
	/// <exception cref="MzLibUtil.MzLibException">An entry is already a decoy, or the motifs' sites cannot be held in place.</exception>
	/// <param name="decoyIdentifier">Prefixed, with an underscore, to each accession, as mzLib's decoy generators name theirs.</param>
	public static List<MslLibraryEntry> CreateDecoys(IReadOnlyList<MslLibraryEntry> targets, List<DigestionMotif> motifs,
		out MslDecoyAugmentationResult result, Dictionary<string, Modification>? knownModifications = null, string decoyIdentifier = "DECOY")
	{
		RefuseDecoys(targets.Any(entry => entry.IsDecoy));
		var maker = new DecoyMaker(targets.Select(entry => entry.BaseSequence), motifs, knownModifications, decoyIdentifier);
		var decoys = new List<MslLibraryEntry>();
		foreach (var target in targets)
		{
			if (maker.Decoy(target) is { } decoy)
				decoys.Add(decoy);
		}
		result = maker.Result;
		return decoys;
	}

	/// <summary>
	/// Writes <paramref name="outputPath"/>: every entry of the decoy-free library at <paramref name="inputPath"/>, each followed
	/// by its decoy (<see cref="CreateDecoys"/>). The input is read index-only, one entry's fragments at a time, and the output
	/// is written streaming, so a library of any size fits.
	/// </summary>
	/// <exception cref="MzLibUtil.MzLibException">The input already holds decoys, or the motifs' sites cannot be held in place.</exception>
	public static MslDecoyAugmentationResult AddDecoys(string inputPath, string outputPath, List<DigestionMotif> motifs,
		Dictionary<string, Modification>? knownModifications = null, string decoyIdentifier = "DECOY")
	{
		using var library = MslLibrary.LoadIndexOnly(inputPath);
		RefuseDecoys(library.DecoyCount > 0);
		var sequences = Enumerable.Range(0, library.PrecursorCount)
			.Select(i => library.GetEntry(i, includeFragments: false)!.BaseSequence);
		var maker = new DecoyMaker(sequences, motifs, knownModifications, decoyIdentifier);
		MslWriter.WriteStreaming(outputPath, WithDecoys(library.GetAllEntries(), maker));
		return maker.Result;
	}

	private static IEnumerable<MslLibraryEntry> WithDecoys(IEnumerable<MslLibraryEntry> targets, DecoyMaker maker)
	{
		foreach (var target in targets)
		{
			yield return target;
			if (maker.Decoy(target) is { } decoy)
				yield return decoy;
		}
	}

	private static void RefuseDecoys(bool libraryHasDecoys)
	{
		if (libraryHasDecoys)
			throw new MzLibException("The library already holds decoys. Decoys are added only to a library that has none, "
				+ "so that no library mixes two decoy recipes.");
	}

	/// <summary>Makes one entry's decoy at a time, against the whole library's sequences, and counts what it made.</summary>
	private sealed class DecoyMaker
	{
		private readonly List<DigestionMotif> _motifs;
		private readonly Dictionary<string, Modification> _knownModifications;
		private readonly string _decoyIdentifier;
		private readonly HashSet<string> _librarySequences;
		private readonly HashSet<string> _libraryLeucineSequences;
		// Each base sequence's partner, so every charge and modified form of a peptide gets the same rearrangement
		private readonly Dictionary<string, EntrapmentPeptide> _partners = new(StringComparer.Ordinal);
		private int _read, _made, _without, _leftOut;

		public DecoyMaker(IEnumerable<string> librarySequences, List<DigestionMotif> motifs,
			Dictionary<string, Modification>? knownModifications, string decoyIdentifier)
		{
			_motifs = motifs;
			_knownModifications = knownModifications ?? Mods.AllKnownProteinModsDictionary;
			_decoyIdentifier = decoyIdentifier;
			_librarySequences = librarySequences.ToHashSet(StringComparer.Ordinal);
			_libraryLeucineSequences = _librarySequences.Select(IBioPolymerWithSetMods.LeucineEquivalent).ToHashSet(StringComparer.Ordinal);
		}

		public MslDecoyAugmentationResult Result => new(_read, _made, _without, _leftOut);

		public MslLibraryEntry? Decoy(MslLibraryEntry target)
		{
			_read++;
			if (!_partners.TryGetValue(target.BaseSequence, out var partner))
			{
				// Fold 1 of 3 and seed 1, as the predicted library's permuted decoys; a candidate equal to a library sequence up to
				// I/L is passed over for the next, so a decoy is never a target's spectrum under another name
				partner = EntrapmentPeptideGenerator.Create(target.BaseSequence, _motifs, _librarySequences, fold: 1, foldCount: 3, seed: 1,
					rejectInContext: candidate => _libraryLeucineSequences.Contains(IBioPolymerWithSetMods.LeucineEquivalent(candidate)));
				_partners.Add(target.BaseSequence, partner);
			}

			Dictionary<int, Modification>? modifications = null;
			try
			{
				modifications = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(target.FullSequence, _knownModifications);
			}
			catch (MzLibException)
			{
				// a modification the catalog does not know: the entry is counted as without a decoy
			}
			if (!partner.Succeeded || modifications is null)
			{
				_without++;
				return null;
			}

			string fullSequence = IBioPolymerWithSetMods.DetermineFullSequence(partner.EntrapmentSequence!,
				IBioPolymerWithSetMods.PermuteModifications(modifications, partner.SwappedPositions!));
			var decoy = target.WithFragments(DecoyIons(target, fullSequence));
			decoy.FullSequence = fullSequence;
			decoy.BaseSequence = partner.EntrapmentSequence!;
			decoy.IsDecoy = true;
			decoy.ProteinAccession = string.Join("|", target.ProteinAccession.Split('|', StringSplitOptions.RemoveEmptyEntries)
				.Select(accession => _decoyIdentifier + "_" + accession));
			_made++;
			return decoy;
		}

		/// <summary>Each target ion at the decoy's m/z, with the target's type, number, charge, loss and intensity.</summary>
		private List<MslFragmentIon> DecoyIons(MslLibraryEntry target, string decoyFullSequence)
		{
			var peptide = new PeptideWithSetModifications(decoyFullSequence, _knownModifications);
			// The product types do not depend on the dissociation type's name, and Fragment makes none for Unknown
			var dissociation = target.DissociationType == DissociationType.Unknown ? DissociationType.HCD : target.DissociationType;
			var products = new List<Product>();
			peptide.Fragment(dissociation, FragmentationTerminus.Both, products);
			var internalIons = target.MatchedFragmentIons.Where(ion => ion.IsInternalFragment).ToList();
			if (internalIons.Count > 0)
			{
				// FragmentInternally clears the list it is given, so the internal products are made apart and appended
				int shortest = internalIons.Min(ion => ion.SecondaryFragmentNumber - ion.FragmentNumber);
				var internalProducts = new List<Product>();
				peptide.FragmentInternally(dissociation, Math.Max(1, shortest), internalProducts);
				products.AddRange(internalProducts);
			}

			var mapped = LibrarySpectrum.GetDecoyProducts(target.MatchedFragmentIons, products);
			var ions = new List<MslFragmentIon>(mapped.Length);
			for (int i = 0; i < mapped.Length; i++)
			{
				var ion = target.MatchedFragmentIons[i];
				if (mapped[i] is not { } product)
				{
					_leftOut++;
					continue;
				}
				ions.Add(new MslFragmentIon
				{
					ProductType = ion.ProductType,
					SecondaryProductType = ion.SecondaryProductType,
					FragmentNumber = ion.FragmentNumber,
					SecondaryFragmentNumber = ion.SecondaryFragmentNumber,
					ResiduePosition = product.ResiduePosition,
					Charge = ion.Charge,
					NeutralLoss = ion.NeutralLoss,
					Intensity = ion.Intensity,
					ExcludeFromQuant = ion.ExcludeFromQuant,
					Mz = (float)product.NeutralMass.ToMz(ion.Charge),
				});
			}
			return ions;
		}
	}
}

/// <summary>What <see cref="MslLibraryAugmenter"/> made.</summary>
/// <param name="TargetsRead">Library entries read (targets and entrapments).</param>
/// <param name="DecoysMade">Decoys made, one per entry that has a possible partner.</param>
/// <param name="TargetsWithoutDecoy">Entries with no decoy: every rearrangement of their free residues equals a library
/// sequence up to I/L (proven by enumeration, not given up on), or their full sequence names a modification that is not known.</param>
/// <param name="FragmentsLeftOut">Target ions with no decoy counterpart (an internal ion, for example, when the decoy has no
/// such internal product), left out of their decoys.</param>
public sealed record MslDecoyAugmentationResult(int TargetsRead, int DecoysMade, int TargetsWithoutDecoy, int FragmentsLeftOut);
