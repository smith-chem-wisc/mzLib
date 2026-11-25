using Omics.Digestion;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;

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
	public static List<MslLibraryEntry> CreateDecoys(IReadOnlyList<MslLibraryEntry> targets, List<DigestionMotif> motifs,
		out MslDecoyAugmentationResult result, Dictionary<string, Modification>? knownModifications = null)
	{
		result = new MslDecoyAugmentationResult(0, 0, 0, 0); // red: not written yet
		return new List<MslLibraryEntry>();
	}

	/// <summary>
	/// Writes <paramref name="outputPath"/>: every entry of the decoy-free library at <paramref name="inputPath"/>, plus its
	/// decoy (<see cref="CreateDecoys"/>). The input is read index-only and the output written streaming.
	/// </summary>
	/// <exception cref="MzLibUtil.MzLibException">The input already holds decoys, or the motifs' sites cannot be held in place.</exception>
	public static MslDecoyAugmentationResult AddDecoys(string inputPath, string outputPath, List<DigestionMotif> motifs,
		Dictionary<string, Modification>? knownModifications = null)
	{
		return new MslDecoyAugmentationResult(0, 0, 0, 0); // red: not written yet
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
