using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Decoys for a decoy-free library, each derived from its target's own entry (<see cref="MslLibraryAugmenter"/>), and the two
/// helpers it is built on: <see cref="IBioPolymerWithSetMods.PermuteModifications"/> and
/// <see cref="LibrarySpectrum.GetDecoyProducts(IReadOnlyList{MslFragmentIon}, IReadOnlyList{Product})"/>.
/// </summary>
[TestFixture]
public class TestMslLibraryAugmenter
{
	private static readonly List<DigestionMotif> Trypsin = ProteaseDictionary.Dictionary["trypsin"].DigestionMotifs;
	private static readonly double Water = ChemicalFormula.ParseFormula("H2O").MonoisotopicMass;

	private static readonly Modification Oxidation = Mod("Oxidation", "M", "Anywhere.", "O");
	private static readonly Modification Acetyl = Mod("Acetyl", "X", "Peptide N-terminal.", "C2H2O");
	private static readonly Modification Amidation = Mod("Amidation", "X", "Peptide C-terminal.", "H1N1O-1");
	private static readonly Dictionary<string, Modification> Known = new[] { Oxidation, Acetyl, Amidation }.ToDictionary(m => m.IdWithMotif);

	// ── the mod-mover ────────────────────────────────────────────────────────

	[Test]
	public void ModificationsMoveWithTheirResiduesAndTerminalOnesStayOnTheirTermini()
	{
		// ABCDEK, residues 0..2 rearranged: 0 -> 2, 1 -> 0, 2 -> 1; D, E and K stay
		var swapped = new[] { 2, 0, 1, 3, 4, 5 };
		var mods = new Dictionary<int, Modification> { [1] = Acetyl, [3] = Oxidation, [5] = Oxidation, [8] = Amidation };

		var moved = IBioPolymerWithSetMods.PermuteModifications(mods, swapped);

		Assert.That(moved.Keys, Is.EquivalentTo(new[] { 1, 2, 5, 8 }), "residue 1 moved to 0 (key 2), residue 3 stayed (key 5)");
		Assert.That(moved[1], Is.SameAs(Acetyl), "the N-terminal modification stays on the N-terminus");
		Assert.That(moved[2], Is.SameAs(Oxidation));
		Assert.That(moved[8], Is.SameAs(Amidation), "the C-terminal modification (key length + 2) stays on the C-terminus");
		Assert.That(mods.Keys, Is.EquivalentTo(new[] { 1, 3, 5, 8 }), "the input is left as it was");
	}

	// ── the I/L fold, moved down from MslLibraryBuilder ─────────────────────

	[Test]
	public void LeucineEquivalentReadsEveryIAsLAndLeavesModificationNamesAlone()
	{
		Assert.That(IBioPolymerWithSetMods.LeucineEquivalent("PEIPTIDE[Test:Imaginary on I]K"),
			Is.EqualTo("PELPTLDE[Test:Imaginary on I]K"));
	}

	// ── the decoy-fragments helper ──────────────────────────────────────────

	[Test]
	public void ATargetIonMapsToTheDecoysProductOfTheSameTypeAndNumber()
	{
		var decoy = Products("EPTPIDEK");
		var target = new List<MslFragmentIon> { Ion(ProductType.b, 3), Ion(ProductType.y, 2) };

		var mapped = LibrarySpectrum.GetDecoyProducts(target, decoy);

		Assert.That(mapped[0]!.NeutralMass, Is.EqualTo(Single(decoy, ProductType.b, 3).NeutralMass).Within(1e-9));
		Assert.That(mapped[1]!.NeutralMass, Is.EqualTo(Single(decoy, ProductType.y, 2).NeutralMass).Within(1e-9));
	}

	[Test]
	public void ATargetIonsNeutralLossIsCarriedOntoItsDecoyIon()
	{
		// The old method paired on type and number alone, so b5-H2O took b5's m/z
		var decoy = Products("EPTPIDEK");
		var target = new List<MslFragmentIon> { Ion(ProductType.b, 5, neutralLoss: Water) };

		var mapped = LibrarySpectrum.GetDecoyProducts(target, decoy);

		Assert.That(mapped[0]!.NeutralMass, Is.EqualTo(Single(decoy, ProductType.b, 5).NeutralMass - Water).Within(1e-9));
		Assert.That(mapped[0]!.NeutralLoss, Is.EqualTo(Water).Within(1e-9));
	}

	[Test]
	public void AnInternalIonMapsOnlyToTheDecoysInternalIonOfTheSameSpan()
	{
		var peptide = new PeptideWithSetModifications("EPTPIDEK", new Dictionary<string, Modification>());
		var internalProducts = new List<Product>();
		peptide.FragmentInternally(DissociationType.HCD, 2, internalProducts);
		var span = internalProducts.First(p => p.SecondaryProductType != null);
		var target = new List<MslFragmentIon>
		{
			new() { ProductType = span.ProductType, SecondaryProductType = span.SecondaryProductType, FragmentNumber = span.FragmentNumber,
				SecondaryFragmentNumber = span.SecondaryFragmentNumber, Charge = 1, Intensity = 1f }
		};

		Assert.That(LibrarySpectrum.GetDecoyProducts(target, Products("EPTPIDEK"))[0], Is.Null, "no internal products were given");
		Assert.That(LibrarySpectrum.GetDecoyProducts(target, internalProducts)[0]!.NeutralMass, Is.EqualTo(span.NeutralMass).Within(1e-9));
	}

	[Test]
	public void LibrarySpectrumIonsMapTheSameWayAsMslIons()
	{
		// The overload LibrarySpectrum.GetDecoyLibrarySpectrumFromTargetByReverse can move onto (mzLib #1630)
		var decoy = Products("EPTPIDEK");
		var b5 = Single(Products("PEPTIDEK"), ProductType.b, 5);
		var lossy = new Product(ProductType.b, FragmentationTerminus.N, b5.NeutralMass - Water, 5, 5, Water);
		var target = new List<MatchedFragmentIon> { new(lossy, lossy.NeutralMass.ToMz(1), 1.0, 1) };

		var mapped = LibrarySpectrum.GetDecoyProducts(target, decoy);

		Assert.That(mapped[0]!.NeutralMass, Is.EqualTo(Single(decoy, ProductType.b, 5).NeutralMass - Water).Within(1e-9));
	}

	// ── the augmenter ───────────────────────────────────────────────────────

	[Test]
	public void EveryTargetGetsOneDecoyWithItsChargePrecursorMzAndRetentionTime()
	{
		var targets = new List<MslLibraryEntry> { Entry("PEPTIDEK", 2, irt: 12.5), Entry("SAMPLERK", 3, irt: 40) };

		var decoys = MslLibraryAugmenter.CreateDecoys(targets, Trypsin, out var result, Known);

		Assert.That(decoys, Has.Count.EqualTo(2), "one decoy per target, in the targets' order");
		Assert.That(result, Is.EqualTo(new MslDecoyAugmentationResult(2, 2, 0, 0)));
		for (int i = 0; i < targets.Count; i++)
		{
			Assert.That(decoys[i].IsDecoy, Is.True);
			Assert.That(decoys[i].ChargeState, Is.EqualTo(targets[i].ChargeState));
			Assert.That(decoys[i].PrecursorMz, Is.EqualTo(targets[i].PrecursorMz).Within(1e-6));
			Assert.That(decoys[i].RetentionTime, Is.EqualTo(targets[i].RetentionTime));
			Assert.That(targets[i].IsDecoy, Is.False, "the targets are left as they were");
		}
	}

	[Test]
	public void TheDecoyIsARearrangementOfItsTargetWithItsCleavageSitesInPlace()
	{
		var target = Entry("PEPKTIDEK"); // K at 3 (a missed cleavage) and K at 8 are trypsin's sites

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out _, Known).Single();

		Assert.That(decoy.BaseSequence, Is.Not.EqualTo(target.BaseSequence));
		Assert.That(decoy.BaseSequence.OrderBy(c => c), Is.EqualTo(target.BaseSequence.OrderBy(c => c)), "the same residues");
		Assert.That(decoy.BaseSequence[3], Is.EqualTo('K'));
		Assert.That(decoy.BaseSequence[8], Is.EqualTo('K'));
	}

	[Test]
	public void NoDecoyEqualsALibrarySequenceUpToIsoleucineAndLeucine()
	{
		// Crowded with I/L twins, so a decoy that ignored I = L would very likely land on one
		var sequences = new[] { "LIAGK", "ILAGK", "LLAGK", "IIAGK", "LAIGK", "IALGK", "AILGK", "ALIGK", "GLIAK", "GILAK" };
		var targets = sequences.Select(s => Entry(s)).ToList();

		var decoys = MslLibraryAugmenter.CreateDecoys(targets, Trypsin, out var result, Known);

		var library = sequences.Select(s => s.Replace('I', 'L')).ToHashSet();
		Assert.That(decoys, Is.Not.Empty);
		Assert.That(decoys.Select(d => d.BaseSequence.Replace('I', 'L')).Where(library.Contains), Is.Empty);
		Assert.That(result.DecoysMade + result.TargetsWithoutDecoy, Is.EqualTo(sequences.Length));
	}

	[Test]
	public void EveryChargeOfAPeptideGetsTheSameDecoy()
	{
		var decoys = MslLibraryAugmenter.CreateDecoys(new[] { Entry("PEPTIDEK", 2), Entry("PEPTIDEK", 3) }, Trypsin, out _, Known);

		Assert.That(decoys, Has.Count.EqualTo(2));
		Assert.That(decoys[1].FullSequence, Is.EqualTo(decoys[0].FullSequence), "peptide-level FDR groups by full sequence");
	}

	[Test]
	public void TheSameLibraryGetsTheSameDecoysInAnyOrder()
	{
		var targets = new List<MslLibraryEntry> { Entry("PEPTIDEK"), Entry("SAMPLERK"), Entry("ELVISLIVESK") };

		var first = MslLibraryAugmenter.CreateDecoys(targets, Trypsin, out _, Known).Select(d => d.FullSequence).ToList();
		var again = MslLibraryAugmenter.CreateDecoys(targets, Trypsin, out _, Known).Select(d => d.FullSequence).ToList();
		var reversed = MslLibraryAugmenter.CreateDecoys(Enumerable.Reverse(targets).ToList(), Trypsin, out _, Known)
			.Select(d => d.FullSequence).Reverse().ToList();

		Assert.That(first, Has.Count.EqualTo(3));
		Assert.That(again, Is.EqualTo(first));
		Assert.That(reversed, Is.EqualTo(first));
	}

	[Test]
	public void AModificationTravelsWithItsResidueSoThePrecursorMassIsUnchanged()
	{
		var target = Entry("PEPM[Test:Oxidation on M]TIDEK");

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out _, Known).Single();

		int m = decoy.BaseSequence.IndexOf('M');
		Assert.That(decoy.FullSequence, Does.Contain("M[Test:Oxidation on M]"));
		Assert.That(IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(decoy.FullSequence, Known).Keys, Is.EqualTo(new[] { m + 2 }));
		Assert.That(new PeptideWithSetModifications(decoy.FullSequence, Known).MonoisotopicMass,
			Is.EqualTo(new PeptideWithSetModifications(target.FullSequence, Known).MonoisotopicMass).Within(1e-6));
	}

	[Test]
	public void TerminalModificationsStayOnTheirTermini()
	{
		// mzLib #1630: GetReverseDecoyFromTarget drops the C-terminal one
		var target = Entry("[Test:Acetyl on X]PEPTIDEK-[Test:Amidation on X]");

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out _, Known).Single();

		Assert.That(decoy.FullSequence, Does.StartWith("[Test:Acetyl on X]"));
		Assert.That(decoy.FullSequence, Does.EndWith("-[Test:Amidation on X]"));
		Assert.That(decoy.PrecursorMz, Is.EqualTo(target.PrecursorMz).Within(1e-6));
	}

	[Test]
	public void FragmentsKeepTheirIntensitiesAtTheDecoysMz()
	{
		var target = Entry("PEPTIDEK");

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out _, Known).Single();

		var decoyProducts = Products(decoy.FullSequence);
		Assert.That(decoy.MatchedFragmentIons, Has.Count.EqualTo(target.MatchedFragmentIons.Count));
		for (int i = 0; i < target.MatchedFragmentIons.Count; i++)
		{
			MslFragmentIon t = target.MatchedFragmentIons[i], d = decoy.MatchedFragmentIons[i];
			Assert.That((d.ProductType, d.FragmentNumber, d.Charge, d.Intensity), Is.EqualTo((t.ProductType, t.FragmentNumber, t.Charge, t.Intensity)));
			Assert.That(d.Mz, Is.EqualTo(Single(decoyProducts, d.ProductType, d.FragmentNumber).NeutralMass.ToMz(d.Charge)).Within(1e-3));
		}
		Assert.That(decoy.MatchedFragmentIons, Is.Not.SameAs(target.MatchedFragmentIons));
		Assert.That(target.MatchedFragmentIons[0].Mz, Is.EqualTo(Entry("PEPTIDEK").MatchedFragmentIons[0].Mz), "the target's ions are not moved");
	}

	[Test]
	public void AnIonsNeutralLossCarriesIntoTheDecoy()
	{
		var target = Entry("PEPTIDEK");
		var y3 = target.MatchedFragmentIons.First(f => f.ProductType == ProductType.y && f.FragmentNumber == 3);
		y3.NeutralLoss = Water;
		y3.Mz -= (float)Water;

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out _, Known).Single();

		var decoyY3 = decoy.MatchedFragmentIons.Single(f => f.ProductType == ProductType.y && f.FragmentNumber == 3);
		Assert.That(decoyY3.NeutralLoss, Is.EqualTo(Water));
		Assert.That(decoyY3.Mz, Is.EqualTo((Single(Products(decoy.FullSequence), ProductType.y, 3).NeutralMass - Water).ToMz(1)).Within(1e-3));
	}

	[Test]
	public void AnInternalIonGetsTheDecoysInternalIonOfTheSameSpan()
	{
		var target = Entry("PEPTIDEK");
		var targetInternal = new List<Product>();
		new PeptideWithSetModifications("PEPTIDEK", Known).FragmentInternally(DissociationType.HCD, 3, targetInternal);
		var span = targetInternal.First(p => p.SecondaryProductType != null);
		target.MatchedFragmentIons.Add(new MslFragmentIon
		{
			ProductType = span.ProductType, SecondaryProductType = span.SecondaryProductType, FragmentNumber = span.FragmentNumber,
			SecondaryFragmentNumber = span.SecondaryFragmentNumber, Charge = 1, Mz = (float)span.NeutralMass.ToMz(1), Intensity = 0.5f
		});

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out var result, Known).Single();

		var decoyInternal = new List<Product>();
		new PeptideWithSetModifications(decoy.FullSequence, Known).FragmentInternally(DissociationType.HCD, 3, decoyInternal);
		var expected = decoyInternal.Single(p => p.ProductType == span.ProductType && p.SecondaryProductType == span.SecondaryProductType
			&& p.FragmentNumber == span.FragmentNumber && p.SecondaryFragmentNumber == span.SecondaryFragmentNumber);
		var internalIon = decoy.MatchedFragmentIons.Single(f => f.IsInternalFragment);
		Assert.That(internalIon.Mz, Is.EqualTo(expected.NeutralMass.ToMz(1)).Within(1e-3));
		Assert.That(result.FragmentsLeftOut, Is.Zero);
	}

	[Test]
	public void WithNoModificationsGivenMzLibsCatalogIsRead()
	{
		var target = Entry("PEPTIDEK");
		target.FullSequence = "PEPT[Common Biological:Phosphorylation on T]IDEK";
		Assume.That(Mods.AllKnownProteinModsDictionary.ContainsKey("Phosphorylation on T"), "the catalog names phosphorylation on T");

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out var result).Single();

		Assert.That(decoy.FullSequence, Does.Contain("Phosphorylation on T]"));
		Assert.That(result.TargetsWithoutDecoy, Is.Zero);
	}

	[Test]
	public void AnEntryOfUnknownDissociationTypeStillGetsAllItsIons()
	{
		// Libraries converted from .msp carry Unknown today, and PWSM.Fragment yields nothing for it
		var target = Entry("PEPTIDEK");
		target.DissociationType = DissociationType.Unknown;

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out var result, Known).Single();

		Assert.That(decoy.MatchedFragmentIons, Has.Count.EqualTo(target.MatchedFragmentIons.Count));
		Assert.That(result.FragmentsLeftOut, Is.Zero);
	}

	[Test]
	public void AnIonWithNoDecoyCounterpartIsLeftOutAndCounted()
	{
		var target = Entry("PEPTIDEK");
		target.MatchedFragmentIons.Add(Ion(ProductType.c, 3)); // HCD makes no c ions

		var decoy = MslLibraryAugmenter.CreateDecoys(new[] { target }, Trypsin, out var result, Known).Single();

		Assert.That(decoy.MatchedFragmentIons, Has.Count.EqualTo(target.MatchedFragmentIons.Count - 1));
		Assert.That(decoy.MatchedFragmentIons.Select(f => f.ProductType), Has.No.Member(ProductType.c));
		Assert.That(result.FragmentsLeftOut, Is.EqualTo(1));
	}

	[Test]
	public void TheEntrapmentFlagAndTheAccessionsCarryOver()
	{
		var targets = new List<MslLibraryEntry> { Entry("PEPTIDEK", accession: "P1|P2", entrapment: true), Entry("SAMPLERK", accession: "") };

		var decoys = MslLibraryAugmenter.CreateDecoys(targets, Trypsin, out _, Known);

		Assert.That(decoys[0].ProteinAccession, Is.EqualTo("DECOY_P1|DECOY_P2"));
		Assert.That(decoys[0].IsEntrapment, Is.True);
		Assert.That(decoys[1].ProteinAccession, Is.EqualTo(""), "a library without accessions (a DDA library's) stays without");
		Assert.That(decoys[1].IsEntrapment, Is.False);
	}

	[Test]
	public void ALibraryThatAlreadyHasDecoysIsRefused()
	{
		var decoy = Entry("EPTPIDEK");
		decoy.IsDecoy = true;

		var e = Assert.Throws<MzLibException>(() => MslLibraryAugmenter.CreateDecoys(new[] { Entry("PEPTIDEK"), decoy }, Trypsin, out _, Known));

		Assert.That(e!.Message, Does.Contain("decoy"));
	}

	[Test]
	public void AnEntryWithNoPossibleDecoyGetsNoneAndIsCounted()
	{
		// AAAAK's free residues have one arrangement, itself
		var decoys = MslLibraryAugmenter.CreateDecoys(new[] { Entry("AAAAK"), Entry("PEPTIDEK") }, Trypsin, out var result, Known);

		Assert.That(decoys.Select(d => d.BaseSequence.Length), Is.EqualTo(new[] { 8 }));
		Assert.That(result, Is.EqualTo(new MslDecoyAugmentationResult(2, 1, 1, 0)));
	}

	[Test]
	public void AnEntryNamingAnUnknownModificationGetsNoDecoyAndIsCounted()
	{
		var unknown = Entry("PEPTIDEK");
		unknown.FullSequence = "PEPT[Test:Mystery on T]IDEK";

		var decoys = MslLibraryAugmenter.CreateDecoys(new[] { unknown, Entry("SAMPLERK") }, Trypsin, out var result, Known);

		Assert.That(decoys, Has.Count.EqualTo(1));
		Assert.That(result.TargetsWithoutDecoy, Is.EqualTo(1));
	}

	[Test]
	public void ADigestionWhoseCleavageSitesCannotBeHeldIsRefused()
	{
		var trypsinP = ProteaseDictionary.Dictionary["trypsin|P"].DigestionMotifs;

		Assert.Throws<MzLibException>(() => MslLibraryAugmenter.CreateDecoys(new[] { Entry("PEPTIDEK") }, trypsinP, out _, Known));
	}

	[Test]
	public void TheFileOverloadWritesEveryTargetAndItsDecoy()
	{
		string folder = Path.Combine(TestContext.CurrentContext.WorkDirectory, "MslAugmenter_" + Path.GetRandomFileName());
		Directory.CreateDirectory(folder);
		try
		{
			string input = Path.Combine(folder, "targets.msl"), output = Path.Combine(folder, "complete.msl");
			var targets = new List<MslLibraryEntry> { Entry("PEPTIDEK", 2), Entry("PEPTIDEK", 3), Entry("SAMPLERK", 2), Entry("AAAAK", 2) };
			MslLibrary.Save(input, targets);

			var result = MslLibraryAugmenter.AddDecoys(input, output, Trypsin, Known);

			Assert.That(result, Is.EqualTo(new MslDecoyAugmentationResult(4, 3, 1, 0)));
			using var written = MslLibrary.Load(output);
			var entries = written.GetAllEntries().ToList();
			Assert.That(written.DecoyCount, Is.EqualTo(3));
			Assert.That(entries.Where(e => !e.IsDecoy).Select(e => e.Name), Is.EquivalentTo(targets.Select(t => t.Name)));
			Assert.That(entries.Where(e => e.IsDecoy).All(e => e.MatchedFragmentIons.Count > 0), Is.True);

			var refused = Assert.Throws<MzLibException>(() => MslLibraryAugmenter.AddDecoys(output, Path.Combine(folder, "again.msl"), Trypsin, Known));
			Assert.That(refused!.Message, Does.Contain("decoy"));
		}
		finally
		{
			Directory.Delete(folder, true);
		}
	}

	// ── helpers ─────────────────────────────────────────────────────────────

	private static Modification Mod(string id, string residue, string location, string formula)
	{
		ModificationMotif.TryGetMotif(residue, out var motif);
		return new Modification(_originalId: id, _modificationType: "Test", _target: motif, _locationRestriction: location,
			_chemicalFormula: ChemicalFormula.ParseFormula(formula));
	}

	/// <summary>An empirical HCD entry: every b and y ion at charge 1, intensities falling with the ion's order.</summary>
	private static MslLibraryEntry Entry(string fullSequence, int charge = 2, double irt = 35, string accession = "P1", bool entrapment = false)
	{
		var peptide = new PeptideWithSetModifications(fullSequence, Known);
		var ions = Products(fullSequence).Where(p => p.NeutralLoss == 0)
			.Select((p, i) => new MslFragmentIon
			{
				ProductType = p.ProductType, FragmentNumber = p.FragmentNumber, ResiduePosition = p.ResiduePosition, Charge = 1,
				Mz = (float)p.NeutralMass.ToMz(1), Intensity = 1f / (i + 1), NeutralLoss = 0
			}).ToList();
		return new MslLibraryEntry
		{
			FullSequence = fullSequence, BaseSequence = peptide.BaseSequence, PrecursorMz = peptide.MonoisotopicMass.ToMz(charge),
			ChargeState = charge, RetentionTime = irt, DissociationType = DissociationType.HCD, Nce = 28,
			MoleculeType = MslFormat.MoleculeType.Peptide, Source = MslFormat.SourceType.Empirical, ProteinAccession = accession,
			ProteinName = string.Empty, GeneName = string.Empty, QValue = 0.01f, IsEntrapment = entrapment, MatchedFragmentIons = ions
		};
	}

	private static List<Product> Products(string fullSequence)
	{
		var products = new List<Product>();
		new PeptideWithSetModifications(fullSequence, Known).Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
		return products;
	}

	private static Product Single(IEnumerable<Product> products, ProductType type, int number) =>
		products.Single(p => p.ProductType == type && p.FragmentNumber == number && p.NeutralLoss == 0 && p.SecondaryProductType == null);

	private static MslFragmentIon Ion(ProductType type, int number, double neutralLoss = 0) =>
		new() { ProductType = type, FragmentNumber = number, Charge = 1, Intensity = 1f, NeutralLoss = neutralLoss };
}
