using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics.Occupancy;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ModificationOccupancyCalculatorTests
{
    #region CalculateProteinLevelOccupancy Tests

    [Test]
    public void ProteinLevelWithSingleModOnSinglePeptide()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result.ContainsKey(4), Is.True);
        Assert.That(result[4].Count, Is.EqualTo(1));
        Assert.That(result[4][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[4][0].TotalCount, Is.EqualTo(1));
        Assert.That(result[4][0].CountBasedOccupancy, Is.EqualTo(1.0));
    }

    [Test]
    public void ProteinLevelWithModifiedAndUnmodifiedPeptides()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2]);

        Assert.That(result.ContainsKey(4), Is.True);
        Assert.That(result[4][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[4][0].TotalCount, Is.EqualTo(2));
        Assert.That(result[4][0].CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void ProteinLevelModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var commonMod = new Modification("Oxidation", null, "Common Variable", null, motif, "Anywhere.", null, 15.995);

        var mods = new Dictionary<int, Modification> { { 4, commonMod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Oxidation]EF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "ACD[Oxidation]EF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ProteinLevelPeptideTerminalModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("A", out var motif);
        var pepNMod = new Modification("Acetylation", null, "Biological", null, motif, "NPep", null, 42.011);

        var mods = new Dictionary<int, Modification> { { 1, pepNMod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "[Acetylation]ACDEF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "[Acetylation]ACDEF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ProteinLevelWithIntensities()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [1_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [3_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2]);

        var site = result[4][0];
        Assert.That(site.ModifiedIntensity, Is.EqualTo(1_000_000));
        Assert.That(site.TotalIntensity, Is.EqualTo(4_000_000));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.25));
    }

    [Test]
    public void ProteinLevelWithOverlappingPeptidesCoveringPosition()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        // Peptide 1: ACDEF (positions 1-5), modified at D (position 3 in protein)
        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);

        // Peptide 2: CDEFG (positions 2-6), unmodified but covers position 3
        var overlappingPeptide = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("test.raw", "CDEFG", "CDEFG", 1.0, 2, [overlappingPeptide]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm1, psm2]);

        Assert.That(result[4][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[4][0].TotalCount, Is.EqualTo(2));
        Assert.That(result[4][0].CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void ProteinLevelWithNoPeptides()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, Enumerable.Empty<ISpectralMatch>());

        Assert.That(result, Is.Empty);
    }

    /// <summary>
    /// Regression test: a single PSM whose <see cref="ISpectralMatch.GetIdentifiedBioPolymersWithSetMods"/>
    /// returns two ambiguous forms must not inflate TotalCount.
    /// Only the form whose FullSequence matches <see cref="ISpectralMatch.FullSequence"/> is counted;
    /// the alternative form is ignored entirely.
    /// </summary>
    [Test]
    public void ProteinLevel_SinglePsmTwoAmbiguousInterpretations_OccupancyIsNotInflated()
    {
        var protein = new MockBioPolymer("IVENGSEQGSYDADK", "Q6PI26");
        ModificationMotif.TryGetMotif("N", out var motif);
        var deamidation = new Modification("Deamidation on N", null, "Biological", null, motif, "Anywhere.", null, 0.984);
        var deamidatedAsp = new Modification("Deamidated asparagine on N", null, "Biological", null, motif, "Anywhere.", null, 0.984);

        var form1 = new MockBioPolymerWithSetMods(
            "IVEN", "IVEN[Deamidation on N]", protein, 1, 4,
            new Dictionary<int, Modification> { { 5, deamidation } });
        var form2 = new MockBioPolymerWithSetMods(
            "IVEN", "IVEN[Deamidated asparagine on N]", protein, 1, 4,
            new Dictionary<int, Modification> { { 5, deamidatedAsp } });

        // Single PSM: FullSequence matches form1. form2 is an alternative returned by
        // GetIdentifiedBioPolymersWithSetMods() but must not be counted.
        var psm = new MockSpectralMatch("test.raw", "IVEN[Deamidation on N]", "IVEN", 1.0, 1, [form1, form2]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(protein, [psm]);

        Assert.That(result.ContainsKey(5), Is.True);
        var modsAtSite = result[5];

        // Only the form matching PSM.FullSequence should be discovered.
        var deamSite = modsAtSite.FirstOrDefault(s => s.ModificationIdWithMotif == "Deamidation on N");
        Assert.That(deamSite, Is.Not.Null);
        Assert.That(deamSite!.TotalCount, Is.EqualTo(1), "TotalCount must be 1 (one PSM), not 2 (two forms)");
        Assert.That(deamSite.ModifiedCount, Is.EqualTo(1));
        Assert.That(deamSite.CountBasedOccupancy, Is.EqualTo(1.0));

        // The unmatched alternative form's mod must not appear.
        var deamAspSite = modsAtSite.FirstOrDefault(s => s.ModificationIdWithMotif == "Deamidated asparagine on N");
        Assert.That(deamAspSite, Is.Null, "Alternative form not matching PSM.FullSequence must be excluded");
    }

    /// <summary>
    /// Regression test: when a PSM maps to two proteins and protein B's form appears first in
    /// <see cref="ISpectralMatch.GetIdentifiedBioPolymersWithSetMods"/>, TotalCount for protein A
    /// must still be 1 — the Accession filter inside the calculator ensures the correct form is found.
    /// </summary>
    [Test]
    public void ProteinLevel_SharedPeptideTwoProteins_TotalCountIsNotUnderCounted()
    {
        var proteinA = new MockBioPolymer("ACDEFGHIK", "P00001");
        var proteinB = new MockBioPolymer("ACDEFKLMN", "P00002");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var formB = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", proteinB, 1, 5, mods);
        var formA = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", proteinA, 1, 5, mods);

        // PSM returns protein B's form first; the calculator must still resolve protein A's form.
        var psm = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [formB, formA]);

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(proteinA, [psm]);

        Assert.That(result.ContainsKey(4), Is.True);
        var site = result[4][0];
        Assert.That(site.TotalCount, Is.EqualTo(1), "TotalCount must be 1 — protein-A form must be found even when protein B's form comes first");
        Assert.That(site.ModifiedCount, Is.EqualTo(1));
        Assert.That(site.CountBasedOccupancy, Is.LessThanOrEqualTo(1.0));
    }

    #endregion

    #region CalculatePeptideLevelOccupancy Tests

    [Test]
    public void PeptideLevelOccupancyReturnedPerGroup()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm1, psm2]);

        Assert.That(result.ContainsKey(4), Is.True);
        var site = result[4][0];
        Assert.That(site.ModifiedCount, Is.EqualTo(1));
        Assert.That(site.TotalCount, Is.EqualTo(2));
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void PeptideLevelWithIntensities()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);
        var unmodifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [modifiedPeptide]);
        psm1.Intensities = [2_000_000.0];
        var psm2 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 1.0, 2, [unmodifiedPeptide]);
        psm2.Intensities = [8_000_000.0];

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm1, psm2]);

        var site = result[4][0];
        Assert.That(site.ModifiedIntensity, Is.EqualTo(2_000_000));
        Assert.That(site.TotalIntensity, Is.EqualTo(10_000_000));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0.2));
    }

    [Test]
    public void PeptideLevelCommonFixedModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("C", out var motif);
        var fixedMod = new Modification("Carbamidomethyl", null, "Common Fixed", null, motif, "Anywhere.", null, 57.021);

        var mods = new Dictionary<int, Modification> { { 3, fixedMod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "AC[Carbamidomethyl]DEF", protein, 1, 5, mods);
        var psm = new MockSpectralMatch("test.raw", "AC[Carbamidomethyl]DEF", "ACDEF", 1.0, 1, [peptide]);

        var result = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm]);

        Assert.That(result, Is.Empty);
    }

    /// <summary>
    /// The method requires all PSMs to share the same BaseSequence.
    /// Passing PSMs with different base sequences must throw <see cref="ArgumentException"/>.
    /// To calculate occupancy across multiple peptides, call the method separately per base sequence.
    /// </summary>
    [Test]
    public void PeptideLevelWithMultipleBaseSequences()
    {
        var protein = new MockBioPolymer("ACDEFGHIKLMNPQR", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var mod = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

        var mods1 = new Dictionary<int, Modification> { { 4, mod } };
        var peptide1 = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods1);

        var mods2 = new Dictionary<int, Modification> { { 3, mod } };
        var peptide2 = new MockBioPolymerWithSetMods("GHIKLM", "GH[Phosphorylation]IKLM", protein, 6, 11, mods2);

        var psm1 = new MockSpectralMatch("test.raw", "ACD[Phosphorylation]EF", "ACDEF", 1.0, 1, [peptide1]);
        var psm2 = new MockSpectralMatch("test.raw", "GH[Phosphorylation]IKLM", "GHIKLM", 1.0, 2, [peptide2]);

        // PSMs with different base sequences must not be mixed in a single call.
        Assert.That(
            () => ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([psm1, psm2]),
            Throws.ArgumentException);
    }

    #endregion

    #region SiteSpecificModificationOccupancy Tests

    [Test]
    public void ToSpectralCountModInfoStringMatchesExpectedFormat()
    {
        var site = new SiteSpecificModificationOccupancy(5, "Phosphorylation on S")
        {
            ModifiedCount = 3,
            TotalCount = 10
        };

        string expected = "pos4[Phosphorylation on S,info:fraction=0.30(3/10)]";
        Assert.That(site.ToModInfoString(intensityBased: false), Is.EqualTo(expected));
    }

    [Test]
    public void IntensityBasedStoichiometryZeroTotalIntensityDoesNotThrowDivByZero()
    {
        var site = new SiteSpecificModificationOccupancy(1, "TestMod")
        {
            ModifiedIntensity = 0,
            TotalIntensity = 0
        };

        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0));
    }

    [Test]
    public void CountBasedOccupancyZeroTotalIntensityDoesNotThrowDivByZero()
    {
        var site = new SiteSpecificModificationOccupancy(1, "TestMod")
        {
            ModifiedCount = 0,
            TotalCount = 0
        };

        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0));
    }

    #endregion
}
