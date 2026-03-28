using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics;

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

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, new[] { peptide });

        Assert.That(result.ContainsKey(3), Is.True);
        Assert.That(result[3].Count, Is.EqualTo(1));
        Assert.That(result[3][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[3][0].TotalCount, Is.EqualTo(1));
        Assert.That(result[3][0].CountBasedOccupancy, Is.EqualTo(1.0));
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

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, new IBioPolymerWithSetMods[] { modifiedPeptide, unmodifiedPeptide });

        Assert.That(result.ContainsKey(3), Is.True);
        Assert.That(result[3][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[3][0].TotalCount, Is.EqualTo(2));
        Assert.That(result[3][0].CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void ProteinLevelModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        var commonMod = new Modification("Oxidation", null, "Common Variable", null, motif, "Anywhere.", null, 15.995);

        var mods = new Dictionary<int, Modification> { { 4, commonMod } };
        var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Oxidation]EF", protein, 1, 5, mods);

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, new[] { peptide });

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

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, new[] { peptide });

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

        var intensities = new Dictionary<string, double>
        {
            ["ACD[Phosphorylation]EF"] = 1_000_000,
            ["ACDEF"] = 3_000_000
        };

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, new IBioPolymerWithSetMods[] { modifiedPeptide, unmodifiedPeptide }, null, intensities);

        var site = result[3][0];
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

        // Peptide 1: ACDEF (positions 1-5), modified at D (position 3)
        var mods = new Dictionary<int, Modification> { { 4, mod } };
        var modifiedPeptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5, mods);

        // Peptide 2: CDEFG (positions 2-6), unmodified but covers position 3
        var overlappingPeptide = new MockBioPolymerWithSetMods("CDEFG", "CDEFG", protein, 2, 6);

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, new IBioPolymerWithSetMods[] { modifiedPeptide, overlappingPeptide });

        Assert.That(result[3][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[3][0].TotalCount, Is.EqualTo(2));
        Assert.That(result[3][0].CountBasedOccupancy, Is.EqualTo(0.5));
    }

    [Test]
    public void ProteinLevelWithNoPeptides()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, Enumerable.Empty<IBioPolymerWithSetMods>());

        Assert.That(result, Is.Empty);
    }

    /// <summary>
    /// Regression test for the "1 PSM → occupancy 1/2" bug.
    ///
    /// When a single PSM has two ambiguous interpretations of the same peptide
    /// (e.g. "Deamidation on N" vs "Deamidated asparagine on N" at the same site),
    /// PopulateOccupancy previously expanded them via SelectMany, causing TotalCount = 2
    /// and occupancy = 0.50 (1/2) instead of 1.0 (1/1).
    ///
    /// The fix: callers pass a deduplicated <paramref name="sequencesForTotalCount"/> list
    /// (one entry per PSM) separately from the full <paramref name="localizedSequences"/>
    /// list (all forms, for ModifiedCount).  Both modifications should show occupancy = 1.0.
    /// </summary>
    [Test]
    public void ProteinLevel_SinglePsmTwoAmbiguousInterpretations_OccupancyIsNotInflated()
    {
        var protein = new MockBioPolymer("IVENGSEQGSYDADK", "Q6PI26");
        ModificationMotif.TryGetMotif("N", out var motif);
        var deamidation = new Modification("Deamidation on N", null, "Biological", null, motif, "Anywhere.", null, 0.984);
        var deamidatedAsp = new Modification("Deamidated asparagine on N", null, "Biological", null, motif, "Anywhere.", null, 0.984);

        // Two interpretation forms from a single PSM: same base sequence, different mod identity
        var form1 = new MockBioPolymerWithSetMods(
            "IVEN", "IVEN[Deamidation on N]", protein, 1, 4,
            new Dictionary<int, Modification> { { 5, deamidation } });
        var form2 = new MockBioPolymerWithSetMods(
            "IVEN", "IVEN[Deamidated asparagine on N]", protein, 1, 4,
            new Dictionary<int, Modification> { { 5, deamidatedAsp } });

        // allSequences = both forms (used for ModifiedCount numerator)
        IBioPolymerWithSetMods[] allSequences = [form1, form2];
        // coverageSequences = one representative per PSM (used for TotalCount denominator)
        IBioPolymerWithSetMods[] coverageSequences = [form1];

        var result = ModificationOccupancyCalculator.CalculateProteinLevelOccupancy(
            protein, allSequences, coverageSequences);

        Assert.That(result.ContainsKey(4), Is.True, "Expected occupancy data at protein position 4 (N)");
        var modsAtSite = result[4];

        var deamSite = modsAtSite.FirstOrDefault(s => s.ModificationIdWithMotif == "Deamidation on N");
        var deamAspSite = modsAtSite.FirstOrDefault(s => s.ModificationIdWithMotif == "Deamidated asparagine on N");

        Assert.That(deamSite, Is.Not.Null, "Deamidation on N should be present");
        Assert.That(deamSite!.TotalCount, Is.EqualTo(1), "TotalCount must be 1 (one PSM), not 2");
        Assert.That(deamSite.ModifiedCount, Is.EqualTo(1));
        Assert.That(deamSite.CountBasedOccupancy, Is.EqualTo(1.0), "Occupancy must be 1/1 = 100%, not 1/2 = 50%");

        Assert.That(deamAspSite, Is.Not.Null, "Deamidated asparagine on N should be present");
        Assert.That(deamAspSite!.TotalCount, Is.EqualTo(1), "TotalCount must be 1 (one PSM), not 2");
        Assert.That(deamAspSite.ModifiedCount, Is.EqualTo(1));
        Assert.That(deamAspSite.CountBasedOccupancy, Is.EqualTo(1.0), "Occupancy must be 1/1 = 100%, not 1/2 = 50%");
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

        var result = ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy(
            new IBioPolymerWithSetMods[] { modifiedPeptide, unmodifiedPeptide });

        Assert.That(result.ContainsKey(4), Is.True); // peptide-local position (AllModsOneIsNterminus key)
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

        var intensities = new Dictionary<string, double>
        {
            ["ACD[Phosphorylation]EF"] = 2_000_000,
            ["ACDEF"] = 8_000_000
        };

        var result = ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy(
            new IBioPolymerWithSetMods[] { modifiedPeptide, unmodifiedPeptide }, intensities);

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

        var result = ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy(new[] { peptide });

        Assert.That(result, Is.Empty);
    }

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

        // New signature operates on a single base sequence group; test each separately
        var result1 = ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy(new[] { peptide1 });
        var result2 = ModificationOccupancyCalculator.CalculatePeptideLevelOccupancy(new[] { peptide2 });

        Assert.That(result1.Count, Is.EqualTo(1));
        Assert.That(result2.Count, Is.EqualTo(1));
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

        string expected = "#aa5[Phosphorylation on S,info:occupancy=0.30(3/10)]";
        Assert.That(site.ToSpectralCountModInfoString(), Is.EqualTo(expected));
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
