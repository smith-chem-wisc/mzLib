using Easy.Common.Extensions;
using Newtonsoft.Json.Bson;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics;

/// <summary>
/// Educational unit tests for understanding how PTM occupancy is calculated
/// in ModificationOccupancyCalculator.
///
/// KEY CONCEPTS:
/// =============
/// PTM occupancy answers the question: "At a given amino acid position, what fraction
/// of the observed peptides carry a specific modification?"
///
/// Two metrics are computed:
///   1. Count-Based Occupancy = ModifiedCount / TotalCount
///      - ModifiedCount: number of PSMs carrying this mod at this position
///      - TotalCount: total PSMs covering this position (modified + unmodified)
///
///   2. Intensity-Based Stoichiometry = ModifiedIntensity / TotalIntensity
///      - ModifiedIntensity: sum of intensities from PSMs with the mod at this position
///      - TotalIntensity: sum of intensities from ALL PSMs covering this position
///
/// POSITION MAPPING (AllModsOneIsNterminus convention):
///   - Key 1 = N-terminal modification slot           → result position 1
///   - Key 2 = first amino acid residue               → result position 2 (for peptide at protein pos 1)
///   - Key (n+1) = nth amino acid residue
///   - For "Anywhere." mods, result position = OneBasedStartResidue + key - 1
///   - For "N-terminal." mods, result position = 1 (always)
///   - For "C-terminal." mods, result position = bioPolymerLength + 2 (always)
///
/// IMPORTANT: The calculator only reports positions where a modification EXISTS.
/// Unmodified positions produce no entries in the result dictionary.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PtmOccupancyLearningTests
{
    // ========================================================================
    // HELPER: Creates a Modification with "Anywhere." location restriction
    // ========================================================================
    private static Modification CreateMod(string name, string motifChar)
    {
        ModificationMotif.TryGetMotif(motifChar, out var motif);
        return new Modification(name, null, "Biological", null, motif, "Anywhere.", null, 79.966);
    }

    #region Test 1: Single unmodified peptide — no occupancy to report

    /// <summary>
    /// TEST 1: One protein, one peptide (whole protein sequence), completely unmodified, 1 PSM.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK  (9 amino acids)
    ///   Peptide:  ACDEFGHIK  (spans the entire protein, positions 1–9)
    ///   PSMs:     1 unmodified PSM
    ///
    /// EXPECTED RESULT:
    ///   The result dictionary is EMPTY. The calculator only creates entries at positions
    ///   where a modification is observed. Since this peptide has no modifications,
    ///   there is nothing to report — the occupancy of any hypothetical PTM at any
    ///   position is implicitly 0/1 = 0%, but this is not explicitly stored.
    ///
    ///   This is a fundamental design choice: the calculator answers "what is the
    ///   occupancy of modifications that WERE observed?" not "what is the occupancy
    ///   of modifications that COULD exist?"
    /// </summary>
    [Test]
    public void Test1_SingleUnmodifiedPeptide_NoOccupancyReported()
    {
        // Arrange: one protein, one unmodified peptide spanning the whole protein
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var unmodifiedPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK",    // base sequence
            "ACDEFGHIK",    // full sequence (no modification brackets)
            protein,        // parent protein
            1, 9);          // spans positions 1 through 9

        // The API expects spectral matches as input, so we create a mock PSM that identifies this unmodified peptide.
        var psm = new MockSpectralMatch("test.mz", unmodifiedPeptide.FullSequence, unmodifiedPeptide.BaseSequence, 1, 1, [unmodifiedPeptide]);

        // Set the intensity for this PSM (arbitrary value since it won't affect the result)
        psm.Intensities = new double[] { 1000.0 };

        // Act: calculate protein-level occupancy
        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, new[] { psm });

        // Assert: result is empty because no modifications were observed
        // The calculator does not create entries for unmodified positions.
        // Even though this PSM "covers" all 9 positions, there are no modifications
        // at any position, so there is nothing to report.
        Assert.That(result, Is.Empty,
            "No modifications exist, so the result dictionary should be empty. " +
            "PTM occupancy is only reported at positions where a modification was actually observed.");
    }

    #endregion

    #region Test 2: One modified + one unmodified PSM at a single position

    /// <summary>
    /// TEST 2: One protein, one peptide (whole protein), sometimes modified, sometimes not.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK
    ///   PSM 1:    ACDEFGHIK         (unmodified, intensity = 1)
    ///   PSM 2:    ACD[Phospho]EFGHIK (Phosphorylation on D at protein position 4, intensity = 2)
    ///
    /// This tests the core occupancy calculation: of the 2 PSMs covering position 4,
    /// only 1 carries the modification.
    ///
    /// OCCUPANCY AT THE MODIFIED POSITION (D, protein position 4):
    ///   Count-Based:     ModifiedCount=1, TotalCount=2  → 1/2 = 0.50 (50%)
    ///   Intensity-Based: ModifiedIntensity=2, TotalIntensity=3 → 2/3 ≈ 0.667 (66.7%)
    ///
    ///   Notice the two metrics give DIFFERENT answers because the modified PSM
    ///   has higher intensity (2) than the unmodified (1). Intensity-based stoichiometry
    ///   weights each PSM by its signal strength.
    ///
    /// OCCUPANCY AT ANY OTHER POSITION (e.g., A at position 2):
    ///   Not reported — the calculator only tracks positions where mods exist.
    /// </summary>
    [Test]
    public void Test2_OneModOneUnmod_OccupancyAtModifiedAndUnmodifiedSites()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoOnD = CreateMod("Phosphorylation", "D");

        // PSM 1: unmodified peptide, intensity = 1
        var unmodifiedPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein, 1, 9);
        var unmodifiedPsm = new MockSpectralMatch("test.mz", unmodifiedPeptide.FullSequence, unmodifiedPeptide.BaseSequence, 1, 1, [unmodifiedPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        // PSM 2: Phosphorylation on D (3rd residue → AllModsOneIsNterminus key = 4), intensity = 2
        // Key 4 maps to protein position: 1 + 4 - 1 = 4
        var modifiedPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoOnD } });
        var modifiedPsm = new MockSpectralMatch("test.mz", modifiedPeptide.FullSequence, modifiedPeptide.BaseSequence, 1, 1, [modifiedPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { unmodifiedPsm, modifiedPsm });

        // --- Occupancy at the MODIFIED position (D, protein position 4) ---
        // Both PSMs cover position 4, but only 1 carries the phosphorylation.
        Assert.That(result.ContainsKey(4), Is.True,
            "Position 4 (D) should have occupancy data because a modification was observed there.");

        var siteD = result[4][0];

        // Count-based: 1 modified out of 2 total = 50%
        Assert.That(siteD.ModifiedCount, Is.EqualTo(1),
            "Only 1 of the 2 PSMs carries Phosphorylation at position 4.");
        Assert.That(siteD.TotalCount, Is.EqualTo(2),
            "Both PSMs (modified and unmodified) cover position 3, so TotalCount = 2.");
        Assert.That(siteD.CountBasedOccupancy, Is.EqualTo(0.5),
            "Count-based occupancy = 1/2 = 0.50. Half the PSMs are modified at this site.");

        // Intensity-based: modified intensity = 2, total intensity = 1 + 2 = 3
        Assert.That(siteD.ModifiedIntensity, Is.EqualTo(2.0),
            "The modified PSM has intensity 2.");
        Assert.That(siteD.TotalIntensity, Is.EqualTo(3.0),
            "Total intensity = 1 (unmodified) + 2 (modified) = 3.");
        Assert.That(siteD.IntensityBasedStoichiometry, Is.EqualTo(2.0 / 3.0).Within(1e-10),
            "Intensity-based stoichiometry = 2/3 ≈ 0.667. Higher than count-based because " +
            "the modified PSM has higher intensity than the unmodified one.");

        // --- Occupancy at an UNMODIFIED position (e.g., position 2, A) ---
        // No modification was observed at position 2, so the calculator does not report it.
        Assert.That(result.ContainsKey(2), Is.False,
            "Position 2 (A) has no modification, so it does not appear in the result. " +
            "The calculator only tracks positions where modifications were observed.");
    }

    #endregion

    #region Test 3: Modifications at two different positions

    /// <summary>
    /// TEST 3: One protein, one peptide, modifications at two separate positions.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK
    ///   PSM 1:    ACDEFGHIK            (unmodified, intensity = 1)
    ///   PSM 2:    ACD[Phospho]EFGHIK   (Phospho on D at position 4, intensity = 2)
    ///   PSM 3:    ACDEFG[Phospho]HIK   (Phospho on G at position 7, intensity = 3)
    ///
    /// Each PSM represents a different observation from mass spec. All 3 PSMs cover
    /// ALL positions in the protein because they all span the full sequence.
    ///
    /// AT POSITION 4 (D, Phosphorylation):
    ///   Count:     1 modified / 3 total = 0.333 (33.3%)
    ///   Intensity: 2 / (1+2+3) = 2/6 = 0.333 (33.3%)
    ///
    /// AT POSITION 7 (G, Phosphorylation):
    ///   Count:     1 modified / 3 total = 0.333 (33.3%)
    ///   Intensity: 3 / (1+2+3) = 3/6 = 0.500 (50.0%)
    ///
    /// KEY INSIGHT: The count-based occupancy is the same at both sites (1/3),
    /// but intensity-based stoichiometry differs because the PSM modified at G
    /// has higher intensity (3) than the PSM modified at D (2). This shows how
    /// intensity weighting can reveal that one modification site may be more
    /// abundantly occupied than another, even when the same number of PSMs
    /// carry each modification.
    ///
    /// AT AN UNMODIFIED POSITION (e.g., position 2):
    ///   Not reported — no modification was observed there.
    /// </summary>
    [Test]
    public void Test3_TwoModificationsAtDifferentPositions()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoD = CreateMod("Phosphorylation", "D");
        var phosphoG = CreateMod("Phosphorylation", "G");

        // PSM 1: unmodified, intensity = 1
        var unmodPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein, 1, 9);
        var unmodPsm = new MockSpectralMatch("test.mz", unmodPeptide.FullSequence, unmodPeptide.BaseSequence, 1, 1, [unmodPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        // PSM 2: Phospho on D (key 4 → protein position 1+4-1=4), intensity = 2
        var modDPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var modDPsm = new MockSpectralMatch("test.mz", modDPeptide.FullSequence, modDPeptide.BaseSequence, 1, 1, [modDPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        // PSM 3: Phospho on G (key 7 → protein position 1+7-1=7), intensity = 3
        var modGPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFG[Phosphorylation]HIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 7, phosphoG } });
        var modGPsm = new MockSpectralMatch("test.mz", modGPeptide.FullSequence, modGPeptide.BaseSequence, 1, 1, [modGPeptide])
        {
            Intensities = new double[] { 3.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { unmodPsm, modDPsm, modGPsm });

        // --- Position 4 (D): Phosphorylation ---
        // All 3 PSMs cover this position. Only 1 carries Phospho here.
        Assert.That(result.ContainsKey(4), Is.True);
        var siteD = result[4][0];
        Assert.That(siteD.ModifiedCount, Is.EqualTo(1),
            "Only 1 PSM has Phosphorylation at position 4 (D).");
        Assert.That(siteD.TotalCount, Is.EqualTo(3),
            "All 3 PSMs span the full protein and cover position 4.");
        Assert.That(siteD.CountBasedOccupancy, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Count occupancy at D = 1/3 ≈ 33.3%.");
        Assert.That(siteD.ModifiedIntensity, Is.EqualTo(2.0),
            "The PSM modified at D has intensity 2.");
        Assert.That(siteD.TotalIntensity, Is.EqualTo(6.0),
            "Total intensity = 1 + 2 + 3 = 6 (all PSMs covering this position).");
        Assert.That(siteD.IntensityBasedStoichiometry, Is.EqualTo(2.0 / 6.0).Within(1e-10),
            "Intensity stoichiometry at D = 2/6 ≈ 33.3%. Same as count here by coincidence.");

        // --- Position 7 (G): Phosphorylation ---
        // All 3 PSMs cover this position. Only 1 carries Phospho here.
        Assert.That(result.ContainsKey(7), Is.True);
        var siteG = result[7][0];
        Assert.That(siteG.ModifiedCount, Is.EqualTo(1),
            "Only 1 PSM has Phosphorylation at position 7 (G).");
        Assert.That(siteG.TotalCount, Is.EqualTo(3),
            "All 3 PSMs cover position 7.");
        Assert.That(siteG.CountBasedOccupancy, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Count occupancy at G = 1/3. Same as D because each site has exactly 1 modified PSM out of 3.");
        Assert.That(siteG.ModifiedIntensity, Is.EqualTo(3.0),
            "The PSM modified at G has intensity 3.");
        Assert.That(siteG.TotalIntensity, Is.EqualTo(6.0),
            "Total intensity is 6 (same denominator as D — all PSMs cover all positions).");
        Assert.That(siteG.IntensityBasedStoichiometry, Is.EqualTo(3.0 / 6.0).Within(1e-10),
            "Intensity stoichiometry at G = 3/6 = 50%. HIGHER than D's 33.3% because the " +
            "PSM modified at G has higher intensity (3 vs 2). This demonstrates how " +
            "intensity-based stoichiometry can differentiate site occupancy even when " +
            "count-based occupancy is the same.");

        // --- Unmodified position (e.g., position 2, A) ---
        Assert.That(result.ContainsKey(2), Is.False,
            "Position 2 (A) has no modification observed, so it's not in the result.");

        // Only 2 positions are in the result: 3 and 6 (the two modified sites)
        Assert.That(result.Count, Is.EqualTo(2),
            "Only the 2 modified positions appear in the result dictionary.");
    }

    #endregion

    #region Test 4: Two peptides (full + half length), both unmodified

    /// <summary>
    /// TEST 4: Two peptides of different lengths, both unmodified.
    ///
    /// SCENARIO:
    ///   Protein:       ACDEFGHIK  (positions 1–9)
    ///   Long peptide:  ACDEFGHIK  (positions 1–9, intensity = 1)
    ///   Short peptide: ACDEF      (positions 1–5, intensity = 2)
    ///
    ///   Shared positions: 1–5 (covered by BOTH peptides)
    ///   Non-shared positions: 6–9 (covered ONLY by the long peptide)
    ///
    /// EXPECTED RESULT:
    ///   Empty — neither peptide has modifications, so there is nothing to report.
    ///
    ///   Even though position 3 is covered by 2 PSMs and position 7 by only 1 PSM,
    ///   the calculator does not create entries for unmodified positions. The "coverage"
    ///   only matters as a denominator when there IS a modification to report.
    ///
    /// This test establishes the baseline for Tests 5a and 5b, which add modifications
    /// to these same peptides.
    /// </summary>
    [Test]
    public void Test4_TwoPeptidesBothUnmodified_NoOccupancyReported()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");

        // Long peptide: full protein, positions 1–9, intensity = 1
        var longPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein, 1, 9);
        var longPsm = new MockSpectralMatch("test.mz", longPeptide.FullSequence, longPeptide.BaseSequence, 1, 1, [longPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        // Short peptide: first half, positions 1–5, intensity = 2
        var shortPeptide = new MockBioPolymerWithSetMods(
            "ACDEF", "ACDEF", protein, 1, 5);
        var shortPsm = new MockSpectralMatch("test.mz", shortPeptide.FullSequence, shortPeptide.BaseSequence, 1, 1, [shortPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { longPsm, shortPsm });

        // No modifications on either peptide → empty result
        Assert.That(result, Is.Empty,
            "Both peptides are unmodified. The calculator only reports positions with " +
            "observed modifications. Coverage information (2 PSMs at positions 1-5, " +
            "1 PSM at positions 6-9) is not stored unless a modification triggers it.");
    }

    #endregion

    #region Test 5a: Overlapping peptides, modification at a SHARED position

    /// <summary>
    /// TEST 5a: Long peptide modified at a position covered by both peptides.
    ///
    /// SCENARIO:
    ///   Protein:       ACDEFGHIK  (positions 1–9)
    ///   Long peptide:  ACD[Phospho]EFGHIK  (positions 1–9, mod at D = position 4, intensity = 1)
    ///   Short peptide: ACDEF               (positions 1–5, unmodified, intensity = 2)
    ///
    ///   Position 4 (D) is SHARED — both peptides cover it.
    ///   The short peptide does not carry the modification at position 4.
    ///
    /// AT THE MODIFIED POSITION (D, position 4 — shared by both peptides):
    ///   TotalCount = 2 (both peptides cover position 4)
    ///   ModifiedCount = 1 (only the long peptide has Phospho at D)
    ///   Count Occupancy = 1/2 = 0.50
    ///
    ///   TotalIntensity = 1 + 2 = 3 (intensities of ALL peptides covering position 4)
    ///   ModifiedIntensity = 1 (only the long peptide's intensity counts as modified)
    ///   Intensity Stoichiometry = 1/3 ≈ 0.333
    ///
    ///   KEY INSIGHT: The short peptide acts as evidence AGAINST the modification.
    ///   It covers position 4 but does NOT carry Phospho, so it increases the denominator
    ///   (TotalCount and TotalIntensity) without increasing the numerator. This pulls
    ///   the occupancy DOWN from what it would be if only the long peptide were observed.
    ///
    /// AT AN UNMODIFIED SHARED POSITION (e.g., position 2):
    ///   Not reported — no modification observed there.
    ///
    /// AT A NON-SHARED POSITION (e.g., position 7 — only long peptide covers it):
    ///   Not reported — no modification observed there.
    /// </summary>
    [Test]
    public void Test5a_ModificationAtSharedPosition_BothPeptidesContributeToDenominator()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoD = CreateMod("Phosphorylation", "D");

        // Long peptide: full protein, Phospho at D (key=4 → protein pos 1+4-1=4), intensity=1
        var longPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var longPsm = new MockSpectralMatch("test.mz", longPeptide.FullSequence, longPeptide.BaseSequence, 1, 1, [longPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        // Short peptide: first half, unmodified, intensity=2
        var shortPeptide = new MockBioPolymerWithSetMods(
            "ACDEF", "ACDEF", protein, 1, 5);
        var shortPsm = new MockSpectralMatch("test.mz", shortPeptide.FullSequence, shortPeptide.BaseSequence, 1, 1, [shortPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { longPsm, shortPsm });

        // --- Modified position (D, position 4) — SHARED by both peptides ---
        Assert.That(result.ContainsKey(4), Is.True);
        var site = result[4][0];

        // Both peptides cover position 4, so TotalCount = 2
        Assert.That(site.TotalCount, Is.EqualTo(2),
            "Both the long peptide (1-9) and short peptide (1-5) cover position 4, so TotalCount = 2.");
        Assert.That(site.ModifiedCount, Is.EqualTo(1),
            "Only the long peptide carries Phosphorylation at position 4.");
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0.5),
            "Count occupancy = 1/2 = 50%. The short unmodified peptide dilutes the occupancy.");

        // Intensity: total = 1 (long) + 2 (short) = 3
        Assert.That(site.TotalIntensity, Is.EqualTo(3.0),
            "Both peptides contribute intensity to the denominator: 1 + 2 = 3.");
        Assert.That(site.ModifiedIntensity, Is.EqualTo(1.0),
            "Only the long peptide's intensity (1) counts as modified.");
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Intensity stoichiometry = 1/3 ≈ 33.3%. Lower than count-based 50% because " +
            "the unmodified short peptide has higher intensity (2) than the modified long peptide (1).");

        // --- Unmodified shared position (e.g., position 2) ---
        Assert.That(result.ContainsKey(2), Is.False,
            "Position 2 has no modification, so it's not reported.");

        // --- Non-shared position (e.g., position 7) ---
        Assert.That(result.ContainsKey(7), Is.False,
            "Position 7 is only covered by the long peptide, but it has no modification there.");
    }

    #endregion

    #region Test 5b: Overlapping peptides, modification at a NON-SHARED position

    /// <summary>
    /// TEST 5b: Long peptide modified at a position NOT covered by the short peptide.
    ///
    /// SCENARIO:
    ///   Protein:       ACDEFGHIK  (positions 1–9)
    ///   Long peptide:  ACDEFGH[Phospho]IK  (positions 1–9, mod at H = position 8, intensity = 1)
    ///   Short peptide: ACDEF                (positions 1–5, unmodified, intensity = 2)
    ///
    ///   Position 8 (H) is NOT SHARED — only the long peptide covers it.
    ///   The short peptide ends at position 5 and cannot contribute evidence at position 8.
    ///
    /// AT THE MODIFIED POSITION (H, position 8 — NOT shared):
    ///   TotalCount = 1 (only long peptide covers position 8)
    ///   ModifiedCount = 1
    ///   Count Occupancy = 1/1 = 1.00 (100%)
    ///
    ///   TotalIntensity = 1 (only long peptide's intensity)
    ///   ModifiedIntensity = 1
    ///   Intensity Stoichiometry = 1/1 = 1.00 (100%)
    ///
    ///   KEY INSIGHT: The short peptide cannot dilute the occupancy here because it
    ///   doesn't cover position 8. Contrast this with Test 5a where the short peptide
    ///   DID cover the modified position and reduced occupancy to 50%. This shows
    ///   how peptide coverage geometry affects occupancy calculations.
    ///
    /// AT SHARED UNMODIFIED POSITIONS (1–5):
    ///   Not reported — no modification there.
    /// </summary>
    [Test]
    public void Test5b_ModificationAtNonSharedPosition_OnlyLongPeptideContributes()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoH = CreateMod("Phosphorylation", "H");

        // Long peptide: Phospho at H (key=8 → protein pos 1+8-1=8), intensity=1
        var longPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGH[Phosphorylation]IK", protein, 1, 9,
            new Dictionary<int, Modification> { { 8, phosphoH } });
        var longPsm = new MockSpectralMatch("test.mz", longPeptide.FullSequence, longPeptide.BaseSequence, 1, 1, [longPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        // Short peptide: positions 1–5, unmodified, intensity=2
        var shortPeptide = new MockBioPolymerWithSetMods(
            "ACDEF", "ACDEF", protein, 1, 5);
        var shortPsm = new MockSpectralMatch("test.mz", shortPeptide.FullSequence, shortPeptide.BaseSequence, 1, 1, [shortPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { longPsm, shortPsm });

        // --- Modified position (H, position 8) — only long peptide covers it ---
        Assert.That(result.ContainsKey(8), Is.True);
        var site = result[8][0];

        Assert.That(site.TotalCount, Is.EqualTo(1),
            "Only the long peptide covers position 8. The short peptide (1-5) does NOT reach position 8.");
        Assert.That(site.ModifiedCount, Is.EqualTo(1),
            "The long peptide carries Phospho at position 8.");
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(1.0),
            "Count occupancy = 1/1 = 100%. Compare to Test 5a where sharing diluted it to 50%.");

        Assert.That(site.TotalIntensity, Is.EqualTo(1.0),
            "Only the long peptide's intensity counts — short peptide doesn't cover this position.");
        Assert.That(site.ModifiedIntensity, Is.EqualTo(1.0));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(1.0),
            "Intensity stoichiometry = 1/1 = 100%. The short peptide's intensity (2) " +
            "is NOT included because it doesn't cover position 8.");

        // --- Shared unmodified positions (1–5): not reported ---
        Assert.That(result.ContainsKey(4), Is.False,
            "Position 4 is covered by both peptides but has no modification.");

        Assert.That(result.Count, Is.EqualTo(1),
            "Only the modified position appears in the result.");
    }

    #endregion

    #region Test 6: Two proteins, identical sequences, shared unmodified peptide

    /// <summary>
    /// TEST 6: Two proteins with identical sequences share one unmodified peptide.
    ///
    /// SCENARIO:
    ///   Protein 1: ACDEFGHIK  (accession P1)
    ///   Protein 2: ACDEFGHIK  (accession P2)
    ///   1 PSM:     ACDEFGHIK  (unmodified, intensity = 1)
    ///
    ///   The peptide maps to both proteins (shared/ambiguous peptide).
    ///   In the real software, PopulateOccupancy filters peptides by Parent.Accession,
    ///   so each protein's occupancy is calculated independently with only its own peptides.
    ///
    /// HOW OCCUPANCY IS DISTRIBUTED:
    ///   Both proteins get empty results — the peptide is unmodified.
    ///
    ///   The key point about shared peptides is that each protein gets its own copy
    ///   of the peptide in the occupancy calculation. But since there are no modifications,
    ///   there's nothing to distribute.
    /// </summary>
    [Test]
    public void Test6_TwoProteinsSharedUnmodifiedPeptide_BothEmpty()
    {
        var protein1 = new MockBioPolymer("ACDEFGHIK", "P00001");
        var protein2 = new MockBioPolymer("ACDEFGHIK", "P00002");

        // The same PSM maps to both proteins → create one peptide form per protein
        var peptideForP1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein1, 1, 9);
        var peptideForP2 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein2, 1, 9);

        var psm = new MockSpectralMatch("test.mz", peptideForP1.FullSequence, peptideForP1.BaseSequence, 1, 1, [peptideForP1, peptideForP2])
        {
            Intensities = new double[] { 1.0 }
        };

        // Calculate occupancy for each protein independently
        // (mimicking how PopulateOccupancy filters by Parent.Accession)
        var resultP1 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein1, new[] { psm });

        var resultP2 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein2, new[] { psm });

        // Both proteins get empty results — no modifications to report
        Assert.That(resultP1, Is.Empty,
            "Protein 1 has no modifications from this unmodified shared peptide.");
        Assert.That(resultP2, Is.Empty,
            "Protein 2 has no modifications from this unmodified shared peptide. " +
            "For shared peptides, occupancy is calculated per-protein but since the " +
            "peptide is unmodified, both proteins show nothing.");
    }

    #endregion

    #region Test 7: Two proteins, identical sequences, shared modified peptide

    /// <summary>
    /// TEST 7: Two proteins with identical sequences share one MODIFIED peptide.
    ///
    /// SCENARIO:
    ///   Protein 1: ACDEFGHIK  (accession P1)
    ///   Protein 2: ACDEFGHIK  (accession P2)
    ///   1 PSM:     ACD[Phospho]EFGHIK  (modified at D, position 4, intensity = 1)
    ///
    ///   The PSM maps to both proteins. Each protein gets its own copy of the
    ///   modified peptide for its occupancy calculation.
    ///
    /// HOW OCCUPANCY IS DISTRIBUTED:
    ///   Each protein independently shows:
    ///     Count Occupancy = 1/1 = 100%
    ///     Intensity Stoichiometry = 1/1 = 100%
    ///
    ///   Both proteins show identical, full occupancy. This is because from each
    ///   protein's perspective, the ONLY peptide covering it is modified. The occupancy
    ///   is not "split" between proteins — each protein gets the full 100%.
    ///
    ///   This makes biological sense: if the only evidence you have for a protein
    ///   is a modified peptide, then 100% of the observed evidence is modified.
    /// </summary>
    [Test]
    public void Test7_TwoProteinsSharedModifiedPeptide_BothShow100Percent()
    {
        var protein1 = new MockBioPolymer("ACDEFGHIK", "P00001");
        var protein2 = new MockBioPolymer("ACDEFGHIK", "P00002");
        var phosphoD = CreateMod("Phosphorylation", "D");

        // One modified PSM → one peptide form per protein
        var modPeptideP1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein1, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var modPeptideP2 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein2, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });

        var psm = new MockSpectralMatch("test.mz", modPeptideP1.FullSequence, modPeptideP1.BaseSequence, 1, 1, [modPeptideP1, modPeptideP2])
        {
            Intensities = new double[] { 1.0 }
        };

        // Calculate independently for each protein
        var resultP1 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein1, new[] { psm });
        var resultP2 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein2, new[] { psm });

        // Protein 1 at position 4
        Assert.That(resultP1.ContainsKey(4), Is.True);
        Assert.That(resultP1[4][0].CountBasedOccupancy, Is.EqualTo(1.0),
            "Protein 1: 1 modified PSM / 1 total PSM = 100% occupancy.");
        Assert.That(resultP1[4][0].IntensityBasedStoichiometry, Is.EqualTo(1.0),
            "Protein 1: intensity 1 / total intensity 1 = 100%.");

        // Protein 2 at position 4 — SAME result
        Assert.That(resultP2.ContainsKey(4), Is.True);
        Assert.That(resultP2[4][0].CountBasedOccupancy, Is.EqualTo(1.0),
            "Protein 2: also 100%. Occupancy is NOT split between proteins. " +
            "Each protein independently sees 100% of its evidence as modified.");
        Assert.That(resultP2[4][0].IntensityBasedStoichiometry, Is.EqualTo(1.0),
            "Protein 2: intensity stoichiometry also 100%.");

        // CONCERN: Both proteins report 100% occupancy from a single shared PSM, but the
        // modification physically exists on only ONE protein molecule. The calculator does
        // not apportion shared peptide evidence between proteins — it duplicates it. A consumer
        // summing occupancy across proteins in a group could overcount the total modification
        // burden. For example, if Protein 1 and Protein 2 are in the same protein group, a
        // naive sum would suggest 200% total modification, which is physically impossible.
        // Whether this is a problem depends on how downstream code consumes these values.
    }

    #endregion

    #region Test 8: Two proteins, shared peptide, modified + unmodified PSMs

    /// <summary>
    /// TEST 8: Two proteins with identical sequences. Both a modified and unmodified PSM are observed.
    ///
    /// SCENARIO:
    ///   Protein 1: ACDEFGHIK  (accession P1)
    ///   Protein 2: ACDEFGHIK  (accession P2)
    ///   PSM 1:     ACD[Phospho]EFGHIK  (modified at D, intensity = 1)
    ///   PSM 2:     ACDEFGHIK           (unmodified, intensity = 2)
    ///
    ///   Both PSMs map to both proteins (shared peptides).
    ///
    /// HOW OCCUPANCY IS DISTRIBUTED:
    ///   Each protein receives BOTH PSMs for its calculation. The result is identical
    ///   for both proteins:
    ///     Count:     1 modified / 2 total = 50%
    ///     Intensity: 1 / (1+2) = 1/3 ≈ 33.3%
    ///
    ///   The occupancy is NOT split or halved between proteins. Each protein independently
    ///   sees the same 2 PSMs and computes the same occupancy. This means if a peptide
    ///   is shared between N proteins, all N proteins report the same occupancy values.
    /// </summary>
    [Test]
    public void Test8_TwoProteinsSharedPeptide_ModifiedAndUnmodified()
    {
        var protein1 = new MockBioPolymer("ACDEFGHIK", "P00001");
        var protein2 = new MockBioPolymer("ACDEFGHIK", "P00002");
        var phosphoD = CreateMod("Phosphorylation", "D");

        // Modified PSM → one peptide form per protein
        var modP1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein1, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var modP2 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein2, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var modPsm = new MockSpectralMatch("test.mz", modP1.FullSequence, modP1.BaseSequence, 1, 1, [modP1, modP2])
        {
            Intensities = new double[] { 1.0 }
        };

        // Unmodified PSM → one peptide form per protein
        var unmodP1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein1, 1, 9);
        var unmodP2 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein2, 1, 9);
        var unmodPsm = new MockSpectralMatch("test.mz", unmodP1.FullSequence, unmodP1.BaseSequence, 1, 1, [unmodP1, unmodP2])
        {
            Intensities = new double[] { 2.0 }
        };

        // Protein 1: receives modP1 + unmodP1
        var resultP1 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein1,
            new[] { modPsm, unmodPsm});

        // Protein 2: receives modP2 + unmodP2
        var resultP2 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein2,
            new[] { modPsm, unmodPsm });

        // --- Protein 1 ---
        var siteP1 = resultP1[4][0];
        Assert.That(siteP1.ModifiedCount, Is.EqualTo(1));
        Assert.That(siteP1.TotalCount, Is.EqualTo(2));
        Assert.That(siteP1.CountBasedOccupancy, Is.EqualTo(0.5),
            "Protein 1: 1 modified / 2 total = 50%.");
        Assert.That(siteP1.IntensityBasedStoichiometry, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Protein 1: 1/(1+2) ≈ 33.3%. Lower than 50% because the unmodified PSM " +
            "has higher intensity.");

        // --- Protein 2 — results are IDENTICAL ---
        var siteP2 = resultP2[4][0];
        Assert.That(siteP2.CountBasedOccupancy, Is.EqualTo(0.5),
            "Protein 2: same 50% as Protein 1.");
        Assert.That(siteP2.IntensityBasedStoichiometry, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Protein 2: same 33.3% as Protein 1. Both proteins see the same shared " +
            "peptide data, so occupancy is identical. The occupancy is NOT split — " +
            "it is DUPLICATED across both proteins.");
    }

    #endregion

    #region Test 9: Two proteins with shared + unique regions (missed cleavage), all unmodified

    /// <summary>
    /// TEST 9: Two proteins each with a shared peptide and a unique peptide,
    ///         observed as missed cleavage (both peptides joined), all unmodified.
    ///
    /// SCENARIO:
    ///   Protein 1: ACDEFGHIK   (accession P1)
    ///   Protein 2: ACDEFLMNPQ  (accession P2)
    ///
    ///   Shared region:  ACDEF (positions 1–5 in both proteins)
    ///   P1 unique:      GHIK  (positions 6–9 in Protein 1)
    ///   P2 unique:      LMNPQ (positions 6–10 in Protein 2)
    ///
    ///   Missed cleavage PSM for P1: ACDEFGHIK   (spans full P1, intensity = 1)
    ///   Missed cleavage PSM for P2: ACDEFLMNPQ  (spans full P2, intensity = 2)
    ///
    ///   The missed cleavage sequences are DIFFERENT (ACDEFGHIK vs ACDEFLMNPQ),
    ///   so they map unambiguously to their respective proteins.
    ///
    /// EXPECTED RESULT:
    ///   Empty for both proteins — no modifications observed.
    ///
    ///   KEY INSIGHT: Even though positions 1–5 contain the same amino acids in both
    ///   proteins, the missed cleavage PSMs have different full sequences. Each PSM
    ///   maps to exactly one protein. The shared peptide region is only "shared" in
    ///   the biological sense — the PSMs themselves are unambiguous.
    /// </summary>
    [Test]
    public void Test9_TwoProteinsMissedCleavageUnmodified_BothEmpty()
    {
        var protein1 = new MockBioPolymer("ACDEFGHIK", "P00001");
        var protein2 = new MockBioPolymer("ACDEFLMNPQ", "P00002");

        // Missed cleavage for P1: full protein sequence, unmodified, intensity = 1
        var peptideP1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein1, 1, 9);
        var psmP1 = new MockSpectralMatch("test.mz", peptideP1.FullSequence, peptideP1.BaseSequence, 1, 1, [peptideP1])
        {
            Intensities = new double[] { 1.0 }
        };

        // Missed cleavage for P2: full protein sequence, unmodified, intensity = 2
        var peptideP2 = new MockBioPolymerWithSetMods(
            "ACDEFLMNPQ", "ACDEFLMNPQ", protein2, 1, 10);
        var psmP2 = new MockSpectralMatch("test.mz", peptideP2.FullSequence, peptideP2.BaseSequence, 1, 1, [peptideP2])
        {
            Intensities = new double[] { 2.0 }
        };

        // Each protein gets only its own PSM
        var resultP1 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein1, new[] { psmP1 });
        var resultP2 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein2, new[] { psmP2 });

        Assert.That(resultP1, Is.Empty,
            "Protein 1: no modifications on the missed cleavage PSM → empty.");
        Assert.That(resultP2, Is.Empty,
            "Protein 2: no modifications on the missed cleavage PSM → empty.");
    }

    #endregion

    #region Test 10: Two proteins with missed cleavage, modification in UNSHARED region

    /// <summary>
    /// TEST 10: Same as Test 9, but Protein 1's PSM is modified in its UNIQUE region.
    ///
    /// SCENARIO:
    ///   Protein 1: ACDEFGHIK   (accession P1)
    ///   Protein 2: ACDEFLMNPQ  (accession P2)
    ///
    ///   PSM for P1: ACDEFG[Phospho]HIK  (modified at G, position 6 — UNIQUE to P1, intensity = 1)
    ///   PSM for P2: ACDEFLMNPQ          (unmodified, intensity = 2)
    ///
    /// FOR PROTEIN 1 (modified position G at position 6 — unique region):
    ///   ModifiedCount = 1, TotalCount = 1 → Count Occupancy = 100%
    ///   ModifiedIntensity = 1, TotalIntensity = 1 → Intensity Stoichiometry = 100%
    ///
    ///   The modification is in P1's unique region, so only P1's PSM covers it.
    ///   With only one PSM and it being modified, occupancy is 100%.
    ///
    /// FOR PROTEIN 2 (unmodified):
    ///   Empty — no modifications on P2's PSM.
    ///
    /// FOR A SHARED POSITION (e.g., position 3):
    ///   Protein 1: not reported (position 3 has no modification on P1's PSM)
    ///   Protein 2: not reported (P2's PSM is entirely unmodified)
    /// </summary>
    [Test]
    public void Test10_ModificationInUnsharedRegion_OnlyAffectsOneProtein()
    {
        var protein1 = new MockBioPolymer("ACDEFGHIK", "P00001");
        var protein2 = new MockBioPolymer("ACDEFLMNPQ", "P00002");
        var phosphoG = CreateMod("Phosphorylation", "G");

        // P1's PSM: missed cleavage with Phospho at G (key=7 → pos 1+7-2=6), intensity=1
        var peptideP1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFG[Phosphorylation]HIK", protein1, 1, 9,
            new Dictionary<int, Modification> { { 7, phosphoG } });
        var psmP1 = new MockSpectralMatch("test.mz", peptideP1.FullSequence, peptideP1.BaseSequence, 1, 1, [peptideP1])
        {
            Intensities = new double[] { 1.0 }
        };

        // P2's PSM: missed cleavage, unmodified, intensity=2
        var peptideP2 = new MockBioPolymerWithSetMods(
            "ACDEFLMNPQ", "ACDEFLMNPQ", protein2, 1, 10);
        var psmP2 = new MockSpectralMatch("test.mz", peptideP2.FullSequence, peptideP2.BaseSequence, 1, 1, [peptideP2])
        {
            Intensities = new double[] { 2.0 }
        };

        // Protein 1: gets only its own PSM (which is modified)
        var resultP1 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein1, new[] { psmP1 });

        // Protein 2: gets only its own PSM (which is unmodified)
        var resultP2 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein2, new[] { psmP2 });

        // --- Protein 1: modified position (G, position 7, unique region) ---
        Assert.That(resultP1.ContainsKey(7), Is.True,
            "Protein 1 has a modification at position 7 (G).");
        var siteP1 = resultP1[7][0];
        Assert.That(siteP1.ModifiedCount, Is.EqualTo(1));
        Assert.That(siteP1.TotalCount, Is.EqualTo(1),
            "Only P1's own PSM covers position 7. P2's PSM is not involved.");
        Assert.That(siteP1.CountBasedOccupancy, Is.EqualTo(1.0),
            "Count occupancy = 1/1 = 100%. The sole PSM is modified.");
        Assert.That(siteP1.IntensityBasedStoichiometry, Is.EqualTo(1.0),
            "Intensity stoichiometry = 1/1 = 100%.");

        // --- Protein 2: unmodified → empty ---
        Assert.That(resultP2, Is.Empty,
            "Protein 2's PSM is unmodified, so no occupancy is reported for P2.");

        // --- Shared position (e.g., position 4): not reported for either protein ---
        Assert.That(resultP1.ContainsKey(4), Is.False,
            "Shared position 4 has no modification on P1's PSM.");
    }

    #endregion

    #region Test 11: Two proteins with missed cleavage, modification in SHARED region

    /// <summary>
    /// TEST 11: Same as Test 9, but Protein 1's PSM is modified in the SHARED region.
    ///
    /// SCENARIO:
    ///   Protein 1: ACDEFGHIK   (accession P1)
    ///   Protein 2: ACDEFLMNPQ  (accession P2)
    ///
    ///   PSM for P1: ACD[Phospho]EFGHIK  (modified at D, position 4 — SHARED region, intensity = 1)
    ///   PSM for P2: ACDEFLMNPQ          (unmodified, intensity = 2)
    ///
    ///   Position 4 (D) exists in BOTH proteins, but the modification is only on P1's PSM.
    ///   Because the missed cleavage sequences are different, each PSM maps unambiguously
    ///   to its own protein.
    ///
    /// FOR PROTEIN 1 (modified at shared position D, position 4):
    ///   ModifiedCount = 1, TotalCount = 1 → Count Occupancy = 100%
    ///   ModifiedIntensity = 1, TotalIntensity = 1 → Intensity Stoichiometry = 100%
    ///
    ///   Even though position 4 is biologically "shared," P1's occupancy is calculated
    ///   using only P1's own PSM. The fact that P2's PSM also covers the same amino acid
    ///   is irrelevant — P2's PSM maps to a different protein.
    ///
    /// FOR PROTEIN 2 (unmodified):
    ///   Empty — no modifications on P2's PSM, even at position 4 (D).
    ///
    ///   KEY INSIGHT: Protein 2 does NOT get occupancy information for position 4 (D),
    ///   even though it has the same amino acid there, because Protein 2's PSM is
    ///   unmodified. Each protein's occupancy is completely independent.
    ///
    /// FOR AN UNMODIFIED POSITION ON PROTEIN 1 (e.g., position 7, G):
    ///   Not reported — no modification at position 7 on P1's PSM.
    /// </summary>
    [Test]
    public void Test11_ModificationInSharedRegion_OnlyAffectsProteinWithModifiedPsm()
    {
        var protein1 = new MockBioPolymer("ACDEFGHIK", "P00001");
        var protein2 = new MockBioPolymer("ACDEFLMNPQ", "P00002");
        var phosphoD = CreateMod("Phosphorylation", "D");

        // P1's PSM: missed cleavage with Phospho at D (key=4 → pos 1+4-1=4), intensity=1
        var peptideP1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein1, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var psmP1 = new MockSpectralMatch("test.mz", peptideP1.FullSequence, peptideP1.BaseSequence, 1, 1, [peptideP1])
        {
            Intensities = new double[] { 1.0 }
        };

        // P2's PSM: missed cleavage, unmodified, intensity=2
        var peptideP2 = new MockBioPolymerWithSetMods(
            "ACDEFLMNPQ", "ACDEFLMNPQ", protein2, 1, 10);
        var psmP2 = new MockSpectralMatch("test.mz", peptideP2.FullSequence, peptideP2.BaseSequence, 1, 1, [peptideP2])
        {
            Intensities = new double[] { 2.0 }
        };

        // Protein 1: gets its own PSM (modified at shared position)
        var resultP1 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein1, new[] { psmP1 });

        // Protein 2: gets its own PSM (unmodified)
        var resultP2 = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein2, new[] { psmP2 });

        // --- Protein 1: modified at position 4 (D, in shared region) ---
        Assert.That(resultP1.ContainsKey(4), Is.True,
            "Protein 1 has Phospho at position 4 (D), which is in the shared region.");
        var siteP1 = resultP1[4][0];
        Assert.That(siteP1.ModifiedCount, Is.EqualTo(1));
        Assert.That(siteP1.TotalCount, Is.EqualTo(1),
            "Only P1's PSM is considered for P1's occupancy. P2's PSM (even though it " +
            "covers the same amino acid sequence at position 4) belongs to a different protein.");
        Assert.That(siteP1.CountBasedOccupancy, Is.EqualTo(1.0),
            "Count occupancy = 1/1 = 100% for Protein 1.");
        Assert.That(siteP1.IntensityBasedStoichiometry, Is.EqualTo(1.0),
            "Intensity stoichiometry = 1/1 = 100% for Protein 1.");

        // --- Protein 2: no modifications → empty ---
        Assert.That(resultP2, Is.Empty,
            "Protein 2's PSM is unmodified. Even though position 4 has the SAME amino acid (D) " +
            "as Protein 1, Protein 2 shows no occupancy because its own PSM has no modifications. " +
            "Occupancy is computed per-protein, not per-amino-acid-across-proteins.");

        // --- Protein 1 at an unmodified position (e.g., position 7, G) ---
        Assert.That(resultP1.ContainsKey(7), Is.False,
            "Position 7 on Protein 1 has no modification, so it's not reported.");

        // Only 1 position reported for Protein 1 (position 4)
        Assert.That(resultP1.Count, Is.EqualTo(1),
            "Only the modified position appears in Protein 1's result.");
    }

    #endregion

    // ========================================================================
    // GAP-FILLING TESTS: Scenarios not covered by the original 11 test prompts
    // ========================================================================

    #region Gap A: Competing modifications at the SAME position

    /// <summary>
    /// GAP TEST A: Two different modification types at the same amino acid position.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK
    ///   PSM 1:    ACD[Phospho]EFGHIK    (Phospho on D at position 3, intensity = 1)
    ///   PSM 2:    ACD[Acetyl]EFGHIK     (Acetyl on D at position 3, intensity = 2)
    ///   PSM 3:    ACDEFGHIK             (unmodified, intensity = 3)
    ///
    /// WHY THIS MATTERS:
    ///   When two different modifications compete for the same site, the calculator
    ///   uses a "positionTotals" cache to ensure they SHARE the same denominator.
    ///   Without this cache, each mod would independently count total coverage,
    ///   and their occupancies could (incorrectly) sum to more than 100%.
    ///
    ///   With the cache, both mods share TotalCount=3, so:
    ///     Phospho occupancy = 1/3 ≈ 33.3%
    ///     Acetyl occupancy  = 1/3 ≈ 33.3%
    ///     Sum               = 2/3 ≈ 66.7%  (leaves room for the unmodified 1/3)
    ///
    ///   This correctly reflects that 1/3 of observations are Phospho, 1/3 are Acetyl,
    ///   and 1/3 are unmodified. The occupancies are coherent.
    /// </summary>
    [Test]
    public void GapA_CompetingModsAtSamePosition_ShareDenominator()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoD = CreateMod("Phosphorylation", "D");

        // Acetyl on D — different modification at the same amino acid
        ModificationMotif.TryGetMotif("D", out var motifD);
        var acetylD = new Modification("Acetylation", null, "Biological", null, motifD, "Anywhere.", null, 42.011);

        // PSM 1: Phospho at D (key=4 → protein position 3), intensity = 1
        var phosphoPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var phosphoPsm = new MockSpectralMatch("test.mz", phosphoPeptide.FullSequence, phosphoPeptide.BaseSequence, 1, 1, [phosphoPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        // PSM 2: Acetyl at D (key=4 → protein position 3), intensity = 2
        var acetylPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Acetylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, acetylD } });
        var acetylPsm = new MockSpectralMatch("test.mz", acetylPeptide.FullSequence, acetylPeptide.BaseSequence, 1, 1, [acetylPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        // PSM 3: unmodified, intensity = 3
        var unmodPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein, 1, 9);
        var unmodPsm = new MockSpectralMatch("test.mz", unmodPeptide.FullSequence, unmodPeptide.BaseSequence, 1, 1, [unmodPeptide])
        {
            Intensities = new double[] { 3.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { phosphoPsm, acetylPsm, unmodPsm });

        // Position 4 should have TWO entries: one for Phospho, one for Acetyl
        Assert.That(result.ContainsKey(4), Is.True);
        Assert.That(result[4].Count, Is.EqualTo(2),
            "Two different mods at position 4 → two SiteSpecificModificationOccupancy entries.");

        var phosphoSite = result[4].First(s => s.ModificationIdWithMotif == "Phosphorylation on D");
        var acetylSite = result[4].First(s => s.ModificationIdWithMotif == "Acetylation on D");

        // Both mods share the SAME denominator (TotalCount = 3, TotalIntensity = 6)
        // This is the key behavior of the positionTotals cache.
        Assert.That(phosphoSite.TotalCount, Is.EqualTo(3),
            "Phospho shares denominator: all 3 PSMs cover position 4.");
        Assert.That(acetylSite.TotalCount, Is.EqualTo(3),
            "Acetyl shares the SAME denominator as Phospho. The positionTotals cache " +
            "ensures that the denominator is computed once per position, not once per mod type.");

        Assert.That(phosphoSite.TotalIntensity, Is.EqualTo(6.0),
            "Shared total intensity = 1 + 2 + 3 = 6.");
        Assert.That(acetylSite.TotalIntensity, Is.EqualTo(6.0),
            "Same shared total intensity for Acetyl.");

        // Each mod has its own numerator
        Assert.That(phosphoSite.ModifiedCount, Is.EqualTo(1));
        Assert.That(phosphoSite.CountBasedOccupancy, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Phospho count occupancy = 1/3 ≈ 33.3%.");
        Assert.That(phosphoSite.ModifiedIntensity, Is.EqualTo(1.0));
        Assert.That(phosphoSite.IntensityBasedStoichiometry, Is.EqualTo(1.0 / 6.0).Within(1e-10),
            "Phospho intensity stoichiometry = 1/6 ≈ 16.7%.");

        Assert.That(acetylSite.ModifiedCount, Is.EqualTo(1));
        Assert.That(acetylSite.CountBasedOccupancy, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Acetyl count occupancy = 1/3 ≈ 33.3%.");
        Assert.That(acetylSite.ModifiedIntensity, Is.EqualTo(2.0));
        Assert.That(acetylSite.IntensityBasedStoichiometry, Is.EqualTo(2.0 / 6.0).Within(1e-10),
            "Acetyl intensity stoichiometry = 2/6 ≈ 33.3%.");

        // The sum of count-based occupancies = 2/3, leaving 1/3 for unmodified. Coherent!
        double sumCountOccupancy = phosphoSite.CountBasedOccupancy + acetylSite.CountBasedOccupancy;
        Assert.That(sumCountOccupancy, Is.EqualTo(2.0 / 3.0).Within(1e-10),
            "Sum of occupancies = 2/3 ≈ 66.7%. The remaining 1/3 is the unmodified fraction. " +
            "The shared denominator guarantees occupancies are coherent and sum to ≤ 1.0.");
    }

    #endregion

    #region Gap B: Ambiguous PSM interpretations (sequencesForTotalCount deduplication)

    /// <summary>
    /// GAP TEST B: One PSM with two ambiguous modification localizations.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK
    ///   1 PSM, but the search engine reports two possible interpretations:
    ///     Form 1: ACD[Phospho]EFGHIK  (Phospho on D, position 3)
    ///     Form 2: ACDE[Phospho]FGHIK  (Phospho on E, position 4)
    ///
    ///   This is a SINGLE observation from the mass spec — the PSM could be either form,
    ///   but we don't know which.
    ///
    /// THE BUG THIS PREVENTS:
    ///   If we naively pass both forms as the full peptide list, the calculator sees
    ///   2 "peptides" covering each position → TotalCount = 2. But only 1 PSM exists!
    ///   This would give occupancy = 1/2 = 50% at each site, implying the modification
    ///   is absent half the time — which is wrong, because ALL observations show the mod.
    ///
    /// THE FIX:
    ///   Pass the full list as `localizedSequences` (for the numerator — both forms count),
    ///   but pass a DEDUPLICATED list (one entry per PSM) as `sequencesForTotalCount`
    ///   (for the denominator). This gives TotalCount = 1 and occupancy = 1/1 = 100%.
    ///
    /// This test demonstrates the fix both with and without deduplication so you can
    /// see the difference.
    /// </summary>
    [Test]
    public void GapB_AmbiguousPsm_WithoutDeduplication_InflatesDenominator()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoD = CreateMod("Phosphorylation", "D");
        var phosphoE = CreateMod("Phosphorylation", "E");

        // Two interpretations of the SAME PSM
        var form1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var form2 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDE[Phosphorylation]FGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 5, phosphoE } });
        var psm1 = new MockSpectralMatch("test.mz", form1.FullSequence, form1.BaseSequence, 1, 1, [form1])
        {
            Intensities = new double[] { 1.0 }
        };
        var psm2 = new MockSpectralMatch("test.mz", form2.FullSequence, form2.BaseSequence, 1, 1, [form2])
        {
            Intensities = new double[] { 1.0 }
        };

        // WITHOUT deduplication: pass both forms as both localizedSequences AND coverage
        // (sequencesForTotalCount = null means localizedSequences is reused for denominator)
        var resultBuggy = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { psm1, psm2 }); 

        // The denominator is inflated: TotalCount = 2 (both forms counted)
        var siteD_buggy = resultBuggy[4][0];
        Assert.That(siteD_buggy.TotalCount, Is.EqualTo(2),
            "WITHOUT deduplication: TotalCount = 2 because both interpretations are counted " +
            "as separate observations. But there was really only 1 PSM!");
        Assert.That(siteD_buggy.CountBasedOccupancy, Is.EqualTo(0.5),
            "WITHOUT deduplication: occupancy = 1/2 = 50%. This is MISLEADING — " +
            "it suggests the mod is absent half the time, but every observation has the mod.");
    }

    [Test]
    public void GapBTemp_AmbiguousPsm_Unreported()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoD = CreateMod("Phosphorylation", "D");
        var phosphoE = CreateMod("Phosphorylation", "E");

        var form1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var form2 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDE[Phosphorylation]FGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 5, phosphoE } });

        // This represents a SINGLE PSM with two ambiguous forms. In this instance, since the 
        // code is currently designed to only report unambiguous PSMs (full sequence != null), 
        // the psm will NOT contribute towards occupancy. This is the correct conservative behaviour
        // for now. 
        var psm = new MockSpectralMatch("test.mz", null, form1.BaseSequence, 1, 1, [form1, form2])
        {
            Intensities = new double[] { 1.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] {psm}); // <-- 1 per PSM

        Assert.That(result, Is.Empty);
    }

    #endregion

    #region Gap C: N-terminal and C-terminal modifications

    /// <summary>
    /// GAP TEST C: Modifications with N-terminal and C-terminal location restrictions.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK  (length 9)
    ///   PSM 1:    [Acetyl]ACDEFGHIK  (N-terminal Acetylation, intensity = 1)
    ///   PSM 2:    ACDEFGHIK          (unmodified, intensity = 2)
    ///
    /// POSITION MAPPING:
    ///   N-terminal mods (LocationRestriction = "N-terminal.") ALWAYS map to protein position 1,
    ///   regardless of where the peptide starts in the protein. This is a special case in
    ///   TryGetProteinPosition.
    ///
    ///   Similarly, C-terminal mods ("C-terminal.") ALWAYS map to position
    ///   (bioPolymerLength + 2) in the result dictionary.
    ///
    ///   This differs from "Anywhere." mods which use the formula:
    ///     proteinPosition = OneBasedStartResidue + key - 1    
    ///
    /// AT PROTEIN POSITION 1 (N-terminal Acetylation):
    ///   Count: 1/2 = 50%
    ///   Intensity: 1/3 ≈ 33.3%
    /// </summary>
    [Test]
    public void GapC_NTerminalModification_MapsToProteinPosition1()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");

        // N-terminal mod uses LocationRestriction = "N-terminal." (note the period)
        ModificationMotif.TryGetMotif("A", out var motif);
        var nTermAcetyl = new Modification("Acetylation", null, "Biological", null, motif,
            "N-terminal.", null, 42.011);

        // PSM 1: N-terminal acetylation (key=1 in AllModsOneIsNterminus = N-terminal slot)
        var modPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "[Acetylation]ACDEFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 1, nTermAcetyl } });
        var modPsm = new MockSpectralMatch("test.mz", modPeptide.FullSequence, modPeptide.BaseSequence, 1, 1, [modPeptide])
        {
            Intensities = new double[] { 1.0 }
        };
        // PSM 2: unmodified
        var unmodPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein, 1, 9);
        var unmodPsm = new MockSpectralMatch("test.mz", unmodPeptide.FullSequence, unmodPeptide.BaseSequence, 1, 1, [unmodPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { modPsm, unmodPsm });

        // N-terminal mods always map to protein position 1
        Assert.That(result.ContainsKey(1), Is.True,
            "N-terminal mods map to protein position 1, regardless of AllModsOneIsNterminus key. " +
            "The TryGetProteinPosition method has special handling: if LocationRestriction is " +
            "'N-terminal.', it sets indexInProtein = 1.");

        var site = result[1][0];
        Assert.That(site.ModifiedCount, Is.EqualTo(1));
        Assert.That(site.TotalCount, Is.EqualTo(2));
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0.5),
            "Count-based occupancy = 1/2 = 50%.");
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "Intensity-based stoichiometry = 1/3 ≈ 33.3%.");

        // --- C-terminal modification (applied directly to the protein) ---
        var cTermAcetyl = new Modification("Acetylation", null, "Biological", null, motif,
            "C-terminal.", null, 42.011);

        // PSM with C-terminal acetylation
        var cTermPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK[Acetylation]", protein, 1, 9,
            new Dictionary<int, Modification> { { 11, cTermAcetyl } });
        var cTermPsm = new MockSpectralMatch("test.mz", cTermPeptide.FullSequence, cTermPeptide.BaseSequence, 1, 1, [cTermPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        var resultCterm = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { cTermPsm });

        // C-terminal mods always map to bioPolymerLength + 2 (position 11 here: 9 + 2 = 11)
        Assert.That(resultCterm.ContainsKey(11), Is.True,
            "C-terminal mods map to bioPolymerLength + 2 = 11 in the result dictionary. " +
            "TryGetProteinPosition sets indexInProtein = bioPolymerLength + 2 for 'C-terminal.' mods.");

        Assert.That(resultCterm[11][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(resultCterm[11][0].TotalCount, Is.EqualTo(1));
        Assert.That(resultCterm[11][0].CountBasedOccupancy, Is.EqualTo(1.0));
    }

    #endregion

    #region Gap D: Peptide starting in the MIDDLE of the protein

    /// <summary>
    /// GAP TEST D: A peptide that does not start at position 1 in the protein.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK  (positions 1–9)
    ///   Peptide:  FGHIK      (positions 5–9, a tryptic peptide from the C-terminal half)
    ///   Modification: Phospho on G (2nd residue of peptide, key=3 in AllModsOneIsNterminus)
    ///
    /// POSITION MAPPING:
    ///   For "Anywhere." mods: proteinPosition = OneBasedStartResidue + key - 1
    ///   Here: proteinPosition = 5 + 3 - 1 = 7
    ///   So key=3 in the peptide maps to protein position 7 (G). Correct!
    ///
    ///   This test verifies the position mapping formula when OneBasedStartResidue ≠ 1.
    ///   In Tests 1–11, all peptides started at position 1, so the formula simplified to
    ///   proteinPosition = key - 1. Here we confirm the general formula works.
    ///
    /// WHY THIS MATTERS:
    ///   In real experiments, proteins are digested into peptides by trypsin. Most peptides
    ///   do NOT start at position 1 of the protein. The position mapping formula must
    ///   correctly translate peptide-local modification positions to absolute protein coordinates.
    /// </summary>
    [Test]
    public void GapD_MidProteinPeptide_PositionMappingUsesStartResidue()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoG = CreateMod("Phosphorylation", "G");

        // Peptide FGHIK starts at position 5 in the protein
        // G is the 2nd residue of the peptide → AllModsOneIsNterminus key = 3
        // (key 1 = N-term, key 2 = F, key 3 = G, key 4 = H, ...)
        // Protein position = 5 + 3 - 1 = 7 → G is at protein position 7 ✓
        var peptide = new MockBioPolymerWithSetMods(
            "FGHIK", "FG[Phosphorylation]HIK", protein, 5, 9,
            new Dictionary<int, Modification> { { 3, phosphoG } });
        var psm = new MockSpectralMatch("test.mz", peptide.FullSequence, peptide.BaseSequence, 1, 1, [peptide])
        {
            Intensities = new double[] { 1.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, new[] { psm });

        // The modification should appear at protein position 7 (G), NOT position 3
        Assert.That(result.ContainsKey(7), Is.True,
            "Key=3 in a peptide starting at position 5 maps to protein position 5+3-1=7. " +
            "The formula accounts for the peptide's offset within the protein.");
        Assert.That(result.ContainsKey(3), Is.False,
            "Position 3 would be wrong — that would be the result if OneBasedStartResidue were ignored.");

        Assert.That(result[7][0].ModifiedCount, Is.EqualTo(1));
        Assert.That(result[7][0].TotalCount, Is.EqualTo(1));
        Assert.That(result[7][0].CountBasedOccupancy, Is.EqualTo(1.0));
    }

    #endregion

    #region Gap E: Peptide-level vs Protein-level coordinate systems

    /// <summary>
    /// GAP TEST E: Same modification analyzed at both protein-level and peptide-level.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK  (positions 1–9)
    ///   Peptide:  FGHIK     (positions 5–9 in protein)
    ///   Mod:      Phospho on G, AllModsOneIsNterminus key = 3
    ///
    /// PROTEIN-LEVEL result:
    ///   Position key = 6 (mapped to protein coordinates: 5 + 3 - 2 = 6)
    ///
    /// PEPTIDE-LEVEL result:
    ///   Position key = 3 (the raw AllModsOneIsNterminus key, NOT mapped to protein)
    ///   This means: key 1 = N-terminal, key 2 = 1st residue (F), key 3 = 2nd residue (G)
    ///
    /// WHY THIS MATTERS:
    ///   The two calculators use DIFFERENT coordinate systems:
    ///   - Protein-level: absolute position in the protein (1-based)
    ///   - Peptide-level: position within the peptide using AllModsOneIsNterminus convention
    ///
    ///   When reviewing results, you must know which calculator produced them to interpret
    ///   the position numbers correctly.
    /// </summary>
    [Test]
    public void GapE_PeptideLevelVsProteinLevel_DifferentCoordinates()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoG = CreateMod("Phosphorylation", "G");

        var modPeptide = new MockBioPolymerWithSetMods(
            "FGHIK", "FG[Phosphorylation]HIK", protein, 5, 9,
            new Dictionary<int, Modification> { { 3, phosphoG } });
        var unmodPeptide = new MockBioPolymerWithSetMods(
            "FGHIK", "FGHIK", protein, 5, 9);

        var modPsm = new MockSpectralMatch("test.mz", modPeptide.FullSequence, modPeptide.BaseSequence, 1, 1, [modPeptide])
        {
            Intensities = new double[] { 1.0 }
        };
        var unmodPsm = new MockSpectralMatch("test.mz", unmodPeptide.FullSequence, unmodPeptide.BaseSequence, 1, 1, [unmodPeptide])
        {
            Intensities = new double[] { 2.0 }
        };

        // --- Protein-level: maps to PROTEIN coordinates ---
        var proteinResult = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { modPsm , unmodPsm });

        Assert.That(proteinResult.ContainsKey(7), Is.True,
            "PROTEIN-level uses absolute protein coordinates. Key=3 → position 5+3-1=7.");
        Assert.That(proteinResult.ContainsKey(4), Is.False,
            "Position 4 would be wrong for protein-level — that's where D is, not G.");

        // --- Peptide-level: uses raw AllModsOneIsNterminus keys ---
        var peptideResult = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy(
            new[] { modPsm, unmodPsm });

        Assert.That(peptideResult.ContainsKey(3), Is.True,
            "PEPTIDE-level uses the raw AllModsOneIsNterminus key directly. " +
            "Key=3 means '2nd residue' in the peptide (key 1=N-term, 2=1st residue, 3=2nd residue). " +
            "No mapping to protein coordinates is performed.");
        Assert.That(peptideResult.ContainsKey(6), Is.False,
            "Position 6 would be wrong for peptide-level — peptide-level doesn't know about " +
            "protein coordinates.");

        // Both calculators report the same occupancy values — only the position keys differ
        Assert.That(proteinResult[7][0].CountBasedOccupancy, Is.EqualTo(0.5));
        Assert.That(peptideResult[3][0].CountBasedOccupancy, Is.EqualTo(0.5));
        Assert.That(proteinResult[7][0].CountBasedOccupancy,
            Is.EqualTo(peptideResult[3][0].CountBasedOccupancy),
            "Same modification, same PSMs → same occupancy. Only the position key differs " +
            "between protein-level (7) and peptide-level (3).");
    }

    #endregion

    #region Gap F: Excluded modification types are silently filtered

    /// <summary>
    /// GAP TEST F: Certain modification types are automatically excluded from occupancy.
    ///
    /// The calculator filters out:
    ///   1. "Common Variable" mods (e.g., Oxidation of M) — these are search artifacts,
    ///      not biologically meaningful PTMs
    ///   2. "Common Fixed" mods (e.g., Carbamidomethylation of C) — applied uniformly
    ///      during sample preparation, always present, occupancy is always 100% and meaningless
    ///   3. Peptide-terminal mods with LocationRestriction "NPep" or "PepC" — these are
    ///      peptide-level artifacts from digestion, not true protein modifications
    ///
    /// Protein-terminal mods ("N-terminal." and "C-terminal.") are NOT excluded — they
    /// represent real protein modifications.
    ///
    /// WHY THIS MATTERS:
    ///   If you expect to see occupancy for a modification and the result is empty,
    ///   check whether the modification type is in the excluded list. This is a common
    ///   source of confusion when analyzing results.
    /// </summary>
    [Test]
    public void GapF_CommonVariableModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);

        // "Common Variable" type → excluded from occupancy
        var oxidation = new Modification("Oxidation", null, "Common Variable", null, motif,
            "Anywhere.", null, 15.995);

        var peptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Oxidation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, oxidation } });
        var psm = new MockSpectralMatch("test.mz", peptide.FullSequence, peptide.BaseSequence, 1, 1, [peptide])
        {
            Intensities = new double[] { 1.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, new[] { psm });

        Assert.That(result, Is.Empty,
            "'Common Variable' modifications are excluded from occupancy calculations. " +
            "These are typically search engine artifacts (like oxidation) that don't represent " +
            "biologically meaningful PTMs. The modification is silently skipped.");
    }

    [Test]
    public void GapF_CommonFixedModIsExcluded()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("C", out var motif);

        // "Common Fixed" type → excluded from occupancy
        var carbamido = new Modification("Carbamidomethyl", null, "Common Fixed", null, motif,
            "Anywhere.", null, 57.021);

        var peptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "AC[Carbamidomethyl]DEFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 3, carbamido } });
        var psm = new MockSpectralMatch("test.mz", peptide.FullSequence, peptide.BaseSequence, 1, 1, [peptide])
        {
            Intensities = new double[] { 1.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, new[] { psm });

        Assert.That(result, Is.Empty,
            "'Common Fixed' modifications are excluded. These mods (like carbamidomethylation " +
            "of cysteine) are applied during sample preparation and present on every peptide — " +
            "their occupancy would always be 100% and carry no biological information.");
    }

    [Test]
    public void GapF_PeptideTerminalModIsExcluded_ButProteinTerminalIsKept()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("A", out var motif);

        // Peptide N-terminal mod (LocationRestriction = "NPep") → excluded
        var pepNterm = new Modification("PyroGlu", null, "Biological", null, motif,
            "NPep", null, -17.027);
        var pepNtermPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "[PyroGlu]ACDEFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 1, pepNterm } });
        var pepNtermPsm = new MockSpectralMatch("test.mz", pepNtermPeptide.FullSequence, pepNtermPeptide.BaseSequence, 1, 1, [pepNtermPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        var resultPepTerm = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, new[] { pepNtermPsm });

        Assert.That(resultPepTerm, Is.Empty,
            "'NPep' (peptide N-terminal) mods are excluded. These are artifacts of enzymatic " +
            "digestion, not true protein modifications.");

        // Protein N-terminal mod (LocationRestriction = "N-terminal.") → KEPT
        var protNterm = new Modification("Acetylation", null, "Biological", null, motif,
            "N-terminal.", null, 42.011);
        var protNtermPeptide = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "[Acetylation]ACDEFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 1, protNterm } });
        var protNtermPsm = new MockSpectralMatch("test.mz", protNtermPeptide.FullSequence, protNtermPeptide.BaseSequence, 1, 1, [protNtermPeptide])
        {
            Intensities = new double[] { 1.0 }
        };

        var resultProtTerm = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein, new[] { protNtermPsm });

        Assert.That(resultProtTerm, Is.Not.Empty,
            "'N-terminal.' (protein N-terminal) mods are NOT excluded. These represent real " +
            "biological modifications of the protein's N-terminus. The subtle difference in " +
            "LocationRestriction strings ('NPep' vs 'N-terminal.') determines the behavior.");
    }

    #endregion

    #region Gap G: Multiple PSMs of the same modified form

    /// <summary>
    /// GAP TEST G: Multiple PSMs all carrying the same modification.
    ///
    /// SCENARIO:
    ///   Protein:  ACDEFGHIK
    ///   PSM 1:    ACD[Phospho]EFGHIK  (intensity = 1)
    ///   PSM 2:    ACD[Phospho]EFGHIK  (intensity = 3)
    ///   PSM 3:    ACD[Phospho]EFGHIK  (intensity = 5)
    ///   PSM 4:    ACDEFGHIK           (unmodified, intensity = 1)
    ///
    ///   Three PSMs carry the modification, one does not.
    ///
    /// AT POSITION 4 (D):
    ///   Count:     3 modified / 4 total = 75%  (CORRECT)
    ///   Intensity: expected 9/10 = 90%
    ///
    /// WHY THIS MATTERS:
    ///   In real experiments, you often see the same modification in many PSMs.
    ///   The ModifiedCount accumulates — it's not just 0 or 1.
    ///
    ///   Note on intensity: the intensity dictionary is keyed by FullSequence, not by PSM.
    ///   All 3 modified PSMs share the same FullSequence ("ACD[Phosphorylation]EFGHIK"),
    ///   so they map to a single entry. The caller (PopulateOccupancy) sums their
    ///   intensities into one dictionary entry before calling the calculator.
    ///   For this test, we sum them as the caller would: 1 + 3 + 5 = 9.
    /// </summary>
    [Test]
    public void GapG_MultiplePsmsWithSameModification()
    {
        var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        var phosphoD = CreateMod("Phosphorylation", "D");

        // 3 PSMs with the same modification (same FullSequence)
        var mod1 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var mod2 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });
        var mod3 = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACD[Phosphorylation]EFGHIK", protein, 1, 9,
            new Dictionary<int, Modification> { { 4, phosphoD } });

        var psm1 = new MockSpectralMatch("test.mz", mod1.FullSequence, mod1.BaseSequence, 1, 1, [mod1])
        {
            Intensities = new double[] { 1.0 }
        };
        var psm2 = new MockSpectralMatch("test.mz", mod2.FullSequence, mod2.BaseSequence, 1, 1, [mod2])
        {
            Intensities = new double[] { 3.0 }
        };
        var psm3 = new MockSpectralMatch("test.mz", mod3.FullSequence, mod3.BaseSequence, 1, 1, [mod3])
        {
            Intensities = new double[] { 5.0 }
        };

        // 1 unmodified PSM
        var unmod = new MockBioPolymerWithSetMods(
            "ACDEFGHIK", "ACDEFGHIK", protein, 1, 9);
        var unmodPsm = new MockSpectralMatch("test.mz", unmod.FullSequence, unmod.BaseSequence, 1, 1, [unmod])
        {
            Intensities = new double[] { 1.0 }
        };

        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            protein,
            new[] { psm1, psm2, psm3, unmodPsm });

        var site = result[4][0];

        Assert.That(site.ModifiedCount, Is.EqualTo(3),
            "ModifiedCount = 3 because three separate PSM forms carry Phospho at this site. " +
            "Each peptide in the localizedSequences list that has this mod increments the count.");
        Assert.That(site.TotalCount, Is.EqualTo(4),
            "TotalCount = 4: all 4 peptide forms cover position 4.");
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0.75),
            "Count occupancy = 3/4 = 75%. Three-quarters of observations are modified.");

        // Intensity: the calculator looks up each peptide's FullSequence in the intensity dict.
        // All 3 modified peptides have FullSequence "ACD[Phosphorylation]EFGHIK" → intensity 9.
        // The unmodified peptide has FullSequence "ACDEFGHIK" → intensity 1.
        // So the calculator sums: ModifiedIntensity = 9, TotalIntensity = 9 + 1 = 10
        Assert.That(site.ModifiedIntensity, Is.EqualTo(9.0),
            "ModifiedIntensity = 1 + 3 + 5 = 9. Each PSM contributes its own intensity " +
            "individually, so the three modified PSMs are summed correctly.");
        Assert.That(site.TotalIntensity, Is.EqualTo(10.0),
            "TotalIntensity = 9 (modified) + 1 (unmodified) = 10.");
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(9.0 / 10.0).Within(1e-10),
            "Intensity stoichiometry = 9/10 = 90%.");

        // NOTE: Each PSM contributes its own intensity value independently.
        // When 3 PSMs share FullSequence "ACD[Phospho]EFGHIK" with intensities 1, 3, and 5,
        // the calculator accumulates 1 + 3 + 5 = 9 for ModifiedIntensity (not 9+9+9 = 27).
        // The denominator is likewise correct: 9 + 1 = 10.
        // Intensity-based stoichiometry therefore equals the true signal ratio: 9/10 = 90%.
        double expectedStoichiometry = 9.0 / 10.0;
        Assert.That(site.IntensityBasedStoichiometry,
            Is.EqualTo(expectedStoichiometry).Within(1e-10),
            "The computed stoichiometry (90%) matches the expected true stoichiometry. " +
            "Each PSM's intensity is counted once, regardless of how many PSMs share " +
            "the same FullSequence.");
    }

    #endregion
}
