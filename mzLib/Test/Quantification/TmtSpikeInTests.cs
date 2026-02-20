using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Quantification.Strategies;
using Test.Quantification.TestHelpers;

namespace Test.Quantification;

[TestFixture]
[ExcludeFromCodeCoverage]
public class TmtSpikeInTests
{
    #region Private Helper Classes

    private class TmtTestExperimentalDesign : IExperimentalDesign
    {
        public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

        public TmtTestExperimentalDesign(Dictionary<string, ISampleInfo[]> dict)
        {
            FileNameSampleInfoDictionary = dict;
        }
    }

    #endregion

    #region Helper Methods

    /// <summary>
    /// Creates a TMT experimental design where all files share the same channel-to-condition mapping.
    /// </summary>
    /// <param name="fileNames">File names (with extension, e.g. "file1.raw")</param>
    /// <param name="channelLabels">TMT channel labels, e.g. { "126", "127N", "127C" }</param>
    /// <param name="conditions">Condition per channel (same length as channelLabels)</param>
    /// <param name="bioReps">Biological replicate per channel (same length as channelLabels)</param>
    /// <param name="techReps">Technical replicate per file (same length as fileNames)</param>
    private static IExperimentalDesign CreateTmtExperimentalDesign(
        string[] fileNames,
        string[] channelLabels,
        string[] conditions,
        int[] bioReps,
        int[] techReps)
    {
        var dict = new Dictionary<string, ISampleInfo[]>();
        for (int f = 0; f < fileNames.Length; f++)
        {
            var samples = new ISampleInfo[channelLabels.Length];
            for (int c = 0; c < channelLabels.Length; c++)
            {
                bool isRef = conditions[c] == "Reference";
                samples[c] = new IsobaricQuantSampleInfo(
                    fileNames[f], conditions[c], bioReps[c], techReps[f], 0, 0,
                    channelLabels[c], 126.0 + c * 0.1, isRef);
            }
            dict[fileNames[f]] = samples;
        }
        return new TmtTestExperimentalDesign(dict);
    }

    #endregion

    /// <summary>
    /// Tests that PivotByFile creates one matrix per file, with correct dimensions and values.
    /// </summary>
    [Test]
    public void PivotByFile_CreatesCorrectPerFileMatrices()
    {
        // Setup: 2 files, 3 channels each
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2 }, channels, conditions, bioReps, techReps);

        // PSMs do not need identifiedBioPolymers for PivotByFile — only QuantValues matter
        // PSMs are ordered by FullSequence within each file: "PEPTIDEK" < "SEQUENCEK"
        var file1Psm1 = new BaseSpectralMatch(file1, 1, 100.0, "PEPTIDEK", "PEPTIDEK")
            { QuantValues = new double[] { 100, 200, 300 } };
        var file1Psm2 = new BaseSpectralMatch(file1, 2, 95.0, "SEQUENCEK", "SEQUENCEK")
            { QuantValues = new double[] { 400, 500, 600 } };
        var file2Psm1 = new BaseSpectralMatch(file2, 3, 90.0, "PEPTIDEK", "PEPTIDEK")
            { QuantValues = new double[] { 700, 800, 900 } };
        var file2Psm2 = new BaseSpectralMatch(file2, 4, 85.0, "SEQUENCEK", "SEQUENCEK")
            { QuantValues = new double[] { 1000, 1100, 1200 } };

        var allPsms = new List<ISpectralMatch> { file1Psm1, file1Psm2, file2Psm1, file2Psm2 };

        // Act
        var perFileMatrices = QuantificationEngine.PivotByFile(allPsms, design);

        // Assert: 2 entries (one per file)
        Assert.That(perFileMatrices.Count, Is.EqualTo(2));
        Assert.That(perFileMatrices.ContainsKey(file1), Is.True);
        Assert.That(perFileMatrices.ContainsKey(file2), Is.True);

        // File1 matrix: 2 rows (PSMs), 3 columns (channels)
        var m1 = perFileMatrices[file1];
        Assert.That(m1.RowKeys.Count, Is.EqualTo(2));
        Assert.That(m1.ColumnKeys.Count, Is.EqualTo(3));

        // Row order is by FullSequence: "PEPTIDEK" (row 0) < "SEQUENCEK" (row 1)
        var f1Row0 = m1.GetRow(0); // PEPTIDEK
        Assert.That(f1Row0[0], Is.EqualTo(100.0));
        Assert.That(f1Row0[1], Is.EqualTo(200.0));
        Assert.That(f1Row0[2], Is.EqualTo(300.0));

        var f1Row1 = m1.GetRow(1); // SEQUENCEK
        Assert.That(f1Row1[0], Is.EqualTo(400.0));
        Assert.That(f1Row1[1], Is.EqualTo(500.0));
        Assert.That(f1Row1[2], Is.EqualTo(600.0));

        // File2 matrix: 2 rows, 3 columns
        var m2 = perFileMatrices[file2];
        Assert.That(m2.RowKeys.Count, Is.EqualTo(2));
        Assert.That(m2.ColumnKeys.Count, Is.EqualTo(3));

        var f2Row0 = m2.GetRow(0); // PEPTIDEK
        Assert.That(f2Row0[0], Is.EqualTo(700.0));
        Assert.That(f2Row0[1], Is.EqualTo(800.0));
        Assert.That(f2Row0[2], Is.EqualTo(900.0));

        var f2Row1 = m2.GetRow(1); // SEQUENCEK
        Assert.That(f2Row1[0], Is.EqualTo(1000.0));
        Assert.That(f2Row1[1], Is.EqualTo(1100.0));
        Assert.That(f2Row1[2], Is.EqualTo(1200.0));
    }

    /// <summary>
    /// Tests that CombinePeptideMatrices correctly merges per-file matrices.
    /// Shared peptides get values from both files; absent peptides get zeros.
    /// </summary>
    [Test]
    public void CombinePeptideMatrices_MergesCorrectly()
    {
        // Create 3 peptides from 3 proteins, with sequences that sort alphabetically:
        // AAAAK < LLLLK < SSSSSK
        var protA = new Protein("AAAAK", "PA");
        var protB = new Protein("LLLLK", "PB");
        var protC = new Protein("SSSSSK", "PC");

        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 4);
        var pepA = protA.Digest(digParams, new List<Modification>(), new List<Modification>()).First();
        var pepB = protB.Digest(digParams, new List<Modification>(), new List<Modification>()).First();
        var pepC = protC.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        // Experimental design: 2 files, 3 channels each
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2 }, channels, conditions, bioReps, techReps);

        var file1Samples = design.FileNameSampleInfoDictionary[file1];
        var file2Samples = design.FileNameSampleInfoDictionary[file2];

        // File1: PepA[100, 200, 300], PepB[400, 500, 600]
        var matrix1 = new PeptideMatrix(
            new List<IBioPolymerWithSetMods> { pepA, pepB },
            file1Samples,
            design);
        matrix1.SetRow(pepA, new double[] { 100, 200, 300 });
        matrix1.SetRow(pepB, new double[] { 400, 500, 600 });

        // File2: PepA[700, 800, 900], PepC[1000, 1100, 1200]
        var matrix2 = new PeptideMatrix(
            new List<IBioPolymerWithSetMods> { pepA, pepC },
            file2Samples,
            design);
        matrix2.SetRow(pepA, new double[] { 700, 800, 900 });
        matrix2.SetRow(pepC, new double[] { 1000, 1100, 1200 });

        var perFilePeptideMatrices = new Dictionary<string, QuantMatrix<IBioPolymerWithSetMods>>
        {
            { file1, matrix1 },
            { file2, matrix2 }
        };

        // Act
        var combined = QuantificationEngine.CombinePeptideMatrices(perFilePeptideMatrices, design);

        // Assert structure: 3 peptides, 6 columns (3 from file1 + 3 from file2)
        Assert.That(combined.RowKeys.Count, Is.EqualTo(3));
        Assert.That(combined.ColumnKeys.Count, Is.EqualTo(6));

        // Peptides are sorted by FullSequence: AAAAK, LLLLK, SSSSSK
        // Columns: file1[0-2], file2[3-5] (alphabetical file order: file1 < file2)

        // PepA (AAAAK): present in both files
        var pepARow = combined.GetRow(pepA);
        Assert.That(pepARow[0], Is.EqualTo(100.0)); // file1, ch0
        Assert.That(pepARow[1], Is.EqualTo(200.0)); // file1, ch1
        Assert.That(pepARow[2], Is.EqualTo(300.0)); // file1, ch2
        Assert.That(pepARow[3], Is.EqualTo(700.0)); // file2, ch0
        Assert.That(pepARow[4], Is.EqualTo(800.0)); // file2, ch1
        Assert.That(pepARow[5], Is.EqualTo(900.0)); // file2, ch2

        // PepB (LLLLK): present only in file1, zeros in file2
        var pepBRow = combined.GetRow(pepB);
        Assert.That(pepBRow[0], Is.EqualTo(400.0));
        Assert.That(pepBRow[1], Is.EqualTo(500.0));
        Assert.That(pepBRow[2], Is.EqualTo(600.0));
        Assert.That(pepBRow[3], Is.EqualTo(0.0));
        Assert.That(pepBRow[4], Is.EqualTo(0.0));
        Assert.That(pepBRow[5], Is.EqualTo(0.0));

        // PepC (SSSSSK): absent in file1, zeros there; present in file2
        var pepCRow = combined.GetRow(pepC);
        Assert.That(pepCRow[0], Is.EqualTo(0.0));
        Assert.That(pepCRow[1], Is.EqualTo(0.0));
        Assert.That(pepCRow[2], Is.EqualTo(0.0));
        Assert.That(pepCRow[3], Is.EqualTo(1000.0));
        Assert.That(pepCRow[4], Is.EqualTo(1100.0));
        Assert.That(pepCRow[5], Is.EqualTo(1200.0));
    }

    /// <summary>
    /// Tests the full TMT pipeline end-to-end with known values.
    /// Uses NoNormalization, SumRollUp, NoCollapse so final values are predictable sums.
    /// </summary>
    [Test]
    public void RunTmt_FullPipeline_ProducesCorrectProteinMatrix()
    {
        // Setup: 2 files, 3 channels, 2 proteins with 1 unique peptide each, 1 PSM per peptide per file
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2 }, channels, conditions, bioReps, techReps);

        // Create proteins and peptides
        var protein1 = new Protein("PEPTIDEK", "P001");
        var protein2 = new Protein("SEQUENCEK", "P002");
        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
        var pep1 = protein1.Digest(digParams, new List<Modification>(), new List<Modification>()).First();
        var pep2 = protein2.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        var allPeptides = new List<IBioPolymerWithSetMods> { pep1, pep2 };

        var group1 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 });
        var group2 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein2 },
            new HashSet<IBioPolymerWithSetMods> { pep2 },
            new HashSet<IBioPolymerWithSetMods> { pep2 });
        var proteinGroups = new List<IBioPolymerGroup> { group1, group2 };

        // PSMs: 1 per peptide per file, each with known QuantValues
        // Protein1/pep1: file1=[100, 200, 300], file2=[1000, 2000, 3000]
        // Protein2/pep2: file1=[400, 500, 600], file2=[4000, 5000, 6000]
        var sm1 = new BaseSpectralMatch(file1, 1, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 100, 200, 300 } };
        var sm2 = new BaseSpectralMatch(file1, 2, 100.0, pep2.FullSequence, pep2.BaseSequence, new[] { pep2 })
            { QuantValues = new double[] { 400, 500, 600 } };
        var sm3 = new BaseSpectralMatch(file2, 3, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 1000, 2000, 3000 } };
        var sm4 = new BaseSpectralMatch(file2, 4, 100.0, pep2.FullSequence, pep2.BaseSequence, new[] { pep2 })
            { QuantValues = new double[] { 4000, 5000, 6000 } };

        var spectralMatches = new List<ISpectralMatch> { sm1, sm2, sm3, sm4 };

        var parameters = QuantificationParameters.GetSimpleParameters();
        parameters.WriteRawInformation = false;
        parameters.WritePeptideInformation = false;
        parameters.WriteProteinInformation = false;

        var engine = new QuantificationEngine(parameters, design, spectralMatches, allPeptides, proteinGroups);

        // Act
        var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

        // Assert dimensions: 2 proteins × 6 columns (3 channels × 2 files)
        Assert.That(proteinMatrix.RowKeys.Count, Is.EqualTo(2));
        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(6));

        // Column order: file1 channels (0-2), file2 channels (3-5) — alphabetical file order
        // With NoNormalization + SumRollUp + NoCollapse: protein value = single peptide value = single PSM value

        // Protein1: file1=[100, 200, 300], file2=[1000, 2000, 3000]
        var p1Row = proteinMatrix.GetRow(group1);
        Assert.That(p1Row[0], Is.EqualTo(100.0));    // file1, ch126
        Assert.That(p1Row[1], Is.EqualTo(200.0));    // file1, ch127N
        Assert.That(p1Row[2], Is.EqualTo(300.0));    // file1, ch127C
        Assert.That(p1Row[3], Is.EqualTo(1000.0));   // file2, ch126
        Assert.That(p1Row[4], Is.EqualTo(2000.0));   // file2, ch127N
        Assert.That(p1Row[5], Is.EqualTo(3000.0));   // file2, ch127C

        // Protein2: file1=[400, 500, 600], file2=[4000, 5000, 6000]
        var p2Row = proteinMatrix.GetRow(group2);
        Assert.That(p2Row[0], Is.EqualTo(400.0));
        Assert.That(p2Row[1], Is.EqualTo(500.0));
        Assert.That(p2Row[2], Is.EqualTo(600.0));
        Assert.That(p2Row[3], Is.EqualTo(4000.0));
        Assert.That(p2Row[4], Is.EqualTo(5000.0));
        Assert.That(p2Row[5], Is.EqualTo(6000.0));
    }

    /// <summary>
    /// Tests that SumCollapse combines technical replicates correctly.
    /// 3 files (tech reps) × 3 channels → 9 columns before collapse → 3 columns after.
    /// Each collapsed column = sum of that condition across all 3 technical replicates.
    /// </summary>
    [Test]
    public void RunTmt_WithSumCollapse_CombinesTechnicalReplicates()
    {
        // Setup: 3 files (technical replicates 1, 2, 3), same conditions, 3 channels
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string file3 = "file3.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2, 3 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2, file3 }, channels, conditions, bioReps, techReps);

        // Create protein and peptide
        var protein1 = new Protein("PEPTIDEK", "P001");
        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
        var pep1 = protein1.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        var allPeptides = new List<IBioPolymerWithSetMods> { pep1 };
        var group1 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 });
        var proteinGroups = new List<IBioPolymerGroup> { group1 };

        // Each tech rep gives the same intensities: [100, 200, 300]
        // After SumCollapse across 3 tech reps: [300, 600, 900]
        var sm1 = new BaseSpectralMatch(file1, 1, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 100, 200, 300 } };
        var sm2 = new BaseSpectralMatch(file2, 2, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 100, 200, 300 } };
        var sm3 = new BaseSpectralMatch(file3, 3, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 100, 200, 300 } };

        var spectralMatches = new List<ISpectralMatch> { sm1, sm2, sm3 };

        // Use SumCollapse to combine the 3 technical replicates
        var parameters = new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = new NoNormalization(),
            SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
            PeptideNormalizationStrategy = new NoNormalization(),
            CollapseStrategy = new SumCollapse(),
            PeptideToProteinRollUpStrategy = new SumRollUp(),
            ProteinNormalizationStrategy = new NoNormalization(),
            OutputDirectory = string.Empty,
            UseSharedPeptidesForProteinQuant = false,
            WriteRawInformation = false,
            WritePeptideInformation = false,
            WriteProteinInformation = false
        };

        var engine = new QuantificationEngine(parameters, design, spectralMatches, allPeptides, proteinGroups);

        // Act
        var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

        // Assert: columns collapsed from 9 (3 files × 3 channels) to 3 (unique conditions)
        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(3));

        // SumCollapse groups by Condition + '_' + BiologicalReplicate:
        //   Reference_1, CondA_1, CondB_1
        // Each group sums 3 tech reps: 100×3=300, 200×3=600, 300×3=900
        var p1Row = proteinMatrix.GetRow(group1);
        Assert.That(p1Row[0], Is.EqualTo(300.0)); // Reference_1: 100+100+100
        Assert.That(p1Row[1], Is.EqualTo(600.0)); // CondA_1:     200+200+200
        Assert.That(p1Row[2], Is.EqualTo(900.0)); // CondB_1:     300+300+300
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Real SpikeIn-5mix-MS3 integration tests
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Returns the path to the solution root directory (the mzLib/ folder containing
    /// TMT_Spike-In_Info/) by walking up from the NUnit test output directory until
    /// the directory containing TMT_Spike-In_Info/ is found. This is robust to
    /// varying output path depths (Debug vs Release, with or without framework suffix).
    /// </summary>
    private static string GetSolutionDir()
    {
        var dir = new System.IO.DirectoryInfo(TestContext.CurrentContext.TestDirectory);
        while (dir != null)
        {
            if (System.IO.Directory.Exists(System.IO.Path.Combine(dir.FullName, "TMT_Spike-In_Info")))
                return dir.FullName;
            dir = dir.Parent;
        }
        throw new System.IO.DirectoryNotFoundException(
            $"Could not find TMT_Spike-In_Info directory starting from: {TestContext.CurrentContext.TestDirectory}");
    }

    /// <summary>
    /// Loads the real UPS spike-in PSM file, runs the full TMT quantification pipeline,
    /// and performs basic sanity checks on the resulting protein matrix.
    /// </summary>
    [Test]
    public void LoadAndRunSpikeInData_BasicPipeline()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(upsPsmPath), Is.True, $"UPS psmtsv not found at: {upsPsmPath}");

        // Build all quantification inputs from the PSM file
        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputs(upsPsmPath);

        Assert.That(spectralMatches.Count, Is.GreaterThan(0), "No spectral matches loaded");
        Assert.That(peptides.Count, Is.GreaterThan(0), "No peptides created");
        Assert.That(proteinGroups.Count, Is.GreaterThan(0), "No protein groups created");

        var design = new SpikeInExperimentalDesign();
        var parameters = QuantificationParameters.GetSimpleParameters();
        parameters.WriteRawInformation = false;
        parameters.WritePeptideInformation = false;
        parameters.WriteProteinInformation = false;

        var engine = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups);
        var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

        // Basic structure checks
        Assert.That(proteinMatrix, Is.Not.Null);
        Assert.That(proteinMatrix.RowKeys.Count, Is.GreaterThan(0), "Protein matrix has no rows");

        // 14 files × 10 channels = 140 columns (with NoCollapse)
        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(140),
            $"Expected 140 columns (14 files × 10 channels), got {proteinMatrix.ColumnKeys.Count}");

        // At least some non-zero values
        bool anyNonZero = proteinMatrix.RowKeys
            .Any(p => proteinMatrix.GetRow(p).Any(v => v > 0));
        Assert.That(anyNonZero, Is.True, "All protein matrix values are zero");
    }

    /// <summary>
    /// Verifies that UPS proteins show fold changes in the correct direction:
    /// higher-concentration conditions should have higher quantified intensities.
    /// Uses NoNormalization + SumRollUp, so absolute fold changes may deviate from
    /// the true values (1.5 / 2.0 / 8.0 etc.), but direction should be preserved.
    /// </summary>
    [Test]
    public void EvaluateFoldChanges_UPSProteins()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(upsPsmPath), Is.True, $"UPS psmtsv not found at: {upsPsmPath}");

        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputs(upsPsmPath);

        var design = new SpikeInExperimentalDesign();
        var parameters = QuantificationParameters.GetSimpleParameters();
        parameters.WriteRawInformation = false;
        parameters.WritePeptideInformation = false;
        parameters.WriteProteinInformation = false;

        var engine = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups);
        var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

        // Compute fold changes across all quantifiable proteins: condition "1" vs "0.125" (expected 8×)
        var foldChanges_1_vs_0125 = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.125")
            .Select(x => x.foldChange)
            .ToList();

        // Compute fold changes: condition "1" vs "0.5" (expected 2×)
        var foldChanges_1_vs_0500 = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.5")
            .Select(x => x.foldChange)
            .ToList();

        Assert.That(foldChanges_1_vs_0125.Count, Is.GreaterThan(0),
            "No proteins had quantifiable intensities in both '1' and '0.125' conditions");

        // Median fold change should be > 1.5 (direction correct; true value is 8.0)
        // This is a relaxed sanity check — without normalization, ratio compression is expected
        double median_1_vs_0125 = QuantificationEvaluator.Median(foldChanges_1_vs_0125);
        Assert.That(median_1_vs_0125, Is.GreaterThan(1.5),
            $"Median fold change '1' vs '0.125' was {median_1_vs_0125:F2}, expected > 1.5 (true value is 8.0)");

        if (foldChanges_1_vs_0500.Count > 0)
        {
            double median_1_vs_0500 = QuantificationEvaluator.Median(foldChanges_1_vs_0500);
            Assert.That(median_1_vs_0500, Is.GreaterThan(1.0),
                $"Median fold change '1' vs '0.5' was {median_1_vs_0500:F2}, expected > 1.0 (true value is 2.0)");
        }
    }

    /// <summary>
    /// Verifies that HeLa background proteins show fold changes close to 1.0 across conditions.
    /// HeLa proteins are present at constant amounts across all channels and should not
    /// appear differentially abundant. Marked [Explicit] because the HeLa search file is large.
    /// </summary>
    [Test, Explicit("Processes large HeLa psmtsv file (~1.35M lines); run explicitly to evaluate background protein stability.")]
    public void EvaluateFoldChanges_HelaBackground()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string helaPsmPath = Path.Combine(dataDir, "TMT3_HeLa_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(helaPsmPath), Is.True, $"HeLa psmtsv not found at: {helaPsmPath}");

        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputs(helaPsmPath);

        Assert.That(spectralMatches.Count, Is.GreaterThan(0), "No HeLa spectral matches loaded");

        var design = new SpikeInExperimentalDesign();
        var parameters = QuantificationParameters.GetSimpleParameters();
        parameters.WriteRawInformation = false;
        parameters.WritePeptideInformation = false;
        parameters.WriteProteinInformation = false;

        var engine = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups);
        var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

        Assert.That(proteinMatrix.RowKeys.Count, Is.GreaterThan(0), "HeLa protein matrix is empty");

        // For HeLa background: fold change '1' vs '0.5' should be close to 1.0
        var foldChanges = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.5")
            .Select(x => x.foldChange)
            .ToList();

        Assert.That(foldChanges.Count, Is.GreaterThan(0),
            "No HeLa proteins had quantifiable intensities in both '1' and '0.5' conditions");

        double medianFc = QuantificationEvaluator.Median(foldChanges);

        // Without normalization, expect the median HeLa fold change to be between 0.3 and 3.0
        // (a very loose bound — normalization would tighten this toward 1.0)
        Assert.That(medianFc, Is.InRange(0.3, 3.0),
            $"Median HeLa fold change '1' vs '0.5' was {medianFc:F3}, expected in [0.3, 3.0]");
    }
}
