using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Quantification.Interfaces;
using Quantification.Strategies;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
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
        string helaPsmPath = Path.Combine(dataDir, "HeLa_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
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

    // ──────────────────────────────────────────────────────────────────────────
    // Task 10: GlobalMedianNormalization unit test
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Verifies that GlobalMedianNormalization equalizes column medians in log2 space
    /// and that zero values are preserved unchanged.
    /// </summary>
    [Test]
    public void GlobalMedianNormalization_EqualizesColumnMedians()
    {
        // Setup: 2 rows × 3 columns. Before normalization the column medians (in log2 space) differ:
        //   Col 0: [100, 100] → log2 median = log2(100) ≈ 6.644
        //   Col 1: [200, 200] → log2 median = log2(200) ≈ 7.644  (global median)
        //   Col 2: [400, 400] → log2 median = log2(400) ≈ 8.644
        // Global median of [6.644, 7.644, 8.644] = 7.644
        // After normalization every column median should equal log2(200).
        var design = CreateTmtExperimentalDesign(
            new[] { "file1.raw" },
            new[] { "126", "127N", "127C" },
            new[] { "CondA", "CondB", "CondC" },
            new[] { 1, 1, 1 }, new[] { 1 });
        var columns = design.FileNameSampleInfoDictionary["file1.raw"];

        var matrix = new QuantMatrix<string>(new[] { "Entity1", "Entity2" }, columns, design);
        matrix.SetRow("Entity1", new double[] { 100, 200, 400 });
        matrix.SetRow("Entity2", new double[] { 100, 200, 400 });

        var norm = new GlobalMedianNormalization();
        var result = norm.NormalizeIntensities(matrix);

        // After normalization all columns should contain [200, 200].
        Assert.That(result.Matrix[0, 0], Is.EqualTo(200.0).Within(1e-6), "Col0 row0 should be 200");
        Assert.That(result.Matrix[1, 0], Is.EqualTo(200.0).Within(1e-6), "Col0 row1 should be 200");
        Assert.That(result.Matrix[0, 1], Is.EqualTo(200.0).Within(1e-6), "Col1 row0 should be unchanged");
        Assert.That(result.Matrix[1, 1], Is.EqualTo(200.0).Within(1e-6), "Col1 row1 should be unchanged");
        Assert.That(result.Matrix[0, 2], Is.EqualTo(200.0).Within(1e-6), "Col2 row0 should be 200");
        Assert.That(result.Matrix[1, 2], Is.EqualTo(200.0).Within(1e-6), "Col2 row1 should be 200");

        // Zeros must remain zero after normalization.
        var matrixWithZero = new QuantMatrix<string>(new[] { "E1", "E2", "E3" }, columns, design);
        matrixWithZero.SetRow("E1", new double[] { 100, 0, 400 });
        matrixWithZero.SetRow("E2", new double[] { 200, 500, 600 });
        matrixWithZero.SetRow("E3", new double[] { 400, 800, 1600 });
        var resultWithZero = norm.NormalizeIntensities(matrixWithZero);
        Assert.That(resultWithZero.Matrix[0, 1], Is.EqualTo(0.0), "Zero should remain zero after normalization");
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Task 11: ReferenceChannelNormalization unit test
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Verifies that ReferenceChannelNormalization converts intensities to ratios relative to
    /// the mean of the reference channels in the same file.
    /// </summary>
    [Test]
    public void ReferenceChannelNormalization_ConvertsToRatios()
    {
        // Setup: 1 file, 3 channels (126=reference, 127N=non-reference, 131N=reference).
        // Using "Reference" condition name so CreateTmtExperimentalDesign sets IsReferenceChannel=true.
        var design = CreateTmtExperimentalDesign(
            new[] { "file1.raw" },
            new[] { "126", "127N", "131N" },
            new[] { "Reference", "CondA", "Reference" },
            new[] { 1, 1, 1 }, new[] { 1 });
        var columns = design.FileNameSampleInfoDictionary["file1.raw"];

        // Row 0: ref1=100, nonRef=200, ref2=100  → refMean=100, ratio=2.0
        // Row 1: ref1=200, nonRef=400, ref2=200  → refMean=200, ratio=2.0
        var matrix = new QuantMatrix<string>(new[] { "Entity1", "Entity2" }, columns, design);
        matrix.SetRow("Entity1", new double[] { 100, 200, 100 });
        matrix.SetRow("Entity2", new double[] { 200, 400, 200 });

        var norm = new ReferenceChannelNormalization();
        var result = norm.NormalizeIntensities(matrix);

        // Reference channels (col 0 and col 2) → 1.0
        Assert.That(result.Matrix[0, 0], Is.EqualTo(1.0).Within(1e-9), "Entity1 ref ch 126 should be 1.0");
        Assert.That(result.Matrix[0, 2], Is.EqualTo(1.0).Within(1e-9), "Entity1 ref ch 131N should be 1.0");
        Assert.That(result.Matrix[1, 0], Is.EqualTo(1.0).Within(1e-9), "Entity2 ref ch 126 should be 1.0");
        Assert.That(result.Matrix[1, 2], Is.EqualTo(1.0).Within(1e-9), "Entity2 ref ch 131N should be 1.0");

        // Non-reference (col 1): 200/100=2.0, 400/200=2.0
        Assert.That(result.Matrix[0, 1], Is.EqualTo(2.0).Within(1e-9), "Entity1 non-ref ch 127N should be 2.0");
        Assert.That(result.Matrix[1, 1], Is.EqualTo(2.0).Within(1e-9), "Entity2 non-ref ch 127N should be 2.0");
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Task 12: MedianRollUp unit test
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Verifies that MedianRollUp computes the per-column median across rolled-up rows,
    /// correctly ignoring zero values (missing data).
    /// </summary>
    [Test]
    public void MedianRollUp_ComputesMedianPerColumn()
    {
        var design = CreateTmtExperimentalDesign(
            new[] { "file1.raw" },
            new[] { "126", "127N" },
            new[] { "CondA", "CondB" },
            new[] { 1, 1 }, new[] { 1 });
        var columns = design.FileNameSampleInfoDictionary["file1.raw"];

        // 3 PSM rows → 1 protein, 2 columns
        var inputMatrix = new QuantMatrix<string>(
            new[] { "PSM1", "PSM2", "PSM3" }, columns, design);
        inputMatrix.SetRow("PSM1", new double[] { 100, 200 });
        inputMatrix.SetRow("PSM2", new double[] { 300, 400 });
        inputMatrix.SetRow("PSM3", new double[] { 500, 100 });

        var map = new Dictionary<string, List<int>>
        {
            { "Protein1", new List<int> { 0, 1, 2 } }
        };

        var rollUp = new MedianRollUp();
        var result = rollUp.RollUp(inputMatrix, map);

        // Col 0: median([100, 300, 500]) = 300
        Assert.That(result.Matrix[0, 0], Is.EqualTo(300.0).Within(1e-6), "Median of [100,300,500] should be 300");
        // Col 1: median([200, 400, 100]) = sorted [100,200,400] → middle = 200
        Assert.That(result.Matrix[0, 1], Is.EqualTo(200.0).Within(1e-6), "Median of [200,400,100] should be 200");

        // Test with zeros: zeros are excluded from the median.
        var inputMatrixWithZeros = new QuantMatrix<string>(
            new[] { "PSM1", "PSM2", "PSM3" }, columns, design);
        inputMatrixWithZeros.SetRow("PSM1", new double[] { 100, 0 });
        inputMatrixWithZeros.SetRow("PSM2", new double[] { 300, 400 });
        inputMatrixWithZeros.SetRow("PSM3", new double[] { 0, 100 });

        var resultWithZeros = rollUp.RollUp(inputMatrixWithZeros, map);

        // Col 0: non-zero values [100, 300] → median = (100+300)/2 = 200
        Assert.That(resultWithZeros.Matrix[0, 0], Is.EqualTo(200.0).Within(1e-6), "Median of [100,300] (ignoring zeros) should be 200");
        // Col 1: non-zero values [400, 100] → sorted [100,400] → median = (100+400)/2 = 250
        Assert.That(resultWithZeros.Matrix[0, 1], Is.EqualTo(250.0).Within(1e-6), "Median of [400,100] (ignoring zeros) should be 250");
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Task 13: MeanCollapse unit test
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Verifies that MeanCollapse averages technical replicates, producing values on the
    /// same scale as individual measurements (unlike SumCollapse which produces inflated sums).
    /// </summary>
    [Test]
    public void MeanCollapse_AveragesTechnicalReplicates()
    {
        // Same setup as RunTmt_WithSumCollapse_CombinesTechnicalReplicates:
        // 3 files (tech reps 1, 2, 3), 3 channels, 1 protein with values [100, 200, 300] per file.
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string file3 = "file3.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2, 3 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2, file3 }, channels, conditions, bioReps, techReps);

        var protein1 = new Protein("PEPTIDEK", "P001");
        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
        var pep1 = protein1.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        var allPeptides = new List<IBioPolymerWithSetMods> { pep1 };
        var group1 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 });
        var proteinGroups = new List<IBioPolymerGroup> { group1 };

        // Each tech rep contributes identical intensities [100, 200, 300].
        var sm1 = new BaseSpectralMatch(file1, 1, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 100, 200, 300 } };
        var sm2 = new BaseSpectralMatch(file2, 2, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 100, 200, 300 } };
        var sm3 = new BaseSpectralMatch(file3, 3, 100.0, pep1.FullSequence, pep1.BaseSequence, new[] { pep1 })
            { QuantValues = new double[] { 100, 200, 300 } };

        var spectralMatches = new List<ISpectralMatch> { sm1, sm2, sm3 };

        var parameters = new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = new NoNormalization(),
            SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
            PeptideNormalizationStrategy = new NoNormalization(),
            CollapseStrategy = new MeanCollapse(),
            PeptideToProteinRollUpStrategy = new SumRollUp(),
            ProteinNormalizationStrategy = new NoNormalization(),
            OutputDirectory = string.Empty,
            UseSharedPeptidesForProteinQuant = false,
            WriteRawInformation = false,
            WritePeptideInformation = false,
            WriteProteinInformation = false
        };

        var engine = new QuantificationEngine(parameters, design, spectralMatches, allPeptides, proteinGroups);
        var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

        // After MeanCollapse, 9 columns (3 files × 3 channels) collapse to 3 groups.
        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(3),
            "Expected 3 columns after MeanCollapse (Reference_1, CondA_1, CondB_1)");

        // mean(100,100,100)/3 would be sum/3 — but group size=3 so mean = 300/3 = 100
        var p1Row = proteinMatrix.GetRow(group1);
        Assert.That(p1Row[0], Is.EqualTo(100.0).Within(1e-6), "Reference_1: mean(100,100,100)=100");
        Assert.That(p1Row[1], Is.EqualTo(200.0).Within(1e-6), "CondA_1: mean(200,200,200)=200");
        Assert.That(p1Row[2], Is.EqualTo(300.0).Within(1e-6), "CondB_1: mean(300,300,300)=300");
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Task 9: Pipeline helper + baseline metrics tests
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Loads PSMs from a psmtsv file, configures a QuantificationEngine with the given strategies,
    /// runs the TMT pipeline, and returns the protein matrix.
    /// </summary>
    private static QuantMatrix<IBioPolymerGroup> RunPipelineWithStrategies(
        string psmtsvPath,
        INormalizationStrategy psmNorm,
        IRollUpStrategy psmToPeptideRollUp,
        INormalizationStrategy peptideNorm,
        ICollapseStrategy collapse,
        IRollUpStrategy peptideToProteinRollUp,
        INormalizationStrategy proteinNorm)
    {
        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputs(psmtsvPath);
        return RunPipelineWithStrategies(
            spectralMatches, peptides, proteinGroups,
            psmNorm, psmToPeptideRollUp, peptideNorm, collapse, peptideToProteinRollUp, proteinNorm);
    }

    /// <summary>
    /// Configures a QuantificationEngine with the given strategies using pre-loaded inputs,
    /// runs the TMT pipeline, and returns the protein matrix. Use this overload to avoid
    /// re-loading the PSM file for every combination in the strategy evaluation test.
    /// </summary>
    private static QuantMatrix<IBioPolymerGroup> RunPipelineWithStrategies(
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> peptides,
        List<IBioPolymerGroup> proteinGroups,
        INormalizationStrategy psmNorm,
        IRollUpStrategy psmToPeptideRollUp,
        INormalizationStrategy peptideNorm,
        ICollapseStrategy collapse,
        IRollUpStrategy peptideToProteinRollUp,
        INormalizationStrategy proteinNorm)
    {
        var design = new SpikeInExperimentalDesign();
        var parameters = new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = psmNorm,
            SpectralMatchToPeptideRollUpStrategy = psmToPeptideRollUp,
            PeptideNormalizationStrategy = peptideNorm,
            CollapseStrategy = collapse,
            PeptideToProteinRollUpStrategy = peptideToProteinRollUp,
            ProteinNormalizationStrategy = proteinNorm,
            OutputDirectory = string.Empty,
            WriteRawInformation = false,
            WritePeptideInformation = false,
            WriteProteinInformation = false
        };
        var engine = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups);
        return engine.RunTmtAndReturnProteinMatrix();
    }

    /// <summary>
    /// Runs the NoNormalization + SumRollUp + NoCollapse baseline pipeline on UPS spike-in data
    /// and prints accuracy metrics for all three expected fold change comparisons.
    /// No assertions are made — this test documents the baseline for strategy comparison.
    /// </summary>
    [Test]
    public void BaselineMetrics_NoNormalization_UPS()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(upsPsmPath), Is.True, $"UPS psmtsv not found at: {upsPsmPath}");

        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputs(upsPsmPath);

        var proteinMatrix = RunPipelineWithStrategies(
            spectralMatches, peptides, proteinGroups,
            new NoNormalization(), new SumRollUp(),
            new NoNormalization(), new NoCollapse(),
            new SumRollUp(), new NoNormalization());

        var comparisons = new[]
        {
            ("1 vs 0.125", "1", "0.125", 8.0),
            ("1 vs 0.5",   "1", "0.5",   2.0),
            ("1 vs 0.667", "1", "0.667", 1.5),
        };

        TestContext.WriteLine("=== Baseline Metrics (NoNormalization + SumRollUp + NoCollapse) ===");

        foreach (var (label, numerator, denominator, expectedFc) in comparisons)
        {
            var foldChanges = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, numerator, denominator)
                .Select(x => x.foldChange)
                .ToList();

            if (foldChanges.Count == 0)
            {
                TestContext.WriteLine($"\nComparison: {label} (expected FC = {expectedFc}) — no quantifiable proteins");
                continue;
            }

            double medianFc  = QuantificationEvaluator.Median(foldChanges);
            double meanFc    = foldChanges.Average();
            double mae       = QuantificationEvaluator.MeanAbsoluteError(foldChanges, expectedFc);
            double maeLog2   = QuantificationEvaluator.MeanAbsoluteLog2Error(foldChanges, expectedFc);
            double frac2x    = QuantificationEvaluator.FractionWithinFactor(foldChanges, expectedFc, 2.0);

            TestContext.WriteLine($"\nComparison: {label} (expected FC = {expectedFc})");
            TestContext.WriteLine($"  Proteins quantified: {foldChanges.Count}");
            TestContext.WriteLine($"  Median FC:           {medianFc:F3}");
            TestContext.WriteLine($"  Mean FC:             {meanFc:F3}");
            TestContext.WriteLine($"  MAE:                 {mae:F3}");
            TestContext.WriteLine($"  MAE(log2):           {maeLog2:F3}");
            TestContext.WriteLine($"  Fraction within 2x:  {frac2x:F3}");
        }
    }

    /// <summary>
    /// Baseline metrics for HeLa background proteins (expected fold change = 1.0).
    /// Marked explicit because the HeLa file is large (~1.35M lines, ~90s loading time).
    /// </summary>
    [Test, Explicit("Large HeLa file; run explicitly to evaluate background protein stability")]
    public void BaselineMetrics_NoNormalization_HeLa()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string helaPsmPath = Path.Combine(dataDir, "HeLa_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(helaPsmPath), Is.True, $"HeLa psmtsv not found at: {helaPsmPath}");

        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputs(helaPsmPath);

        var proteinMatrix = RunPipelineWithStrategies(
            spectralMatches, peptides, proteinGroups,
            new NoNormalization(), new SumRollUp(),
            new NoNormalization(), new NoCollapse(),
            new SumRollUp(), new NoNormalization());

        var foldChanges = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.5")
            .Select(x => x.foldChange)
            .ToList();

        TestContext.WriteLine("=== Baseline Metrics HeLa (NoNormalization + SumRollUp + NoCollapse) ===");
        TestContext.WriteLine($"Comparison: 1 vs 0.5 (expected FC = 1.0)");
        TestContext.WriteLine($"  Proteins quantified:   {foldChanges.Count}");

        if (foldChanges.Count > 0)
        {
            double medianFc  = QuantificationEvaluator.Median(foldChanges);
            double mae       = QuantificationEvaluator.MeanAbsoluteError(foldChanges, 1.0);
            double frac15x   = QuantificationEvaluator.FractionWithinFactor(foldChanges, 1.0, 1.5);
            double frac2x    = QuantificationEvaluator.FractionWithinFactor(foldChanges, 1.0, 2.0);

            TestContext.WriteLine($"  Median FC:             {medianFc:F3}");
            TestContext.WriteLine($"  MAE:                   {mae:F3}");
            TestContext.WriteLine($"  Fraction within 1.5x:  {frac15x:F3}");
            TestContext.WriteLine($"  Fraction within 2.0x:  {frac2x:F3}");
        }
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Protein overlap analysis: HeLa vs UPS
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Reads both the HeLa and UPS AllPSMs files, extracts unique protein accessions from each,
    /// and reports which proteins are shared vs distinct to each search.
    /// </summary>
    [Test]
    public void ProteinOverlap_HelaVsUPS()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        string helaPsmPath = Path.Combine(dataDir, "HeLa_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(upsPsmPath), Is.True, $"UPS psmtsv not found at: {upsPsmPath}");
        Assert.That(File.Exists(helaPsmPath), Is.True, $"HeLa psmtsv not found at: {helaPsmPath}");

        // Read target PSMs (q <= 0.01) from each file and collect unique accessions
        // Use the first non-decoy accession from pipe-delimited accession strings
        var upsRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(upsPsmPath, out _)
            .Where(r => r.QValue <= 0.01 && r.DecoyContamTarget.Contains('T')
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

        var helaRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(helaPsmPath, out _)
            .Where(r => r.QValue <= 0.01 && r.DecoyContamTarget.Contains('T')
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

        // Extract unique accessions (first non-decoy accession from each pipe-delimited string)
        var upsAccessions = new HashSet<string>(
            upsRecords.Select(r => GetFirstNonDecoyAccession(r.Accession))
                      .Where(a => a != null));

        var helaAccessions = new HashSet<string>(
            helaRecords.Select(r => GetFirstNonDecoyAccession(r.Accession))
                       .Where(a => a != null));

        var shared = new HashSet<string>(upsAccessions);
        shared.IntersectWith(helaAccessions);

        var upsOnly = new HashSet<string>(upsAccessions);
        upsOnly.ExceptWith(helaAccessions);

        var helaOnly = new HashSet<string>(helaAccessions);
        helaOnly.ExceptWith(upsAccessions);

        // Print summary
        TestContext.WriteLine("=== Protein Overlap: HeLa vs UPS ===");
        TestContext.WriteLine($"UPS total accessions:  {upsAccessions.Count}");
        TestContext.WriteLine($"HeLa total accessions: {helaAccessions.Count}");
        TestContext.WriteLine($"Shared:                {shared.Count}");
        TestContext.WriteLine($"UPS-only:              {upsOnly.Count}");
        TestContext.WriteLine($"HeLa-only:             {helaOnly.Count}");

        // Print shared proteins
        TestContext.WriteLine($"\n--- Shared proteins ({shared.Count}) ---");
        foreach (var acc in shared.OrderBy(a => a))
        {
            TestContext.WriteLine($"  {acc}");
        }

        // Print UPS-only proteins
        TestContext.WriteLine($"\n--- UPS-only proteins ({upsOnly.Count}) ---");
        foreach (var acc in upsOnly.OrderBy(a => a))
        {
            TestContext.WriteLine($"  {acc}");
        }

        // Print HeLa-only proteins (just count + first 50 to avoid massive output)
        TestContext.WriteLine($"\n--- HeLa-only proteins ({helaOnly.Count}) ---");
        foreach (var acc in helaOnly.OrderBy(a => a).Take(50))
        {
            TestContext.WriteLine($"  {acc}");
        }
        if (helaOnly.Count > 50)
            TestContext.WriteLine($"  ... and {helaOnly.Count - 50} more");

        // Basic sanity assertions
        Assert.That(upsAccessions.Count, Is.GreaterThan(0), "No UPS accessions found");
        Assert.That(helaAccessions.Count, Is.GreaterThan(0), "No HeLa accessions found");
    }

    public static string[] UPSOnlyAccessions =
    {
        "O00762",
        "P00167",
        "P00709",
        "P00915",
        "P01008",
        "P01031",
        "P01112",
        "P01127",
        "P01133",
        "P01344",
        "P01375",
        "P01579",
        "P02144",
        "P02741",
        "P02753",
        "P02768",
        "P02787",
        "P02788",
        "P05413",
        "P08263",
        "P10145",
        "P10636-8",
        "P61626",
        "P62988",
        "P63165",
        "P68871",
        "P69905"
    };

    public static string[] SharedAccessions =
    {
        "O76070",
        "P00441",
        "P00918",
        "P04040",
        "P06396",
        "P06732",
        "P08758",
        "P09211",
        "P10599",
        "P12081",
        "P15559",
        "P16083",
        "P41159",
        "P51965",
        "P55957",
        "P62937",
        "P63279",
        "P99999",
        "Q06830",
        "Q15843"
    };

    private static string GetFirstNonDecoyAccession(string rawAccession)
    {
        var parts = rawAccession.Split('|');
        return parts.FirstOrDefault(a => !a.StartsWith("DECOY_")) ?? parts[0];
    }

    [Test]
    public void CombinedQuantHelaUPS()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string upsDir = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask");
        string helaDir = Path.Combine(dataDir, "HeLa_TMT3_Search", "Task1-SearchTask");
        Assert.That(Directory.Exists(upsDir), Is.True, $"UPS directory not found at: {upsDir}");
        Assert.That(Directory.Exists(helaDir), Is.True, $"HeLa directory not found at: {helaDir}");

        var sw = new Stopwatch();
        sw.Start();

        // Load combined inputs from both search directories
        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputsFromMultipleDirectories(
                new[] { upsDir, helaDir });

           
        sw.Stop();
        TestContext.WriteLine("All data loaded in :" + sw.Elapsed);
            // Write out time it took to load in data


        TestContext.WriteLine($"Combined: {spectralMatches.Count} PSMs, {peptides.Count} peptides, {proteinGroups.Count} protein groups");

        // Run quantification with a single strategy set (no combinatorial sweep)
        var design = new SpikeInExperimentalDesign();
        var parameters = new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = new NoNormalization(),
            SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
            PeptideNormalizationStrategy = new NoNormalization(),
            CollapseStrategy = new NoCollapse(),
            PeptideToProteinRollUpStrategy = new SumRollUp(),
            ProteinNormalizationStrategy = new NoNormalization(),
            OutputDirectory = string.Empty,
            UseSharedPeptidesForProteinQuant = false,
            WriteRawInformation = false,
            WritePeptideInformation = false,
            WriteProteinInformation = false
        };

        var engine = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups);
        var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

        Assert.That(proteinMatrix.RowKeys.Count, Is.GreaterThan(0), "Combined protein matrix is empty");

        // Classify proteins by accession
        var upsOnlySet = new HashSet<string>(UPSOnlyAccessions);
        var sharedSet = new HashSet<string>(SharedAccessions);

        // All protein accessions in the result matrix
        var allAccessions = proteinMatrix.RowKeys
            .Select(pg => pg.BioPolymerGroupName)
            .ToHashSet();

        // HeLa-only = everything that's not UPS-only and not shared
        var helaOnlyAccessions = allAccessions
            .Where(a => !upsOnlySet.Contains(a) && !sharedSet.Contains(a))
            .ToHashSet();

        TestContext.WriteLine($"\nProtein classification in combined matrix:");
        TestContext.WriteLine($"  Total proteins:    {allAccessions.Count}");
        TestContext.WriteLine($"  UPS-only found:    {allAccessions.Count(a => upsOnlySet.Contains(a))}");
        TestContext.WriteLine($"  Shared found:      {allAccessions.Count(a => sharedSet.Contains(a))}");
        TestContext.WriteLine($"  HeLa-only found:   {helaOnlyAccessions.Count}");

        // Evaluate fold changes for each protein group
        var comparisons = new[]
        {
            ("1 vs 0.125", "1", "0.125", 8.0),
            ("1 vs 0.5",   "1", "0.5",   2.0),
            ("1 vs 0.667", "1", "0.667", 1.5),
        };

        foreach (var (label, numerator, denominator, expectedFc) in comparisons)
        {
            var allFoldChanges = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, numerator, denominator);

            var upsOnlyFCs = allFoldChanges.Where(x => upsOnlySet.Contains(x.accession)).Select(x => x.foldChange).ToList();
            var sharedFCs = allFoldChanges.Where(x => sharedSet.Contains(x.accession)).Select(x => x.foldChange).ToList();
            var helaOnlyFCs = allFoldChanges.Where(x => helaOnlyAccessions.Contains(x.accession)).Select(x => x.foldChange).ToList();

            TestContext.WriteLine($"\n=== {label} (expected UPS FC = {expectedFc}, expected HeLa FC = 1.0) ===");

            PrintFoldChangeGroup("UPS-only", upsOnlyFCs, expectedFc);
            PrintFoldChangeGroup("Shared (UPS+HeLa)", sharedFCs, expectedFc);
            PrintFoldChangeGroup("HeLa-only", helaOnlyFCs, 1.0);
        }

        // Basic assertions: UPS-only proteins should show directional fold change for 8x comparison
        var fc8x_upsOnly = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.125")
            .Where(x => upsOnlySet.Contains(x.accession))
            .Select(x => x.foldChange).ToList();

        if (fc8x_upsOnly.Count > 0)
        {
            double medianFc = QuantificationEvaluator.Median(fc8x_upsOnly);
            Assert.That(medianFc, Is.GreaterThan(1.5),
                $"UPS-only median FC '1' vs '0.125' was {medianFc:F2}, expected > 1.5 (true value = 8.0)");
        }

        // HeLa-only proteins should be close to 1.0 for any comparison
        var fc8x_helaOnly = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.125")
            .Where(x => helaOnlyAccessions.Contains(x.accession))
            .Select(x => x.foldChange).ToList();

        if (fc8x_helaOnly.Count > 0)
        {
            double medianFc = QuantificationEvaluator.Median(fc8x_helaOnly);
            Assert.That(medianFc, Is.InRange(0.3, 3.0),
                $"HeLa-only median FC '1' vs '0.125' was {medianFc:F3}, expected near 1.0");
        }
    }

    private static void PrintFoldChangeGroup(string groupName, List<double> foldChanges, double expectedFc)
    {
        TestContext.WriteLine($"\n  {groupName} ({foldChanges.Count} proteins):");
        if (foldChanges.Count == 0)
        {
            TestContext.WriteLine("    No quantifiable proteins in this group");
            return;
        }

        double medianFc = QuantificationEvaluator.Median(foldChanges);
        double meanFc = foldChanges.Average();
        double maeLog2 = QuantificationEvaluator.MeanAbsoluteLog2Error(foldChanges, expectedFc);
        double frac2x = QuantificationEvaluator.FractionWithinFactor(foldChanges, expectedFc, 2.0);

        TestContext.WriteLine($"    Median FC:           {medianFc:F3}");
        TestContext.WriteLine($"    Mean FC:             {meanFc:F3}");
        TestContext.WriteLine($"    MAE(log2):           {maeLog2:F3}");
        TestContext.WriteLine($"    Fraction within 2x:  {frac2x:F3}");
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Task 14: Combinatorial strategy evaluation
    // ──────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Evaluates all meaningful combinations of normalization, roll-up, and collapse strategies
    /// against UPS spike-in data and prints a ranked table sorted by combined MAE(log2).
    /// PSMs are loaded once to avoid repeated file I/O across the ~48 combinations.
    /// </summary>
    [Test]
    public void StrategyEvaluation_AllCombinations_UPS()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        //string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        //Assert.That(File.Exists(upsPsmPath), Is.True, $"UPS psmtsv not found at: {upsPsmPath}");

        //// Load PSMs once; reuse across all strategy combinations.
        //var (spectralMatches, peptides, proteinGroups) =
        //    PsmTsvQuantAdapter.BuildQuantificationInputs(upsPsmPath);

        var (spectralMatches, peptides, proteinGroups) = PsmTsvQuantAdapter.BuildQuantificationInputsFromDirectory(
            Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info", "UPS_TMT3_Search", "Task1-SearchTask"));

        // Strategy options
        var psmNorms = new INormalizationStrategy[]
        {
            new NoNormalization(),
            new GlobalMedianNormalization(),
            new ReferenceChannelNormalization()
        };
        var rollUps = new IRollUpStrategy[]
        {
            new SumRollUp(),
            new MedianRollUp()
        };
        var peptideNorms = new INormalizationStrategy[]
        {
            new NoNormalization(),
            new GlobalMedianNormalization()
        };
        var collapses = new ICollapseStrategy[]
        {
            new NoCollapse(),
            new MeanCollapse()
        };
        var proteinNorms = new INormalizationStrategy[]
        {
            new NoNormalization(),
            new GlobalMedianNormalization()
        };

        var results = new List<(string config, double mae_8x, double mae_2x, double medianFc_8x, double medianFc_2x)>();

        foreach (var psmNorm in psmNorms)
        foreach (var rollUp in rollUps)
        foreach (var pepNorm in peptideNorms)
        foreach (var collapse in collapses)
        foreach (var protNorm in proteinNorms)
        {
            string config = $"{psmNorm.Name} | {rollUp.Name} | {pepNorm.Name} | {collapse.Name} | {rollUp.Name} | {protNorm.Name}";

            try
            {
                var proteinMatrix = RunPipelineWithStrategies(
                    spectralMatches, peptides, proteinGroups,
                    psmNorm, rollUp, pepNorm, collapse, rollUp, protNorm);

                var fc_8x = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.125")
                    .Select(x => x.foldChange).ToList();
                var fc_2x = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.5")
                    .Select(x => x.foldChange).ToList();

                if (fc_8x.Count > 0 && fc_2x.Count > 0)
                {
                    double mae_8x    = QuantificationEvaluator.MeanAbsoluteLog2Error(fc_8x, 8.0);
                    double mae_2x    = QuantificationEvaluator.MeanAbsoluteLog2Error(fc_2x, 2.0);
                    double medFc_8x  = QuantificationEvaluator.Median(fc_8x);
                    double medFc_2x  = QuantificationEvaluator.Median(fc_2x);
                    results.Add((config, mae_8x, mae_2x, medFc_8x, medFc_2x));
                }
            }
            catch (Exception ex)
            {
                TestContext.WriteLine($"FAILED: {config} → {ex.Message}");
            }
        }

        // Sort by combined MAE(log2).
        results = results.OrderBy(r => r.mae_8x + r.mae_2x).ToList();

        TestContext.WriteLine("=== Strategy Evaluation Results (sorted by combined MAE_log2) ===");
        TestContext.WriteLine($"{"Rank",-5} {"MAE_log2(8x)",-14} {"MAE_log2(2x)",-14} {"MedFC(8x)",-12} {"MedFC(2x)",-12} Config");
        TestContext.WriteLine(new string('-', 120));

        for (int i = 0; i < results.Count; i++)
        {
            var r = results[i];
            TestContext.WriteLine($"{i + 1,-5} {r.mae_8x,-14:F3} {r.mae_2x,-14:F3} {r.medianFc_8x,-12:F3} {r.medianFc_2x,-12:F3} {r.config}");
        }

        Assert.That(results.Count, Is.GreaterThan(0), "No strategy combinations could be evaluated");

        if (results.Count > 1)
        {
            double bestScore  = results.First().mae_8x + results.First().mae_2x;
            double worstScore = results.Last().mae_8x  + results.Last().mae_2x;
            TestContext.WriteLine($"\nBest combined MAE_log2:  {bestScore:F3}");
            TestContext.WriteLine($"Worst combined MAE_log2: {worstScore:F3}");
        }
    }
}
