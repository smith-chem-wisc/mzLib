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
using Development.QuantificationDevelopment.TestHelpers;

namespace Development.QuantificationDevelopment;

/// <summary>
/// Integration tests that use the real TMT spike-in data files.
/// These tests reference local files (TMT_Spike-In_Info/) and should NOT be pushed to CI.
/// They exist to enable evaluation of normalization/roll-up/collapse strategies.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TmtSpikeInDevelopmentTests
{
    #region Helper Methods

    /// <summary>
    /// Returns the path to the solution root directory (the mzLib/ folder containing
    /// TMT_Spike-In_Info/) by walking up from the NUnit test output directory until
    /// the directory containing TMT_Spike-In_Info/ is found.
    /// </summary>
    private static string GetSolutionDir()
    {
        var dir = new DirectoryInfo(TestContext.CurrentContext.TestDirectory);
        while (dir != null)
        {
            if (Directory.Exists(Path.Combine(dir.FullName, "TMT_Spike-In_Info")))
                return dir.FullName;
            dir = dir.Parent;
        }
        throw new DirectoryNotFoundException(
            $"Could not find TMT_Spike-In_Info directory starting from: {TestContext.CurrentContext.TestDirectory}");
    }

    /// <summary>
    /// Runs the TMT pipeline with the given strategies using pre-loaded inputs.
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

    private static string GetFirstNonDecoyAccession(string rawAccession)
    {
        var parts = rawAccession.Split('|');
        return parts.FirstOrDefault(a => !a.StartsWith("DECOY_")) ?? parts[0];
    }

    #endregion

    #region Known Accessions

    public static string[] UPSOnlyAccessions =
    {
        "O00762", "P00167", "P00709", "P00915", "P01008", "P01031", "P01112", "P01127",
        "P01133", "P01344", "P01375", "P01579", "P02144", "P02741", "P02753", "P02768",
        "P02787", "P02788", "P05413", "P08263", "P10145", "P10636-8", "P61626", "P62988",
        "P63165", "P68871", "P69905"
    };

    public static string[] SharedAccessions =
    {
        "O76070", "P00441", "P00918", "P04040", "P06396", "P06732", "P08758", "P09211",
        "P10599", "P12081", "P15559", "P16083", "P41159", "P51965", "P55957", "P62937",
        "P63279", "P99999", "Q06830", "Q15843"
    };

    #endregion

    // ──────────────────────────────────────────────────────────────────────────
    // Real SpikeIn-5mix-MS3 integration tests
    // ──────────────────────────────────────────────────────────────────────────

    [Test]
    public void LoadAndRunSpikeInData_BasicPipeline()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(upsPsmPath), Is.True, $"UPS psmtsv not found at: {upsPsmPath}");

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

        Assert.That(proteinMatrix, Is.Not.Null);
        Assert.That(proteinMatrix.RowKeys.Count, Is.GreaterThan(0), "Protein matrix has no rows");
        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(140),
            $"Expected 140 columns (14 files × 10 channels), got {proteinMatrix.ColumnKeys.Count}");

        bool anyNonZero = proteinMatrix.RowKeys
            .Any(p => proteinMatrix.GetRow(p).Any(v => v > 0));
        Assert.That(anyNonZero, Is.True, "All protein matrix values are zero");
    }

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

        var foldChanges_1_vs_0125 = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.125")
            .Select(x => x.foldChange).ToList();
        var foldChanges_1_vs_0500 = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.5")
            .Select(x => x.foldChange).ToList();

        Assert.That(foldChanges_1_vs_0125.Count, Is.GreaterThan(0),
            "No proteins had quantifiable intensities in both '1' and '0.125' conditions");

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

        var foldChanges = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.5")
            .Select(x => x.foldChange).ToList();

        Assert.That(foldChanges.Count, Is.GreaterThan(0),
            "No HeLa proteins had quantifiable intensities in both '1' and '0.5' conditions");

        double medianFc = QuantificationEvaluator.Median(foldChanges);
        Assert.That(medianFc, Is.InRange(0.3, 3.0),
            $"Median HeLa fold change '1' vs '0.5' was {medianFc:F3}, expected in [0.3, 3.0]");
    }

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
                .Select(x => x.foldChange).ToList();

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
            .Select(x => x.foldChange).ToList();

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

    [Test]
    public void ProteinOverlap_HelaVsUPS()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
        string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        string helaPsmPath = Path.Combine(dataDir, "HeLa_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
        Assert.That(File.Exists(upsPsmPath), Is.True, $"UPS psmtsv not found at: {upsPsmPath}");
        Assert.That(File.Exists(helaPsmPath), Is.True, $"HeLa psmtsv not found at: {helaPsmPath}");

        var upsRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(upsPsmPath, out _)
            .Where(r => r.QValue <= 0.01 && r.DecoyContamTarget.Contains('T')
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

        var helaRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(helaPsmPath, out _)
            .Where(r => r.QValue <= 0.01 && r.DecoyContamTarget.Contains('T')
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

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

        TestContext.WriteLine("=== Protein Overlap: HeLa vs UPS ===");
        TestContext.WriteLine($"UPS total accessions:  {upsAccessions.Count}");
        TestContext.WriteLine($"HeLa total accessions: {helaAccessions.Count}");
        TestContext.WriteLine($"Shared:                {shared.Count}");
        TestContext.WriteLine($"UPS-only:              {upsOnly.Count}");
        TestContext.WriteLine($"HeLa-only:             {helaOnly.Count}");

        TestContext.WriteLine($"\n--- Shared proteins ({shared.Count}) ---");
        foreach (var acc in shared.OrderBy(a => a))
            TestContext.WriteLine($"  {acc}");

        TestContext.WriteLine($"\n--- UPS-only proteins ({upsOnly.Count}) ---");
        foreach (var acc in upsOnly.OrderBy(a => a))
            TestContext.WriteLine($"  {acc}");

        TestContext.WriteLine($"\n--- HeLa-only proteins ({helaOnly.Count}) ---");
        foreach (var acc in helaOnly.OrderBy(a => a).Take(50))
            TestContext.WriteLine($"  {acc}");
        if (helaOnly.Count > 50)
            TestContext.WriteLine($"  ... and {helaOnly.Count - 50} more");

        Assert.That(upsAccessions.Count, Is.GreaterThan(0), "No UPS accessions found");
        Assert.That(helaAccessions.Count, Is.GreaterThan(0), "No HeLa accessions found");
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

        var (spectralMatches, peptides, proteinGroups) =
            PsmTsvQuantAdapter.BuildQuantificationInputsFromMultipleDirectories(
                new[] { upsDir, helaDir });


        sw.Stop();
        TestContext.WriteLine("All data loaded in :" + sw.Elapsed);

        TestContext.WriteLine($"Combined: {spectralMatches.Count} PSMs, {peptides.Count} peptides, {proteinGroups.Count} protein groups");

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

        var upsOnlySet = new HashSet<string>(UPSOnlyAccessions);
        var sharedSet = new HashSet<string>(SharedAccessions);

        var allAccessions = proteinMatrix.RowKeys
            .Select(pg => pg.BioPolymerGroupName)
            .ToHashSet();

        var helaOnlyAccessions = allAccessions
            .Where(a => !upsOnlySet.Contains(a) && !sharedSet.Contains(a))
            .ToHashSet();

        TestContext.WriteLine($"\nProtein classification in combined matrix:");
        TestContext.WriteLine($"  Total proteins:    {allAccessions.Count}");
        TestContext.WriteLine($"  UPS-only found:    {allAccessions.Count(a => upsOnlySet.Contains(a))}");
        TestContext.WriteLine($"  Shared found:      {allAccessions.Count(a => sharedSet.Contains(a))}");
        TestContext.WriteLine($"  HeLa-only found:   {helaOnlyAccessions.Count}");

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

        var fc8x_upsOnly = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.125")
            .Where(x => upsOnlySet.Contains(x.accession))
            .Select(x => x.foldChange).ToList();

        if (fc8x_upsOnly.Count > 0)
        {
            double medianFc = QuantificationEvaluator.Median(fc8x_upsOnly);
            Assert.That(medianFc, Is.GreaterThan(1.5),
                $"UPS-only median FC '1' vs '0.125' was {medianFc:F2}, expected > 1.5 (true value = 8.0)");
        }

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

    [Test]
    public void StrategyEvaluation_AllCombinations_UPS()
    {
        string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");

        var (spectralMatches, peptides, proteinGroups) = PsmTsvQuantAdapter.BuildQuantificationInputsFromDirectory(
            Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info", "UPS_TMT3_Search", "Task1-SearchTask"));

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
