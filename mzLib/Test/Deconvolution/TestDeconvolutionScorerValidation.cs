using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test
{
    /// <summary>
    /// Validation harness for <see cref="DeconvolutionScorer"/> weight calibration
    /// using deconvolution-level target/decoy labelling.
    ///
    /// Labels are assigned at the deconvolution level:
    ///   TARGET — envelope found with the real isotope spacing (C13-C12 = 1.003355 Da)
    ///   DECOY  — envelope found on the same spectrum with the shifted spacing (0.9444 Da)
    ///
    /// Three test cases:
    ///   BottomUp  — A549 cell line, Velos Orbitrap, Classic z=1-12
    ///               Known result: AUC = 0.50 (Classic pre-filters by Averagine)
    ///   TopDown   — Yeast top-down, Classic z=1-60
    ///               Known result: AUC = 0.50 (same reason)
    ///   IsoDec    — A549 bottom-up, IsoDec
    ///               Expected: AUC > 0.50 (DecoyIsotopeDistance reaches the DLL)
    ///               Requires isodeclib.dll at runtime — [Category("RequiresIsoDec")]
    ///
    /// Output TSV columns:
    ///   Label, ScanNumber, Charge, MonoMass, TotalIntensity,
    ///   Cosine, PpmError, Completeness, RatioConsistency, GenericScore
    ///
    /// Running all:
    ///   dotnet test --filter "Category=DeconvolutionValidation"
    /// Running IsoDec only (requires DLL):
    ///   dotnet test --filter "IsoDec_ExtractFeaturesUsingDeconvolutionDecoys"
    /// </summary>
    [TestFixture]
    [Category("DeconvolutionValidation")]
    public sealed class TestDeconvolutionScorerValidation
    {
        // ── Paths ─────────────────────────────────────────────────────────────

        private const string BottomUpMzmlPath =
            @"E:\Projects\Mann_11cell_lines\A549\A549_1\20100604_Velos1_TaGe_SA_A549_3.mzML";

        private const string BottomUpOutputTsv =
            @"E:\Projects\Mann_11cell_lines\A549\A549_1\MetaMorpheus_Output\DeconvolutionScorerFeatures_BottomUp.tsv";

        private const string TopDownMzmlPath =
            @"E:\Projects\LVS_TD_Yeast\05-26-17_B7A_yeast_td_fract8_rep2.mzML";

        private const string TopDownOutputTsv =
            @"E:\Projects\LVS_TD_Yeast\DeconvolutionScorerFeatures_TopDown.tsv";

        private const string IsoDecOutputTsv =
            @"E:\Projects\LVS_TD_Yeast\DeconvolutionScorerFeatures_IsoDec_TopDown.tsv";

        // ── Parameters ────────────────────────────────────────────────────────

        private static ClassicDeconvolutionParameters BottomUpParams =>
            new ClassicDeconvolutionParameters(
                minCharge: 1,
                maxCharge: 12,
                deconPpm: 4.0,
                intensityRatio: 3.0,
                polarity: Polarity.Positive);

        private static ClassicDeconvolutionParameters TopDownParams =>
            new ClassicDeconvolutionParameters(
                minCharge: 1,
                maxCharge: 60,
                deconPpm: 10.0,
                intensityRatio: 3.0,
                polarity: Polarity.Positive);

        private static IsoDecDeconvolutionParameters IsoDecParams =>
            new IsoDecDeconvolutionParameters(polarity: Polarity.Positive);

        private static readonly AverageResidue Model = new Averagine();

        // Process at most this many MS1 scans to keep runtime reasonable.
        // Set to int.MaxValue to process the full file.
        private const int MaxScansToProcess = 2000;

        // ══════════════════════════════════════════════════════════════════════
        // Bottom-up Classic z=1-12
        // Known result: AUC = 0.50
        // Classic's intensity ratio filter pre-enforces Averagine shape on all
        // recruited envelopes, so the spectrum-shift decoy produces features
        // that are indistinguishable from targets. No score ordering assertion.
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void BottomUp_ExtractFeaturesUsingDeconvolutionDecoys()
        {
            Assert.That(File.Exists(BottomUpMzmlPath), Is.True,
                $"mzML not found: {BottomUpMzmlPath}");

            var (tRows, dRows) = RunExtraction(
                BottomUpMzmlPath, BottomUpParams, BottomUpOutputTsv);

            Console.WriteLine("\n--- Bottom-up Classic: feature means target vs decoy ---");
            PrintFeatureMeans(tRows, dRows);

            Assert.That(tRows.Count, Is.GreaterThan(0));
            Assert.That(dRows.Count, Is.GreaterThan(0));
            foreach (var r in tRows.Concat(dRows))
                Assert.That(r.GenericScore, Is.InRange(0.0, 1.0));

            // No score ordering assertion — AUC = 0.50 is the known result.
            // See class summary for explanation.
        }

        // ══════════════════════════════════════════════════════════════════════
        // Top-down Classic z=1-60
        // Known result: AUC = 0.50 for the same structural reason as bottom-up.
        // Classic's Peak2satisfiesRatio filter is independent of charge range.
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void TopDown_ExtractFeaturesUsingDeconvolutionDecoys()
        {
            Assert.That(File.Exists(TopDownMzmlPath), Is.True,
                $"mzML not found: {TopDownMzmlPath}");

            var (tRows, dRows) = RunExtraction(
                TopDownMzmlPath, TopDownParams, TopDownOutputTsv);

            Console.WriteLine("\n--- Top-down Classic: feature means target vs decoy ---");
            PrintFeatureMeans(tRows, dRows);

            Assert.That(tRows.Count, Is.GreaterThan(0));
            Assert.That(dRows.Count, Is.GreaterThan(0));
            foreach (var r in tRows.Concat(dRows))
                Assert.That(r.GenericScore, Is.InRange(0.0, 1.0));

            double meanT = tRows.Average(r => r.GenericScore);
            double meanD = dRows.Average(r => r.GenericScore);
            Console.WriteLine($"\nMean score — targets: {meanT:F4}  decoys: {meanD:F4}  " +
                              $"delta: {meanT - meanD:+0.0000;-0.0000}");
        }

        // ══════════════════════════════════════════════════════════════════════
        // IsoDec bottom-up
        //
        // IsoDec passes DecoyIsotopeDistance directly to the C DLL as mass_diff_c.
        // The decoy pass therefore searches for isotope series at 0.9444 Da spacing,
        // producing candidate masses at genuinely different positions. The post-hoc
        // scorer evaluates those masses against Averagine at the correct spacing,
        // so decoy envelopes should score worse than targets.
        //
        // This is the first algorithm where a non-trivial AUC is expected.
        //
        // Requires isodeclib.dll at runtime.
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        [Category("RequiresIsoDec")]
        public void IsoDec_ExtractFeaturesUsingDeconvolutionDecoys()
        {
            Assert.That(File.Exists(TopDownMzmlPath), Is.True,
                $"mzML not found: {TopDownMzmlPath}");

            var (tRows, dRows) = RunExtraction(
                TopDownMzmlPath, IsoDecParams, IsoDecOutputTsv);
            Console.WriteLine("\n--- IsoDec top-down: feature means target vs decoy ---");
            PrintFeatureMeans(tRows, dRows);

            Assert.That(tRows.Count, Is.GreaterThan(0));
            Assert.That(dRows.Count, Is.GreaterThan(0),
                "No decoy envelopes found — verify isodeclib.dll is present and " +
                "that DecoyIsotopeDistance is reaching the DLL as mass_diff_c. " +
                "Run TestShallowClone.IsoDec_ShallowClone_NullsIsoSettingsCache " +
                "to confirm the cache invalidation is working.");

            foreach (var r in tRows.Concat(dRows))
                Assert.That(r.GenericScore, Is.InRange(0.0, 1.0));

            double meanT = tRows.Average(r => r.GenericScore);
            double meanD = dRows.Average(r => r.GenericScore);
            Console.WriteLine($"\nMean score — targets: {meanT:F4}  decoys: {meanD:F4}  " +
                              $"delta: {meanT - meanD:+0.0000;-0.0000}");

            // Hard assertion — unlike Classic, IsoDec should discriminate.
            Assert.That(meanT, Is.GreaterThan(meanD),
                "IsoDec target envelopes should score higher than decoys. " +
                "The decoy pass uses mass_diff_c = 0.9444 in the DLL, producing " +
                "candidate masses at genuinely different positions.");
        }

        // ══════════════════════════════════════════════════════════════════════
        // Shared helpers
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Runs DeconvoluteWithDecoys on all MS1 scans in the mzML, computes
        /// features for every target and decoy envelope, writes a TSV, and
        /// returns the two row lists for assertion.
        /// </summary>
        private static (List<FeatureRow> Targets, List<FeatureRow> Decoys) RunExtraction(
            string mzmlPath,
            DeconvolutionParameters deconParams,   // base type — accepts Classic or IsoDec
            string outputTsvPath)
        {
            Console.WriteLine($"Loading: {Path.GetFileName(mzmlPath)}");
            MsDataFile mzmlFile = MsDataFileReader.GetDataFile(mzmlPath).LoadAllStaticData();

            var ms1Scans = mzmlFile.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .Take(MaxScansToProcess)
                .ToList();
            Console.WriteLine($"  Processing {ms1Scans.Count} MS1 scans");

            var rows = new List<FeatureRow>();
            int scansWithTargets = 0;
            int scansWithDecoys = 0;

            foreach (var scan in ms1Scans)
            {
                var (targets, decoys) = Deconvoluter.DeconvoluteWithDecoys(scan, deconParams);

                if (targets.Count > 0) scansWithTargets++;
                if (decoys.Count > 0) scansWithDecoys++;

                foreach (var env in targets)
                {
                    var f = DeconvolutionScorer.ComputeFeatures(env, Model);
                    rows.Add(MakeRow("T", scan.OneBasedScanNumber, env, f));
                }
                foreach (var env in decoys)
                {
                    var f = DeconvolutionScorer.ComputeFeatures(env, Model);
                    rows.Add(MakeRow("D", scan.OneBasedScanNumber, env, f));
                }
            }

            var tRows = rows.Where(r => r.Label == "T").ToList();
            var dRows = rows.Where(r => r.Label == "D").ToList();

            Console.WriteLine($"  Scans with targets: {scansWithTargets}  " +
                              $"with decoys: {scansWithDecoys}");
            Console.WriteLine($"  Target rows: {tRows.Count}  Decoy rows: {dRows.Count}");

            Directory.CreateDirectory(Path.GetDirectoryName(outputTsvPath)!);
            using var writer = new StreamWriter(outputTsvPath);
            writer.WriteLine(
                "Label\tScanNumber\tCharge\tMonoMass\tTotalIntensity\t" +
                "Cosine\tPpmError\tCompleteness\tRatioConsistency\tGenericScore");

            foreach (var r in rows)
                writer.WriteLine(string.Join("\t",
                    r.Label,
                    r.ScanNumber,
                    r.Charge,
                    r.MonoMass.ToString("F4"),
                    r.TotalIntensity.ToString("F0"),
                    r.Cosine.ToString("F6"),
                    r.PpmError.ToString("F4"),
                    r.Completeness.ToString("F6"),
                    r.RatioConsistency.ToString("F6"),
                    r.GenericScore.ToString("F6")));

            Console.WriteLine($"  Written: {outputTsvPath}");
            return (tRows, dRows);
        }

        private static void PrintFeatureMeans(List<FeatureRow> tRows, List<FeatureRow> dRows)
        {
            if (tRows.Count == 0 || dRows.Count == 0) return;
            Console.WriteLine($"  Cosine:           {tRows.Average(r => r.Cosine):F4}  vs  {dRows.Average(r => r.Cosine):F4}");
            Console.WriteLine($"  PpmError:         {tRows.Average(r => r.PpmError):F4}  vs  {dRows.Average(r => r.PpmError):F4}");
            Console.WriteLine($"  Completeness:     {tRows.Average(r => r.Completeness):F4}  vs  {dRows.Average(r => r.Completeness):F4}");
            Console.WriteLine($"  RatioConsistency: {tRows.Average(r => r.RatioConsistency):F4}  vs  {dRows.Average(r => r.RatioConsistency):F4}");
            Console.WriteLine($"  GenericScore:     {tRows.Average(r => r.GenericScore):F4}  vs  {dRows.Average(r => r.GenericScore):F4}");
        }

        private static FeatureRow MakeRow(
            string label, int scanNumber,
            IsotopicEnvelope env, EnvelopeScoreFeatures features)
            => new FeatureRow
            {
                Label = label,
                ScanNumber = scanNumber,
                Charge = Math.Abs(env.Charge),
                MonoMass = env.MonoisotopicMass,
                TotalIntensity = env.TotalIntensity,
                Cosine = features.AveragineCosineSimilarity,
                PpmError = features.AvgPpmError,
                Completeness = features.PeakCompleteness,
                RatioConsistency = features.IntensityRatioConsistency,
                GenericScore = DeconvolutionScorer.ComputeScore(features)
            };

        private sealed class FeatureRow
        {
            public string Label { get; set; } = "";
            public int ScanNumber { get; set; }
            public int Charge { get; set; }
            public double MonoMass { get; set; }
            public double TotalIntensity { get; set; }
            public double Cosine { get; set; }
            public double PpmError { get; set; }
            public double Completeness { get; set; }
            public double RatioConsistency { get; set; }
            public double GenericScore { get; set; }
        }
    }
}