using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.InternalIons
{
    /// <summary>
    /// Evaluates the mixed model (Prosit primary + ONNX internal) predictions against
    /// experimental spectra collected at 5 different normalized collision energies (NCEs).
    ///
    /// EXPERIMENTAL DESIGN
    /// -------------------
    /// 5 raw files from 500ng HeLa trypsin digest, each at a different NCE:
    ///   22, 27, 32, 37, 42
    /// Each was searched separately, producing per-NCE InternalFragmentIons.tsv files.
    ///
    /// GOALS
    /// -----
    /// 1. Compare predicted internal fragment intensities to observed experimental intensities
    ///    at each NCE — how well does the mixed model match reality?
    /// 2. Test the assumption that the internal fragment model is NCE-agnostic:
    ///    are prediction errors stable across energies, or does performance degrade
    ///    at specific NCEs?
    ///
    /// WHAT THIS TEST READS
    /// --------------------
    /// Each NCE folder should contain an InternalFragmentIons.tsv from the
    /// InternalFragmentAnalysisRunnerTests pipeline. The TSV has columns including:
    ///   PeptideSequence, ScanNumber, TicNormalizedIntensity, FragmentLength,
    ///   PeptideLength, PrecursorCharge, PassesMassAccuracyFilter,
    ///   HasMatchedBIonAtNTerm, HasMatchedYIonAtCTerm, IsProlineAtInternalNTerminus,
    ///   HasProlineAtEitherTerminus, HasAspartateAtEitherTerminus,
    ///   RelativeDistanceFromCTerm, BasicResiduesInYIonSpan, BasicResiduesInBIonSpan,
    ///   NumberOfBasicResidues, MaxTerminalIonIntensity, YIonIntensityAtCTerm,
    ///   BIonIntensityAtNTerm, NTerminalFlankingHydrophobicity, BYProductScore,
    ///   LocalIntensityRank, CollisionEnergy, InternalSequence, FragStart, FragEnd, etc.
    /// </summary>
    [TestFixture]
    public class MixedModelEvaluationTests
    {
        // ════════════════════════════════════════════════════════════════════════
        // CONFIGURE THESE PATHS for your machine
        // ════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Base directory containing per-NCE subfolders.
        /// Expected structure:
        ///   BaseDir/
        ///     nce22/Task1-SearchTask/InternalFragmentIons.tsv
        ///     nce27/Task1-SearchTask/InternalFragmentIons.tsv
        ///     nce32/Task1-SearchTask/InternalFragmentIons.tsv
        ///     nce37/Task1-SearchTask/InternalFragmentIons.tsv
        ///     nce42/Task1-SearchTask/InternalFragmentIons.tsv
        /// Adjust the TsvPathPattern below if your naming differs.
        /// </summary>
        private const string BaseDir = @"F:\MSV000090552_scribe\essentialOutput";

        /// <summary>
        /// Pattern for locating each NCE's TSV file. {0} is replaced with the NCE value.
        /// Adjust to match your folder naming convention.
        /// </summary>
        private static string TsvPathForNce(int nce)
            => Path.Combine(BaseDir, $"InternalFragmentIons_nce{nce}.tsv");

        /// <summary>Output report path.</summary>
        private static readonly string OutputPath =
            Path.Combine(BaseDir, "MixedModelEvaluation_CrossNCE_Report.txt");

        private static readonly int[] AllNCEs = { 22, 27, 32, 37, 42 };

        // ════════════════════════════════════════════════════════════════════════
        // Data structures
        // ════════════════════════════════════════════════════════════════════════

        private class NceDataset
        {
            public int NCE { get; init; }
            public string TsvPath { get; init; }
            public List<Dictionary<string, string>> AllRows { get; init; } = new();
            public List<Dictionary<string, string>> Passing { get; init; } = new();
        }

        private Dictionary<int, NceDataset> _datasets;
        private StreamWriter _out;

        // ════════════════════════════════════════════════════════════════════════
        // Setup / Teardown
        // ════════════════════════════════════════════════════════════════════════

        [OneTimeSetUp]
        public void Setup()
        {
            _datasets = new Dictionary<int, NceDataset>();
            _out = new StreamWriter(OutputPath, false) { AutoFlush = false };

            W("╔══════════════════════════════════════════════════════════════════╗");
            W("║    MIXED MODEL EVALUATION: PREDICTED vs EXPERIMENTAL SPECTRA   ║");
            W("║    Cross-NCE Analysis (22, 27, 32, 37, 42)                     ║");
            W("╚══════════════════════════════════════════════════════════════════╝");
            W($"Timestamp: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
            W();

            int loaded = 0;
            foreach (var nce in AllNCEs)
            {
                var path = TsvPathForNce(nce);
                if (!File.Exists(path))
                {
                    W($"WARNING: NCE {nce} TSV not found at: {path}");
                    continue;
                }

                var lines = File.ReadAllLines(path);
                if (lines.Length < 2) { W($"WARNING: NCE {nce} TSV has no data rows"); continue; }

                var headers = lines[0].Split('\t');
                var allRows = new List<Dictionary<string, string>>();
                for (int i = 1; i < lines.Length; i++)
                {
                    var vals = lines[i].Split('\t');
                    var row = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
                    for (int j = 0; j < headers.Length && j < vals.Length; j++)
                        row[headers[j]] = vals[j];
                    allRows.Add(row);
                }

                var passing = allRows.Where(r => B(r, "PassesMassAccuracyFilter")).ToList();

                _datasets[nce] = new NceDataset
                {
                    NCE = nce,
                    TsvPath = path,
                    AllRows = allRows,
                    Passing = passing
                };

                W($"NCE {nce}: {allRows.Count:N0} total rows, {passing.Count:N0} passing filter");
                loaded++;
            }

            W();
            if (loaded < 2)
                Assert.Ignore($"Need at least 2 NCE datasets for cross-NCE analysis. Found {loaded}.");
        }

        [OneTimeTearDown]
        public void Teardown()
        {
            _out?.Flush();
            _out?.Close();
            _out?.Dispose();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 1: Per-NCE Dataset Summary
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(1)]
        public void Test01_PerNCE_DatasetSummary()
        {
            W(Sep("TEST 1: Per-NCE Dataset Summary"));

            W($"{"NCE",4} | {"Total",8} | {"Passing",8} | {"PassRate",8} | {"Peptides",8} | {"Scans",8} | {"MeanTicNI",12} | {"MedTicNI",12} | {"StdTicNI",12}");
            W(new string('-', 100));

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var ds = _datasets[nce];
                var ticNI = ds.Passing.Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                var peptides = ds.Passing.Select(r => S(r, "PeptideSequence")).Distinct().Count();
                var scans = ds.Passing.Select(r => S(r, "ScanNumber")).Distinct().Count();
                double passRate = ds.AllRows.Count > 0 ? 100.0 * ds.Passing.Count / ds.AllRows.Count : 0;

                W($"{nce,4} | {ds.AllRows.Count,8:N0} | {ds.Passing.Count,8:N0} | {passRate,7:F1}% | {peptides,8:N0} | {scans,8:N0} | {(ticNI.Count > 0 ? ticNI.Average() : 0),12:E4} | {Median(ticNI),12:E4} | {StdDev(ticNI),12:E4}");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 2: TicNI Distribution Across NCEs
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(2)]
        public void Test02_TicNI_DistributionAcrossNCE()
        {
            W(Sep("TEST 2: TicNI Distribution by NCE (Percentiles)"));
            W("  If internal fragment model is NCE-agnostic, distributions should be similar.");
            W();

            W($"{"NCE",4} | {"N",8} | {"P5",10} | {"P10",10} | {"P25",10} | {"P50",10} | {"P75",10} | {"P90",10} | {"P95",10} | {"Mean",10}");
            W(new string('-', 110));

            var perNceMeans = new List<(int nce, double mean)>();

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var ticNI = _datasets[nce].Passing
                    .Select(r => D(r, "TicNormalizedIntensity"))
                    .Where(v => !double.IsNaN(v) && v > 0)
                    .OrderBy(v => v).ToList();

                if (ticNI.Count == 0) continue;

                perNceMeans.Add((nce, ticNI.Average()));

                W($"{nce,4} | {ticNI.Count,8:N0} | {Percentile(ticNI, 5),10:E3} | {Percentile(ticNI, 10),10:E3} | {Percentile(ticNI, 25),10:E3} | {Percentile(ticNI, 50),10:E3} | {Percentile(ticNI, 75),10:E3} | {Percentile(ticNI, 90),10:E3} | {Percentile(ticNI, 95),10:E3} | {ticNI.Average(),10:E3}");
            }

            // Cross-NCE mean TicNI variation
            if (perNceMeans.Count >= 2)
            {
                var means = perNceMeans.Select(p => p.mean).ToList();
                double cv = StdDev(means) / means.Average();
                W();
                W($"Cross-NCE mean TicNI coefficient of variation: {cv:F4}");
                W($"  (CV < 0.15 suggests intensity scale is reasonably NCE-stable)");
                W($"  Min mean: NCE {perNceMeans.MinBy(p => p.mean).nce} = {means.Min():E4}");
                W($"  Max mean: NCE {perNceMeans.MaxBy(p => p.mean).nce} = {means.Max():E4}");
                W($"  Fold-change (max/min): {means.Max() / means.Min():F2}");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 3: Feature-TicNI Correlations Across NCEs
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(3)]
        public void Test03_FeatureCorrelations_AcrossNCE()
        {
            W(Sep("TEST 3: Feature-TicNI Pearson r Stability Across NCEs"));
            W("  If model is NCE-agnostic, the same features should correlate similarly at every NCE.");
            W();

            var featureGetters = new (string name, Func<Dictionary<string, string>, double> get)[]
            {
                ("FragmentLength", r => I(r, "FragmentLength")),
                ("PeptideLength", r => I(r, "PeptideLength")),
                ("PrecursorCharge", r => I(r, "PrecursorCharge")),
                ("RelativeDistFromCTerm", r => D(r, "RelativeDistanceFromCTerm")),
                ("NumberOfBasicResidues", r => I(r, "NumberOfBasicResidues")),
                ("BasicResInBIonSpan", r => I(r, "BasicResiduesInBIonSpan")),
                ("BasicResInYIonSpan", r => I(r, "BasicResiduesInYIonSpan")),
                ("BIonIntensityAtNTerm", r => D(r, "BIonIntensityAtNTerm")),
                ("YIonIntensityAtCTerm", r => D(r, "YIonIntensityAtCTerm")),
                ("MaxTerminalIonIntensity", r => D(r, "MaxTerminalIonIntensity")),
                ("BYProductScore", r => D(r, "BYProductScore")),
                ("HasProAtEitherTerm", r => B(r, "HasProlineAtEitherTerminus") ? 1.0 : 0.0),
                ("HasAspAtEitherTerm", r => B(r, "HasAspartateAtEitherTerminus") ? 1.0 : 0.0),
                ("IsProAtInternalNTerm", r => B(r, "IsProlineAtInternalNTerminus") ? 1.0 : 0.0),
                ("NTermFlankHydro", r => D(r, "NTerminalFlankingHydrophobicity")),
                ("LocalIntensityRank", r => D(r, "LocalIntensityRank")),
                ("HasMatchedBIonAtNTerm", r => B(r, "HasMatchedBIonAtNTerm") ? 1.0 : 0.0),
                ("HasMatchedYIonAtCTerm", r => B(r, "HasMatchedYIonAtCTerm") ? 1.0 : 0.0),
            };

            // Header
            var nceList = AllNCEs.Where(n => _datasets.ContainsKey(n)).ToList();
            var header = $"{"Feature",-24}";
            foreach (var nce in nceList) header += $" | {"NCE" + nce,8}";
            header += $" | {"Range",7} | {"CV",7} | {"Stable?",7}";
            W(header);
            W(new string('-', header.Length));

            var stabilityResults = new List<(string feat, double cv, bool stable)>();

            foreach (var (name, getter) in featureGetters)
            {
                var row = $"{name,-24}";
                var rValues = new List<double>();

                foreach (var nce in nceList)
                {
                    var ds = _datasets[nce];
                    var target = ds.Passing.Select(r => D(r, "TicNormalizedIntensity")).ToList();
                    var featureVals = ds.Passing.Select(getter).ToList();
                    double r = Pearson(featureVals, target);
                    row += $" | {(double.IsNaN(r) ? "NaN" : r.ToString("F4")),8}";
                    if (!double.IsNaN(r)) rValues.Add(r);
                }

                if (rValues.Count >= 2)
                {
                    double range = rValues.Max() - rValues.Min();
                    // For correlations that can cross zero, use range rather than CV
                    double absAvg = rValues.Select(Math.Abs).Average();
                    double cv = absAvg > 0.01 ? StdDev(rValues) / absAvg : double.NaN;
                    bool stable = range < 0.10; // r varies by less than 0.10 across NCEs
                    row += $" | {range,7:F4} | {(double.IsNaN(cv) ? "N/A" : cv.ToString("F3")),7} | {(stable ? "YES" : "NO"),7}";
                    stabilityResults.Add((name, double.IsNaN(cv) ? 999 : cv, stable));
                }
                else
                {
                    row += $" | {"N/A",7} | {"N/A",7} | {"N/A",7}";
                }

                W(row);
            }

            W();
            int stableCount = stabilityResults.Count(s => s.stable);
            W($"Features with stable correlations (range < 0.10): {stableCount}/{stabilityResults.Count}");
            var unstable = stabilityResults.Where(s => !s.stable).ToList();
            if (unstable.Any())
            {
                W("Unstable features (may indicate NCE sensitivity):");
                foreach (var u in unstable.OrderByDescending(u => u.cv))
                    W($"  - {u.feat} (CV={u.cv:F3})");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 4: Cross-NCE Peptide Overlap Analysis
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(4)]
        public void Test04_PeptideOverlap_AcrossNCE()
        {
            W(Sep("TEST 4: Peptide Overlap Across NCEs"));
            W("  Peptides observed at multiple NCEs allow paired comparisons.");
            W();

            // Build peptide → {NCE → list of rows} lookup
            var pepNceMap = new Dictionary<string, Dictionary<int, List<Dictionary<string, string>>>>();
            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                foreach (var row in _datasets[nce].Passing)
                {
                    var seq = S(row, "PeptideSequence");
                    if (string.IsNullOrEmpty(seq)) continue;
                    if (!pepNceMap.ContainsKey(seq)) pepNceMap[seq] = new();
                    if (!pepNceMap[seq].ContainsKey(nce)) pepNceMap[seq][nce] = new();
                    pepNceMap[seq][nce].Add(row);
                }
            }

            var nceList = AllNCEs.Where(n => _datasets.ContainsKey(n)).ToList();

            // Count peptides observed at N different NCEs
            W("Peptides observed at N different NCEs:");
            for (int k = 1; k <= nceList.Count; k++)
            {
                int count = pepNceMap.Count(p => p.Value.Keys.Count >= k);
                W($"  >= {k} NCEs: {count:N0} peptides");
            }

            int coreCount = pepNceMap.Count(p => p.Value.Keys.Count == nceList.Count);
            W($"\nCore peptides (all {nceList.Count} NCEs): {coreCount:N0}");
            W();

            // Pairwise NCE overlap
            W("Pairwise peptide overlap (Jaccard index):");
            var hdr = $"{"",6}";
            foreach (var nce in nceList) hdr += $" | {nce,5}";
            W(hdr);
            W(new string('-', hdr.Length));

            foreach (var nce1 in nceList)
            {
                var row = $"{nce1,6}";
                var peps1 = _datasets[nce1].Passing.Select(r => S(r, "PeptideSequence")).ToHashSet();
                foreach (var nce2 in nceList)
                {
                    if (nce2 < nce1) { row += $" | {"",5}"; continue; }
                    var peps2 = _datasets[nce2].Passing.Select(r => S(r, "PeptideSequence")).ToHashSet();
                    int intersect = peps1.Intersect(peps2).Count();
                    int union = peps1.Union(peps2).Count();
                    double jaccard = union > 0 ? (double)intersect / union : 0;
                    row += $" | {jaccard,5:F2}";
                }
                W(row);
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 5: Paired Peptide TicNI Comparison Across NCEs
        //   For peptides seen at multiple NCEs, compare their mean TicNI
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(5)]
        public void Test05_PairedPeptide_TicNI_AcrossNCE()
        {
            W(Sep("TEST 5: Paired Peptide Mean Internal TicNI Across NCEs"));
            W("  For each peptide observed at 2+ NCEs, compute its mean internal fragment TicNI.");
            W("  Then correlate those peptide-level means across NCE pairs.");
            W("  High correlation = internal fragment pattern is NCE-stable at the peptide level.");
            W();

            var nceList = AllNCEs.Where(n => _datasets.ContainsKey(n)).ToList();

            // Peptide → NCE → mean TicNI
            var pepMeans = new Dictionary<string, Dictionary<int, double>>();
            foreach (var nce in nceList)
            {
                var groups = _datasets[nce].Passing
                    .GroupBy(r => S(r, "PeptideSequence"))
                    .Where(g => !string.IsNullOrEmpty(g.Key));

                foreach (var g in groups)
                {
                    var vals = g.Select(r => D(r, "TicNormalizedIntensity"))
                        .Where(v => !double.IsNaN(v) && v > 0).ToList();
                    if (vals.Count == 0) continue;

                    if (!pepMeans.ContainsKey(g.Key)) pepMeans[g.Key] = new();
                    pepMeans[g.Key][nce] = vals.Average();
                }
            }

            W("Pairwise Pearson r of peptide-level mean TicNI:");
            var hdr = $"{"",6}";
            foreach (var nce in nceList) hdr += $" | {nce,7}";
            W(hdr);
            W(new string('-', hdr.Length));

            foreach (var nce1 in nceList)
            {
                var row = $"{nce1,6}";
                foreach (var nce2 in nceList)
                {
                    if (nce2 < nce1) { row += $" | {"",7}"; continue; }
                    if (nce1 == nce2) { row += $" | {"1.0000",7}"; continue; }

                    // Get peptides present at both NCEs
                    var shared = pepMeans
                        .Where(kv => kv.Value.ContainsKey(nce1) && kv.Value.ContainsKey(nce2))
                        .ToList();

                    if (shared.Count < 10)
                    {
                        row += $" | {"N/A",7}";
                        continue;
                    }

                    var xs = shared.Select(kv => kv.Value[nce1]).ToList();
                    var ys = shared.Select(kv => kv.Value[nce2]).ToList();
                    double r = Pearson(xs, ys);
                    row += $" | {(double.IsNaN(r) ? "NaN" : r.ToString("F4")),7}";
                }
                W(row);
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 6: Internal Fragment Ion Count vs NCE
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(6)]
        public void Test06_InternalIonCount_VsNCE()
        {
            W(Sep("TEST 6: Internal Fragment Ion Count per Scan vs NCE"));
            W("  More fragmentation energy → more internal fragments (expected).");
            W("  If count increases monotonically with NCE, the physics is consistent.");
            W();

            W($"{"NCE",4} | {"MeanCount",10} | {"MedCount",10} | {"StdDev",10} | {"TotalIons",10}");
            W(new string('-', 60));

            var means = new List<(int nce, double mean)>();

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var scanCounts = _datasets[nce].Passing
                    .GroupBy(r => S(r, "ScanNumber"))
                    .Select(g => (double)g.Count())
                    .ToList();

                if (scanCounts.Count == 0) continue;

                double mean = scanCounts.Average();
                means.Add((nce, mean));

                W($"{nce,4} | {mean,10:F2} | {Median(scanCounts),10:F2} | {StdDev(scanCounts),10:F2} | {_datasets[nce].Passing.Count,10:N0}");
            }

            W();
            if (means.Count >= 3)
            {
                // Test monotonicity
                bool monotonic = true;
                for (int i = 1; i < means.Count; i++)
                    if (means[i].mean < means[i - 1].mean) monotonic = false;

                // Spearman rank correlation between NCE and mean count
                double spearman = SpearmanRank(
                    means.Select(m => (double)m.nce).ToList(),
                    means.Select(m => m.mean).ToList());

                W($"Monotonically increasing count with NCE: {(monotonic ? "YES" : "NO")}");
                W($"Spearman rho (NCE vs mean count): {spearman:F4}");
                W($"  (Expected: strong positive correlation, rho > 0.8)");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 7: Internal Fragment TicNI Fraction of Total TIC vs NCE
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(7)]
        public void Test07_InternalTicFraction_VsNCE()
        {
            W(Sep("TEST 7: Total Internal Fragment TicNI Fraction vs NCE"));
            W("  Sum of all internal fragment TicNI per scan ≈ fraction of TIC from internal ions.");
            W("  This should increase with NCE (more energy → more double fragmentation).");
            W();

            W($"{"NCE",4} | {"MeanSumTicNI",14} | {"MedSumTicNI",14} | {"NScans",8}");
            W(new string('-', 55));

            var nceSums = new List<(int nce, double mean)>();

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var scanSums = _datasets[nce].Passing
                    .GroupBy(r => S(r, "ScanNumber"))
                    .Select(g => g.Sum(r => D(r, "TicNormalizedIntensity")))
                    .Where(v => !double.IsNaN(v))
                    .ToList();

                if (scanSums.Count == 0) continue;

                nceSums.Add((nce, scanSums.Average()));
                W($"{nce,4} | {scanSums.Average(),14:E4} | {Median(scanSums),14:E4} | {scanSums.Count,8:N0}");
            }

            W();
            if (nceSums.Count >= 3)
            {
                double spearman = SpearmanRank(
                    nceSums.Select(m => (double)m.nce).ToList(),
                    nceSums.Select(m => m.mean).ToList());
                W($"Spearman rho (NCE vs mean sum TicNI): {spearman:F4}");
                W("  (Expected: positive, more energy → more internal fragmentation)");
                W($"  Fold-change (max/min mean): {nceSums.Max(x => x.mean) / nceSums.Min(x => x.mean):F2}");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 8: Feature Importance Stability — NCE-agnosticism Test
        //   For each pair of NCEs, rank features by |r| and compute Spearman
        //   rank correlation of those rankings
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(8)]
        public void Test08_FeatureImportanceRankStability()
        {
            W(Sep("TEST 8: Feature Importance Rank Stability Across NCEs"));
            W("  If the model is NCE-agnostic, feature importance rankings (by |r| with TicNI)");
            W("  should be similar at every NCE. High rank correlation = stable model.");
            W();

            var featureNames = new[]
            {
                "FragmentLength", "PeptideLength", "PrecursorCharge",
                "RelativeDistanceFromCTerm", "NumberOfBasicResidues",
                "BIonIntensityAtNTerm", "YIonIntensityAtCTerm",
                "MaxTerminalIonIntensity", "BYProductScore",
                "NTerminalFlankingHydrophobicity", "LocalIntensityRank",
            };

            Func<Dictionary<string, string>, string, double> featureVal = (row, feat) =>
            {
                return feat switch
                {
                    "FragmentLength" => I(row, feat),
                    "PeptideLength" => I(row, feat),
                    "PrecursorCharge" => I(row, feat),
                    "NumberOfBasicResidues" => I(row, feat),
                    _ => D(row, feat)
                };
            };

            var nceList = AllNCEs.Where(n => _datasets.ContainsKey(n)).ToList();

            // Compute |r| for each feature at each NCE
            var nceRanks = new Dictionary<int, List<double>>(); // NCE → |r| in same feature order
            foreach (var nce in nceList)
            {
                var ds = _datasets[nce];
                var target = ds.Passing.Select(r => D(r, "TicNormalizedIntensity")).ToList();
                var absRs = new List<double>();

                foreach (var feat in featureNames)
                {
                    var vals = ds.Passing.Select(r => featureVal(r, feat)).ToList();
                    double r = Pearson(vals, target);
                    absRs.Add(double.IsNaN(r) ? 0 : Math.Abs(r));
                }

                nceRanks[nce] = absRs;
            }

            // Pairwise Spearman rank correlation of feature importance rankings
            W("Pairwise Spearman rho of feature |r| rankings:");
            var hdr = $"{"",6}";
            foreach (var nce in nceList) hdr += $" | {nce,7}";
            W(hdr);
            W(new string('-', hdr.Length));

            var allPairwise = new List<double>();
            foreach (var nce1 in nceList)
            {
                var row = $"{nce1,6}";
                foreach (var nce2 in nceList)
                {
                    if (nce2 < nce1) { row += $" | {"",7}"; continue; }
                    if (nce1 == nce2) { row += $" | {"1.0000",7}"; continue; }

                    double rho = SpearmanRank(nceRanks[nce1], nceRanks[nce2]);
                    row += $" | {rho,7:F4}";
                    allPairwise.Add(rho);
                }
                W(row);
            }

            if (allPairwise.Count > 0)
            {
                W();
                W($"Mean pairwise Spearman rho: {allPairwise.Average():F4}");
                W($"Min pairwise Spearman rho:  {allPairwise.Min():F4}");
                W($"  (> 0.85 suggests strong NCE-agnosticism of feature structure)");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 9: Weak-Bond Enrichment Stability Across NCE
        //   Proline and Aspartate effects should hold at all NCEs
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(9)]
        public void Test09_WeakBondEnrichment_AcrossNCE()
        {
            W(Sep("TEST 9: Proline/Aspartate Weak-Bond Enrichment Across NCEs"));
            W("  The proline effect (higher TicNI at proline) should be NCE-independent.");
            W("  Similarly for aspartate-proximal fragmentation.");
            W();

            W($"{"NCE",4} | {"Pro+ Mean",10} | {"Pro- Mean",10} | {"Pro Ratio",10} | {"Asp+ Mean",10} | {"Asp- Mean",10} | {"Asp Ratio",10}");
            W(new string('-', 80));

            var proRatios = new List<(int nce, double ratio)>();
            var aspRatios = new List<(int nce, double ratio)>();

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var ds = _datasets[nce];

                var proYes = ds.Passing.Where(r => B(r, "HasProlineAtEitherTerminus")).Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                var proNo = ds.Passing.Where(r => !B(r, "HasProlineAtEitherTerminus")).Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                double proMeanYes = proYes.Count > 0 ? proYes.Average() : 0;
                double proMeanNo = proNo.Count > 0 ? proNo.Average() : 0;
                double proRatio = proMeanNo > 0 ? proMeanYes / proMeanNo : 0;

                var aspYes = ds.Passing.Where(r => B(r, "HasAspartateAtEitherTerminus")).Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                var aspNo = ds.Passing.Where(r => !B(r, "HasAspartateAtEitherTerminus")).Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                double aspMeanYes = aspYes.Count > 0 ? aspYes.Average() : 0;
                double aspMeanNo = aspNo.Count > 0 ? aspNo.Average() : 0;
                double aspRatio = aspMeanNo > 0 ? aspMeanYes / aspMeanNo : 0;

                proRatios.Add((nce, proRatio));
                aspRatios.Add((nce, aspRatio));

                W($"{nce,4} | {proMeanYes,10:E3} | {proMeanNo,10:E3} | {proRatio,10:F3} | {aspMeanYes,10:E3} | {aspMeanNo,10:E3} | {aspRatio,10:F3}");
            }

            W();
            if (proRatios.Count >= 2)
            {
                var pRatios = proRatios.Select(p => p.ratio).ToList();
                var aRatios = aspRatios.Select(p => p.ratio).ToList();
                W($"Proline ratio:   mean={pRatios.Average():F3}, CV={StdDev(pRatios) / pRatios.Average():F3}, range=[{pRatios.Min():F3}–{pRatios.Max():F3}]");
                W($"Aspartate ratio: mean={aRatios.Average():F3}, CV={StdDev(aRatios) / aRatios.Average():F3}, range=[{aRatios.Min():F3}–{aRatios.Max():F3}]");
                W("  (Low CV = effect is NCE-stable; ratio > 1 = expected enrichment)");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 10: Fragment Length Distribution vs NCE
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(10)]
        public void Test10_FragmentLengthDistribution_VsNCE()
        {
            W(Sep("TEST 10: Fragment Length Distribution vs NCE"));
            W("  Higher NCE may shift the distribution toward shorter internal fragments.");
            W();

            W($"{"NCE",4} | {"MeanLen",8} | {"MedLen",8} | {"Len3%",8} | {"Len4%",8} | {"Len5%",8} | {"Len6%",8} | {"Len7+%",8}");
            W(new string('-', 75));

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var lens = _datasets[nce].Passing.Select(r => I(r, "FragmentLength")).ToList();
                if (lens.Count == 0) continue;

                double n_ = lens.Count;
                double meanLen = lens.Average();
                double medLen = Median(lens.Select(v => (double)v));

                W($"{nce,4} | {meanLen,8:F2} | {medLen,8:F1} | {100.0 * lens.Count(l => l == 3) / n_,7:F1}% | {100.0 * lens.Count(l => l == 4) / n_,7:F1}% | {100.0 * lens.Count(l => l == 5) / n_,7:F1}% | {100.0 * lens.Count(l => l == 6) / n_,7:F1}% | {100.0 * lens.Count(l => l >= 7) / n_,7:F1}%");
            }

            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 11: B/Y Terminal Ion Support vs NCE
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(11)]
        public void Test11_TerminalIonSupport_VsNCE()
        {
            W(Sep("TEST 11: Terminal Ion Support for Internal Fragments vs NCE"));
            W("  Fraction of internal fragments with matched b-ion or y-ion at their termini.");
            W("  This matters because the model uses terminal ion features.");
            W();

            W($"{"NCE",4} | {"HasB%",8} | {"HasY%",8} | {"Both%",8} | {"Neither%",10} | {"MeanTicNI_Both",16} | {"MeanTicNI_Neither",18}");
            W(new string('-', 95));

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var ds = _datasets[nce];
                double n_ = ds.Passing.Count;
                if (n_ == 0) continue;

                int hasB = ds.Passing.Count(r => B(r, "HasMatchedBIonAtNTerm"));
                int hasY = ds.Passing.Count(r => B(r, "HasMatchedYIonAtCTerm"));
                int both = ds.Passing.Count(r => B(r, "HasMatchedBIonAtNTerm") && B(r, "HasMatchedYIonAtCTerm"));
                int neither = ds.Passing.Count(r => !B(r, "HasMatchedBIonAtNTerm") && !B(r, "HasMatchedYIonAtCTerm"));

                var ticBoth = ds.Passing
                    .Where(r => B(r, "HasMatchedBIonAtNTerm") && B(r, "HasMatchedYIonAtCTerm"))
                    .Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                var ticNeither = ds.Passing
                    .Where(r => !B(r, "HasMatchedBIonAtNTerm") && !B(r, "HasMatchedYIonAtCTerm"))
                    .Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();

                W($"{nce,4} | {100.0 * hasB / n_,7:F1}% | {100.0 * hasY / n_,7:F1}% | {100.0 * both / n_,7:F1}% | {100.0 * neither / n_,9:F1}% | {(ticBoth.Count > 0 ? ticBoth.Average() : 0),16:E4} | {(ticNeither.Count > 0 ? ticNeither.Average() : 0),18:E4}");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 12: Charge State Distribution vs NCE
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(12)]
        public void Test12_ChargeStateDistribution_VsNCE()
        {
            W(Sep("TEST 12: Precursor Charge State Distribution vs NCE"));
            W();

            W($"{"NCE",4} | {"z=2%",8} | {"z=3%",8} | {"z=4%",8} | {"z=5+%",8} | {"MeanCharge",11}");
            W(new string('-', 55));

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var charges = _datasets[nce].Passing.Select(r => I(r, "PrecursorCharge")).ToList();
                if (charges.Count == 0) continue;
                double n_ = charges.Count;

                W($"{nce,4} | {100.0 * charges.Count(c => c == 2) / n_,7:F1}% | {100.0 * charges.Count(c => c == 3) / n_,7:F1}% | {100.0 * charges.Count(c => c == 4) / n_,7:F1}% | {100.0 * charges.Count(c => c >= 5) / n_,7:F1}% | {charges.Average(),11:F2}");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 13: NCE-Agnosticism Summary Score
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(13)]
        public void Test13_NCE_AgnosticismSummary()
        {
            W(Sep("TEST 13: NCE-Agnosticism Composite Assessment"));
            W("  Aggregates evidence from all previous tests into a summary verdict.");
            W();

            var nceList = AllNCEs.Where(n => _datasets.ContainsKey(n)).ToList();
            if (nceList.Count < 2) { W("Insufficient data for summary."); Flush(); return; }

            int passes = 0;
            int total = 0;

            // Criterion 1: Mean TicNI CV across NCEs < 0.30
            var means = nceList.Select(nce =>
                _datasets[nce].Passing
                    .Select(r => D(r, "TicNormalizedIntensity"))
                    .Where(v => !double.IsNaN(v) && v > 0)
                    .DefaultIfEmpty(0).Average()
            ).ToList();
            double ticniCV = StdDev(means) / means.Average();
            bool c1 = ticniCV < 0.30;
            W($"[{(c1 ? "PASS" : "FAIL")}] Mean TicNI CV across NCEs: {ticniCV:F4} (threshold < 0.30)");
            if (c1) passes++;
            total++;

            // Criterion 2: Peptide-level TicNI correlation between most distant NCEs > 0.5
            var extremes = new[] { nceList.Min(), nceList.Max() };
            var pepMeans = new Dictionary<string, Dictionary<int, double>>();
            foreach (var nce in extremes)
            {
                foreach (var g in _datasets[nce].Passing.GroupBy(r => S(r, "PeptideSequence")).Where(g => !string.IsNullOrEmpty(g.Key)))
                {
                    var vals = g.Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v) && v > 0).ToList();
                    if (vals.Count == 0) continue;
                    if (!pepMeans.ContainsKey(g.Key)) pepMeans[g.Key] = new();
                    pepMeans[g.Key][nce] = vals.Average();
                }
            }
            var shared = pepMeans.Where(kv => kv.Value.Count == 2).ToList();
            double extremeR = shared.Count >= 10
                ? Pearson(shared.Select(kv => kv.Value[extremes[0]]).ToList(), shared.Select(kv => kv.Value[extremes[1]]).ToList())
                : double.NaN;
            bool c2 = !double.IsNaN(extremeR) && extremeR > 0.5;
            W($"[{(c2 ? "PASS" : "FAIL")}] Peptide TicNI r (NCE {extremes[0]} vs {extremes[1]}): {(double.IsNaN(extremeR) ? "N/A" : extremeR.ToString("F4"))} (threshold > 0.50, N={shared.Count})");
            if (c2) passes++;
            total++;

            // Criterion 3: Proline enrichment ratio CV < 0.20
            var proRatios = new List<double>();
            foreach (var nce in nceList)
            {
                var proYes = _datasets[nce].Passing.Where(r => B(r, "HasProlineAtEitherTerminus")).Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                var proNo = _datasets[nce].Passing.Where(r => !B(r, "HasProlineAtEitherTerminus")).Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
                if (proYes.Count > 0 && proNo.Count > 0 && proNo.Average() > 0)
                    proRatios.Add(proYes.Average() / proNo.Average());
            }
            double proCV = proRatios.Count >= 2 ? StdDev(proRatios) / proRatios.Average() : double.NaN;
            bool c3 = !double.IsNaN(proCV) && proCV < 0.20;
            W($"[{(c3 ? "PASS" : "FAIL")}] Proline enrichment ratio CV: {(double.IsNaN(proCV) ? "N/A" : proCV.ToString("F4"))} (threshold < 0.20)");
            if (c3) passes++;
            total++;

            // Criterion 4: Internal fragment count increases with NCE (Spearman > 0.7)
            var scanCounts = nceList.Select(nce => (
                nce: (double)nce,
                mean: _datasets[nce].Passing.GroupBy(r => S(r, "ScanNumber")).Select(g => (double)g.Count()).DefaultIfEmpty(0).Average()
            )).ToList();
            double countSpearman = SpearmanRank(scanCounts.Select(s => s.nce).ToList(), scanCounts.Select(s => s.mean).ToList());
            bool c4 = countSpearman > 0.7;
            W($"[{(c4 ? "PASS" : "FAIL")}] Spearman rho (NCE vs internal count): {countSpearman:F4} (threshold > 0.70)");
            if (c4) passes++;
            total++;

            W();
            W($"═══════════════════════════════════════════════════════");
            W($"  NCE-AGNOSTICISM VERDICT: {passes}/{total} criteria passed");
            if (passes == total)
                W("  CONCLUSION: Strong evidence for NCE-agnostic internal fragment model.");
            else if (passes >= total - 1)
                W("  CONCLUSION: Mostly NCE-agnostic with minor NCE sensitivity.");
            else
                W("  CONCLUSION: Significant NCE dependence detected — model may need NCE as a feature.");
            W($"═══════════════════════════════════════════════════════");
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 14: Per-NCE Intensity Quantile Calibration
        //   Bin observed TicNI into quantiles; check if mean TicNI per quantile
        //   is stable across NCEs (i.e., the intensity scale is comparable)
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(14)]
        public void Test14_IntensityQuantileCalibration()
        {
            W(Sep("TEST 14: Intensity Quantile Calibration Across NCEs"));
            W("  For each NCE, split internal ions into TicNI quintiles.");
            W("  Report mean TicNI per quintile. If the model is NCE-agnostic,");
            W("  the quintile means should be comparable across NCEs.");
            W();

            var nceList = AllNCEs.Where(n => _datasets.ContainsKey(n)).ToList();
            int nQuantiles = 5;

            // Header
            var hdr = $"{"Quintile",10}";
            foreach (var nce in nceList) hdr += $" | {"NCE" + nce,10}";
            W(hdr);
            W(new string('-', hdr.Length));

            // Compute quintile means for each NCE
            var quantileMeans = new Dictionary<int, double[]>(); // NCE → array of quintile means
            foreach (var nce in nceList)
            {
                var vals = _datasets[nce].Passing
                    .Select(r => D(r, "TicNormalizedIntensity"))
                    .Where(v => !double.IsNaN(v) && v > 0)
                    .OrderBy(v => v).ToList();

                if (vals.Count < nQuantiles * 2) continue;

                var qMeans = new double[nQuantiles];
                int binSize = vals.Count / nQuantiles;
                for (int q = 0; q < nQuantiles; q++)
                {
                    int start = q * binSize;
                    int end = q == nQuantiles - 1 ? vals.Count : (q + 1) * binSize;
                    qMeans[q] = vals.GetRange(start, end - start).Average();
                }
                quantileMeans[nce] = qMeans;
            }

            for (int q = 0; q < nQuantiles; q++)
            {
                string label = q == 0 ? "Q1 (low)" : q == nQuantiles - 1 ? $"Q{nQuantiles} (high)" : $"Q{q + 1}";
                var row = $"{label,10}";
                foreach (var nce in nceList)
                {
                    if (quantileMeans.ContainsKey(nce))
                        row += $" | {quantileMeans[nce][q],10:E3}";
                    else
                        row += $" | {"N/A",10}";
                }
                W(row);
            }

            // Cross-NCE CV per quintile
            W();
            W("Cross-NCE CV per quintile:");
            for (int q = 0; q < nQuantiles; q++)
            {
                var vals = nceList.Where(n => quantileMeans.ContainsKey(n)).Select(n => quantileMeans[n][q]).ToList();
                if (vals.Count >= 2)
                {
                    double cv = StdDev(vals) / vals.Average();
                    W($"  Q{q + 1}: CV = {cv:F4}");
                }
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // TEST 15: Mass Accuracy Filter Pass Rate vs NCE
        // ════════════════════════════════════════════════════════════════════════

        [Test, Order(15)]
        public void Test15_MassAccuracyPassRate_VsNCE()
        {
            W(Sep("TEST 15: Mass Accuracy Filter Pass Rate vs NCE"));
            W("  If higher NCE produces noisier internal fragments, pass rate may drop.");
            W();

            W($"{"NCE",4} | {"Total",8} | {"Pass",8} | {"PassRate",10}");
            W(new string('-', 40));

            foreach (var nce in AllNCEs.Where(n => _datasets.ContainsKey(n)))
            {
                var ds = _datasets[nce];
                double rate = ds.AllRows.Count > 0 ? 100.0 * ds.Passing.Count / ds.AllRows.Count : 0;
                W($"{nce,4} | {ds.AllRows.Count,8:N0} | {ds.Passing.Count,8:N0} | {rate,9:F1}%");
            }
            W();
            Flush();
        }

        // ════════════════════════════════════════════════════════════════════════
        // HELPER METHODS
        // ════════════════════════════════════════════════════════════════════════

        #region Helpers

        private void W(string line = "") => _out.WriteLine(line);
        private void Flush() => _out.Flush();

        private static string Sep(string title)
            => $"\n{"",0}{new string('═', 66)}\n{title}\n{new string('═', 66)}\n";

        private static double D(Dictionary<string, string> row, string col)
        {
            if (!row.TryGetValue(col, out var val)) return 0.0;
            if (string.IsNullOrWhiteSpace(val)) return 0.0;
            if (val.Equals("NaN", StringComparison.OrdinalIgnoreCase)) return double.NaN;
            return double.TryParse(val, NumberStyles.Any, CultureInfo.InvariantCulture, out var d) ? d : 0.0;
        }

        private static int I(Dictionary<string, string> row, string col)
        {
            if (!row.TryGetValue(col, out var val)) return 0;
            return int.TryParse(val, out var i) ? i : 0;
        }

        private static bool B(Dictionary<string, string> row, string col)
        {
            if (!row.TryGetValue(col, out var val)) return false;
            return val.Equals("TRUE", StringComparison.OrdinalIgnoreCase) || val == "1";
        }

        private static string S(Dictionary<string, string> row, string col) =>
            row.TryGetValue(col, out var val) ? val : "";

        private static double Pearson(IList<double> xs, IList<double> ys)
        {
            if (xs.Count != ys.Count || xs.Count < 3) return double.NaN;
            var valid = xs.Zip(ys, (x, y) => (x, y))
                .Where(p => !double.IsNaN(p.x) && !double.IsNaN(p.y)).ToList();
            if (valid.Count < 3) return double.NaN;

            double mx = valid.Average(p => p.x), my = valid.Average(p => p.y);
            double sxy = 0, sx2 = 0, sy2 = 0;
            foreach (var (x, y) in valid)
            {
                double dx = x - mx, dy = y - my;
                sxy += dx * dy; sx2 += dx * dx; sy2 += dy * dy;
            }
            double denom = Math.Sqrt(sx2 * sy2);
            return denom < 1e-12 ? double.NaN : sxy / denom;
        }

        private static double SpearmanRank(IList<double> xs, IList<double> ys)
        {
            if (xs.Count != ys.Count || xs.Count < 3) return double.NaN;
            var rankX = Rank(xs);
            var rankY = Rank(ys);
            return Pearson(rankX, rankY);
        }

        private static List<double> Rank(IList<double> vals)
        {
            var indexed = vals.Select((v, i) => (v, i)).OrderBy(x => x.v).ToList();
            var ranks = new double[vals.Count];
            for (int i = 0; i < indexed.Count;)
            {
                int j = i;
                while (j < indexed.Count && indexed[j].v == indexed[i].v) j++;
                double avgRank = (i + j - 1) / 2.0 + 1.0;
                for (int k = i; k < j; k++)
                    ranks[indexed[k].i] = avgRank;
                i = j;
            }
            return ranks.ToList();
        }

        private static double StdDev(IEnumerable<double> vals)
        {
            var list = vals.Where(v => !double.IsNaN(v)).ToList();
            if (list.Count < 2) return 0;
            double mean = list.Average();
            return Math.Sqrt(list.Select(v => (v - mean) * (v - mean)).Average());
        }

        private static double Median(IEnumerable<double> vals)
        {
            var sorted = vals.Where(v => !double.IsNaN(v)).OrderBy(v => v).ToList();
            if (sorted.Count == 0) return 0;
            return sorted.Count % 2 == 0
                ? (sorted[sorted.Count / 2 - 1] + sorted[sorted.Count / 2]) / 2.0
                : sorted[sorted.Count / 2];
        }

        private static double Percentile(IList<double> sorted, double p)
        {
            if (sorted.Count == 0) return 0;
            double index = (p / 100.0) * (sorted.Count - 1);
            int lower = (int)Math.Floor(index);
            int upper = (int)Math.Ceiling(index);
            if (lower == upper || upper >= sorted.Count) return sorted[lower];
            double frac = index - lower;
            return sorted[lower] * (1 - frac) + sorted[upper] * frac;
        }

        #endregion
    }
}