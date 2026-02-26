// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Development.Dia
{
    /// <summary>
    /// Benchmarks the DIA extraction engine against DIA-NN ground truth.
    /// 
    /// Workflow:
    ///   1. Load a real DIA mzML file (or slice)
    ///   2. Build DiaScanIndex — validate window count, scan count, RT range
    ///   3. Load DIA-NN ground truth TSV → LibraryPrecursorInput[]
    ///   4. Filter to precursors whose RT falls within the mzML's RT range
    ///   5. Generate queries, extract, score
    ///   6. Report correctness metrics:
    ///      - Fragment detection rate
    ///      - Score distributions
    ///      - Precursor recovery rate
    ///      - Timing
    /// 
    /// Ground truth TSV format (produced by convert_diann_lib_to_benchmark_tsv.py):
    ///   Sequence  PrecursorMz  Charge  RT  IsDecoy  ProteinGroup  Gene  QValue  FragmentMzs  FragmentIntensities  FragmentAnnotations
    ///   Fields are tab-delimited. FragmentMzs/Intensities are semicolon-delimited floats.
    /// 
    /// Usage:
    ///   RealDataBenchmark.Run(
    ///       mzmlPath: @"Test\Dia\DiaTestData\Fig2HeLa-0-5h_MHRM_R01_T0_Scans0to500.mzML",
    ///       groundTruthTsvPath: @"Test\Dia\DiaTestData\diann_ground_truth.tsv"
    ///   );
    /// </summary>
    public static class RealDataBenchmark
    {
        public static void Run(string mzmlPath, string groundTruthTsvPath,
            float ppmTolerance = 20f, float rtToleranceMinutes = 2.0f, int minFragments = 3)
        {
            Console.WriteLine("=== DIA Real Data Benchmark (vs DIA-NN Ground Truth) ===");
            Console.WriteLine();

            // ─── Step 1: Load mzML ──────────────────────────────────────────
            Console.WriteLine($"Loading: {Path.GetFileName(mzmlPath)}");
            var sw = Stopwatch.StartNew();

            var msDataFile = MsDataFileReader.GetDataFile(mzmlPath);
            msDataFile.LoadAllStaticData();
            var allScans = msDataFile.GetAllScansList().ToArray();

            var loadTime = sw.Elapsed;
            Console.WriteLine($"  Load time:    {loadTime.TotalSeconds:F2}s");
            Console.WriteLine($"  Total scans:  {allScans.Length}");

            int ms1Count = allScans.Count(s => s.MsnOrder == 1);
            int ms2Count = allScans.Count(s => s.MsnOrder == 2);
            Console.WriteLine($"  MS1 scans:    {ms1Count}");
            Console.WriteLine($"  MS2 scans:    {ms2Count}");

            float firstRt = (float)allScans.First().RetentionTime;
            float lastRt = (float)allScans.Last().RetentionTime;
            Console.WriteLine($"  RT range:     {firstRt:F2} – {lastRt:F2} min");
            Console.WriteLine();

            // ─── Step 2: Build SoA index ────────────────────────────────────
            sw.Restart();
            using var index = DiaScanIndexBuilder.Build(allScans);
            var buildTime = sw.Elapsed;

            Console.WriteLine("SoA Index built:");
            Console.WriteLine($"  Build time:   {buildTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Windows:      {index.WindowCount}");
            Console.WriteLine($"  MS2 scans:    {index.ScanCount}");
            Console.WriteLine($"  Total peaks:  {index.TotalPeakCount:N0}");
            Console.WriteLine($"  Peaks/scan:   {(index.ScanCount > 0 ? index.TotalPeakCount / index.ScanCount : 0)}");
            Console.WriteLine($"  Build rate:   {index.TotalPeakCount / buildTime.TotalSeconds:N0} peaks/sec");

            float indexRtMin = index.GetGlobalRtMin();
            float indexRtMax = index.GetGlobalRtMax();
            Console.WriteLine($"  Index RT:     {indexRtMin:F2} – {indexRtMax:F2} min");

            // Print window boundaries
            Console.WriteLine($"  Window boundaries:");
            foreach (int wid in index.GetWindowIds().OrderBy(x => x))
            {
                var (lo, hi) = index.GetWindowBounds(wid);
                index.TryGetScanRangeForWindow(wid, out _, out int count);
                Console.WriteLine($"    Window {wid,2}: {lo:F1} – {hi:F1} m/z ({count} scans)");
            }
            Console.WriteLine();

            // ─── Step 3: Load ground truth ──────────────────────────────────
            Console.WriteLine($"Loading ground truth: {Path.GetFileName(groundTruthTsvPath)}");
            var allPrecursors = LoadGroundTruthTsv(groundTruthTsvPath);
            Console.WriteLine($"  Total precursors: {allPrecursors.Count:N0}");
            Console.WriteLine($"  RT range:         {allPrecursors.Min(p => p.RetentionTime):F2} – " +
                              $"{allPrecursors.Max(p => p.RetentionTime):F2} min");

            // Filter to precursors whose RT falls within the mzML slice (with tolerance)
            float rtFilterMin = indexRtMin - rtToleranceMinutes;
            float rtFilterMax = indexRtMax + rtToleranceMinutes;
            var precursorsInRange = allPrecursors
                .Where(p => p.RetentionTime >= rtFilterMin && p.RetentionTime <= rtFilterMax)
                .ToList();
            Console.WriteLine($"  In RT range:      {precursorsInRange.Count:N0} " +
                              $"(filter: {rtFilterMin:F2} – {rtFilterMax:F2} min)");

            // Also filter: precursor m/z must map to a valid window
            var precursorsWithWindow = precursorsInRange
                .Where(p => index.FindWindowForPrecursorMz(p.PrecursorMz) >= 0)
                .ToList();
            Console.WriteLine($"  With valid window: {precursorsWithWindow.Count:N0}");
            Console.WriteLine();

            if (precursorsWithWindow.Count == 0)
            {
                Console.WriteLine("ERROR: No precursors fall within the mzML slice's RT range + windows.");
                Console.WriteLine("This likely means the 500-scan slice covers very early gradient time");
                Console.WriteLine("where few peptides elute. Try using a larger slice or the full file.");
                return;
            }

            // ─── Step 4: Generate queries ───────────────────────────────────
            var parameters = new DiaSearchParameters
            {
                PpmTolerance = ppmTolerance,
                RtToleranceMinutes = rtToleranceMinutes,
                MinFragmentsRequired = minFragments,
                MinScoreThreshold = 0f, // Keep everything for analysis
            };

            sw.Restart();
            var genResult = DiaLibraryQueryGenerator.Generate(precursorsWithWindow, index, parameters);
            var genTime = sw.Elapsed;

            Console.WriteLine("Query generation:");
            Console.WriteLine($"  Time:         {genTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Queries:      {genResult.Queries.Length:N0}");
            Console.WriteLine($"  Groups:       {genResult.PrecursorGroups.Length:N0}");
            Console.WriteLine($"  Skipped (no window): {genResult.SkippedNoWindow}");
            Console.WriteLine($"  Skipped (no frags):  {genResult.SkippedNoFragments}");
            Console.WriteLine();

            // ─── Step 5: Extract ────────────────────────────────────────────
            using var extractor = new CpuFragmentExtractor(index);
            var results = new FragmentResult[genResult.Queries.Length];
            int bufferSize = Math.Max(genResult.Queries.Length * 200, 1024);
            var rtBuf = new float[bufferSize];
            var intBuf = new float[bufferSize];

            sw.Restart();
            int totalPoints = extractor.ExtractBatch(genResult.Queries, results, rtBuf, intBuf);
            var extractTime = sw.Elapsed;

            Console.WriteLine("Extraction:");
            Console.WriteLine($"  Time:         {extractTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Data points:  {totalPoints:N0}");
            Console.WriteLine($"  Queries/sec:  {genResult.Queries.Length / extractTime.TotalSeconds:N0}");
            Console.WriteLine();

            // ─── Step 6: Score ──────────────────────────────────────────────
            var dotScorer = new NormalizedDotProductScorer();
            var saScorer = new SpectralAngleScorer();

            sw.Restart();
            var extractionResult = new ExtractionResult(results, rtBuf, intBuf, 0);
            var searchResults = DiaLibraryQueryGenerator.AssembleResults(
                precursorsWithWindow, genResult, extractionResult, parameters, dotScorer, saScorer);
            var scoreTime = sw.Elapsed;

            Console.WriteLine("Scoring + Assembly:");
            Console.WriteLine($"  Time:         {scoreTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Results:      {searchResults.Count:N0} (passing MinFragments={minFragments})");
            Console.WriteLine();

            // ─── Step 7: Evaluate against ground truth ──────────────────────
            EvaluateResults(searchResults, precursorsWithWindow, genResult);

            // ─── Summary ────────────────────────────────────────────────────
            Console.WriteLine("=== Timing Summary ===");
            Console.WriteLine($"  mzML load:       {loadTime.TotalSeconds:F2}s");
            Console.WriteLine($"  Index build:     {buildTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Query gen:       {genTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Extraction:      {extractTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Scoring:         {scoreTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Total pipeline:  {(loadTime + buildTime + genTime + extractTime + scoreTime).TotalSeconds:F2}s");
        }

        // ─── Evaluation ─────────────────────────────────────────────────────

        private static void EvaluateResults(
            List<DiaSearchResult> searchResults,
            List<LibraryPrecursorInput> groundTruth,
            DiaLibraryQueryGenerator.GenerationResult genResult)
        {
            Console.WriteLine("=== Correctness Metrics (vs DIA-NN Ground Truth) ===");
            Console.WriteLine();

            // Metric 1: Precursor recovery rate
            float recoveryRate = (float)searchResults.Count / groundTruth.Count;
            Console.WriteLine($"Precursor Recovery:");
            Console.WriteLine($"  Ground truth in range:  {groundTruth.Count}");
            Console.WriteLine($"  Recovered (>= min frags): {searchResults.Count}");
            Console.WriteLine($"  Recovery rate:          {recoveryRate:P1}");
            Console.WriteLine();

            if (searchResults.Count == 0)
            {
                Console.WriteLine("No results to evaluate. The mzML slice may not cover enough RT range.");
                return;
            }

            // Metric 2: Fragment detection rate
            var detectionRates = searchResults.Select(r => r.FragmentDetectionRate).ToArray();
            Console.WriteLine($"Fragment Detection Rate:");
            Console.WriteLine($"  Mean:   {detectionRates.Average():P1}");
            Console.WriteLine($"  Median: {Median(detectionRates):P1}");
            Console.WriteLine($"  Min:    {detectionRates.Min():P1}");
            Console.WriteLine($"  Max:    {detectionRates.Max():P1}");
            Console.WriteLine();

            // Metric 3: Score distributions
            var dotScores = searchResults
                .Where(r => !float.IsNaN(r.DotProductScore))
                .Select(r => r.DotProductScore).ToArray();
            var saScores = searchResults
                .Where(r => !float.IsNaN(r.SpectralAngleScore))
                .Select(r => r.SpectralAngleScore).ToArray();

            if (dotScores.Length > 0)
            {
                Console.WriteLine($"Dot Product Score (library vs extracted intensities):");
                Console.WriteLine($"  Mean:   {dotScores.Average():F4}");
                Console.WriteLine($"  Median: {Median(dotScores):F4}");
                Console.WriteLine($"  Q25:    {Percentile(dotScores, 25):F4}");
                Console.WriteLine($"  Q75:    {Percentile(dotScores, 75):F4}");
                Console.WriteLine();
            }

            if (saScores.Length > 0)
            {
                Console.WriteLine($"Spectral Angle Score:");
                Console.WriteLine($"  Mean:   {saScores.Average():F4}");
                Console.WriteLine($"  Median: {Median(saScores):F4}");
                Console.WriteLine($"  Q25:    {Percentile(saScores, 25):F4}");
                Console.WriteLine($"  Q75:    {Percentile(saScores, 75):F4}");
                Console.WriteLine();
            }

            // Metric 4: Fragments detected distribution
            var fragsDetected = searchResults.Select(r => (float)r.FragmentsDetected).ToArray();
            var fragsQueried = searchResults.Select(r => (float)r.FragmentsQueried).ToArray();
            Console.WriteLine($"Fragments per Precursor:");
            Console.WriteLine($"  Queried:  mean={fragsQueried.Average():F1}, median={Median(fragsQueried):F0}");
            Console.WriteLine($"  Detected: mean={fragsDetected.Average():F1}, median={Median(fragsDetected):F0}");
            Console.WriteLine();

            // Metric 5: Top hits (sanity check — print best-scoring precursors)
            Console.WriteLine("Top 10 precursors by dot product score:");
            var top10 = searchResults
                .Where(r => !float.IsNaN(r.DotProductScore))
                .OrderByDescending(r => r.DotProductScore)
                .Take(10);
            foreach (var r in top10)
            {
                Console.WriteLine($"  {r.Sequence,-30} z={r.ChargeState} " +
                    $"DP={r.DotProductScore:F3} SA={r.SpectralAngleScore:F3} " +
                    $"frags={r.FragmentsDetected}/{r.FragmentsQueried}");
            }
            Console.WriteLine();

            // Metric 6: Bottom hits (what's failing?)
            Console.WriteLine("Bottom 10 precursors by dot product score (still passing):");
            var bottom10 = searchResults
                .Where(r => !float.IsNaN(r.DotProductScore))
                .OrderBy(r => r.DotProductScore)
                .Take(10);
            foreach (var r in bottom10)
            {
                Console.WriteLine($"  {r.Sequence,-30} z={r.ChargeState} " +
                    $"DP={r.DotProductScore:F3} SA={r.SpectralAngleScore:F3} " +
                    $"frags={r.FragmentsDetected}/{r.FragmentsQueried}");
            }
            Console.WriteLine();
        }

        // ─── Ground Truth Loader ────────────────────────────────────────────

        /// <summary>
        /// Loads the ground truth TSV produced by convert_diann_lib_to_benchmark_tsv.py.
        /// One row per precursor, semicolon-delimited fragment arrays.
        /// </summary>
        public static List<LibraryPrecursorInput> LoadGroundTruthTsv(string path)
        {
            var precursors = new List<LibraryPrecursorInput>();
            bool headerRead = false;
            int colSequence = -1, colPrecursorMz = -1, colCharge = -1, colRt = -1;
            int colIsDecoy = -1, colFragMzs = -1, colFragInts = -1;

            foreach (var line in File.ReadLines(path))
            {
                if (!headerRead)
                {
                    var headers = line.Split('\t');
                    for (int i = 0; i < headers.Length; i++)
                    {
                        switch (headers[i].Trim())
                        {
                            case "Sequence": colSequence = i; break;
                            case "PrecursorMz": colPrecursorMz = i; break;
                            case "Charge": colCharge = i; break;
                            case "RT": colRt = i; break;
                            case "IsDecoy": colIsDecoy = i; break;
                            case "FragmentMzs": colFragMzs = i; break;
                            case "FragmentIntensities": colFragInts = i; break;
                        }
                    }
                    headerRead = true;
                    continue;
                }

                var cols = line.Split('\t');
                if (cols.Length <= colFragInts) continue;

                string sequence = cols[colSequence];
                double precursorMz = double.Parse(cols[colPrecursorMz], CultureInfo.InvariantCulture);
                int charge = int.Parse(cols[colCharge]);
                double rt = double.Parse(cols[colRt], CultureInfo.InvariantCulture);
                bool isDecoy = colIsDecoy >= 0 && cols[colIsDecoy] == "1";

                float[] fragMzs = cols[colFragMzs]
                    .Split(';', StringSplitOptions.RemoveEmptyEntries)
                    .Select(s => float.Parse(s, CultureInfo.InvariantCulture))
                    .ToArray();
                float[] fragInts = cols[colFragInts]
                    .Split(';', StringSplitOptions.RemoveEmptyEntries)
                    .Select(s => float.Parse(s, CultureInfo.InvariantCulture))
                    .ToArray();

                if (fragMzs.Length == 0 || fragMzs.Length != fragInts.Length)
                    continue;

                precursors.Add(new LibraryPrecursorInput(
                    sequence, precursorMz, charge, rt, isDecoy, fragMzs, fragInts));
            }

            return precursors;
        }

        // ─── Stats Helpers ──────────────────────────────────────────────────

        private static float Median(float[] values)
        {
            if (values.Length == 0) return float.NaN;
            var sorted = values.OrderBy(x => x).ToArray();
            int mid = sorted.Length / 2;
            return sorted.Length % 2 == 0
                ? (sorted[mid - 1] + sorted[mid]) / 2f
                : sorted[mid];
        }

        private static float Percentile(float[] values, int percentile)
        {
            if (values.Length == 0) return float.NaN;
            var sorted = values.OrderBy(x => x).ToArray();
            int idx = (int)(sorted.Length * percentile / 100.0);
            idx = Math.Clamp(idx, 0, sorted.Length - 1);
            return sorted[idx];
        }
    }
}