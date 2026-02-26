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
    ///   5. Pass 1: Broad extraction with wide RT windows → anchor selection → RT calibration fit
    ///   6. Pass 2: Narrow extraction with calibrated RT windows → apex intensity scoring
    ///   7. Report correctness metrics:
    ///      - Fragment detection rate
    ///      - Score distributions (summed vs apex, broad vs calibrated)
    ///      - Precursor recovery rate
    ///      - Timing
    /// 
    /// Ground truth TSV format (produced by convert_diann_lib_to_benchmark_tsv.py):
    ///   Sequence  PrecursorMz  Charge  RT  IsDecoy  ProteinGroup  Gene  QValue  FragmentMzs  FragmentIntensities  FragmentAnnotations
    ///   Fields are tab-delimited. FragmentMzs/Intensities are semicolon-delimited floats.
    /// 
    /// Usage:
    ///   RealDataBenchmark.Run(
    ///       mzmlPath: @"Test\Dia\DiaTestData\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
    ///       groundTruthTsvPath: @"Test\Dia\DiaTestData\diann_ground_truth.tsv"
    ///   );
    /// </summary>
    public static class RealDataBenchmark
    {
        public static void Run(string mzmlPath, string groundTruthTsvPath,
            float ppmTolerance = 20f, float rtToleranceMinutes = 5.0f, int minFragments = 3)
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

            var dotScorer = new NormalizedDotProductScorer();
            var saScorer = new SpectralAngleScorer();
            int threads = Math.Min(Environment.ProcessorCount, 16);

            // ─── Step 4: Pass 1 — Broad extraction (for RT calibration anchors) ─
            Console.WriteLine("=== Pass 1: Broad Extraction (for RT calibration) ===");

            var broadParams = new DiaSearchParameters
            {
                PpmTolerance = ppmTolerance,
                RtToleranceMinutes = rtToleranceMinutes,
                MinFragmentsRequired = minFragments,
                MinScoreThreshold = 0f,
            };

            sw.Restart();
            var broadGenResult = DiaLibraryQueryGenerator.Generate(precursorsWithWindow, index, broadParams);
            var broadGenTime = sw.Elapsed;

            Console.WriteLine($"  Query generation:  {broadGenTime.TotalMilliseconds:F0} ms  " +
                              $"({broadGenResult.Queries.Length:N0} queries, " +
                              $"{broadGenResult.PrecursorGroups.Length:N0} groups)");

            // Use orchestrator for parallel extraction
            sw.Restart();
            ExtractionResult broadExtraction;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
            {
                broadExtraction = orchestrator.ExtractAll(broadGenResult.Queries,
                    maxDegreeOfParallelism: threads);
            }
            var broadExtractTime = sw.Elapsed;

            Console.WriteLine($"  Extraction:        {broadExtractTime.TotalMilliseconds:F0} ms  " +
                              $"({broadExtraction.TotalDataPoints:N0} data points, " +
                              $"{broadGenResult.Queries.Length / broadExtractTime.TotalSeconds:N0} q/sec)");

            // Score Pass 1 with apex intensity
            sw.Restart();
            var broadResults = DiaLibraryQueryGenerator.AssembleResults(
                precursorsWithWindow, broadGenResult, broadExtraction.Results, broadParams,
                dotScorer, saScorer);
            var broadScoreTime = sw.Elapsed;

            var broadDotScores = broadResults
                .Where(r => !float.IsNaN(r.DotProductScore))
                .Select(r => r.DotProductScore).ToArray();

            Console.WriteLine($"  Scoring:           {broadScoreTime.TotalMilliseconds:F0} ms  " +
                              $"({broadResults.Count:N0} results)");
            if (broadDotScores.Length > 0)
            {
                Console.WriteLine($"  DP scores (summed): median={Median(broadDotScores):F4}, " +
                                  $"mean={broadDotScores.Average():F4}, " +
                                  $"Q25={Percentile(broadDotScores, 25):F4}, " +
                                  $"Q75={Percentile(broadDotScores, 75):F4}");
            }
            Console.WriteLine();

            // ─── Step 5: RT Calibration ─────────────────────────────────────
            Console.WriteLine("=== RT Calibration ===");

            // Select high-confidence anchors: DP > 0.6, with library RT
            float anchorDpThreshold = 0.6f;
            var anchors = new List<(double LibraryRt, double ObservedRt)>();

            for (int g = 0; g < broadGenResult.PrecursorGroups.Length; g++)
            {
                // Find the matching search result for this group
                var group = broadGenResult.PrecursorGroups[g];
                var input = precursorsWithWindow[group.InputIndex];

                // Find the result by sequence+charge (results are in same order minus filtered ones)
                var matchingResult = broadResults.FirstOrDefault(r =>
                    r.Sequence == input.Sequence && r.ChargeState == input.ChargeState);

                if (matchingResult == null || float.IsNaN(matchingResult.DotProductScore))
                    continue;
                if (matchingResult.DotProductScore < anchorDpThreshold)
                    continue;
                if (!input.RetentionTime.HasValue)
                    continue;

                // Find the apex RT: scan with highest total fragment intensity
                float bestIntensity = 0f;
                float bestRt = (float)input.RetentionTime.Value;

                // Walk through all fragments of this precursor to find the overall apex RT
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int qi = group.QueryOffset + f;
                    var fr = broadExtraction.Results[qi];
                    if (fr.DataPointCount == 0) continue;

                    for (int p = 0; p < fr.DataPointCount; p++)
                    {
                        float intensity = broadExtraction.IntensityBuffer[fr.IntensityBufferOffset + p];
                        if (intensity > bestIntensity)
                        {
                            bestIntensity = intensity;
                            bestRt = broadExtraction.RtBuffer[fr.RtBufferOffset + p];
                        }
                    }
                }

                anchors.Add((input.RetentionTime.Value, bestRt));
            }

            Console.WriteLine($"  Anchor candidates (DP > {anchorDpThreshold}): {anchors.Count}");

            RtCalibrationModel calibration = null;
            if (anchors.Count >= RtCalibrationModel.MinReliableAnchors)
            {
                var anchorLibRts = anchors.Select(a => a.LibraryRt).ToArray();
                var anchorObsRts = anchors.Select(a => a.ObservedRt).ToArray();

                sw.Restart();
                calibration = RtCalibrationFitter.Fit(
                    (ReadOnlySpan<double>)anchorLibRts,
                    (ReadOnlySpan<double>)anchorObsRts);
                var fitTime = sw.Elapsed;

                if (calibration != null)
                {
                    Console.WriteLine($"  Fit time:          {fitTime.TotalMilliseconds:F1} ms");
                    Console.WriteLine($"  Model:             RT = {calibration.Slope:F4} * libRT + {calibration.Intercept:F4}");
                    Console.WriteLine($"  σ:                 {calibration.SigmaMinutes:F4} min");
                    Console.WriteLine($"  R²:                {calibration.RSquared:F4}");
                    Console.WriteLine($"  Anchors used:      {calibration.AnchorCount}");
                    Console.WriteLine($"  Reliable:          {(calibration.IsReliable ? "YES" : "NO")}");

                    double broadHalfWidth = rtToleranceMinutes;
                    double narrowHalfWidth = calibration.GetMinutesWindowHalfWidth(3.0);
                    Console.WriteLine($"  Window reduction:  ±{broadHalfWidth:F2} min → ±{narrowHalfWidth:F2} min " +
                                      $"({broadHalfWidth / narrowHalfWidth:F1}× narrower)");
                }
                else
                {
                    Console.WriteLine("  Calibration fit returned null — falling back to broad windows.");
                }
            }
            else
            {
                Console.WriteLine($"  Too few anchors ({anchors.Count}) for calibration — " +
                                  $"need >= {RtCalibrationModel.MinReliableAnchors}. Using broad windows.");
            }
            Console.WriteLine();

            // ─── Step 6: Pass 2 — Narrow extraction with calibrated RT windows ─
            Console.WriteLine("=== Pass 2: Calibrated Extraction ===");

            DiaLibraryQueryGenerator.GenerationResult narrowGenResult;
            var narrowParams = new DiaSearchParameters
            {
                PpmTolerance = ppmTolerance,
                RtToleranceMinutes = rtToleranceMinutes, // fallback if calibration failed
                MinFragmentsRequired = minFragments,
                MinScoreThreshold = 0f,
                CalibratedWindowSigmaMultiplier = 3.0,
            };

            sw.Restart();
            if (calibration != null && calibration.IsReliable)
            {
                // Use calibrated generation — maps library RT → predicted RT via linear model
                // Since our ground truth has run-specific RT (not iRT), we treat the library RT
                // as if it were iRT for the calibration model. The model maps libRT → observed RT.
                narrowGenResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                    precursorsWithWindow, index, narrowParams, calibration);
            }
            else
            {
                // Fall back to fixed windows
                narrowGenResult = DiaLibraryQueryGenerator.Generate(
                    precursorsWithWindow, index, narrowParams);
            }
            var narrowGenTime = sw.Elapsed;

            Console.WriteLine($"  Query generation:  {narrowGenTime.TotalMilliseconds:F0} ms  " +
                              $"({narrowGenResult.Queries.Length:N0} queries)");

            sw.Restart();
            ExtractionResult narrowExtraction;
            using (var orchestrator = new DiaExtractionOrchestrator(index))
            {
                narrowExtraction = orchestrator.ExtractAll(narrowGenResult.Queries,
                    maxDegreeOfParallelism: threads);
            }
            var narrowExtractTime = sw.Elapsed;

            Console.WriteLine($"  Extraction:        {narrowExtractTime.TotalMilliseconds:F0} ms  " +
                              $"({narrowExtraction.TotalDataPoints:N0} data points, " +
                              $"{narrowGenResult.Queries.Length / narrowExtractTime.TotalSeconds:N0} q/sec)");

            // Score Pass 2 with apex intensity
            sw.Restart();
            var narrowResults = DiaLibraryQueryGenerator.AssembleResults(
                precursorsWithWindow, narrowGenResult, narrowExtraction.Results, narrowParams,
                dotScorer, saScorer);
            var narrowScoreTime = sw.Elapsed;

            Console.WriteLine($"  Scoring:           {narrowScoreTime.TotalMilliseconds:F0} ms  " +
                              $"({narrowResults.Count:N0} results)");
            Console.WriteLine();

            // ─── Step 7: Evaluate Pass 2 against ground truth ───────────────
            EvaluateResults(narrowResults, precursorsWithWindow, narrowGenResult);

            // ─── Step 8: Compare Pass 1 vs Pass 2 ──────────────────────────
            Console.WriteLine("=== Pass 1 vs Pass 2 Comparison ===");

            var narrowDotScores = narrowResults
                .Where(r => !float.IsNaN(r.DotProductScore))
                .Select(r => r.DotProductScore).ToArray();

            Console.WriteLine($"  {"Metric",-30} {"Pass 1 (broad)",-20} {"Pass 2 (calibrated)",-20}");
            Console.WriteLine($"  {"------",-30} {"------",-20} {"------",-20}");
            Console.WriteLine($"  {"RT window",-30} {"±" + rtToleranceMinutes.ToString("F1") + " min",-20} " +
                              $"{"±" + (calibration?.GetMinutesWindowHalfWidth(3.0).ToString("F2") ?? "N/A") + " min",-20}");
            Console.WriteLine($"  {"Extraction time",-30} {broadExtractTime.TotalMilliseconds.ToString("F0") + " ms",-20} " +
                              $"{narrowExtractTime.TotalMilliseconds.ToString("F0") + " ms",-20}");
            Console.WriteLine($"  {"Data points",-30} {broadExtraction.TotalDataPoints.ToString("N0"),-20} " +
                              $"{narrowExtraction.TotalDataPoints.ToString("N0"),-20}");
            Console.WriteLine($"  {"Results",-30} {broadResults.Count.ToString("N0"),-20} " +
                              $"{narrowResults.Count.ToString("N0"),-20}");

            if (broadDotScores.Length > 0 && narrowDotScores.Length > 0)
            {
                Console.WriteLine($"  {"Median DP (summed)",-30} {Median(broadDotScores).ToString("F4"),-20} " +
                                  $"{Median(narrowDotScores).ToString("F4"),-20}");
                Console.WriteLine($"  {"Mean DP (summed)",-30} {broadDotScores.Average().ToString("F4"),-20} " +
                                  $"{narrowDotScores.Average().ToString("F4"),-20}");
                Console.WriteLine($"  {"Q25 DP",-30} {Percentile(broadDotScores, 25).ToString("F4"),-20} " +
                                  $"{Percentile(narrowDotScores, 25).ToString("F4"),-20}");
                Console.WriteLine($"  {"Q75 DP",-30} {Percentile(broadDotScores, 75).ToString("F4"),-20} " +
                                  $"{Percentile(narrowDotScores, 75).ToString("F4"),-20}");
                Console.WriteLine($"  {"DP > 0.7",-30} " +
                                  $"{broadDotScores.Count(s => s > 0.7f).ToString("N0") + " (" + (100.0 * broadDotScores.Count(s => s > 0.7f) / broadDotScores.Length).ToString("F1") + "%)",-20} " +
                                  $"{narrowDotScores.Count(s => s > 0.7f).ToString("N0") + " (" + (100.0 * narrowDotScores.Count(s => s > 0.7f) / narrowDotScores.Length).ToString("F1") + "%)",-20}");
                Console.WriteLine($"  {"DP > 0.8",-30} " +
                                  $"{broadDotScores.Count(s => s > 0.8f).ToString("N0") + " (" + (100.0 * broadDotScores.Count(s => s > 0.8f) / broadDotScores.Length).ToString("F1") + "%)",-20} " +
                                  $"{narrowDotScores.Count(s => s > 0.8f).ToString("N0") + " (" + (100.0 * narrowDotScores.Count(s => s > 0.8f) / narrowDotScores.Length).ToString("F1") + "%)",-20}");
            }
            Console.WriteLine();

            // ─── Summary ────────────────────────────────────────────────────
            var totalTime = loadTime + buildTime + broadGenTime + broadExtractTime + broadScoreTime
                          + narrowGenTime + narrowExtractTime + narrowScoreTime;
            Console.WriteLine("=== Timing Summary ===");
            Console.WriteLine($"  mzML load:         {loadTime.TotalSeconds:F2}s");
            Console.WriteLine($"  Index build:       {buildTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Pass 1 gen:        {broadGenTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Pass 1 extract:    {broadExtractTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Pass 1 score:      {broadScoreTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Calibration fit:   <1 ms");
            Console.WriteLine($"  Pass 2 gen:        {narrowGenTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Pass 2 extract:    {narrowExtractTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Pass 2 score:      {narrowScoreTime.TotalMilliseconds:F0} ms");
            Console.WriteLine($"  Total pipeline:    {totalTime.TotalSeconds:F2}s");
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
                Console.WriteLine($"Dot Product Score (summed intensity, library vs extracted):");
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

            // Metric 5: Top hits
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

            // Metric 6: Bottom hits
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