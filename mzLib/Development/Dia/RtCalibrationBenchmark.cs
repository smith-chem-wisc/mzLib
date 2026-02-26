// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;

namespace Development.Dia
{
    /// <summary>
    /// Benchmarks for iRT calibration components (Phase 8).
    /// 
    /// Measures:
    ///   1. IrtLibraryIndex construction and query throughput
    ///   2. RtCalibrationFitter performance (OLS, RANSAC, iterative refinement)
    ///   3. Calibrated query generation vs fixed-window generation
    ///   4. End-to-end: anchor extraction → calibration fit → re-generation
    /// 
    /// Simulates realistic proteomics scales:
    ///   - 100K+ precursors with iRT values
    ///   - 200–2000 anchor matches for calibration
    ///   - True linear model: observed_RT = 0.5 * iRT + 5.0
    /// </summary>
    public static class RtCalibrationBenchmark
    {
        public static void RunAll()
        {
            Console.WriteLine("========================================");
            Console.WriteLine("  iRT Calibration Benchmarks (Phase 8)");
            Console.WriteLine("========================================\n");

            BenchmarkIrtLibraryIndex();
            BenchmarkIrtLibraryIndexQuery();
            BenchmarkCalibrationFitter();
            BenchmarkCalibratedQueryGeneration();
            BenchmarkEndToEndCalibrationLoop();

            Console.WriteLine("Done.");
        }

        // ── Benchmark 1: IrtLibraryIndex construction ───────────────────────

        public static void BenchmarkIrtLibraryIndex()
        {
            Console.WriteLine("=== IrtLibraryIndex Construction ===");

            foreach (int count in new[] { 10_000, 50_000, 100_000, 500_000 })
            {
                var precursors = GeneratePrecursorsWithIrt(count, seed: 42);

                // Warmup
                _ = new IrtLibraryIndex(precursors);

                GC.Collect();
                GC.WaitForPendingFinalizers();

                var sw = Stopwatch.StartNew();
                const int iterations = 10;
                IrtLibraryIndex index = null;
                for (int i = 0; i < iterations; i++)
                    index = new IrtLibraryIndex(precursors);
                sw.Stop();

                double avgMs = sw.Elapsed.TotalMilliseconds / iterations;
                Console.WriteLine(
                    $"  {count,8:N0} precursors:  {avgMs,8:F2} ms  |  " +
                    $"{count / (avgMs / 1000.0),12:N0} entries/sec  |  " +
                    $"iRT range [{index.MinIrt:F1}, {index.MaxIrt:F1}]");
            }
            Console.WriteLine();
        }

        // ── Benchmark 2: IrtLibraryIndex query throughput ───────────────────

        public static void BenchmarkIrtLibraryIndexQuery()
        {
            Console.WriteLine("=== IrtLibraryIndex Query Throughput ===");

            int precursorCount = 100_000;
            var precursors = GeneratePrecursorsWithIrt(precursorCount, seed: 42);
            var index = new IrtLibraryIndex(precursors);

            var rng = new Random(99);
            int queryCount = 1_000_000;

            // Generate random query windows of varying sizes
            double[] queryCenters = new double[queryCount];
            for (int i = 0; i < queryCount; i++)
                queryCenters[i] = rng.NextDouble() * 120.0;

            foreach (double halfWidth in new[] { 2.0, 5.0, 10.0, 20.0 })
            {
                // Warmup
                for (int i = 0; i < 1000; i++)
                    index.QueryRange(queryCenters[i] - halfWidth, queryCenters[i] + halfWidth,
                        out _, out _);

                GC.Collect();
                var sw = Stopwatch.StartNew();
                long totalCandidates = 0;
                for (int i = 0; i < queryCount; i++)
                {
                    index.QueryRange(
                        queryCenters[i] - halfWidth,
                        queryCenters[i] + halfWidth,
                        out _, out int count);
                    totalCandidates += count;
                }
                sw.Stop();

                double avgCandidates = (double)totalCandidates / queryCount;
                Console.WriteLine(
                    $"  ±{halfWidth,4:F1} iRT:  {sw.ElapsedMilliseconds,6} ms  |  " +
                    $"{queryCount / sw.Elapsed.TotalSeconds,12:N0} queries/sec  |  " +
                    $"avg {avgCandidates:F1} candidates/query");
            }
            Console.WriteLine();
        }

        // ── Benchmark 3: RtCalibrationFitter ────────────────────────────────

        public static void BenchmarkCalibrationFitter()
        {
            Console.WriteLine("=== RtCalibrationFitter Performance ===");

            // True model: RT = 0.5 * iRT + 5.0
            double trueSlope = 0.5, trueIntercept = 5.0;

            foreach (int anchorCount in new[] { 50, 200, 500, 1000, 2000 })
            {
                var rng = new Random(42);
                var xs = new double[anchorCount];
                var ys = new double[anchorCount];

                for (int i = 0; i < anchorCount; i++)
                {
                    xs[i] = rng.NextDouble() * 120.0;
                    ys[i] = trueSlope * xs[i] + trueIntercept + NextGaussian(rng) * 0.3;
                }

                // Add 10% outliers
                int outlierCount = anchorCount / 10;
                for (int i = 0; i < outlierCount; i++)
                {
                    int idx = rng.Next(anchorCount);
                    ys[idx] += (rng.NextDouble() > 0.5 ? 1 : -1) * 15.0;
                }

                // Warmup
                RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys);

                GC.Collect();
                const int iterations = 100;
                var sw = Stopwatch.StartNew();
                RtCalibrationModel model = null;
                for (int i = 0; i < iterations; i++)
                    model = RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys);
                sw.Stop();

                double avgMs = sw.Elapsed.TotalMilliseconds / iterations;
                string reliability = model != null && model.IsReliable ? "RELIABLE" : "UNRELIABLE";
                string slopeStr = model != null ? $"slope={model.Slope:F4}" : "null";
                string anchorStr = model != null ? $"anchors={model.AnchorCount}" : "";

                Console.WriteLine(
                    $"  {anchorCount,5} anchors (10% outlier):  {avgMs,8:F3} ms  |  " +
                    $"{slopeStr}  {anchorStr}  [{reliability}]");
            }

            // Compare RANSAC vs OLS-only
            Console.WriteLine();
            Console.WriteLine("  RANSAC vs OLS-only (500 anchors, 10% outlier):");
            {
                var rng = new Random(42);
                int n = 500;
                var xs = new double[n];
                var ys = new double[n];
                for (int i = 0; i < n; i++)
                {
                    xs[i] = rng.NextDouble() * 120.0;
                    ys[i] = trueSlope * xs[i] + trueIntercept + NextGaussian(rng) * 0.3;
                }
                for (int i = 0; i < 50; i++)
                    ys[rng.Next(n)] += (rng.NextDouble() > 0.5 ? 1 : -1) * 15.0;

                var optRansac = new RtCalibrationFitter.FitOptions { UseRansac = true };
                var optOls = new RtCalibrationFitter.FitOptions { UseRansac = false };

                // Warmup
                RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, optRansac);
                RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, optOls);

                const int iterations = 200;

                var sw1 = Stopwatch.StartNew();
                for (int i = 0; i < iterations; i++)
                    RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, optRansac);
                sw1.Stop();

                var sw2 = Stopwatch.StartNew();
                for (int i = 0; i < iterations; i++)
                    RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, optOls);
                sw2.Stop();

                var modelR = RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, optRansac);
                var modelO = RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, optOls);

                Console.WriteLine(
                    $"    RANSAC:   {sw1.Elapsed.TotalMilliseconds / iterations:F3} ms  |  " +
                    $"slope={modelR?.Slope:F4}  σ={modelR?.SigmaMinutes:F4}  R²={modelR?.RSquared:F4}");
                Console.WriteLine(
                    $"    OLS-only: {sw2.Elapsed.TotalMilliseconds / iterations:F3} ms  |  " +
                    $"slope={modelO?.Slope:F4}  σ={modelO?.SigmaMinutes:F4}  R²={modelO?.RSquared:F4}");
            }
            Console.WriteLine();
        }

        // ── Benchmark 4: Calibrated vs fixed query generation ───────────────

        public static void BenchmarkCalibratedQueryGeneration()
        {
            Console.WriteLine("=== Calibrated vs Fixed Query Generation ===");

            var scanIndex = BuildBenchmarkScanIndex(windowCount: 30, cyclesPerWindow: 200, peaksPerScan: 200);

            foreach (int precursorCount in new[] { 10_000, 50_000, 100_000 })
            {
                var precursors = GeneratePrecursorsWithIrt(precursorCount, seed: 42);
                var parameters = new DiaSearchParameters
                {
                    PpmTolerance = 20f,
                    RtToleranceMinutes = 5f,
                    CalibratedWindowSigmaMultiplier = 3.0
                };

                // Build a calibration model
                var calibration = new RtCalibrationModel(
                    slope: 0.5, intercept: 5.0, sigmaMinutes: 0.3,
                    rSquared: 0.98, anchorCount: 200);

                // Warmup
                DiaLibraryQueryGenerator.Generate(precursors, scanIndex, parameters);
                DiaLibraryQueryGenerator.GenerateCalibrated(precursors, scanIndex, parameters, calibration);

                GC.Collect();
                const int iterations = 5;

                // Fixed window
                var sw1 = Stopwatch.StartNew();
                DiaLibraryQueryGenerator.GenerationResult fixedResult = default;
                for (int i = 0; i < iterations; i++)
                    fixedResult = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, parameters);
                sw1.Stop();

                // Calibrated window
                var sw2 = Stopwatch.StartNew();
                DiaLibraryQueryGenerator.GenerationResult calibResult = default;
                for (int i = 0; i < iterations; i++)
                    calibResult = DiaLibraryQueryGenerator.GenerateCalibrated(precursors, scanIndex, parameters, calibration);
                sw2.Stop();

                double fixedMs = sw1.Elapsed.TotalMilliseconds / iterations;
                double calibMs = sw2.Elapsed.TotalMilliseconds / iterations;

                Console.WriteLine($"  {precursorCount,8:N0} precursors:");
                Console.WriteLine(
                    $"    Fixed ±5min:     {fixedMs,8:F2} ms  |  " +
                    $"{fixedResult.Queries.Length / (fixedMs / 1000.0),12:N0} queries/sec  |  " +
                    $"{fixedResult.Queries.Length:N0} queries");
                Console.WriteLine(
                    $"    Calibrated 3σ:   {calibMs,8:F2} ms  |  " +
                    $"{calibResult.Queries.Length / (calibMs / 1000.0),12:N0} queries/sec  |  " +
                    $"{calibResult.Queries.Length:N0} queries");
            }

            scanIndex.Dispose();
            Console.WriteLine();
        }

        // ── Benchmark 5: End-to-end calibration loop ────────────────────────

        public static void BenchmarkEndToEndCalibrationLoop()
        {
            Console.WriteLine("=== End-to-End Calibration Loop ===");
            Console.WriteLine("  Simulates: broad search → anchor extraction → fit → refined generation");
            Console.WriteLine();

            var scanIndex = BuildBenchmarkScanIndex(windowCount: 30, cyclesPerWindow: 200, peaksPerScan: 200);
            int precursorCount = 50_000;
            var precursors = GeneratePrecursorsWithIrt(precursorCount, seed: 42);

            // True model for generating synthetic match results
            double trueSlope = 0.5, trueIntercept = 5.0;

            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5f
            };

            // iRT calibration settings (used directly, not via DiaSearchParameters,
            // so the benchmark builds even if DiaSearchParameters hasn't been updated yet)
            double initialIrtWindow = 20.0;
            double calibratedWindowK = 3.0;

            var totalSw = Stopwatch.StartNew();

            // Step 1: Build IrtLibraryIndex
            var sw = Stopwatch.StartNew();
            var irtIndex = new IrtLibraryIndex(precursors);
            sw.Stop();
            Console.WriteLine($"  IrtLibraryIndex build:      {sw.ElapsedMilliseconds,6} ms  ({irtIndex.Count:N0} entries)");

            // Step 2: Create provisional calibration
            sw.Restart();
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();
            var provisionalCalibration = RtCalibrationModel.CreateProvisional(
                globalRtMin, globalRtMax, irtIndex.MinIrt, irtIndex.MaxIrt,
                initialIrtWindow);
            sw.Stop();
            Console.WriteLine($"  Provisional calibration:    {sw.ElapsedMilliseconds,6} ms  " +
                              $"(slope={provisionalCalibration.Slope:F4})");

            // Step 3: Phase 1 broad query generation
            sw.Restart();
            var broadResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                precursors, scanIndex, parameters, provisionalCalibration);
            sw.Stop();
            Console.WriteLine($"  Broad query generation:     {sw.ElapsedMilliseconds,6} ms  " +
                              $"({broadResult.Queries.Length:N0} queries)");

            // Step 4: Simulate anchor extraction (synthesize high-confidence matches)
            sw.Restart();
            int anchorCount = 500;
            var rng = new Random(42);
            var anchorIrts = new double[anchorCount];
            var anchorRts = new double[anchorCount];
            for (int i = 0; i < anchorCount; i++)
            {
                anchorIrts[i] = rng.NextDouble() * 120.0;
                anchorRts[i] = trueSlope * anchorIrts[i] + trueIntercept + NextGaussian(rng) * 0.3;
            }
            // Add 5% outliers
            for (int i = 0; i < anchorCount / 20; i++)
                anchorRts[rng.Next(anchorCount)] += (rng.NextDouble() > 0.5 ? 1 : -1) * 10.0;
            sw.Stop();
            Console.WriteLine($"  Anchor generation:          {sw.ElapsedMilliseconds,6} ms  ({anchorCount} anchors)");

            // Step 5: Calibration fitting
            sw.Restart();
            var calibration = RtCalibrationFitter.Fit(
                (ReadOnlySpan<double>)anchorIrts, (ReadOnlySpan<double>)anchorRts);
            sw.Stop();
            Console.WriteLine($"  Calibration fit:            {sw.ElapsedMilliseconds,6} ms  " +
                              $"(slope={calibration?.Slope:F4}, σ={calibration?.SigmaMinutes:F4} min, " +
                              $"R²={calibration?.RSquared:F4}, anchors={calibration?.AnchorCount})");

            // Step 6: Refined query generation with calibrated windows
            sw.Restart();
            var refinedResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                precursors, scanIndex, parameters, calibration);
            sw.Stop();
            Console.WriteLine($"  Refined query generation:   {sw.ElapsedMilliseconds,6} ms  " +
                              $"({refinedResult.Queries.Length:N0} queries)");

            // Step 7: Second iteration (simulate iterative refinement)
            sw.Restart();
            // Regenerate anchors with tighter residuals (simulating better matches)
            var anchorIrts2 = new double[anchorCount];
            var anchorRts2 = new double[anchorCount];
            for (int i = 0; i < anchorCount; i++)
            {
                anchorIrts2[i] = rng.NextDouble() * 120.0;
                anchorRts2[i] = trueSlope * anchorIrts2[i] + trueIntercept + NextGaussian(rng) * 0.15;
            }
            var calibration2 = RtCalibrationFitter.Fit(
                (ReadOnlySpan<double>)anchorIrts2, (ReadOnlySpan<double>)anchorRts2);
            var refinedResult2 = DiaLibraryQueryGenerator.GenerateCalibrated(
                precursors, scanIndex, parameters, calibration2);
            sw.Stop();
            Console.WriteLine($"  Iteration 2 (fit+gen):      {sw.ElapsedMilliseconds,6} ms  " +
                              $"(σ={calibration2?.SigmaMinutes:F4} min, " +
                              $"{refinedResult2.Queries.Length:N0} queries)");

            totalSw.Stop();
            Console.WriteLine();
            Console.WriteLine($"  Total end-to-end:           {totalSw.ElapsedMilliseconds,6} ms");

            // Summary comparison
            double broadWindowMin = provisionalCalibration.GetMinutesWindowHalfWidth(calibratedWindowK) * 2;
            double refinedWindowMin = calibration != null ? calibration.GetMinutesWindowHalfWidth(calibratedWindowK) * 2 : 0;
            Console.WriteLine();
            Console.WriteLine($"  Window comparison:");
            Console.WriteLine($"    Broad (provisional):   ±{broadWindowMin / 2:F2} min ({broadWindowMin:F2} min total)");
            Console.WriteLine($"    Refined (calibrated):  ±{refinedWindowMin / 2:F2} min ({refinedWindowMin:F2} min total)");
            Console.WriteLine($"    Reduction factor:      {broadWindowMin / Math.Max(refinedWindowMin, 0.001):F1}×");

            scanIndex.Dispose();
            Console.WriteLine();
        }

        // ── Helpers ─────────────────────────────────────────────────────────

        /// <summary>
        /// Generates precursors with iRT values uniformly distributed in [0, 120].
        /// Precursor m/z in [400, 1200] to match typical DIA windows.
        /// </summary>
        private static List<LibraryPrecursorInput> GeneratePrecursorsWithIrt(int count, int seed)
        {
            var rng = new Random(seed);
            var precursors = new List<LibraryPrecursorInput>(count);
            int fragmentsPerPrecursor = 6;

            for (int i = 0; i < count; i++)
            {
                double irt = rng.NextDouble() * 120.0;
                double precursorMz = 400.0 + rng.NextDouble() * 800.0;
                double rtMinutes = 0.5 * irt + 5.0; // True model for RT

                var fragMzs = new float[fragmentsPerPrecursor];
                var fragIntensities = new float[fragmentsPerPrecursor];
                for (int f = 0; f < fragmentsPerPrecursor; f++)
                {
                    fragMzs[f] = 100f + (float)(rng.NextDouble() * 1500);
                    fragIntensities[f] = (float)(rng.NextDouble() * 10000);
                }

                precursors.Add(new LibraryPrecursorInput(
                    sequence: $"PEPTIDE{i}",
                    precursorMz: precursorMz,
                    chargeState: rng.Next(2, 5),
                    retentionTime: rtMinutes,
                    isDecoy: false,
                    fragmentMzs: fragMzs,
                    fragmentIntensities: fragIntensities,
                    irtValue: irt));
            }

            return precursors;
        }

        /// <summary>
        /// Builds a realistic DiaScanIndex for benchmarking.
        /// </summary>
        private static DiaScanIndex BuildBenchmarkScanIndex(int windowCount, int cyclesPerWindow, int peaksPerScan)
        {
            var scans = new List<MsDataScan>();
            int scanNumber = 1;
            float startMz = 400f;
            float windowWidth = (1200f - 400f) / windowCount;

            for (int cycle = 0; cycle < cyclesPerWindow; cycle++)
            {
                for (int w = 0; w < windowCount; w++)
                {
                    double rt = cycle * 0.6 + w * 0.01;
                    double center = startMz + w * windowWidth + windowWidth / 2;

                    double[] mzs = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];
                    for (int p = 0; p < peaksPerScan; p++)
                    {
                        mzs[p] = 100.0 + p * 0.5;
                        intensities[p] = 1000.0;
                    }

                    var spectrum = new MzSpectrum(mzs, intensities, false);
                    scans.Add(new MsDataScan(
                        massSpectrum: spectrum,
                        oneBasedScanNumber: scanNumber++,
                        msnOrder: 2,
                        isCentroid: true,
                        polarity: Polarity.Positive,
                        retentionTime: rt,
                        scanWindowRange: new MzRange(100, 1500),
                        scanFilter: "",
                        mzAnalyzer: MZAnalyzerType.Orbitrap,
                        totalIonCurrent: intensities.Sum(),
                        injectionTime: 50,
                        noiseData: null,
                        nativeId: $"scan={scanNumber - 1}",
                        isolationMZ: center,
                        isolationWidth: windowWidth,
                        dissociationType: DissociationType.HCD));
                }
            }

            return DiaScanIndexBuilder.Build(scans.ToArray());
        }

        private static double NextGaussian(Random rng)
        {
            double u1 = 1.0 - rng.NextDouble();
            double u2 = rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
    }
}