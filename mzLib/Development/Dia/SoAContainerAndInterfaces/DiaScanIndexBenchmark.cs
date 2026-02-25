// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using System;
using System.Diagnostics;
using System.Linq;

namespace Development.Dia
{
    /// <summary>
    /// Benchmarks for DiaScanIndex construction and access patterns.
    /// 
    /// This is not a unit test — it generates realistic-scale synthetic DIA data
    /// and measures construction time, memory footprint, and scan access speed.
    /// Run manually via the Development project to evaluate performance.
    /// 
    /// Typical DIA file characteristics being simulated:
    ///   - 20 isolation windows
    ///   - ~50,000 MS2 scans total (~2,500 per window)
    ///   - ~300 peaks per scan
    ///   - ~15 million total peaks
    ///   - Equivalent to a ~3 GB mzML file
    /// </summary>
    public static class DiaScanIndexBenchmark
    {
        /// <summary>
        /// Runs the benchmark suite and prints results to console.
        /// Call from a Development project Main() or test runner.
        /// </summary>
        public static void RunAll()
        {
            Console.WriteLine("=== DiaScanIndex Benchmark ===");
            Console.WriteLine();

            BenchmarkBuild(windowCount: 20, scansPerWindow: 2500, peaksPerScan: 300);
            BenchmarkBuild(windowCount: 40, scansPerWindow: 1500, peaksPerScan: 200);

            BenchmarkSequentialScanAccess(windowCount: 20, scansPerWindow: 2500, peaksPerScan: 300);
            BenchmarkWindowLookup(windowCount: 20, scansPerWindow: 2500, peaksPerScan: 300);

            BenchmarkFragmentExtraction(windowCount: 20, scansPerWindow: 2500, peaksPerScan: 300);
        }

        /// <summary>
        /// Measures time to build a DiaScanIndex from synthetic scans.
        /// This benchmarks the ingestion path: filtering, window discovery, sorting, 
        /// and contiguous array packing.
        /// 
        /// Target: linear in total peak count. For 15M peaks, aim for under 2 seconds.
        /// </summary>
        private static void BenchmarkBuild(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            int totalScans = windowCount * scansPerWindow;
            long totalPeaks = (long)totalScans * peaksPerScan;

            Console.WriteLine($"Build benchmark: {windowCount} windows × {scansPerWindow} scans × {peaksPerScan} peaks");
            Console.WriteLine($"  Total: {totalScans:N0} scans, {totalPeaks:N0} peaks");

            var scans = GenerateSyntheticScans(windowCount, scansPerWindow, peaksPerScan);

            // Warm up (small run to JIT)
            var warmupScans = GenerateSyntheticScans(2, 10, 50);
            using (var _ = DiaScanIndexBuilder.Build(warmupScans)) { }

            // Force GC before measurement
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();

            long memBefore = GC.GetTotalMemory(true);
            var sw = Stopwatch.StartNew();

            using var index = DiaScanIndexBuilder.Build(scans);

            sw.Stop();
            long memAfter = GC.GetTotalMemory(false);

            double mbAllocated = (memAfter - memBefore) / (1024.0 * 1024.0);
            double expectedMb = totalPeaks * 2 * sizeof(float) / (1024.0 * 1024.0); // m/z + intensity

            Console.WriteLine($"  Build time:    {sw.ElapsedMilliseconds} ms");
            Console.WriteLine($"  Memory delta:  {mbAllocated:F1} MB (expected ~{expectedMb:F1} MB for peak data)");
            Console.WriteLine($"  Throughput:    {totalPeaks / sw.Elapsed.TotalSeconds:N0} peaks/sec");
            Console.WriteLine($"  Scan count:    {index.ScanCount}");
            Console.WriteLine($"  Window count:  {index.WindowCount}");
            Console.WriteLine();
        }

        /// <summary>
        /// Measures throughput of sequential scan access (iterating all scans and reading 
        /// their m/z spans). This benchmarks the SoA memory layout benefits.
        /// 
        /// Target: cache-friendly sequential access should saturate memory bandwidth.
        /// </summary>
        private static void BenchmarkSequentialScanAccess(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            Console.WriteLine("Sequential scan access benchmark:");
            var scans = GenerateSyntheticScans(windowCount, scansPerWindow, peaksPerScan);
            using var index = DiaScanIndexBuilder.Build(scans);

            var sw = Stopwatch.StartNew();
            float checksum = 0;
            int iterations = 10;

            for (int iter = 0; iter < iterations; iter++)
            {
                for (int i = 0; i < index.ScanCount; i++)
                {
                    var mz = index.GetScanMzSpan(i);
                    // Touch every element to prevent dead-code elimination
                    if (mz.Length > 0) checksum += mz[0] + mz[mz.Length - 1];
                }
            }

            sw.Stop();
            double scansPerSec = (double)index.ScanCount * iterations / sw.Elapsed.TotalSeconds;
            Console.WriteLine($"  {scansPerSec:N0} scan accesses/sec ({iterations} iterations)");
            Console.WriteLine($"  Checksum (prevent optimization): {checksum}");
            Console.WriteLine();
        }

        /// <summary>
        /// Measures window lookup performance — how fast we can get the scan range for a window.
        /// This is the entry point for window-level parallel processing.
        /// </summary>
        private static void BenchmarkWindowLookup(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            Console.WriteLine("Window lookup benchmark:");
            var scans = GenerateSyntheticScans(windowCount, scansPerWindow, peaksPerScan);
            using var index = DiaScanIndexBuilder.Build(scans);

            var sw = Stopwatch.StartNew();
            int iterations = 100_000;
            int totalScans = 0;

            for (int iter = 0; iter < iterations; iter++)
            {
                for (int w = 0; w < index.WindowCount; w++)
                {
                    if (index.TryGetScanRangeForWindow(w, out _, out int count))
                        totalScans += count;
                }
            }

            sw.Stop();
            double lookupsPerSec = (double)index.WindowCount * iterations / sw.Elapsed.TotalSeconds;
            Console.WriteLine($"  {lookupsPerSec:N0} window lookups/sec");
            Console.WriteLine($"  Total scans found: {totalScans}");
            Console.WriteLine();
        }

        /// <summary>
        /// Measures fragment extraction throughput with realistic targeted queries.
        /// 
        /// In real DIA analysis, fragment queries are not random — they target specific
        /// m/z values that are known to exist (from a spectral library or predicted fragments),
        /// within RT windows centered on predicted elution times. This benchmark simulates
        /// that by:
        ///   1. Building the index from synthetic data
        ///   2. Sampling actual m/z values from actual scans in each window
        ///   3. Centering RT windows on actual scan retention times
        /// 
        /// This ensures a high hit rate (~100%) and realistic data point counts per query
        /// (typically 20–50 XIC points), matching what a real DIA search engine would see.
        /// </summary>
        private static void BenchmarkFragmentExtraction(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            Console.WriteLine("=== Fragment Extraction Benchmarks ===");
            Console.WriteLine();

            var scans = GenerateSyntheticScans(windowCount, scansPerWindow, peaksPerScan);
            using var index = DiaScanIndexBuilder.Build(scans);

            // ── Benchmark 1: Targeted queries (realistic DIA search) ────────
            // Simulates a search with ~6 fragments per precursor, 10K precursors = 60K queries
            // Each query targets an m/z value actually present in a scan, with a ±1 min RT window
            {
                int precursorCount = 10_000;
                int fragmentsPerPrecursor = 6;
                int queryCount = precursorCount * fragmentsPerPrecursor;

                Console.WriteLine($"Targeted extraction: {precursorCount:N0} precursors × {fragmentsPerPrecursor} fragments = {queryCount:N0} queries");

                var queries = GenerateTargetedQueries(index, queryCount, rtHalfWidth: 1.0f, tolerancePpm: 10f, seed: 123);
                RunExtractionBenchmark(index, queries, "  ");
            }

            // ── Benchmark 2: Narrow RT window (well-predicted RT) ───────────
            // ±0.5 min window — fewer scans to check, should be faster per query
            {
                int queryCount = 60_000;
                Console.WriteLine($"Narrow RT window (±0.5 min): {queryCount:N0} queries");
                var queries = GenerateTargetedQueries(index, queryCount, rtHalfWidth: 0.5f, tolerancePpm: 10f, seed: 456);
                RunExtractionBenchmark(index, queries, "  ");
            }

            // ── Benchmark 3: Wide RT window (poor RT prediction) ────────────
            // ±3 min window — more scans per query, stress-tests the RT loop
            {
                int queryCount = 60_000;
                Console.WriteLine($"Wide RT window (±3.0 min): {queryCount:N0} queries");
                var queries = GenerateTargetedQueries(index, queryCount, rtHalfWidth: 3.0f, tolerancePpm: 10f, seed: 789);
                RunExtractionBenchmark(index, queries, "  ");
            }

            // ── Benchmark 4: Scaling test ───────────────────────────────────
            // Run 10K, 50K, 100K, 200K queries to check linearity
            {
                Console.WriteLine("Scaling test (±1.0 min RT window):");
                foreach (int qCount in new[] { 10_000, 50_000, 100_000, 200_000 })
                {
                    var queries = GenerateTargetedQueries(index, qCount, rtHalfWidth: 1.0f, tolerancePpm: 10f, seed: 42);
                    Console.Write($"  {qCount,8:N0} queries → ");
                    RunExtractionBenchmark(index, queries, "", compact: true);
                }
                Console.WriteLine();
            }
        }

        /// <summary>
        /// Generates targeted fragment queries by sampling actual m/z values from the index.
        /// 
        /// For each query:
        ///   - Picks a random window
        ///   - Picks a random scan within that window
        ///   - Picks a random peak m/z from that scan as the target
        ///   - Centers the RT window on that scan's RT
        /// 
        /// This guarantees that every query will find at least one matching peak,
        /// producing realistic hit rates and XIC data point counts.
        /// </summary>
        private static FragmentQuery[] GenerateTargetedQueries(
            DiaScanIndex index, int queryCount, float rtHalfWidth, float tolerancePpm, int seed)
        {
            var rng = new Random(seed);
            var queries = new FragmentQuery[queryCount];
            var windowIds = index.GetWindowIds().ToArray();

            for (int i = 0; i < queryCount; i++)
            {
                // Pick a random window
                int windowId = windowIds[rng.Next(windowIds.Length)];
                index.TryGetScanRangeForWindow(windowId, out int scanStart, out int scanCount);

                // Pick a random scan in this window
                int scanIndex = scanStart + rng.Next(scanCount);

                // Pick a random peak from this scan as the target m/z
                var mzSpan = index.GetScanMzSpan(scanIndex);
                float targetMz = mzSpan[rng.Next(mzSpan.Length)];

                // Center RT window on this scan's RT
                float rt = index.GetScanRt(scanIndex);

                queries[i] = new FragmentQuery(
                    targetMz: targetMz,
                    tolerancePpm: tolerancePpm,
                    rtMin: rt - rtHalfWidth,
                    rtMax: rt + rtHalfWidth,
                    windowId: windowId,
                    queryId: i);
            }

            return queries;
        }

        /// <summary>
        /// Runs the extraction benchmark for a set of queries and prints results.
        /// Handles buffer allocation, warmup, GC, and timing.
        /// </summary>
        private static void RunExtractionBenchmark(
            DiaScanIndex index, FragmentQuery[] queries, string indent, bool compact = false)
        {
            int queryCount = queries.Length;

            using var extractor = new CpuFragmentExtractor(index);

            var results = new FragmentResult[queryCount];
            // Buffer sized for worst case: each query could match many scans in the RT window
            // With 2500 scans/window over 100 min, ±1 min = ~50 scans, ±3 min = ~150 scans
            int bufferSize = queryCount * 200;
            var rtBuf = new float[bufferSize];
            var intBuf = new float[bufferSize];

            // Warm up JIT
            var warmupQ = new FragmentQuery[] { queries[0] };
            var warmupR = new FragmentResult[1];
            extractor.ExtractBatch(warmupQ, warmupR, new float[500], new float[500]);

            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();

            var sw = Stopwatch.StartNew();
            int totalDataPoints = extractor.ExtractBatch(queries, results, rtBuf, intBuf);
            sw.Stop();

            double queriesPerSec = queryCount / sw.Elapsed.TotalSeconds;
            int queriesWithData = 0;
            for (int i = 0; i < results.Length; i++)
            {
                if (results[i].DataPointCount > 0) queriesWithData++;
            }
            double avgPointsPerQuery = queriesWithData > 0 ? (double)totalDataPoints / queriesWithData : 0;

            if (compact)
            {
                Console.WriteLine($"{sw.ElapsedMilliseconds,6} ms | {queriesPerSec,12:N0} q/sec | {totalDataPoints,10:N0} pts | {avgPointsPerQuery,5:F1} pts/q | {queriesWithData * 100.0 / queryCount,5:F1}% hit");
            }
            else
            {
                Console.WriteLine($"{indent}Time:           {sw.ElapsedMilliseconds} ms");
                Console.WriteLine($"{indent}Queries/sec:    {queriesPerSec:N0}");
                Console.WriteLine($"{indent}Data points:    {totalDataPoints:N0}");
                Console.WriteLine($"{indent}Queries w/data: {queriesWithData:N0} / {queryCount:N0} ({queriesWithData * 100.0 / queryCount:F1}% hit rate)");
                Console.WriteLine($"{indent}Avg points/q:   {avgPointsPerQuery:F1}");
                Console.WriteLine();
            }
        }

        /// <summary>
        /// Generates a realistic array of synthetic DIA scans for benchmarking.
        /// 
        /// Simulates a DIA acquisition with:
        ///   - Equally-spaced isolation windows from 400–1200 m/z
        ///   - Each window has a set of "resident" fragment m/z values that appear in
        ///     every scan (simulating real peptide fragments that persist across cycles)
        ///   - Additional random noise peaks fill out each scan to the target peak count
        ///   - Retention times increase linearly (~0.04 min per cycle ≈ 100 min run)
        /// 
        /// This structure ensures that targeted queries (which sample real m/z values)
        /// will find matches in every scan within the RT window, producing realistic
        /// XIC data point counts (25–50 per query with ±1 min RT windows).
        /// </summary>
        private static MsDataScan[] GenerateSyntheticScans(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            var rng = new Random(42);
            int totalMs2 = windowCount * scansPerWindow;
            var scans = new MsDataScan[totalMs2];

            double windowWidth = 25.0;
            double windowStart = 400.0;
            double windowSpacing = (1200.0 - windowStart) / windowCount;

            // Number of persistent fragment m/z values per window.
            // In real DIA, a 25 m/z window might contain 50–200 precursors, each producing
            // 6–10 fragments. We simulate ~100 persistent fragments per window.
            int residentFragmentsPerWindow = Math.Min(100, peaksPerScan / 2);
            int noisePeaksPerScan = peaksPerScan - residentFragmentsPerWindow;

            // Pre-generate resident fragment m/z values for each window.
            // These are the fragments that will appear in every scan of this window.
            double[][] residentMzPerWindow = new double[windowCount][];
            for (int w = 0; w < windowCount; w++)
            {
                residentMzPerWindow[w] = new double[residentFragmentsPerWindow];
                for (int f = 0; f < residentFragmentsPerWindow; f++)
                {
                    // Fragment m/z values spread across typical product ion range (100–1800)
                    residentMzPerWindow[w][f] = 100.0 + rng.NextDouble() * 1700.0;
                }
                Array.Sort(residentMzPerWindow[w]);
            }

            int scanIdx = 0;
            for (int cycle = 0; cycle < scansPerWindow; cycle++)
            {
                double cycleRt = cycle * 0.04;

                for (int w = 0; w < windowCount; w++)
                {
                    double isolationCenter = windowStart + w * windowSpacing + windowSpacing / 2.0;
                    double rt = cycleRt + w * 0.001;

                    // Combine resident fragments + noise peaks
                    double[] mzValues = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];

                    // Copy resident fragments with slight m/z jitter (±0.001 Da, simulating
                    // instrument measurement noise) and varying intensity
                    double[] residents = residentMzPerWindow[w];
                    for (int f = 0; f < residentFragmentsPerWindow; f++)
                    {
                        mzValues[f] = residents[f] + (rng.NextDouble() - 0.5) * 0.002;
                        intensities[f] = 500.0 + rng.NextDouble() * 9500.0;
                    }

                    // Fill remaining slots with random noise peaks
                    for (int p = residentFragmentsPerWindow; p < peaksPerScan; p++)
                    {
                        mzValues[p] = 100.0 + rng.NextDouble() * 1900.0;
                        intensities[p] = rng.NextDouble() * 500.0; // lower intensity noise
                    }

                    Array.Sort(mzValues, intensities);

                    scans[scanIdx] = new MsDataScan(
                        massSpectrum: new MzSpectrum(mzValues, intensities, false),
                        oneBasedScanNumber: scanIdx + 1,
                        msnOrder: 2,
                        isCentroid: true,
                        polarity: Polarity.Positive,
                        retentionTime: rt,
                        scanWindowRange: new MzRange(100, 2000),
                        scanFilter: "FTMS",
                        mzAnalyzer: MZAnalyzerType.Orbitrap,
                        totalIonCurrent: intensities.Sum(),
                        injectionTime: 20.0,
                        noiseData: null,
                        nativeId: $"scan={scanIdx + 1}",
                        isolationMZ: isolationCenter,
                        isolationWidth: windowWidth,
                        dissociationType: DissociationType.HCD);

                    scanIdx++;
                }
            }

            return scans;
        }
    }
}