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

            BenchmarkScoring();

            // Phase 16B, Prompt 3: MS1 SoA storage verification
            VerifyMs1Storage();
            BenchmarkMs1Build(windowCount: 20, scansPerWindow: 2500, ms1ScansPerWindow: 125, peaksPerMs2: 300, peaksPerMs1: 500);
            BenchmarkMs1RtLookup(windowCount: 20, scansPerWindow: 2500, ms1ScansPerWindow: 125, peaksPerMs1: 500);
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
        /// Benchmarks scoring throughput for both NormalizedDotProduct and SpectralAngle scorers.
        /// 
        /// Simulates scoring extracted fragment spectra against library/predicted spectra.
        /// Tests multiple vector lengths to show how SIMD acceleration scales:
        ///   - 6 fragments (typical for a single precursor, too short for SIMD benefit)
        ///   - 20 fragments (transition-level, starts to benefit from SIMD)
        ///   - 100 fragments (full spectral comparison, SIMD should dominate)
        /// 
        /// In a real DIA pipeline, scoring happens after extraction — once per precursor
        /// per candidate. With 50K–200K candidates, scoring must be very fast.
        /// </summary>
        private static void BenchmarkScoring()
        {
            Console.WriteLine("=== Scoring Benchmarks ===");
            Console.WriteLine();

            var dotProduct = new NormalizedDotProductScorer();
            var spectralAngle = new SpectralAngleScorer();

            foreach (int vectorLength in new[] { 6, 20, 100 })
            {
                Console.WriteLine($"Vector length: {vectorLength} fragments");

                // Generate random but consistent test vectors
                var rng = new Random(42);
                int pairCount = 200_000;
                float[][] observed = new float[pairCount][];
                float[][] expected = new float[pairCount][];

                for (int i = 0; i < pairCount; i++)
                {
                    observed[i] = new float[vectorLength];
                    expected[i] = new float[vectorLength];
                    for (int j = 0; j < vectorLength; j++)
                    {
                        // Observed: library-like intensities with noise
                        float baseIntensity = (float)(rng.NextDouble() * 1000.0);
                        observed[i][j] = baseIntensity + (float)(rng.NextDouble() * 200.0 - 100.0);
                        expected[i][j] = baseIntensity;
                        // Clamp to non-negative
                        if (observed[i][j] < 0) observed[i][j] = 0;
                    }
                }

                // Warm up
                for (int i = 0; i < 100; i++)
                {
                    dotProduct.Score(observed[i], expected[i]);
                    spectralAngle.Score(observed[i], expected[i]);
                }

                // Benchmark dot product
                GC.Collect();
                GC.WaitForPendingFinalizers();
                var sw = Stopwatch.StartNew();
                float checksum = 0f;
                for (int i = 0; i < pairCount; i++)
                {
                    checksum += dotProduct.Score(observed[i], expected[i]);
                }
                sw.Stop();
                double dpPerSec = pairCount / sw.Elapsed.TotalSeconds;
                Console.WriteLine($"  Dot product:    {sw.ElapsedMilliseconds,5} ms | {dpPerSec,14:N0} scores/sec | avg score: {checksum / pairCount:F4}");

                // Benchmark spectral angle
                GC.Collect();
                GC.WaitForPendingFinalizers();
                sw.Restart();
                checksum = 0f;
                for (int i = 0; i < pairCount; i++)
                {
                    checksum += spectralAngle.Score(observed[i], expected[i]);
                }
                sw.Stop();
                double saPerSec = pairCount / sw.Elapsed.TotalSeconds;
                Console.WriteLine($"  Spectral angle: {sw.ElapsedMilliseconds,5} ms | {saPerSec,14:N0} scores/sec | avg score: {checksum / pairCount:F4}");

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
            return GenerateSyntheticScansWithMs1(windowCount, scansPerWindow, peaksPerScan,
                ms1ScansPerWindow: 0, peaksPerMs1: 0);
        }

        /// <summary>
        /// Generates synthetic DIA scans including MS1 scans interleaved with MS2.
        /// 
        /// MS1 scans are generated at a rate of 1 per (windowCount / ms1ScansPerWindow)
        /// cycles, interleaved between MS2 cycles at their natural RT positions.
        /// Each MS1 scan covers the full precursor m/z range (300–1800) with
        /// simulated peptide isotope clusters.
        /// 
        /// Phase 16B, Prompt 3: ms1ScansPerWindow parameter added to enable MS1 testing.
        /// </summary>
        private static MsDataScan[] GenerateSyntheticScansWithMs1(
            int windowCount, int scansPerWindow, int peaksPerMs2,
            int ms1ScansPerWindow, int peaksPerMs1)
        {
            var rng = new Random(42);
            int totalMs2 = windowCount * scansPerWindow;

            double windowWidth = 25.0;
            double windowStart = 400.0;
            double windowSpacing = (1200.0 - windowStart) / windowCount;

            int residentFragmentsPerWindow = Math.Min(100, peaksPerMs2 / 2);
            int noisePeaksPerScan = peaksPerMs2 - residentFragmentsPerWindow;

            double[][] residentMzPerWindow = new double[windowCount][];
            for (int w = 0; w < windowCount; w++)
            {
                residentMzPerWindow[w] = new double[residentFragmentsPerWindow];
                for (int f = 0; f < residentFragmentsPerWindow; f++)
                    residentMzPerWindow[w][f] = 100.0 + rng.NextDouble() * 1700.0;
                Array.Sort(residentMzPerWindow[w]);
            }

            // Determine which cycles get an MS1 scan before them.
            // ms1ScansPerWindow total MS1 scans spread across scansPerWindow cycles.
            var ms1CycleSet = new HashSet<int>();
            if (ms1ScansPerWindow > 0 && scansPerWindow > 0)
            {
                int stride = Math.Max(1, scansPerWindow / ms1ScansPerWindow);
                for (int i = 0; i < ms1ScansPerWindow; i++)
                    ms1CycleSet.Add(i * stride);
            }

            // Pre-size: totalMs2 + ms1ScansPerWindow
            var scanList = new System.Collections.Generic.List<MsDataScan>(totalMs2 + ms1ScansPerWindow);
            int oneBasedScanNumber = 1;

            for (int cycle = 0; cycle < scansPerWindow; cycle++)
            {
                double cycleRt = cycle * 0.04;

                // Interleave MS1 before this cycle if it's an MS1 cycle
                if (ms1ScansPerWindow > 0 && ms1CycleSet.Contains(cycle))
                {
                    double ms1Rt = cycleRt - 0.001; // slightly before the MS2 cycle
                    double[] ms1Mz = new double[peaksPerMs1];
                    double[] ms1Int = new double[peaksPerMs1];

                    // Simulate precursor peaks spread across 300–1800 m/z
                    for (int p = 0; p < peaksPerMs1; p++)
                    {
                        ms1Mz[p] = 300.0 + rng.NextDouble() * 1500.0;
                        ms1Int[p] = 1000.0 + rng.NextDouble() * 99000.0;
                    }
                    Array.Sort(ms1Mz, ms1Int);

                    scanList.Add(new MsDataScan(
                        massSpectrum: new MzSpectrum(ms1Mz, ms1Int, false),
                        oneBasedScanNumber: oneBasedScanNumber++,
                        msnOrder: 1,
                        isCentroid: true,
                        polarity: Polarity.Positive,
                        retentionTime: ms1Rt,
                        scanWindowRange: new MzRange(300, 1800),
                        scanFilter: "FTMS",
                        mzAnalyzer: MZAnalyzerType.Orbitrap,
                        totalIonCurrent: ms1Int.Sum(),
                        injectionTime: 50.0,
                        noiseData: null,
                        nativeId: $"scan={oneBasedScanNumber - 1}"));
                }

                for (int w = 0; w < windowCount; w++)
                {
                    double isolationCenter = windowStart + w * windowSpacing + windowSpacing / 2.0;
                    double rt = cycleRt + w * 0.001;

                    double[] mzValues = new double[peaksPerMs2];
                    double[] intensities = new double[peaksPerMs2];

                    double[] residents = residentMzPerWindow[w];
                    for (int f = 0; f < residentFragmentsPerWindow; f++)
                    {
                        mzValues[f] = residents[f] + (rng.NextDouble() - 0.5) * 0.002;
                        intensities[f] = 500.0 + rng.NextDouble() * 9500.0;
                    }
                    for (int p = residentFragmentsPerWindow; p < peaksPerMs2; p++)
                    {
                        mzValues[p] = 100.0 + rng.NextDouble() * 1900.0;
                        intensities[p] = rng.NextDouble() * 500.0;
                    }
                    Array.Sort(mzValues, intensities);

                    scanList.Add(new MsDataScan(
                        massSpectrum: new MzSpectrum(mzValues, intensities, false),
                        oneBasedScanNumber: oneBasedScanNumber++,
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
                        nativeId: $"scan={oneBasedScanNumber - 1}",
                        isolationMZ: isolationCenter,
                        isolationWidth: windowWidth,
                        dissociationType: DissociationType.HCD));
                }
            }

            return scanList.ToArray();
        }

        // ── Phase 16B, Prompt 3: MS1 verification and benchmarks ────────────

        /// <summary>
        /// Correctness verification for MS1 SoA storage.
        /// 
        /// Tests:
        ///   1. Ms1ScanCount matches the number of MS1 scans in the input
        ///   2. GetMs1ScanRt returns correct RT values
        ///   3. GetMs1ScanPeaks returns correct peak counts
        ///   4. FindMs1ScanIndexAtRt returns the correct lower-bound index
        ///   5. Edge cases: query RT before all scans, after all scans, exact match
        ///   6. Zero MS1 scans (MS2-only file) builds correctly
        /// 
        /// All assertions use Debug.Assert so failures are immediate and explicit.
        /// </summary>
        private static void VerifyMs1Storage()
        {
            Console.WriteLine("=== MS1 Storage Verification ===");
            Console.WriteLine();

            // ── Test 1: Basic MS1 ingestion ─────────────────────────────────
            {
                int windowCount = 5, scansPerWindow = 100, ms1Count = 50, peaksPerMs1 = 200;
                var scans = GenerateSyntheticScansWithMs1(
                    windowCount, scansPerWindow, peaksPerMs2: 100,
                    ms1ScansPerWindow: ms1Count, peaksPerMs1: peaksPerMs1);

                int expectedMs1 = scans.Count(s => s.MsnOrder == 1);
                int expectedMs2 = scans.Count(s => s.MsnOrder == 2);

                using var index = DiaScanIndexBuilder.Build(scans);

                Debug.Assert(index.Ms1ScanCount == expectedMs1,
                    $"Ms1ScanCount mismatch: expected {expectedMs1}, got {index.Ms1ScanCount}");
                Debug.Assert(index.ScanCount == expectedMs2,
                    $"ScanCount (MS2) mismatch: expected {expectedMs2}, got {index.ScanCount}");

                Console.WriteLine($"  [PASS] Ms1ScanCount = {index.Ms1ScanCount} (expected {expectedMs1})");
                Console.WriteLine($"  [PASS] ScanCount (MS2) = {index.ScanCount} (expected {expectedMs2})");
            }

            // ── Test 2: MS1 scans sorted ascending by RT ────────────────────
            {
                var scans = GenerateSyntheticScansWithMs1(
                    windowCount: 10, scansPerWindow: 200,
                    peaksPerMs2: 100, ms1ScansPerWindow: 80, peaksPerMs1: 300);

                using var index = DiaScanIndexBuilder.Build(scans);

                float prevRt = float.MinValue;
                bool sorted = true;
                for (int i = 0; i < index.Ms1ScanCount; i++)
                {
                    float rt = index.GetMs1ScanRt(i);
                    if (rt < prevRt) { sorted = false; break; }
                    prevRt = rt;
                }
                Debug.Assert(sorted, "MS1 scans must be sorted ascending by RT.");
                Console.WriteLine($"  [PASS] MS1 scans sorted ascending by RT (n={index.Ms1ScanCount})");
            }

            // ── Test 3: GetMs1ScanPeaks returns correct peak counts ──────────
            {
                int peaksPerMs1 = 150;
                var scans = GenerateSyntheticScansWithMs1(
                    windowCount: 5, scansPerWindow: 50,
                    peaksPerMs2: 100, ms1ScansPerWindow: 20, peaksPerMs1: peaksPerMs1);

                using var index = DiaScanIndexBuilder.Build(scans);

                bool allCorrect = true;
                for (int i = 0; i < index.Ms1ScanCount; i++)
                {
                    index.GetMs1ScanPeaks(i, out var mzs, out var intensities);
                    if (mzs.Length != peaksPerMs1 || intensities.Length != peaksPerMs1)
                    {
                        allCorrect = false;
                        Console.WriteLine($"  [FAIL] MS1 scan {i}: mzs.Length={mzs.Length}, intensities.Length={intensities.Length}, expected {peaksPerMs1}");
                        break;
                    }
                }
                Debug.Assert(allCorrect, "All MS1 scans must have the expected peak count.");
                Console.WriteLine($"  [PASS] All MS1 scans have {peaksPerMs1} peaks");
            }

            // ── Test 4: FindMs1ScanIndexAtRt binary search ──────────────────
            {
                var scans = GenerateSyntheticScansWithMs1(
                    windowCount: 10, scansPerWindow: 100,
                    peaksPerMs2: 50, ms1ScansPerWindow: 50, peaksPerMs1: 100);

                using var index = DiaScanIndexBuilder.Build(scans);

                // Edge case: RT before all scans → should return 0
                {
                    int idx = index.FindMs1ScanIndexAtRt(-999f);
                    Debug.Assert(idx == 0, $"Expected 0 for RT before all scans, got {idx}");
                    Console.WriteLine($"  [PASS] FindMs1ScanIndexAtRt(-999) = {idx} (expected 0)");
                }

                // Edge case: RT after all scans → should return Ms1ScanCount
                {
                    int idx = index.FindMs1ScanIndexAtRt(99999f);
                    Debug.Assert(idx == index.Ms1ScanCount,
                        $"Expected {index.Ms1ScanCount} for RT past end, got {idx}");
                    Console.WriteLine($"  [PASS] FindMs1ScanIndexAtRt(99999) = {idx} (expected {index.Ms1ScanCount})");
                }

                // Mid-run query: verify lower-bound semantics
                if (index.Ms1ScanCount > 2)
                {
                    int mid = index.Ms1ScanCount / 2;
                    float targetRt = index.GetMs1ScanRt(mid);

                    int foundIdx = index.FindMs1ScanIndexAtRt(targetRt);
                    Debug.Assert(foundIdx <= mid,
                        $"FindMs1ScanIndexAtRt(rt[{mid}]) returned {foundIdx} > {mid}");
                    // The returned scan must have RT >= targetRt
                    Debug.Assert(index.GetMs1ScanRt(foundIdx) >= targetRt,
                        $"Scan at foundIdx={foundIdx} has RT={index.GetMs1ScanRt(foundIdx):F4} < targetRt={targetRt:F4}");
                    // The scan before it (if any) must have RT < targetRt
                    if (foundIdx > 0)
                    {
                        Debug.Assert(index.GetMs1ScanRt(foundIdx - 1) < targetRt,
                            $"Scan before foundIdx has RT >= targetRt; lower-bound violated.");
                    }
                    Console.WriteLine($"  [PASS] FindMs1ScanIndexAtRt(rt[{mid}]={targetRt:F4}) = {foundIdx} (lower-bound correct)");
                }

                // RT window iteration test: count scans in [rtMin, rtMax]
                {
                    float rtMin = index.GetMs1ScanRt(index.Ms1ScanCount / 4);
                    float rtMax = index.GetMs1ScanRt(3 * index.Ms1ScanCount / 4);

                    int start = index.FindMs1ScanIndexAtRt(rtMin);
                    int count = 0;
                    for (int i = start; i < index.Ms1ScanCount; i++)
                    {
                        if (index.GetMs1ScanRt(i) > rtMax) break;
                        count++;
                    }

                    // Brute-force count for validation
                    int bruteCount = 0;
                    for (int i = 0; i < index.Ms1ScanCount; i++)
                    {
                        float rt = index.GetMs1ScanRt(i);
                        if (rt >= rtMin && rt <= rtMax) bruteCount++;
                    }

                    Debug.Assert(count == bruteCount,
                        $"RT window iteration mismatch: binary-search count={count}, brute-force={bruteCount}");
                    Console.WriteLine($"  [PASS] RT window [{rtMin:F3}, {rtMax:F3}] → {count} MS1 scans (brute-force verified)");
                }
            }

            // ── Test 5: MS2-only file (no MS1 scans) ────────────────────────
            {
                var scans = GenerateSyntheticScans(windowCount: 5, scansPerWindow: 50, peaksPerScan: 100);
                using var index = DiaScanIndexBuilder.Build(scans);

                Debug.Assert(index.Ms1ScanCount == 0, $"Expected 0 MS1 scans, got {index.Ms1ScanCount}");
                Debug.Assert(index.Ms1TotalPeakCount == 0, "Expected 0 MS1 peaks.");

                // FindMs1ScanIndexAtRt should return 0 on empty index
                int idx = index.FindMs1ScanIndexAtRt(1.0f);
                Debug.Assert(idx == 0, $"FindMs1ScanIndexAtRt on empty MS1 should return 0, got {idx}");

                Console.WriteLine($"  [PASS] MS2-only file: Ms1ScanCount=0, FindMs1ScanIndexAtRt returns 0");
            }

            Console.WriteLine();
            Console.WriteLine("  All MS1 storage verification tests passed.");
            Console.WriteLine();
        }

        /// <summary>
        /// Measures build time for a mixed MS1+MS2 scan array.
        /// 
        /// Reports total build time, MS1 scan count, and memory footprint.
        /// A 45-minute DIA gradient typically produces ~1,480 MS1 scans
        /// (one per ~1.8 seconds). With ~500 peaks per MS1 scan, that's
        /// ~740K additional peaks — negligible overhead on top of 15M MS2 peaks.
        /// </summary>
        private static void BenchmarkMs1Build(
            int windowCount, int scansPerWindow, int ms1ScansPerWindow,
            int peaksPerMs2, int peaksPerMs1)
        {
            int totalMs1 = ms1ScansPerWindow;
            int totalMs2 = windowCount * scansPerWindow;
            long totalMs1Peaks = (long)totalMs1 * peaksPerMs1;
            long totalMs2Peaks = (long)totalMs2 * peaksPerMs2;

            Console.WriteLine($"MS1 build benchmark: {totalMs2:N0} MS2 + {totalMs1:N0} MS1 scans");
            Console.WriteLine($"  MS2 peaks: {totalMs2Peaks:N0} | MS1 peaks: {totalMs1Peaks:N0}");

            var scans = GenerateSyntheticScansWithMs1(
                windowCount, scansPerWindow, peaksPerMs2, ms1ScansPerWindow, peaksPerMs1);

            // Warm up
            var warmupScans = GenerateSyntheticScansWithMs1(2, 10, 50, ms1ScansPerWindow: 5, peaksPerMs1: 50);
            using (var _ = DiaScanIndexBuilder.Build(warmupScans)) { }

            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();

            long memBefore = GC.GetTotalMemory(true);
            var sw = Stopwatch.StartNew();
            using var index = DiaScanIndexBuilder.Build(scans);
            sw.Stop();
            long memAfter = GC.GetTotalMemory(false);

            double mbAllocated = (memAfter - memBefore) / (1024.0 * 1024.0);

            Console.WriteLine($"  Build time:      {sw.ElapsedMilliseconds} ms");
            Console.WriteLine($"  Memory delta:    {mbAllocated:F1} MB");
            Console.WriteLine($"  MS2 scan count:  {index.ScanCount:N0}");
            Console.WriteLine($"  MS1 scan count:  {index.Ms1ScanCount:N0}");
            Console.WriteLine($"  MS1 peak count:  {index.Ms1TotalPeakCount:N0}");
            Console.WriteLine($"  MS1 RT range:    [{index.GetMs1GlobalRtMin():F3}, {index.GetMs1GlobalRtMax():F3}] min");
            Console.WriteLine();
        }

        /// <summary>
        /// Benchmarks FindMs1ScanIndexAtRt (binary search) throughput.
        /// 
        /// In MS1 XIC extraction (Prompt 4), this is called once per precursor per
        /// extraction window. With 77K precursors, it runs ~77K times per search run.
        /// Target: well above 1M calls/sec (trivial for a binary search over ~1,500 elements).
        /// </summary>
        private static void BenchmarkMs1RtLookup(
            int windowCount, int scansPerWindow, int ms1ScansPerWindow, int peaksPerMs1)
        {
            Console.WriteLine($"MS1 RT lookup benchmark: {ms1ScansPerWindow} MS1 scans");

            var scans = GenerateSyntheticScansWithMs1(
                windowCount, scansPerWindow, peaksPerMs2: 100,
                ms1ScansPerWindow: ms1ScansPerWindow, peaksPerMs1: peaksPerMs1);

            using var index = DiaScanIndexBuilder.Build(scans);

            if (index.Ms1ScanCount == 0)
            {
                Console.WriteLine("  Skipped (no MS1 scans).");
                Console.WriteLine();
                return;
            }

            float rtMin = index.GetMs1GlobalRtMin();
            float rtMax = index.GetMs1GlobalRtMax();
            float rtRange = rtMax - rtMin;

            // Generate random RT query values spread across the run
            var rng = new Random(99);
            int queryCount = 500_000;
            var queryRts = new float[queryCount];
            for (int i = 0; i < queryCount; i++)
                queryRts[i] = rtMin + (float)(rng.NextDouble() * rtRange);

            // Warm up
            for (int i = 0; i < 1000; i++)
                index.FindMs1ScanIndexAtRt(queryRts[i]);

            GC.Collect();
            GC.WaitForPendingFinalizers();

            var sw = Stopwatch.StartNew();
            int checksum = 0;
            for (int i = 0; i < queryCount; i++)
                checksum += index.FindMs1ScanIndexAtRt(queryRts[i]);
            sw.Stop();

            double lookupsPerSec = queryCount / sw.Elapsed.TotalSeconds;
            Console.WriteLine($"  {lookupsPerSec:N0} lookups/sec ({queryCount:N0} queries, {sw.ElapsedMilliseconds} ms)");
            Console.WriteLine($"  Checksum (prevent optimization): {checksum}");
            Console.WriteLine();
        }
    }
}