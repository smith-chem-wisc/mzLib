// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using System.Diagnostics;

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
        /// Generates a realistic array of synthetic DIA scans for benchmarking.
        /// 
        /// Simulates a DIA acquisition with:
        ///   - Interleaved MS1 and MS2 scans (MS1 every windowCount scans)
        ///   - Equally-spaced isolation windows from 400-1200 m/z
        ///   - Each MS2 scan has peaksPerScan peaks with random m/z in 100-2000 range
        ///   - Retention times increase linearly
        /// </summary>
        private static MsDataScan[] GenerateSyntheticScans(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            var rng = new Random(42); // Fixed seed for reproducibility
            int totalMs2 = windowCount * scansPerWindow;

            // Insert an MS1 before each cycle of windows
            int totalScans = totalMs2 + scansPerWindow; // approximate MS1 count
            var scans = new MsDataScan[totalMs2]; // Only MS2 for simplicity

            double windowWidth = 25.0;
            double windowStart = 400.0;
            double windowSpacing = (1200.0 - windowStart) / windowCount;

            int scanIdx = 0;
            for (int cycle = 0; cycle < scansPerWindow; cycle++)
            {
                double cycleRt = cycle * 0.04; // ~0.04 min per cycle ≈ 100 min run

                for (int w = 0; w < windowCount; w++)
                {
                    double isolationCenter = windowStart + w * windowSpacing + windowSpacing / 2.0;
                    double rt = cycleRt + w * 0.001; // Slight RT offset within cycle

                    double[] mzValues = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];
                    for (int p = 0; p < peaksPerScan; p++)
                    {
                        mzValues[p] = 100.0 + rng.NextDouble() * 1900.0;
                        intensities[p] = rng.NextDouble() * 10000.0;
                    }
                    Array.Sort(mzValues, intensities); // m/z must be sorted

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
