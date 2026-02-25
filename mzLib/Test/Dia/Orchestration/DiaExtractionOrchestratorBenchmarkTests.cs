// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace Test.Dia
{
    /// <summary>
    /// Performance benchmarks for parallel DIA extraction orchestration.
    /// 
    /// Marked [Explicit] — run manually from Test Explorer to measure speedup.
    /// Compares single-threaded vs multi-threaded extraction across realistic
    /// DIA data scales.
    /// </summary>
    [TestFixture]
    public class DiaExtractionOrchestratorBenchmarkTests
    {
        /// <summary>
        /// Benchmarks parallel extraction speedup at various thread counts.
        /// Uses 20 windows × 500 scans × 200 peaks, with 50 queries per window (1000 total).
        /// 
        /// Expected: near-linear speedup up to the window count, then diminishing returns.
        /// </summary>
        [Test, Explicit("Performance benchmark — run manually from Test Explorer")]
        public void Benchmark_ParallelSpeedup_VariousThreadCounts()
        {
            int windowCount = 20;
            int scansPerWindow = 500;
            int peaksPerScan = 200;
            int queriesPerWindow = 50;

            TestContext.WriteLine($"Setup: {windowCount} windows × {scansPerWindow} scans × {peaksPerScan} peaks");
            TestContext.WriteLine($"Queries: {queriesPerWindow} per window = {windowCount * queriesPerWindow} total");
            TestContext.WriteLine();

            // Build index
            using var index = BuildBenchmarkIndex(windowCount, scansPerWindow, peaksPerScan);
            var queries = GenerateBenchmarkQueries(index, queriesPerWindow);

            // Warmup
            using (var orch = new DiaExtractionOrchestrator(index))
            {
                orch.ExtractAll(queries, maxDegreeOfParallelism: 1);
            }

            int[] threadCounts = { 1, 2, 4, 8, 16 };
            double serialMs = 0;

            foreach (int threads in threadCounts)
            {
                using var orch = new DiaExtractionOrchestrator(index);

                var sw = Stopwatch.StartNew();
                int iterations = 5;
                ExtractionResult result = null;

                for (int i = 0; i < iterations; i++)
                {
                    result = orch.ExtractAll(queries, maxDegreeOfParallelism: threads);
                }

                sw.Stop();
                double avgMs = sw.Elapsed.TotalMilliseconds / iterations;

                if (threads == 1) serialMs = avgMs;
                double speedup = serialMs / avgMs;

                TestContext.WriteLine($"  Threads: {threads,2}  |  Time: {avgMs,8:F2} ms  |  Speedup: {speedup,5:F2}×  |  " +
                    $"Queries/sec: {queries.Length / (avgMs / 1000.0):N0}  |  XIC points: {result.TotalDataPoints:N0}");
            }

            TestContext.WriteLine();
            TestContext.WriteLine($"Available processors: {Environment.ProcessorCount}");
        }

        /// <summary>
        /// Benchmarks scaling with query count at a fixed thread count.
        /// </summary>
        [Test, Explicit("Performance benchmark — run manually from Test Explorer")]
        public void Benchmark_QueryCountScaling()
        {
            int windowCount = 20;
            int scansPerWindow = 500;
            int peaksPerScan = 200;

            using var index = BuildBenchmarkIndex(windowCount, scansPerWindow, peaksPerScan);

            TestContext.WriteLine($"Query count scaling (8 threads, {windowCount} windows × {scansPerWindow} scans):");

            int[] queriesPerWindowValues = { 5, 20, 50, 100, 200 };

            foreach (int qpw in queriesPerWindowValues)
            {
                var queries = GenerateBenchmarkQueries(index, qpw);

                using var orch = new DiaExtractionOrchestrator(index);
                var sw = Stopwatch.StartNew();
                int iterations = 3;
                ExtractionResult result = null;

                for (int i = 0; i < iterations; i++)
                {
                    result = orch.ExtractAll(queries, maxDegreeOfParallelism: 8);
                }

                sw.Stop();
                double avgMs = sw.Elapsed.TotalMilliseconds / iterations;
                double queriesPerSec = queries.Length / (avgMs / 1000.0);

                TestContext.WriteLine($"  Queries: {queries.Length,6}  |  Time: {avgMs,8:F2} ms  |  " +
                    $"Q/sec: {queriesPerSec:N0}  |  XIC pts/query: {(double)result.TotalDataPoints / queries.Length:F1}");
            }
        }

        // ──────────────────────────────────────────────────────────────────────

        private static DiaScanIndex BuildBenchmarkIndex(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            var rng = new Random(42);
            int totalMs2 = windowCount * scansPerWindow;
            var scans = new MsDataScan[totalMs2];

            double windowWidth = 25.0;
            double windowStart = 400.0;
            double windowSpacing = (1200.0 - windowStart) / windowCount;

            // Pre-generate "resident" m/z values per window (fragments that appear in every scan)
            int residentsPerWindow = Math.Min(peaksPerScan / 2, 100);
            var windowResidents = new double[windowCount][];
            for (int w = 0; w < windowCount; w++)
            {
                windowResidents[w] = new double[residentsPerWindow];
                for (int r = 0; r < residentsPerWindow; r++)
                {
                    windowResidents[w][r] = 200.0 + r * 15.0 + rng.NextDouble() * 0.01;
                }
            }

            int scanIdx = 0;
            for (int cycle = 0; cycle < scansPerWindow; cycle++)
            {
                double cycleRt = cycle * 0.04;

                for (int w = 0; w < windowCount; w++)
                {
                    double isolationCenter = windowStart + w * windowSpacing + windowSpacing / 2.0;
                    double rt = cycleRt + w * 0.001;

                    double[] mzValues = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];

                    // Fill resident fragments
                    int residentCount = windowResidents[w].Length;
                    for (int r = 0; r < residentCount; r++)
                    {
                        mzValues[r] = windowResidents[w][r];
                        intensities[r] = 500.0 + rng.NextDouble() * 5000.0;
                    }

                    // Fill remaining with random peaks
                    for (int p = residentCount; p < peaksPerScan; p++)
                    {
                        mzValues[p] = 100.0 + rng.NextDouble() * 1900.0;
                        intensities[p] = rng.NextDouble() * 1000.0;
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

            return DiaScanIndexBuilder.Build(scans);
        }

        private static FragmentQuery[] GenerateBenchmarkQueries(DiaScanIndex index, int queriesPerWindow)
        {
            var windowIds = index.GetWindowIds().ToArray();
            var queries = new List<FragmentQuery>();
            int queryId = 0;

            foreach (int windowId in windowIds)
            {
                index.TryGetScanRangeForWindow(windowId, out int start, out int count);
                var mzSpan = index.GetScanMzSpan(start);
                float minRt = index.GetScanRt(start);
                float maxRt = index.GetScanRt(start + count - 1);

                for (int q = 0; q < queriesPerWindow && q < mzSpan.Length; q++)
                {
                    int peakIdx = (int)((long)q * mzSpan.Length / queriesPerWindow);
                    queries.Add(new FragmentQuery(
                        mzSpan[peakIdx], 20f,
                        minRt + (maxRt - minRt) * 0.2f,  // ±60% RT window
                        maxRt - (maxRt - minRt) * 0.2f,
                        windowId, queryId++));
                }
            }

            return queries.ToArray();
        }
    }
}
