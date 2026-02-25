// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace Development.Dia
{
    /// <summary>
    /// Benchmarks for DIA extraction orchestration.
    /// 
    /// Measures:
    ///   - CPU parallel speedup at various thread counts
    ///   - CPU vs GPU throughput comparison (if GPU available)
    ///   - Query count scaling
    ///   - CPU/GPU result correctness verification
    /// 
    /// Run via Development project: Set as Startup Project → Ctrl+F5.
    /// (Ctrl+F5 = start without debugger = accurate timing)
    /// </summary>
    public static class DiaOrchestrationBenchmark
    {
        public static void RunAll()
        {
            Console.WriteLine("=== DIA Orchestration Benchmark ===");
            Console.WriteLine();
            Console.WriteLine($"  Processors:    {Environment.ProcessorCount}");
            Console.WriteLine($"  OS:            {System.Runtime.InteropServices.RuntimeInformation.OSDescription}");
            Console.WriteLine($"  GPU backend:   {FragmentExtractorFactory.DescribeBackend()}");
            Console.WriteLine();

            // Build shared data
            int windowCount = 20;
            int scansPerWindow = 500;
            int peaksPerScan = 200;

            Console.WriteLine($"  Data:          {windowCount} windows × {scansPerWindow} scans × {peaksPerScan} peaks");
            Console.WriteLine($"                 {windowCount * scansPerWindow:N0} scans, " +
                $"{(long)windowCount * scansPerWindow * peaksPerScan:N0} peaks");

            var sw = Stopwatch.StartNew();
            using var index = BuildRealisticIndex(windowCount, scansPerWindow, peaksPerScan);
            sw.Stop();
            Console.WriteLine($"  Index build:   {sw.ElapsedMilliseconds} ms");
            Console.WriteLine();

            BenchmarkCpuParallelSpeedup(index, queriesPerWindow: 50);
            BenchmarkCpuVsGpu(index, queriesPerWindow: 50);
            BenchmarkQueryScaling(index);
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Benchmark 1: CPU parallel speedup
        // ──────────────────────────────────────────────────────────────────────

        private static void BenchmarkCpuParallelSpeedup(DiaScanIndex index, int queriesPerWindow)
        {
            Console.WriteLine("── CPU Parallel Speedup ──────────────────────────────────────");

            var queries = GenerateQueries(index, queriesPerWindow);
            Console.WriteLine($"  {queries.Length} queries ({queriesPerWindow}/window)");
            Console.WriteLine();

            // Warmup
            var warmupFactory = FragmentExtractorFactory.CreateFactory(preferCpu: true);
            using (var orch = new DiaExtractionOrchestrator(index, warmupFactory))
                orch.ExtractAll(queries, maxDegreeOfParallelism: 1);

            int[] threadCounts = { 1, 2, 4, 8, Math.Min(16, Environment.ProcessorCount) };
            double serialMs = 0;

            foreach (int threads in threadCounts)
            {
                var cpuFactory = FragmentExtractorFactory.CreateFactory(preferCpu: true);
                using var orch = new DiaExtractionOrchestrator(index, cpuFactory);

                const int iterations = 5;
                var timer = Stopwatch.StartNew();
                ExtractionResult result = null;

                for (int i = 0; i < iterations; i++)
                    result = orch.ExtractAll(queries, maxDegreeOfParallelism: threads);

                timer.Stop();
                double avgMs = timer.Elapsed.TotalMilliseconds / iterations;
                if (threads == 1) serialMs = avgMs;

                double speedup = serialMs / avgMs;
                double qPerSec = queries.Length / (avgMs / 1000.0);

                Console.WriteLine(
                    $"  {threads,2} threads:  {avgMs,8:F2} ms  |  " +
                    $"{speedup,5:F2}× speedup  |  " +
                    $"{qPerSec,10:N0} q/sec  |  " +
                    $"{result.TotalDataPoints:N0} XIC pts");
            }
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Benchmark 2: CPU vs GPU
        // ──────────────────────────────────────────────────────────────────────

        private static void BenchmarkCpuVsGpu(DiaScanIndex index, int queriesPerWindow)
        {
            Console.WriteLine("── CPU vs GPU Comparison ─────────────────────────────────────");

            var queries = GenerateQueries(index, queriesPerWindow);
            int cpuThreads = Math.Min(Environment.ProcessorCount, 16);
            const int iterations = 5;

            Console.WriteLine($"  {queries.Length} queries, {iterations} iterations each");
            Console.WriteLine();

            // ── CPU (multi-threaded) ──────────────────────────────────────────
            ExtractionResult cpuResult;
            {
                var cpuFactory = FragmentExtractorFactory.CreateFactory(preferCpu: true);
                using var orch = new DiaExtractionOrchestrator(index, cpuFactory);

                // Warmup
                orch.ExtractAll(queries, maxDegreeOfParallelism: cpuThreads);

                var timer = Stopwatch.StartNew();
                cpuResult = null;
                for (int i = 0; i < iterations; i++)
                    cpuResult = orch.ExtractAll(queries, maxDegreeOfParallelism: cpuThreads);
                timer.Stop();

                double avgMs = timer.Elapsed.TotalMilliseconds / iterations;
                double qPerSec = queries.Length / (avgMs / 1000.0);

                Console.WriteLine(
                    $"  CPU ({cpuThreads} threads): {avgMs,8:F2} ms  |  " +
                    $"{qPerSec,10:N0} q/sec  |  " +
                    $"{cpuResult.TotalDataPoints:N0} XIC pts");
            }

            // ── GPU ───────────────────────────────────────────────────────────
            if (!GpuDeviceDetector.IsGpuAvailable)
            {
                Console.WriteLine($"  GPU:              skipped — {GpuDeviceDetector.Description}");
                Console.WriteLine();
                return;
            }

            try
            {
                // GPU orchestrator uses parallelism=1 because extraction is on GPU
                var gpuFactory = FragmentExtractorFactory.CreateFactory(preferCpu: false);
                using var orch = new DiaExtractionOrchestrator(index, gpuFactory);

                // Warmup (includes kernel JIT compilation + first transfer)
                orch.ExtractAll(queries, maxDegreeOfParallelism: 1);

                var timer = Stopwatch.StartNew();
                ExtractionResult gpuResult = null;
                for (int i = 0; i < iterations; i++)
                    gpuResult = orch.ExtractAll(queries, maxDegreeOfParallelism: 1);
                timer.Stop();

                double avgMs = timer.Elapsed.TotalMilliseconds / iterations;
                double qPerSec = queries.Length / (avgMs / 1000.0);

                Console.WriteLine(
                    $"  GPU:              {avgMs,8:F2} ms  |  " +
                    $"{qPerSec,10:N0} q/sec  |  " +
                    $"{gpuResult.TotalDataPoints:N0} XIC pts");

                // ── Verify CPU and GPU produce matching results ───────────────
                bool match = VerifyResults(cpuResult, gpuResult, queries.Length);
                Console.WriteLine($"  CPU/GPU match:    {(match ? "YES ✓" : "MISMATCH ✗")}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"  GPU:              FAILED — {ex.Message}");
            }

            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Benchmark 3: Query count scaling
        // ──────────────────────────────────────────────────────────────────────

        private static void BenchmarkQueryScaling(DiaScanIndex index)
        {
            Console.WriteLine("── Query Count Scaling ───────────────────────────────────────");

            bool useGpu = GpuDeviceDetector.IsGpuAvailable;
            int threads = useGpu ? 1 : Math.Min(Environment.ProcessorCount, 16);
            string backend = useGpu ? "GPU" : $"CPU ({threads} threads)";
            Console.WriteLine($"  Backend: {backend}");
            Console.WriteLine();

            int[] queriesPerWindowValues = { 5, 20, 50, 100, 200 };

            foreach (int qpw in queriesPerWindowValues)
            {
                var queries = GenerateQueries(index, qpw);
                var factory = FragmentExtractorFactory.CreateFactory(preferCpu: !useGpu);
                using var orch = new DiaExtractionOrchestrator(index, factory);

                // Warmup
                orch.ExtractAll(queries, maxDegreeOfParallelism: threads);

                const int iterations = 3;
                var timer = Stopwatch.StartNew();
                ExtractionResult result = null;
                for (int i = 0; i < iterations; i++)
                    result = orch.ExtractAll(queries, maxDegreeOfParallelism: threads);
                timer.Stop();

                double avgMs = timer.Elapsed.TotalMilliseconds / iterations;
                double qPerSec = queries.Length / (avgMs / 1000.0);
                double ptsPerQ = (double)result.TotalDataPoints / Math.Max(queries.Length, 1);

                Console.WriteLine(
                    $"  {queries.Length,6} queries:  {avgMs,8:F2} ms  |  " +
                    $"{qPerSec,10:N0} q/sec  |  " +
                    $"{ptsPerQ:F1} pts/query");
            }
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Result verification
        // ──────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Verifies that two extraction results match: same data point counts
        /// and same total intensities per query (within float tolerance).
        /// </summary>
        private static bool VerifyResults(ExtractionResult a, ExtractionResult b, int queryCount)
        {
            if (a.TotalDataPoints != b.TotalDataPoints)
            {
                Console.Error.WriteLine(
                    $"  [VERIFY] Total data points differ: {a.TotalDataPoints} vs {b.TotalDataPoints}");
                return false;
            }

            bool allMatch = true;
            int mismatches = 0;
            for (int i = 0; i < queryCount; i++)
            {
                if (a.Results[i].DataPointCount != b.Results[i].DataPointCount)
                {
                    allMatch = false;
                    mismatches++;
                }
                else if (Math.Abs(a.Results[i].TotalIntensity - b.Results[i].TotalIntensity) >
                         a.Results[i].TotalIntensity * 1e-4f) // 0.01% relative tolerance for float
                {
                    allMatch = false;
                    mismatches++;
                }
            }

            if (!allMatch)
            {
                Console.Error.WriteLine(
                    $"  [VERIFY] {mismatches}/{queryCount} queries differ between CPU and GPU");
            }

            return allMatch;
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Data generation (shared across benchmarks)
        // ──────────────────────────────────────────────────────────────────────

        private static DiaScanIndex BuildRealisticIndex(
            int windowCount, int scansPerWindow, int peaksPerScan)
        {
            var rng = new Random(42);
            int totalMs2 = windowCount * scansPerWindow;
            var scans = new MsDataScan[totalMs2];

            double windowWidth = 25.0;
            double windowStart = 400.0;
            double windowSpacing = (1200.0 - windowStart) / windowCount;

            // "Resident" fragments per window — appear in every scan (realistic XIC signal)
            int residentsPerWindow = Math.Min(peaksPerScan / 2, 100);
            var windowResidents = new double[windowCount][];
            for (int w = 0; w < windowCount; w++)
            {
                windowResidents[w] = new double[residentsPerWindow];
                for (int r = 0; r < residentsPerWindow; r++)
                    windowResidents[w][r] = 200.0 + r * 15.0 + rng.NextDouble() * 0.01;
            }

            int scanIdx = 0;
            for (int cycle = 0; cycle < scansPerWindow; cycle++)
            {
                double cycleRt = cycle * 0.04;
                for (int w = 0; w < windowCount; w++)
                {
                    double isoCenter = windowStart + w * windowSpacing + windowSpacing / 2.0;
                    double rt = cycleRt + w * 0.001;

                    double[] mz = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];

                    for (int r = 0; r < residentsPerWindow; r++)
                    {
                        mz[r] = windowResidents[w][r];
                        intensities[r] = 500.0 + rng.NextDouble() * 5000.0;
                    }
                    for (int p = residentsPerWindow; p < peaksPerScan; p++)
                    {
                        mz[p] = 100.0 + rng.NextDouble() * 1900.0;
                        intensities[p] = rng.NextDouble() * 1000.0;
                    }
                    Array.Sort(mz, intensities);

                    scans[scanIdx] = new MsDataScan(
                        new MzSpectrum(mz, intensities, false),
                        scanIdx + 1, 2, true, Polarity.Positive, rt,
                        new MzRange(100, 2000), "FTMS", MZAnalyzerType.Orbitrap,
                        intensities.Sum(), 20.0, null, $"scan={scanIdx + 1}",
                        isolationMZ: isoCenter, isolationWidth: windowWidth,
                        dissociationType: DissociationType.HCD);
                    scanIdx++;
                }
            }

            return DiaScanIndexBuilder.Build(scans);
        }

        private static FragmentQuery[] GenerateQueries(DiaScanIndex index, int queriesPerWindow)
        {
            var queries = new List<FragmentQuery>();
            int qid = 0;
            foreach (int wid in index.GetWindowIds())
            {
                index.TryGetScanRangeForWindow(wid, out int start, out int count);
                var mzSpan = index.GetScanMzSpan(start);
                float rtMin = index.GetScanRt(start);
                float rtMax = index.GetScanRt(start + count - 1);

                for (int q = 0; q < queriesPerWindow && q < mzSpan.Length; q++)
                {
                    int pi = (int)((long)q * mzSpan.Length / queriesPerWindow);
                    queries.Add(new FragmentQuery(
                        mzSpan[pi], 20f,
                        rtMin + (rtMax - rtMin) * 0.2f,
                        rtMax - (rtMax - rtMin) * 0.2f,
                        wid, qid++));
                }
            }
            return queries.ToArray();
        }
    }
}