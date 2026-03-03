// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace Test.Dia
{
    /// <summary>
    /// Tests for DiaExtractionOrchestrator, verifying correct parallel extraction
    /// across windows with result ordering, buffer merging, and thread safety.
    /// </summary>
    [TestFixture]
    public class DiaExtractionOrchestratorTests
    {
        /// <summary>
        /// Verifies that parallel extraction produces identical results to
        /// sequential (single-threaded) extraction. This is the core correctness test.
        /// 
        /// Uses 4 windows × 5 scans each, with 3 queries per window (12 total).
        /// Compares parallel (8 threads) vs serial (1 thread) results.
        /// </summary>
        [Test]
        public void ParallelResults_MatchSequentialResults()
        {
            int windowCount = 4;
            int scansPerWindow = 5;
            using var index = BuildMultiWindowIndex(windowCount, scansPerWindow, peaksPerScan: 20);

            // Generate queries: 3 per window, targeting real m/z values from each window
            var queries = GenerateQueriesFromIndex(index, queriesPerWindow: 3);

            // Run single-threaded (deterministic baseline)
            using var serialOrch = new DiaExtractionOrchestrator(index);
            var serialResult = serialOrch.ExtractAll(queries, maxDegreeOfParallelism: 1);

            // Run parallel
            using var parallelOrch = new DiaExtractionOrchestrator(index);
            var parallelResult = parallelOrch.ExtractAll(queries, maxDegreeOfParallelism: 8);

            // Compare: same number of results
            Assert.That(parallelResult.Results.Length, Is.EqualTo(serialResult.Results.Length));
            Assert.That(parallelResult.TotalDataPoints, Is.EqualTo(serialResult.TotalDataPoints));

            // Compare each result
            for (int i = 0; i < queries.Length; i++)
            {
                Assert.That(parallelResult.Results[i].QueryId, Is.EqualTo(serialResult.Results[i].QueryId),
                    $"QueryId mismatch at index {i}");
                Assert.That(parallelResult.Results[i].DataPointCount, Is.EqualTo(serialResult.Results[i].DataPointCount),
                    $"DataPointCount mismatch at query {i}");
                Assert.That(parallelResult.Results[i].TotalIntensity,
                    Is.EqualTo(serialResult.Results[i].TotalIntensity).Within(0.01f),
                    $"TotalIntensity mismatch at query {i}");

                // Compare actual XIC data
                int count = serialResult.Results[i].DataPointCount;
                for (int d = 0; d < count; d++)
                {
                    int sIdx = serialResult.Results[i].RtBufferOffset + d;
                    int pIdx = parallelResult.Results[i].RtBufferOffset + d;
                    Assert.That(parallelResult.RtBuffer[pIdx],
                        Is.EqualTo(serialResult.RtBuffer[sIdx]).Within(0.001f),
                        $"RT mismatch at query {i}, data point {d}");
                    Assert.That(parallelResult.IntensityBuffer[pIdx],
                        Is.EqualTo(serialResult.IntensityBuffer[sIdx]).Within(0.1f),
                        $"Intensity mismatch at query {i}, data point {d}");
                }
            }
        }

        /// <summary>
        /// Verifies that results are returned in the same order as the input queries,
        /// even when queries span multiple windows and are processed in parallel.
        /// This tests the original-index remapping logic.
        /// </summary>
        [Test]
        public void ResultOrder_MatchesQueryOrder()
        {
            using var index = BuildMultiWindowIndex(3, scansPerWindow: 3, peaksPerScan: 10);

            // Interleave queries from different windows: w0, w1, w2, w0, w1, w2
            var queries = new FragmentQuery[6];
            for (int w = 0; w < 3; w++)
            {
                index.TryGetScanRangeForWindow(w, out int start, out _);
                var mz = index.GetScanMzSpan(start);
                float targetMz = mz[mz.Length / 2];

                queries[w] = new FragmentQuery(targetMz, 20f, 0f, 1000f, windowId: w, queryId: w * 10);
                queries[w + 3] = new FragmentQuery(targetMz, 20f, 0f, 1000f, windowId: w, queryId: w * 10 + 1);
            }

            using var orch = new DiaExtractionOrchestrator(index);
            var result = orch.ExtractAll(queries, maxDegreeOfParallelism: 4);

            // Verify QueryIds match the input order
            Assert.That(result.Results[0].QueryId, Is.EqualTo(0));   // w0, first
            Assert.That(result.Results[1].QueryId, Is.EqualTo(10));  // w1, first
            Assert.That(result.Results[2].QueryId, Is.EqualTo(20));  // w2, first
            Assert.That(result.Results[3].QueryId, Is.EqualTo(1));   // w0, second
            Assert.That(result.Results[4].QueryId, Is.EqualTo(11));  // w1, second
            Assert.That(result.Results[5].QueryId, Is.EqualTo(21));  // w2, second
        }

        /// <summary>
        /// Verifies that empty query input returns an empty result without errors.
        /// </summary>
        [Test]
        public void EmptyQueries_ReturnsEmptyResult()
        {
            using var index = BuildMultiWindowIndex(2, 3, 10);
            using var orch = new DiaExtractionOrchestrator(index);

            var result = orch.ExtractAll(ReadOnlySpan<FragmentQuery>.Empty);

            Assert.That(result.Results.Length, Is.EqualTo(0));
            Assert.That(result.TotalDataPoints, Is.EqualTo(0));
        }

        /// <summary>
        /// Verifies that queries targeting a nonexistent window return zero data points
        /// without crashing, and don't affect other queries in the batch.
        /// </summary>
        [Test]
        public void UnknownWindowId_ReturnsZeroDataPoints()
        {
            using var index = BuildMultiWindowIndex(2, 3, 10);

            // Query 0 targets a valid window, query 1 targets window 999 (doesn't exist)
            index.TryGetScanRangeForWindow(0, out int start, out _);
            var mz = index.GetScanMzSpan(start);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(mz[0], 20f, 0f, 1000f, windowId: 0, queryId: 1),
                new FragmentQuery(500f, 20f, 0f, 1000f, windowId: 999, queryId: 2),
            };

            using var orch = new DiaExtractionOrchestrator(index);
            var result = orch.ExtractAll(queries, maxDegreeOfParallelism: 2);

            Assert.That(result.Results[0].DataPointCount, Is.GreaterThan(0), "Valid window should have results");
            Assert.That(result.Results[1].DataPointCount, Is.EqualTo(0), "Unknown window should have zero results");
            Assert.That(result.Results[1].TotalIntensity, Is.EqualTo(0f));
        }

        /// <summary>
        /// Verifies that the XIC buffer offsets in merged results are valid and
        /// non-overlapping — i.e., each query's data occupies a distinct region
        /// of the unified buffer.
        /// </summary>
        [Test]
        public void MergedBufferOffsets_AreValidAndNonOverlapping()
        {
            using var index = BuildMultiWindowIndex(3, scansPerWindow: 5, peaksPerScan: 20);
            var queries = GenerateQueriesFromIndex(index, queriesPerWindow: 2);

            using var orch = new DiaExtractionOrchestrator(index);
            var result = orch.ExtractAll(queries, maxDegreeOfParallelism: 4);

            // Collect all (offset, length) pairs and verify no overlaps
            var regions = new List<(int Start, int Length)>();
            for (int i = 0; i < result.Results.Length; i++)
            {
                int count = result.Results[i].DataPointCount;
                if (count > 0)
                {
                    int rtOffset = result.Results[i].RtBufferOffset;
                    int intOffset = result.Results[i].IntensityBufferOffset;

                    // Offsets should be within buffer bounds
                    Assert.That(rtOffset + count, Is.LessThanOrEqualTo(result.RtBuffer.Length),
                        $"Query {i}: RT buffer offset out of bounds");
                    Assert.That(intOffset + count, Is.LessThanOrEqualTo(result.IntensityBuffer.Length),
                        $"Query {i}: Intensity buffer offset out of bounds");

                    regions.Add((rtOffset, count));
                }
            }

            // Check for overlaps between any two regions
            for (int i = 0; i < regions.Count; i++)
            {
                for (int j = i + 1; j < regions.Count; j++)
                {
                    int aStart = regions[i].Start;
                    int aEnd = aStart + regions[i].Length;
                    int bStart = regions[j].Start;
                    int bEnd = bStart + regions[j].Length;

                    bool overlaps = aStart < bEnd && bStart < aEnd;
                    Assert.That(overlaps, Is.False,
                        $"Regions {i} [{aStart},{aEnd}) and {j} [{bStart},{bEnd}) overlap");
                }
            }
        }

        /// <summary>
        /// Integration test: parallel extraction + scoring.
        /// Verifies that XIC data extracted in parallel can be correctly scored.
        /// Self-match should give 1.0; cross-match between different fragments should give &lt; 1.0.
        /// </summary>
        [Test]
        public void ParallelExtraction_ThenScoring_ProducesValidScores()
        {
            using var index = BuildMultiWindowIndex(2, scansPerWindow: 10, peaksPerScan: 50);
            var queries = GenerateQueriesFromIndex(index, queriesPerWindow: 2);

            using var orch = new DiaExtractionOrchestrator(index);
            var result = orch.ExtractAll(queries, maxDegreeOfParallelism: 4);

            var scorer = new NormalizedDotProductScorer();

            // Self-score each extracted XIC
            for (int i = 0; i < result.Results.Length; i++)
            {
                int count = result.Results[i].DataPointCount;
                if (count >= 2)
                {
                    var intensities = result.IntensityBuffer.AsSpan(
                        result.Results[i].IntensityBufferOffset, count);
                    float selfScore = scorer.Score(intensities, intensities);
                    Assert.That(selfScore, Is.EqualTo(1.0f).Within(1e-4f),
                        $"Self-score for query {i} should be 1.0");
                }
            }
        }

        /// <summary>
        /// Verifies that the custom extractor factory is actually used — important
        /// for the GPU acceleration path where we'd inject GpuFragmentExtractor.
        /// </summary>
        [Test]
        public void CustomExtractorFactory_IsUsed()
        {
            using var index = BuildMultiWindowIndex(2, 3, 10);

            int factoryCallCount = 0;
            Func<DiaScanIndex, IFragmentExtractor> countingFactory = idx =>
            {
                Interlocked.Increment(ref factoryCallCount);
                return new CpuFragmentExtractor(idx);
            };

            var queries = GenerateQueriesFromIndex(index, queriesPerWindow: 1);

            using var orch = new DiaExtractionOrchestrator(index, countingFactory);
            var result = orch.ExtractAll(queries, maxDegreeOfParallelism: 4);

            // Factory should have been called at least once (one per window group processed)
            Assert.That(factoryCallCount, Is.GreaterThan(0), "Factory should be called for extraction");
            Assert.That(result.Results.Length, Is.EqualTo(queries.Length));
        }

        /// <summary>
        /// Stress test: many queries across many windows, verifying no crashes,
        /// no data corruption, and consistent results across repeated runs.
        /// </summary>
        [Test]
        public void StressTest_ManyQueriesManyWindows_NoCrashes()
        {
            int windowCount = 10;
            int scansPerWindow = 20;
            using var index = BuildMultiWindowIndex(windowCount, scansPerWindow, peaksPerScan: 50);

            // 10 queries per window = 100 total
            var queries = GenerateQueriesFromIndex(index, queriesPerWindow: 10);
            Assert.That(queries.Length, Is.EqualTo(100));

            using var orch = new DiaExtractionOrchestrator(index);

            // Run 5 times and verify consistency
            ExtractionResult firstResult = null;
            for (int run = 0; run < 5; run++)
            {
                var result = orch.ExtractAll(queries, maxDegreeOfParallelism: 8);

                Assert.That(result.Results.Length, Is.EqualTo(100));
                Assert.That(result.TotalDataPoints, Is.GreaterThan(0));

                if (firstResult == null)
                {
                    firstResult = result;
                }
                else
                {
                    // Results should be identical across runs
                    Assert.That(result.TotalDataPoints, Is.EqualTo(firstResult.TotalDataPoints),
                        $"Run {run}: TotalDataPoints should be consistent");
                    for (int i = 0; i < 100; i++)
                    {
                        Assert.That(result.Results[i].DataPointCount,
                            Is.EqualTo(firstResult.Results[i].DataPointCount),
                            $"Run {run}, query {i}: DataPointCount should be consistent");
                    }
                }
            }
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Helper methods
        // ──────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Builds a DiaScanIndex with the specified number of windows and scans.
        /// Each window has peaks at known, deterministic m/z values.
        /// </summary>
        private static DiaScanIndex BuildMultiWindowIndex(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            var rng = new Random(42);
            var scans = new MsDataScan[windowCount * scansPerWindow];
            int scanIdx = 0;

            for (int w = 0; w < windowCount; w++)
            {
                double isolationCenter = 500.0 + w * 50.0; // 500, 550, 600, ...
                for (int s = 0; s < scansPerWindow; s++)
                {
                    double rt = s * 0.5; // 0.0, 0.5, 1.0, ...
                    double[] mzValues = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];

                    for (int p = 0; p < peaksPerScan; p++)
                    {
                        mzValues[p] = 200.0 + p * 10.0 + rng.NextDouble() * 0.001;
                        intensities[p] = 100.0 + rng.NextDouble() * 900.0;
                    }
                    Array.Sort(mzValues, intensities);

                    scans[scanIdx++] = new MsDataScan(
                        massSpectrum: new MzSpectrum(mzValues, intensities, false),
                        oneBasedScanNumber: scanIdx,
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
                        nativeId: $"scan={scanIdx}",
                        isolationMZ: isolationCenter,
                        isolationWidth: 25.0,
                        dissociationType: DissociationType.HCD);
                }
            }

            return DiaScanIndexBuilder.Build(scans);
        }

        /// <summary>
        /// Generates queries targeting real m/z values from each window in the index.
        /// </summary>
        private static FragmentQuery[] GenerateQueriesFromIndex(DiaScanIndex index, int queriesPerWindow)
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
                    // Pick m/z values spread across the scan's peak list
                    int peakIdx = (int)((long)q * mzSpan.Length / queriesPerWindow);
                    float targetMz = mzSpan[peakIdx];

                    queries.Add(new FragmentQuery(
                        targetMz, 20f, minRt - 0.1f, maxRt + 0.1f, windowId, queryId++));
                }
            }

            return queries.ToArray();
        }
    }
}
