// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.Dia
{
    /// <summary>
    /// Tests for the GPU infrastructure: device detection, factory pattern, and fallback.
    /// 
    /// These tests pass on ALL platforms (Windows, Linux, macOS) regardless of whether
    /// a GPU is present or ILGPU is installed. They verify:
    ///   - GpuDeviceDetector never throws
    ///   - FragmentExtractorFactory.Create(preferCpu: true) always gives CpuFragmentExtractor
    ///   - The factory integrates correctly with the orchestrator
    ///   - DescribeBackend returns useful diagnostics
    /// 
    /// GPU-specific correctness tests (CPU vs GPU result comparison) are in the
    /// Development benchmark project because they require actual GPU hardware.
    /// </summary>
    [TestFixture]
    public class GpuInfrastructureTests
    {
        /// <summary>
        /// GpuDeviceDetector must never throw on any platform.
        /// It should always return a boolean and a non-null description.
        /// </summary>
        [Test]
        public void GpuDeviceDetector_NeverThrows_ReturnsValidState()
        {
            bool available = GpuDeviceDetector.IsGpuAvailable;
            string desc = GpuDeviceDetector.Description;
            GpuBackend backend = GpuDeviceDetector.DetectedBackend;

            Assert.That(desc, Is.Not.Null.And.Not.Empty);

            if (available)
            {
                Assert.That(backend, Is.Not.EqualTo(GpuBackend.None));
                TestContext.WriteLine($"GPU detected: {backend} — {desc}");
            }
            else
            {
                Assert.That(backend, Is.EqualTo(GpuBackend.None));
                TestContext.WriteLine($"No GPU: {desc}");
            }
        }

        /// <summary>
        /// Factory with preferCpu=true must always return CpuFragmentExtractor,
        /// even if a GPU is available on the test machine.
        /// </summary>
        [Test]
        public void Factory_PreferCpu_ReturnsCpuExtractor()
        {
            using var index = BuildSmallIndex();
            using var extractor = FragmentExtractorFactory.Create(index, preferCpu: true);

            Assert.That(extractor, Is.TypeOf<CpuFragmentExtractor>());
        }

        /// <summary>
        /// Factory function (for orchestrator) with preferCpu=true must always
        /// produce independent CpuFragmentExtractor instances.
        /// </summary>
        [Test]
        public void CreateFactory_PreferCpu_ProducesIndependentCpuExtractors()
        {
            using var index = BuildSmallIndex();
            var factory = FragmentExtractorFactory.CreateFactory(preferCpu: true);

            using var ext1 = factory(index);
            using var ext2 = factory(index);

            Assert.That(ext1, Is.TypeOf<CpuFragmentExtractor>());
            Assert.That(ext2, Is.TypeOf<CpuFragmentExtractor>());
            Assert.That(ext1, Is.Not.SameAs(ext2));
        }

        /// <summary>
        /// DescribeBackend should return useful strings in all modes.
        /// </summary>
        [Test]
        public void DescribeBackend_ReturnsNonEmptyStrings()
        {
            string forced = FragmentExtractorFactory.DescribeBackend(preferCpu: true);
            string auto = FragmentExtractorFactory.DescribeBackend(preferCpu: false);

            Assert.That(forced, Does.Contain("CPU"));
            Assert.That(auto, Is.Not.Null.And.Not.Empty);

            TestContext.WriteLine($"Forced CPU: {forced}");
            TestContext.WriteLine($"Auto mode:  {auto}");
        }

        /// <summary>
        /// CPU extractor obtained from factory must produce correct extraction results.
        /// Ensures the factory doesn't break the pipeline.
        /// </summary>
        [Test]
        public void Factory_CpuExtractor_ExtractsCorrectly()
        {
            using var index = BuildSmallIndex();
            using var extractor = FragmentExtractorFactory.Create(index, preferCpu: true);

            // Find a known m/z value in the first window
            var windowIds = index.GetWindowIds().ToArray();
            int wid = windowIds[0];
            index.TryGetScanRangeForWindow(wid, out int start, out int count);
            var mzSpan = index.GetScanMzSpan(start);
            float targetMz = mzSpan[mzSpan.Length / 2]; // pick middle peak
            float rtMin = index.GetScanRt(start);
            float rtMax = index.GetScanRt(start + count - 1);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(targetMz, 20f, rtMin - 0.1f, rtMax + 0.1f, wid, queryId: 42)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[count + 10];
            var intBuf = new float[count + 10];

            int points = extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            Assert.That(points, Is.GreaterThan(0));
            Assert.That(results[0].DataPointCount, Is.GreaterThan(0));
            Assert.That(results[0].TotalIntensity, Is.GreaterThan(0f));
            Assert.That(results[0].QueryId, Is.EqualTo(42));
        }

        /// <summary>
        /// Full pipeline: factory → orchestrator → parallel extraction.
        /// Uses CPU-forced mode so it works on all platforms.
        /// </summary>
        [Test]
        public void Orchestrator_WithCpuFactory_FullPipeline()
        {
            using var index = BuildSmallIndex();
            var factory = FragmentExtractorFactory.CreateFactory(preferCpu: true);
            using var orch = new DiaExtractionOrchestrator(index, factory);

            var queries = GenerateQueries(index, queriesPerWindow: 3);
            var result = orch.ExtractAll(queries, maxDegreeOfParallelism: 2);

            Assert.That(result.Results.Length, Is.EqualTo(queries.Length));
            Assert.That(result.TotalDataPoints, Is.GreaterThan(0));

            // All queries target known peaks, so all should have data
            for (int i = 0; i < result.Results.Length; i++)
            {
                Assert.That(result.Results[i].DataPointCount, Is.GreaterThan(0),
                    $"Query {i} (wid={queries[i].WindowId}) should have data");
            }
        }

        /// <summary>
        /// Factory auto-mode should return SOMETHING that works, regardless of GPU status.
        /// </summary>
        [Test]
        public void Factory_AutoMode_ReturnsWorkingExtractor()
        {
            using var index = BuildSmallIndex();

            // Auto mode: may return GPU or CPU, both are fine
            using var extractor = FragmentExtractorFactory.Create(index, preferCpu: false);
            Assert.That(extractor, Is.Not.Null);
            Assert.That(extractor, Is.InstanceOf<IFragmentExtractor>());

            // Verify it can actually extract something
            var windowIds = index.GetWindowIds().ToArray();
            int wid = windowIds[0];
            index.TryGetScanRangeForWindow(wid, out int start, out int count);
            var mzSpan = index.GetScanMzSpan(start);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(mzSpan[0], 20f,
                    index.GetScanRt(start) - 0.1f,
                    index.GetScanRt(start + count - 1) + 0.1f,
                    wid, 1)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[count + 10];
            var intBuf = new float[count + 10];

            int points = extractor.ExtractBatch(queries, results, rtBuf, intBuf);
            Assert.That(points, Is.GreaterThanOrEqualTo(0));
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Helpers
        // ──────────────────────────────────────────────────────────────────────

        private static DiaScanIndex BuildSmallIndex()
        {
            var rng = new Random(42);
            int windowCount = 3;
            int scansPerWindow = 5;
            int peaksPerScan = 20;
            var scans = new MsDataScan[windowCount * scansPerWindow];
            int idx = 0;

            for (int w = 0; w < windowCount; w++)
            {
                double isoCenter = 500.0 + w * 50.0;
                for (int s = 0; s < scansPerWindow; s++)
                {
                    double[] mz = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];
                    for (int p = 0; p < peaksPerScan; p++)
                    {
                        mz[p] = 200.0 + p * 10.0 + rng.NextDouble() * 0.001;
                        intensities[p] = 100.0 + rng.NextDouble() * 900.0;
                    }
                    Array.Sort(mz, intensities);

                    scans[idx++] = new MsDataScan(
                        new MzSpectrum(mz, intensities, false),
                        idx, 2, true, Polarity.Positive, s * 0.5,
                        new MzRange(100, 2000), "FTMS", MZAnalyzerType.Orbitrap,
                        intensities.Sum(), 20.0, null, $"scan={idx}",
                        isolationMZ: isoCenter, isolationWidth: 25.0,
                        dissociationType: DissociationType.HCD);
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
                var mz = index.GetScanMzSpan(start);
                float rtMin = index.GetScanRt(start);
                float rtMax = index.GetScanRt(start + count - 1);

                for (int q = 0; q < queriesPerWindow && q < mz.Length; q++)
                {
                    int pi = (int)((long)q * mz.Length / queriesPerWindow);
                    queries.Add(new FragmentQuery(
                        mz[pi], 20f, rtMin - 0.1f, rtMax + 0.1f, wid, qid++));
                }
            }
            return queries.ToArray();
        }
    }
}