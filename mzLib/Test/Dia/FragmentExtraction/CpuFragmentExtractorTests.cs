// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Linq;

namespace Test.Dia
{
    /// <summary>
    /// Tests for CpuFragmentExtractor, verifying correct fragment ion extraction
    /// from the SoA DiaScanIndex. All tests use synthetic data with known values
    /// so we can verify exact correctness.
    /// </summary>
    [TestFixture]
    public class CpuFragmentExtractorTests
    {
        /// <summary>
        /// Verifies that a single fragment query correctly finds a matching peak
        /// in a scan and returns the correct intensity and RT.
        /// 
        /// Setup: one scan with peaks at m/z 500.0, 600.0, 700.0.
        /// Query: target m/z 600.0, 10 ppm tolerance.
        /// Expected: one XIC data point with the intensity of the 600.0 peak.
        /// </summary>
        [Test]
        public void SingleQuery_ExactMatch_ReturnsCorrectIntensity()
        {
            using var index = BuildIndex(new ScanDef(
                scanNumber: 1, rt: 5.0, isolationMz: 500.0, isolationWidth: 25.0,
                mzValues: new[] { 500.0, 600.0, 700.0 },
                intensities: new[] { 100.0, 250.0, 150.0 }));

            using var extractor = new CpuFragmentExtractor(index);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(targetMz: 600.0f, tolerancePpm: 10f,
                    rtMin: 0f, rtMax: 10f, windowId: 0, queryId: 1)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[100];
            var intBuf = new float[100];

            int totalPoints = extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            Assert.That(totalPoints, Is.EqualTo(1));
            Assert.That(results[0].QueryId, Is.EqualTo(1));
            Assert.That(results[0].DataPointCount, Is.EqualTo(1));
            Assert.That(intBuf[results[0].IntensityBufferOffset], Is.EqualTo(250.0f).Within(0.1f));
            Assert.That(rtBuf[results[0].RtBufferOffset], Is.EqualTo(5.0f).Within(0.01f));
        }

        /// <summary>
        /// Verifies that ppm tolerance correctly includes peaks slightly off-target
        /// and excludes peaks beyond the tolerance.
        /// 
        /// At m/z 1000.0 with 10 ppm tolerance:
        ///   - Match range: [999.99, 1000.01]
        ///   - Peak at 1000.005 should match
        ///   - Peak at 1000.02 should NOT match
        /// </summary>
        [Test]
        public void PpmTolerance_IncludesNearPeaks_ExcludesFarPeaks()
        {
            using var index = BuildIndex(new ScanDef(
                scanNumber: 1, rt: 5.0, isolationMz: 800.0, isolationWidth: 400.0,
                mzValues: new[] { 999.005, 1000.005, 1000.02, 1001.0 },
                intensities: new[] { 50.0, 200.0, 300.0, 400.0 }));

            using var extractor = new CpuFragmentExtractor(index);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(targetMz: 1000.0f, tolerancePpm: 10f,
                    rtMin: 0f, rtMax: 10f, windowId: 0, queryId: 1)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[100];
            var intBuf = new float[100];

            extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            // Only the peak at 1000.005 should match (within 5 ppm of 1000.0)
            // The peak at 1000.02 is 20 ppm away — outside tolerance
            Assert.That(results[0].DataPointCount, Is.EqualTo(1));
            Assert.That(results[0].TotalIntensity, Is.EqualTo(200.0f).Within(0.1f));
        }

        /// <summary>
        /// Verifies that the RT window correctly restricts which scans contribute
        /// to the extracted chromatogram.
        /// 
        /// Setup: 5 scans at RT 1, 2, 3, 4, 5 — all with a peak at m/z 500.0.
        /// Query: RT window [2.0, 4.0].
        /// Expected: 3 XIC data points (from scans at RT 2, 3, 4).
        /// </summary>
        [Test]
        public void RtWindow_RestrictsToCorrectScans()
        {
            var scanDefs = new[]
            {
                new ScanDef(1, 1.0, 600.0, 25.0, new[] { 500.0 }, new[] { 10.0 }),
                new ScanDef(2, 2.0, 600.0, 25.0, new[] { 500.0 }, new[] { 20.0 }),
                new ScanDef(3, 3.0, 600.0, 25.0, new[] { 500.0 }, new[] { 30.0 }),
                new ScanDef(4, 4.0, 600.0, 25.0, new[] { 500.0 }, new[] { 40.0 }),
                new ScanDef(5, 5.0, 600.0, 25.0, new[] { 500.0 }, new[] { 50.0 }),
            };
            using var index = BuildIndex(scanDefs);
            using var extractor = new CpuFragmentExtractor(index);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(targetMz: 500.0f, tolerancePpm: 10f,
                    rtMin: 2.0f, rtMax: 4.0f, windowId: 0, queryId: 1)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[100];
            var intBuf = new float[100];

            int totalPoints = extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            Assert.That(results[0].DataPointCount, Is.EqualTo(3));
            Assert.That(results[0].TotalIntensity, Is.EqualTo(90.0f).Within(0.1f)); // 20 + 30 + 40

            // Verify RT values are correct and ordered
            int offset = results[0].RtBufferOffset;
            Assert.That(rtBuf[offset], Is.EqualTo(2.0f).Within(0.01f));
            Assert.That(rtBuf[offset + 1], Is.EqualTo(3.0f).Within(0.01f));
            Assert.That(rtBuf[offset + 2], Is.EqualTo(4.0f).Within(0.01f));
        }

        /// <summary>
        /// Verifies that multiple queries in a single batch produce independent,
        /// correct results with non-overlapping buffer regions.
        /// </summary>
        [Test]
        public void BatchOfQueries_ProducesIndependentResults()
        {
            var scanDefs = new[]
            {
                // Window 0 (isolation 500): peaks at m/z 300 and 400
                new ScanDef(1, 1.0, 500.0, 25.0, new[] { 300.0, 400.0 }, new[] { 100.0, 200.0 }),
                // Window 1 (isolation 700): peak at m/z 350
                new ScanDef(2, 1.0, 700.0, 25.0, new[] { 350.0 }, new[] { 500.0 }),
            };
            using var index = BuildIndex(scanDefs);
            using var extractor = new CpuFragmentExtractor(index);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(300.0f, 10f, 0f, 10f, windowId: 0, queryId: 10),
                new FragmentQuery(400.0f, 10f, 0f, 10f, windowId: 0, queryId: 20),
                new FragmentQuery(350.0f, 10f, 0f, 10f, windowId: 1, queryId: 30),
            };
            var results = new FragmentResult[3];
            var rtBuf = new float[100];
            var intBuf = new float[100];

            int totalPoints = extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            Assert.That(totalPoints, Is.EqualTo(3)); // One data point per query

            Assert.That(results[0].TotalIntensity, Is.EqualTo(100.0f).Within(0.1f));
            Assert.That(results[0].QueryId, Is.EqualTo(10));

            Assert.That(results[1].TotalIntensity, Is.EqualTo(200.0f).Within(0.1f));
            Assert.That(results[1].QueryId, Is.EqualTo(20));

            Assert.That(results[2].TotalIntensity, Is.EqualTo(500.0f).Within(0.1f));
            Assert.That(results[2].QueryId, Is.EqualTo(30));

            // Buffer offsets should be sequential and non-overlapping
            Assert.That(results[1].RtBufferOffset, Is.EqualTo(results[0].RtBufferOffset + results[0].DataPointCount));
            Assert.That(results[2].RtBufferOffset, Is.EqualTo(results[1].RtBufferOffset + results[1].DataPointCount));
        }

        /// <summary>
        /// Verifies that a query targeting an m/z with no nearby peaks returns
        /// zero data points and zero total intensity.
        /// </summary>
        [Test]
        public void NoMatchingPeaks_ReturnsZeroDataPoints()
        {
            using var index = BuildIndex(new ScanDef(
                1, 5.0, 500.0, 25.0,
                new[] { 300.0, 400.0 }, new[] { 100.0, 200.0 }));

            using var extractor = new CpuFragmentExtractor(index);

            var queries = new FragmentQuery[]
            {
                // Target 800.0 — far from any peaks in the scan
                new FragmentQuery(800.0f, 10f, 0f, 10f, windowId: 0, queryId: 1)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[100];
            var intBuf = new float[100];

            extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            Assert.That(results[0].DataPointCount, Is.EqualTo(0));
            Assert.That(results[0].TotalIntensity, Is.EqualTo(0f));
        }

        /// <summary>
        /// Verifies that when multiple peaks fall within the ppm tolerance of a
        /// single query in one scan, their intensities are summed into one XIC point.
        /// This handles the case of overlapping fragment ions or isotope peaks.
        /// </summary>
        [Test]
        public void MultiplePeaksInTolerance_SumsIntensities()
        {
            // Two peaks very close together: 500.001 and 500.003 
            // At 10 ppm from 500.0, range is [499.995, 500.005] — both should match
            using var index = BuildIndex(new ScanDef(
                1, 5.0, 600.0, 200.0,
                new[] { 499.0, 500.001, 500.003, 501.0 },
                new[] { 10.0, 100.0, 150.0, 20.0 }));

            using var extractor = new CpuFragmentExtractor(index);

            var queries = new FragmentQuery[]
            {
                new FragmentQuery(500.0f, 10f, 0f, 10f, windowId: 0, queryId: 1)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[100];
            var intBuf = new float[100];

            extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            Assert.That(results[0].DataPointCount, Is.EqualTo(1)); // One scan = one XIC point
            Assert.That(results[0].TotalIntensity, Is.EqualTo(250.0f).Within(0.1f)); // 100 + 150
        }

        /// <summary>
        /// Verifies the LowerBound binary search helper returns correct indices
        /// for various target values including edges and out-of-range cases.
        /// </summary>
        [Test]
        public void LowerBound_ReturnsCorrectIndices()
        {
            float[] sorted = { 100f, 200f, 300f, 400f, 500f };

            // Exact match
            Assert.That(CpuFragmentExtractor.LowerBound(sorted, 300f), Is.EqualTo(2));

            // Between values — returns index of next higher
            Assert.That(CpuFragmentExtractor.LowerBound(sorted, 250f), Is.EqualTo(2));

            // Below all values
            Assert.That(CpuFragmentExtractor.LowerBound(sorted, 50f), Is.EqualTo(0));

            // Above all values
            Assert.That(CpuFragmentExtractor.LowerBound(sorted, 600f), Is.EqualTo(5));

            // Empty span
            Assert.That(CpuFragmentExtractor.LowerBound(ReadOnlySpan<float>.Empty, 100f), Is.EqualTo(0));
        }

        // ────────────────────────────────────────────────────────────────────
        //  Helpers for building test indices from compact scan definitions
        // ────────────────────────────────────────────────────────────────────

        /// <summary>Compact representation of a synthetic scan for test setup.</summary>
        private readonly struct ScanDef
        {
            public readonly int ScanNumber;
            public readonly double Rt;
            public readonly double IsolationMz;
            public readonly double IsolationWidth;
            public readonly double[] MzValues;
            public readonly double[] Intensities;

            public ScanDef(int scanNumber, double rt, double isolationMz, double isolationWidth,
                double[] mzValues, double[] intensities)
            {
                ScanNumber = scanNumber;
                Rt = rt;
                IsolationMz = isolationMz;
                IsolationWidth = isolationWidth;
                MzValues = mzValues;
                Intensities = intensities;
            }
        }

        /// <summary>Builds a DiaScanIndex from a single scan definition.</summary>
        private static DiaScanIndex BuildIndex(ScanDef def)
        {
            return BuildIndex(new[] { def });
        }

        /// <summary>Builds a DiaScanIndex from multiple scan definitions.</summary>
        private static DiaScanIndex BuildIndex(ScanDef[] defs)
        {
            var scans = new MsDataScan[defs.Length];
            for (int i = 0; i < defs.Length; i++)
            {
                var d = defs[i];
                var spectrum = new MzSpectrum(d.MzValues, d.Intensities, false);
                scans[i] = new MsDataScan(
                    massSpectrum: spectrum,
                    oneBasedScanNumber: d.ScanNumber,
                    msnOrder: 2,
                    isCentroid: true,
                    polarity: Polarity.Positive,
                    retentionTime: d.Rt,
                    scanWindowRange: new MzRange(100, 2000),
                    scanFilter: "FTMS",
                    mzAnalyzer: MZAnalyzerType.Orbitrap,
                    totalIonCurrent: d.Intensities.Sum(),
                    injectionTime: 50.0,
                    noiseData: null,
                    nativeId: $"scan={d.ScanNumber}",
                    isolationMZ: d.IsolationMz,
                    isolationWidth: d.IsolationWidth,
                    dissociationType: DissociationType.HCD);
            }
            return DiaScanIndexBuilder.Build(scans);
        }
    }
}
