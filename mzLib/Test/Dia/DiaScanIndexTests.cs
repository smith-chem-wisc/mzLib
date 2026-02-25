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
    /// Tests for the DIA SoA container (DiaScanIndex) and its builder (DiaScanIndexBuilder).
    /// All tests use synthetic MsDataScan objects with known values — no file I/O required.
    /// </summary>
    [TestFixture]
    public class DiaScanIndexTests
    {
        /// <summary>
        /// Verifies that the builder correctly packs peak data from multiple MS2 scans
        /// into contiguous float arrays, and that we can retrieve each scan's data 
        /// via offset/length spans without any data corruption or misalignment.
        /// </summary>
        [Test]
        public void Build_ThreeScansInTwoWindows_PacksPeaksContiguously()
        {
            // Arrange: 3 MS2 scans in 2 isolation windows
            //   Window A (center 500.0, width 25): scan 1 (RT 1.0), scan 3 (RT 3.0)
            //   Window B (center 700.0, width 25): scan 2 (RT 2.0)
            var scans = new MsDataScan[]
            {
                CreateMs2Scan(scanNumber: 1, rt: 1.0, isolationMz: 500.0, isolationWidth: 25.0,
                    mzValues: new[] { 200.0, 300.0, 400.0 },
                    intensities: new[] { 100.0, 200.0, 300.0 }),

                CreateMs2Scan(scanNumber: 2, rt: 2.0, isolationMz: 700.0, isolationWidth: 25.0,
                    mzValues: new[] { 250.0, 350.0 },
                    intensities: new[] { 150.0, 250.0 }),

                CreateMs2Scan(scanNumber: 3, rt: 3.0, isolationMz: 500.0, isolationWidth: 25.0,
                    mzValues: new[] { 201.0, 301.0 },
                    intensities: new[] { 110.0, 210.0 }),
            };

            // Act
            using var index = DiaScanIndexBuilder.Build(scans);

            // Assert: basic counts
            Assert.That(index.ScanCount, Is.EqualTo(3));
            Assert.That(index.TotalPeakCount, Is.EqualTo(7)); // 3 + 2 + 2
            Assert.That(index.WindowCount, Is.EqualTo(2));

            // Assert: scans sorted by (windowId, RT), so window A scans come first
            // Window A (500.0): scan1 at RT 1.0, scan3 at RT 3.0
            // Window B (700.0): scan2 at RT 2.0
            // Scan index 0 should be scan1, index 1 should be scan3, index 2 should be scan2

            // Check window A, first scan
            var mzSpan0 = index.GetScanMzSpan(0);
            Assert.That(mzSpan0.Length, Is.EqualTo(3));
            Assert.That(mzSpan0[0], Is.EqualTo(200.0f).Within(0.01f));
            Assert.That(mzSpan0[2], Is.EqualTo(400.0f).Within(0.01f));

            // Check window A, second scan
            var mzSpan1 = index.GetScanMzSpan(1);
            Assert.That(mzSpan1.Length, Is.EqualTo(2));
            Assert.That(mzSpan1[0], Is.EqualTo(201.0f).Within(0.01f));

            // Check window B scan
            var mzSpan2 = index.GetScanMzSpan(2);
            Assert.That(mzSpan2.Length, Is.EqualTo(2));
            Assert.That(mzSpan2[0], Is.EqualTo(250.0f).Within(0.01f));
        }

        /// <summary>
        /// Verifies that the window-to-scan range mapping correctly identifies which
        /// scan indices belong to each window, and that the ranges are non-overlapping.
        /// </summary>
        [Test]
        public void Build_WindowToScanRange_ReturnsCorrectRanges()
        {
            var scans = new MsDataScan[]
            {
                CreateMs2Scan(1, 1.0, 500.0, 25.0, new[] { 200.0 }, new[] { 100.0 }),
                CreateMs2Scan(2, 2.0, 700.0, 25.0, new[] { 300.0 }, new[] { 200.0 }),
                CreateMs2Scan(3, 3.0, 500.0, 25.0, new[] { 400.0 }, new[] { 300.0 }),
                CreateMs2Scan(4, 4.0, 700.0, 25.0, new[] { 500.0 }, new[] { 400.0 }),
            };

            using var index = DiaScanIndexBuilder.Build(scans);

            // Window 0 (500.0) should have 2 scans; Window 1 (700.0) should have 2 scans
            Assert.That(index.TryGetScanRangeForWindow(0, out int start0, out int count0), Is.True);
            Assert.That(count0, Is.EqualTo(2));

            Assert.That(index.TryGetScanRangeForWindow(1, out int start1, out int count1), Is.True);
            Assert.That(count1, Is.EqualTo(2));

            // Ranges should not overlap
            Assert.That(start1, Is.EqualTo(start0 + count0));
        }

        /// <summary>
        /// Verifies that within each window, scans are ordered by retention time (ascending).
        /// This is critical for efficient RT-windowed extraction.
        /// </summary>
        [Test]
        public void Build_ScansWithinWindow_AreSortedByRt()
        {
            // Input scans are NOT in RT order
            var scans = new MsDataScan[]
            {
                CreateMs2Scan(1, 5.0, 500.0, 25.0, new[] { 200.0 }, new[] { 100.0 }),
                CreateMs2Scan(2, 1.0, 500.0, 25.0, new[] { 300.0 }, new[] { 200.0 }),
                CreateMs2Scan(3, 3.0, 500.0, 25.0, new[] { 400.0 }, new[] { 300.0 }),
            };

            using var index = DiaScanIndexBuilder.Build(scans);

            index.TryGetScanRangeForWindow(0, out int start, out int count);
            Assert.That(count, Is.EqualTo(3));

            // RTs should be in ascending order: 1.0, 3.0, 5.0
            Assert.That(index.GetScanRt(start), Is.EqualTo(1.0f).Within(0.01f));
            Assert.That(index.GetScanRt(start + 1), Is.EqualTo(3.0f).Within(0.01f));
            Assert.That(index.GetScanRt(start + 2), Is.EqualTo(5.0f).Within(0.01f));
        }

        /// <summary>
        /// Verifies that MS1 scans, null scans, and MS2 scans without isolation info
        /// are all excluded from the index.
        /// </summary>
        [Test]
        public void Build_FiltersOutNonMs2AndInvalidScans()
        {
            var ms1Scan = CreateMs1Scan(scanNumber: 1, rt: 1.0,
                mzValues: new[] { 500.0, 600.0 }, intensities: new[] { 1000.0, 2000.0 });

            var validMs2 = CreateMs2Scan(2, 2.0, 500.0, 25.0,
                new[] { 200.0 }, new[] { 100.0 });

            // MS2 with no isolation m/z
            var invalidMs2 = CreateMs2Scan(3, 3.0, isolationMz: null, isolationWidth: null,
                new[] { 300.0 }, new[] { 200.0 });

            var scans = new MsDataScan[] { ms1Scan, null, validMs2, invalidMs2 };

            using var index = DiaScanIndexBuilder.Build(scans);

            Assert.That(index.ScanCount, Is.EqualTo(1)); // Only the valid MS2
            Assert.That(index.TotalPeakCount, Is.EqualTo(1));
        }

        /// <summary>
        /// Verifies that the builder returns a valid empty index when given no qualifying scans.
        /// </summary>
        [Test]
        public void Build_NoMs2Scans_ReturnsEmptyIndex()
        {
            var scans = new MsDataScan[]
            {
                CreateMs1Scan(1, 1.0, new[] { 500.0 }, new[] { 1000.0 })
            };

            using var index = DiaScanIndexBuilder.Build(scans);

            Assert.That(index.ScanCount, Is.EqualTo(0));
            Assert.That(index.TotalPeakCount, Is.EqualTo(0));
            Assert.That(index.WindowCount, Is.EqualTo(0));
        }

        /// <summary>
        /// Verifies that isolation window bounds are correctly reported.
        /// </summary>
        [Test]
        public void Build_WindowBounds_CorrectlyComputed()
        {
            var scans = new MsDataScan[]
            {
                CreateMs2Scan(1, 1.0, isolationMz: 500.0, isolationWidth: 26.0,
                    new[] { 200.0 }, new[] { 100.0 }),
            };

            using var index = DiaScanIndexBuilder.Build(scans);

            var bounds = index.GetWindowBounds(0);
            Assert.That(bounds.LowerBound, Is.EqualTo(487.0f).Within(0.01f)); // 500 - 13
            Assert.That(bounds.UpperBound, Is.EqualTo(513.0f).Within(0.01f)); // 500 + 13
        }

        /// <summary>
        /// Verifies that intensity data is correctly packed alongside m/z data,
        /// and that the double-to-float conversion preserves values at float precision.
        /// </summary>
        [Test]
        public void Build_IntensityData_CorrectlyPacked()
        {
            var scans = new MsDataScan[]
            {
                CreateMs2Scan(1, 1.0, 500.0, 25.0,
                    new[] { 200.123456, 300.654321 },
                    new[] { 1234.5, 6789.0 }),
            };

            using var index = DiaScanIndexBuilder.Build(scans);

            var mz = index.GetScanMzSpan(0);
            var intensity = index.GetScanIntensitySpan(0);

            Assert.That(mz[0], Is.EqualTo(200.123456f).Within(0.001f));
            Assert.That(mz[1], Is.EqualTo(300.654321f).Within(0.001f));
            Assert.That(intensity[0], Is.EqualTo(1234.5f).Within(0.1f));
            Assert.That(intensity[1], Is.EqualTo(6789.0f).Within(0.1f));
        }

        /// <summary>
        /// Verifies that the bulk memory accessors (AllMz, AllIntensity) expose the 
        /// full contiguous arrays for GPU transfer scenarios.
        /// </summary>
        [Test]
        public void BulkAccessors_ExposeFullContiguousArrays()
        {
            var scans = new MsDataScan[]
            {
                CreateMs2Scan(1, 1.0, 500.0, 25.0,
                    new[] { 100.0, 200.0 }, new[] { 10.0, 20.0 }),
                CreateMs2Scan(2, 2.0, 500.0, 25.0,
                    new[] { 300.0 }, new[] { 30.0 }),
            };

            using var index = DiaScanIndexBuilder.Build(scans);

            Assert.That(index.AllMz.Length, Is.EqualTo(3));
            Assert.That(index.AllIntensity.Length, Is.EqualTo(3));
        }

        // ────────────────────────────────────────────────────────────────────
        //  Helper methods for creating synthetic scans
        // ────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Creates a synthetic MS2 MsDataScan with the given properties.
        /// Uses MzSpectrum and the standard MsDataScan constructor from mzLib.
        /// </summary>
        private static MsDataScan CreateMs2Scan(
            int scanNumber, double rt,
            double? isolationMz, double? isolationWidth,
            double[] mzValues, double[] intensities)
        {
            var spectrum = new MzSpectrum(mzValues, intensities, false);
            return new MsDataScan(
                massSpectrum: spectrum,
                oneBasedScanNumber: scanNumber,
                msnOrder: 2,
                isCentroid: true,
                polarity: Polarity.Positive,
                retentionTime: rt,
                scanWindowRange: new MzRange(100, 2000),
                scanFilter: "FTMS",
                mzAnalyzer: MZAnalyzerType.Orbitrap,
                totalIonCurrent: intensities.Sum(),
                injectionTime: 50.0,
                noiseData: null,
                nativeId: $"scan={scanNumber}",
                isolationMZ: isolationMz,
                isolationWidth: isolationWidth,
                dissociationType: DissociationType.HCD);
        }

        /// <summary>Convenience overload for MS2 scans with valid isolation info.</summary>
        private static MsDataScan CreateMs2Scan(
            int scanNumber, double rt,
            double isolationMz, double isolationWidth,
            double[] mzValues, double[] intensities)
        {
            return CreateMs2Scan(scanNumber, rt, (double?)isolationMz, (double?)isolationWidth,
                mzValues, intensities);
        }

        /// <summary>Creates a synthetic MS1 scan (for testing that MS1s are filtered out).</summary>
        private static MsDataScan CreateMs1Scan(
            int scanNumber, double rt,
            double[] mzValues, double[] intensities)
        {
            var spectrum = new MzSpectrum(mzValues, intensities, false);
            return new MsDataScan(
                massSpectrum: spectrum,
                oneBasedScanNumber: scanNumber,
                msnOrder: 1,
                isCentroid: true,
                polarity: Polarity.Positive,
                retentionTime: rt,
                scanWindowRange: new MzRange(100, 2000),
                scanFilter: "FTMS",
                mzAnalyzer: MZAnalyzerType.Orbitrap,
                totalIonCurrent: intensities.Sum(),
                injectionTime: 50.0,
                noiseData: null,
                nativeId: $"scan={scanNumber}");
        }
    }
}
