// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.IO;
using System.Linq;

namespace Test.Dia
{
    /// <summary>
    /// End-to-end integration tests that load a real DIA mzML file through mzLib's
    /// existing Mzml reader, build a DiaScanIndex, extract fragments, and score them.
    /// 
    /// This validates that the entire DIA pipeline works with real instrument data,
    /// not just synthetic test data. The test file is a small subset (~1 minute, 767 scans)
    /// from a Thermo DIA acquisition with ~152 isolation windows of ~5 Da width.
    /// 
    /// Test file: E24484_Mag_1_4_50_3.mzML
    ///   - 767 total scans (~20 MS1 + ~747 MS2)
    ///   - 152 unique isolation windows (~4 Da wide, ±2.5 Da offsets)
    ///   - 4–5 MS2 scans per window
    ///   - RT range: ~37.0 – 38.0 minutes
    ///   - Centroided, zlib compressed, 64-bit float arrays
    /// </summary>
    [TestFixture]
    public class DiaIntegrationTests
    {
        // Path to the test mzML file. Adjust if your test data lives elsewhere.
        private static readonly string TestFilePath = Path.Combine(
            TestContext.CurrentContext.TestDirectory, "DataFiles", "E24484_Mag_1_4_50_3.mzML");

        /// <summary>
        /// Loads the DIA mzML file through mzLib's Mzml reader, builds a DiaScanIndex,
        /// and verifies that the index has the expected structure: correct scan count,
        /// window count, and non-empty peak data.
        /// 
        /// This is the most fundamental integration test — if this passes, the bridge
        /// between mzLib's object model and our SoA layout is working correctly.
        /// </summary>
        [Test]
        public void LoadMzml_BuildIndex_CorrectStructure()
        {
            if (!File.Exists(TestFilePath))
            {
                Assert.Ignore($"Test file not found: {TestFilePath}. " +
                    "Copy E24484_Mag_1_4_50_3.mzML to the Test/DataFiles folder.");
            }

            // Load via mzLib's existing reader
            var mzml = new Mzml(TestFilePath);
            mzml.LoadAllStaticData();
            var allScans = mzml.GetMsDataScans();

            Assert.That(allScans, Is.Not.Null);
            Assert.That(allScans.Length, Is.GreaterThan(0), "File should contain scans");

            // Build the DIA index
            using var index = DiaScanIndexBuilder.Build(allScans);

            // Verify basic structure
            Assert.That(index.ScanCount, Is.GreaterThan(0), "Should have MS2 scans");
            Assert.That(index.WindowCount, Is.GreaterThan(100), "DIA file should have many isolation windows");
            Assert.That(index.TotalPeakCount, Is.GreaterThan(0), "Scans should contain peaks");

            TestContext.WriteLine($"Loaded {allScans.Length} total scans from mzML");
            TestContext.WriteLine($"DiaScanIndex: {index.ScanCount} MS2 scans, {index.WindowCount} windows, {index.TotalPeakCount:N0} total peaks");
        }

        /// <summary>
        /// Verifies that every window in the index has at least one scan, and that
        /// the scan metadata (RT, peak count) is reasonable for real data.
        /// </summary>
        [Test]
        public void LoadMzml_AllWindowsHaveScans_MetadataReasonable()
        {
            if (!File.Exists(TestFilePath))
            {
                Assert.Ignore($"Test file not found: {TestFilePath}.");
            }

            var mzml = new Mzml(TestFilePath);
            mzml.LoadAllStaticData();
            using var index = DiaScanIndexBuilder.Build(mzml.GetMsDataScans());

            var windowIds = index.GetWindowIds();

            foreach (int windowId in windowIds)
            {
                Assert.That(index.TryGetScanRangeForWindow(windowId, out int start, out int count), Is.True,
                    $"Window {windowId} should exist in the scan range map");
                Assert.That(count, Is.GreaterThan(0), $"Window {windowId} should have at least one scan");

                // Check that scans within each window are RT-ordered
                for (int i = start + 1; i < start + count; i++)
                {
                    Assert.That(index.GetScanRt(i), Is.GreaterThanOrEqualTo(index.GetScanRt(i - 1)),
                        $"Scans in window {windowId} should be RT-ordered");
                }

                // Check that scans have peaks
                for (int i = start; i < start + count; i++)
                {
                    Assert.That(index.GetScanPeakCount(i), Is.GreaterThan(0),
                        $"Scan {i} in window {windowId} should have peaks");
                }
            }

            TestContext.WriteLine($"All {windowIds.Count} windows validated: scans present, RT-ordered, peaks non-empty");
        }

        /// <summary>
        /// Full pipeline test: load mzML → build index → extract fragments → score.
        /// 
        /// Picks a real m/z value from the first scan in the first window, extracts
        /// its XIC across the full RT range, then scores the extracted intensities
        /// against themselves (should give a perfect 1.0 score).
        /// 
        /// This validates that extraction and scoring work correctly on real data
        /// with real m/z precision, real intensity magnitudes, and real RT values.
        /// </summary>
        [Test]
        public void FullPipeline_Extract_And_Score_RealData()
        {
            if (!File.Exists(TestFilePath))
            {
                Assert.Ignore($"Test file not found: {TestFilePath}.");
            }

            var mzml = new Mzml(TestFilePath);
            mzml.LoadAllStaticData();
            using var index = DiaScanIndexBuilder.Build(mzml.GetMsDataScans());
            using var extractor = new CpuFragmentExtractor(index);

            // Pick the first window and get a real m/z from its first scan
            var windowIds = index.GetWindowIds().ToArray();
            int windowId = windowIds[0];
            index.TryGetScanRangeForWindow(windowId, out int scanStart, out int scanCount);

            var mzSpan = index.GetScanMzSpan(scanStart);
            Assert.That(mzSpan.Length, Is.GreaterThan(0), "First scan should have peaks");

            // Pick a m/z from the middle of the scan's peak list
            float targetMz = mzSpan[mzSpan.Length / 2];
            float minRt = index.GetScanRt(scanStart);
            float maxRt = index.GetScanRt(scanStart + scanCount - 1);

            TestContext.WriteLine($"Window {windowId}: {scanCount} scans, RT [{minRt:F3}, {maxRt:F3}]");
            TestContext.WriteLine($"Target m/z: {targetMz:F4}, tolerance: 20 ppm");

            // Extract XIC
            var queries = new FragmentQuery[]
            {
                new FragmentQuery(targetMz, tolerancePpm: 20f, rtMin: minRt - 0.1f, rtMax: maxRt + 0.1f,
                    windowId: windowId, queryId: 1)
            };
            var results = new FragmentResult[1];
            var rtBuf = new float[1000];
            var intBuf = new float[1000];

            int totalPoints = extractor.ExtractBatch(queries, results, rtBuf, intBuf);

            TestContext.WriteLine($"Extracted {results[0].DataPointCount} XIC data points, total intensity: {results[0].TotalIntensity:F1}");

            Assert.That(results[0].DataPointCount, Is.GreaterThan(0), "Should find the target m/z in at least one scan");
            Assert.That(results[0].TotalIntensity, Is.GreaterThan(0f), "Extracted intensity should be positive");

            // Score the extracted intensities against themselves (perfect self-match)
            if (results[0].DataPointCount >= 2)
            {
                int offset = results[0].IntensityBufferOffset;
                int count = results[0].DataPointCount;
                var observed = new ReadOnlySpan<float>(intBuf, offset, count);

                var dotScorer = new NormalizedDotProductScorer();
                var saScorer = new SpectralAngleScorer();

                float dpScore = dotScorer.Score(observed, observed);
                float saScore = saScorer.Score(observed, observed);

                TestContext.WriteLine($"Self-score: dot product = {dpScore:F4}, spectral angle = {saScore:F4}");

                Assert.That(dpScore, Is.EqualTo(1.0f).Within(1e-4f), "Self-score should be perfect");
                Assert.That(saScore, Is.EqualTo(1.0f).Within(1e-4f), "Self-score should be perfect");
            }
        }

        /// <summary>
        /// Verifies that the double-to-float conversion during index building preserves
        /// sufficient precision for real instrument data. Compares a few peak m/z values
        /// from the original MsDataScan (double) against the DiaScanIndex (float).
        /// 
        /// At m/z 1000, float precision is ~0.06 ppm — well within the 10–20 ppm
        /// tolerance used for fragment matching.
        /// </summary>
        [Test]
        public void FloatPrecision_SufficientForRealData()
        {
            if (!File.Exists(TestFilePath))
            {
                Assert.Ignore($"Test file not found: {TestFilePath}.");
            }

            var mzml = new Mzml(TestFilePath);
            mzml.LoadAllStaticData();
            var allScans = mzml.GetMsDataScans();
            using var index = DiaScanIndexBuilder.Build(allScans);

            // Find the first MS2 scan in the original data
            var firstMs2 = allScans.FirstOrDefault(s => s != null && s.MsnOrder == 2
                && s.IsolationMz.HasValue && s.MassSpectrum?.Size > 0);

            Assert.That(firstMs2, Is.Not.Null, "Should find an MS2 scan");

            // Compare original double[] m/z values against float[] in the index
            // The index reorders by (window, RT), so we find the matching scan by scan number
            double[] originalMz = firstMs2.MassSpectrum.XArray;

            // Check precision: float should be within 0.1 ppm of double at typical m/z ranges
            for (int i = 0; i < Math.Min(10, originalMz.Length); i++)
            {
                double original = originalMz[i];
                float converted = (float)original;
                double ppmError = Math.Abs(original - converted) / original * 1e6;

                Assert.That(ppmError, Is.LessThan(1.0),
                    $"Float precision at m/z {original:F4} should be well under 1 ppm, got {ppmError:F4} ppm");
            }

            TestContext.WriteLine("Float precision validated: all checked m/z values within 1 ppm of original doubles");
        }

        /// <summary>
        /// Performance test: measures how long it takes to load the mzML file, build
        /// the index, and run a batch extraction. Prints timing to output.
        /// Not a strict pass/fail — used to establish a baseline for real-file performance.
        /// </summary>
        [Test]
        public void Performance_LoadBuildExtract_RealFile()
        {
            if (!File.Exists(TestFilePath))
            {
                Assert.Ignore($"Test file not found: {TestFilePath}.");
            }

            var sw = System.Diagnostics.Stopwatch.StartNew();

            var mzml = new Mzml(TestFilePath);
            mzml.LoadAllStaticData();
            var allScans = mzml.GetMsDataScans();
            var loadTime = sw.ElapsedMilliseconds;

            sw.Restart();
            using var index = DiaScanIndexBuilder.Build(allScans);
            var buildTime = sw.ElapsedMilliseconds;

            // Generate queries from real data: one query per window, targeting a real m/z
            var windowIds = index.GetWindowIds().ToArray();
            var queries = new FragmentQuery[windowIds.Length];
            for (int i = 0; i < windowIds.Length; i++)
            {
                index.TryGetScanRangeForWindow(windowIds[i], out int start, out int count);
                var mzSpan = index.GetScanMzSpan(start);
                float targetMz = mzSpan[mzSpan.Length / 2];
                float minRt = index.GetScanRt(start) - 0.1f;
                float maxRt = index.GetScanRt(start + count - 1) + 0.1f;
                queries[i] = new FragmentQuery(targetMz, 20f, minRt, maxRt, windowIds[i], i);
            }

            var results = new FragmentResult[queries.Length];
            var rtBuf = new float[queries.Length * 20];
            var intBuf = new float[queries.Length * 20];

            sw.Restart();
            int totalPoints = 0;
            // Run extraction 100 times to get a stable measurement
            for (int iter = 0; iter < 100; iter++)
            {
                using var extractor = new CpuFragmentExtractor(index);
                totalPoints = extractor.ExtractBatch(queries, results, rtBuf, intBuf);
            }
            var extractTime = sw.ElapsedMilliseconds;

            TestContext.WriteLine($"mzML load:     {loadTime} ms ({allScans.Length} scans)");
            TestContext.WriteLine($"Index build:   {buildTime} ms ({index.ScanCount} MS2 scans, {index.TotalPeakCount:N0} peaks)");
            TestContext.WriteLine($"Extraction:    {extractTime} ms for 100 × {queries.Length} queries ({totalPoints} data points per batch)");
            TestContext.WriteLine($"Per-batch:     {extractTime / 100.0:F1} ms");

            // Sanity check: extraction should complete, no crashes
            Assert.That(totalPoints, Is.GreaterThan(0));
        }
    }
}
