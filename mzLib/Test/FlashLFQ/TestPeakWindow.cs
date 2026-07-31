using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;

namespace Test.FlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestPeakWindow
    {
        private static string SlicedMzmlPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "FlashLFQ", "TestData", @"sliced-mzml.mzml");

        /// <summary>
        /// Five MS1 scans, each holding the same four peaks: three clustered between 500 and 501 m/z and one far
        /// away at 600 m/z. Retention times run 1.0 to 1.4 minutes.
        /// </summary>
        private static MsDataScan[] BuildTestScans()
        {
            double[] mzs = { 500.0, 500.5, 501.0, 600.0 };
            var scans = new MsDataScan[5];

            for (int s = 0; s < scans.Length; s++)
            {
                double[] intensities = mzs.Select((_, i) => 1e5 * (i + 1) * (s + 1)).ToArray();
                scans[s] = new MsDataScan(
                    massSpectrum: new MzSpectrum(mzs.ToArray(), intensities, false), oneBasedScanNumber: s + 1,
                    msnOrder: 1, isCentroid: true, polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0,
                    scanWindowRange: new MzRange(400, 1600), scanFilter: "f", mzAnalyzer: MZAnalyzerType.Orbitrap,
                    totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            return scans;
        }

        private static ChromatographicPeak QuantifySlicedMzml(SpectraFileInfo mzml, out FlashLfqResults results)
        {
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2,
                new List<ProteinGroup> { pg });

            results = new FlashLfqEngine(new FlashLfqParameters { MaxThreads = 1, Silent = true },
                new List<Identification> { id }).Run();

            return results.Peaks[mzml].Single();
        }

        #region IndexingEngine range queries

        [Test]
        public static void TestGetPeaksInRangeByScanIndex()
        {
            PeakIndexingEngine engine = PeakIndexingEngine.InitializeIndexingEngine(BuildTestScans());

            // the three clustered peaks, in the middle three scans
            var peaks = engine.GetPeaksInRange(500.0, 501.0, 1, 3);
            Assert.AreEqual(9, peaks.Count);
            Assert.IsTrue(peaks.All(p => p.M >= 500.0 && p.M <= 501.0));
            Assert.IsTrue(peaks.All(p => p.ZeroBasedScanIndex >= 1 && p.ZeroBasedScanIndex <= 3));

            // peaks come back ordered by scan index, then m/z
            for (int i = 1; i < peaks.Count; i++)
            {
                Assert.IsTrue(peaks[i].ZeroBasedScanIndex > peaks[i - 1].ZeroBasedScanIndex
                    || (peaks[i].ZeroBasedScanIndex == peaks[i - 1].ZeroBasedScanIndex && peaks[i].M >= peaks[i - 1].M));
            }

            // bounds are inclusive on both ends
            Assert.AreEqual(5, engine.GetPeaksInRange(500.0, 500.0, 0, 4).Count);
            Assert.AreEqual(20, engine.GetPeaksInRange(400.0, 1600.0, 0, 4).Count);

            // an empty region, an inverted m/z range, and an inverted scan range all come back empty rather than throwing
            Assert.AreEqual(0, engine.GetPeaksInRange(550.0, 560.0, 0, 4).Count);
            Assert.AreEqual(0, engine.GetPeaksInRange(501.0, 500.0, 0, 4).Count);
            Assert.AreEqual(0, engine.GetPeaksInRange(500.0, 501.0, 3, 1).Count);

            // scan indices outside the file don't run off the end of the index
            Assert.AreEqual(20, engine.GetPeaksInRange(0.0, 5000.0, -10, 500).Count);
        }

        [Test]
        public static void TestGetPeaksInRangeByRetentionTime()
        {
            PeakIndexingEngine engine = PeakIndexingEngine.InitializeIndexingEngine(BuildTestScans());

            var peaks = engine.GetPeaksInRange(500.0, 501.0, 1.1, 1.3);
            Assert.AreEqual(9, peaks.Count);
            Assert.IsTrue(peaks.All(p => p.RetentionTime >= 1.1f && p.RetentionTime <= 1.3f));

            // a retention time range containing no scans
            Assert.AreEqual(0, engine.GetPeaksInRange(500.0, 501.0, 5.0, 6.0).Count);
        }

        [Test]
        public static void TestGetPeaksInRangeThrowsBeforeIndexing()
        {
            PeakIndexingEngine engine = PeakIndexingEngine.InitializeIndexingEngine(BuildTestScans());
            engine.ClearIndex();

            Assert.Throws<MzLibException>(() => engine.GetPeaksInRange(500.0, 501.0, 0, 4));
        }

        #endregion

        #region Scan number mapping

        /// <summary>
        /// The whole reason the map exists: indexed peaks are numbered against an MS1-only array, so a zero based
        /// index is not a scan number and addressing the data file with one lands on the wrong scan.
        /// </summary>
        [Test]
        public static void TestMs1ScanNumberMapIsNotTheZeroBasedIndex()
        {
            MsDataFile dataFile = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            dataFile.LoadAllStaticData();

            int[] map = PeakWindow.BuildMs1ScanNumberMap(dataFile);
            Assert.IsTrue(map.Length > 0);

            // every mapped scan really is an MS1 scan, and the numbers ascend
            for (int i = 0; i < map.Length; i++)
            {
                Assert.AreEqual(1, dataFile.GetOneBasedScan(map[i]).MsnOrder);
                if (i > 0)
                {
                    Assert.IsTrue(map[i] > map[i - 1]);
                }
            }

            // the file interleaves MS2 scans, so the map is not simply index + 1. Were it identity here, this
            // fixture would not be able to catch a regression back to zero based addressing.
            Assert.IsTrue(map.Where((scanNumber, i) => scanNumber != i + 1).Any());

            // and the map covers exactly the MS1 scans of the file
            Assert.AreEqual(dataFile.GetAllScansList().Count(s => s.MsnOrder == 1), map.Length);
        }

        #endregion

        #region Window extraction

        /// <summary>
        /// Runs FlashLFQ over a sliced mzML, then checks that the extracted window spans the scans the peak was
        /// observed in and the expected m/z range, and that it picks up peaks that are not part of the quantified
        /// precursor.
        /// </summary>
        [Test]
        public static void TestPeakWindowExtraction()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            ChromatographicPeak peak = QuantifySlicedMzml(mzml, out _);
            Assert.IsTrue(peak.IsotopicEnvelopes.Any());

            MsDataFile dataFile = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            dataFile.LoadAllStaticData();
            int[] map = PeakWindow.BuildMs1ScanNumberMap(dataFile);

            PeakWindow window = PeakWindow.Create(peak, dataFile, map, mzExpansion: 2.0);

            // the window is bounded by the peak's observed m/z values, expanded by mzExpansion on either side
            Assert.AreEqual(peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.M) - 2.0, window.MinMz, 1e-6);
            Assert.AreEqual(peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.M) + 2.0, window.MaxMz, 1e-6);

            // and it covers exactly the scans the peak was observed in
            int firstMs1Index = peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex);
            int lastMs1Index = peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex);
            Assert.AreEqual(lastMs1Index - firstMs1Index + 1, window.Scans.Length);
            Assert.AreEqual(map[firstMs1Index], window.Scans.First().OneBasedScanNumber);
            Assert.AreEqual(map[lastMs1Index], window.Scans.Last().OneBasedScanNumber);

            // retention times are the scans' own, taken from the file
            Assert.AreEqual(dataFile.GetOneBasedScan(map[firstMs1Index]).RetentionTime, window.MinRetentionTime);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[lastMs1Index]).RetentionTime, window.MaxRetentionTime);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[peak.Apex.IndexedPeak.ZeroBasedScanIndex]).RetentionTime,
                window.ApexRetentionTime);

            // each scan's slice is the m/z window of that scan's spectrum, parallel arrays and all
            Assert.AreEqual(window.Scans.Sum(s => s.Mz.Length), window.PeakCount);
            foreach (PeakWindowScan scan in window.Scans)
            {
                MsDataScan sourceScan = dataFile.GetOneBasedScan(scan.OneBasedScanNumber);
                Assert.AreEqual(sourceScan.RetentionTime, scan.RetentionTime);
                Assert.AreEqual(scan.Mz.Length, scan.Intensity.Length);
                Assert.AreEqual(scan.Mz.Length, scan.IsPeakfindingPeak.Length);
                Assert.IsTrue(scan.Mz.All(mz => mz >= window.MinMz && mz <= window.MaxMz));

                // the slice matches a brute force filter of the whole spectrum, so the binary search found it all
                var expected = sourceScan.MassSpectrum.XArray
                    .Where(mz => mz >= window.MinMz && mz <= window.MaxMz).ToArray();
                CollectionAssert.AreEqual(expected, scan.Mz);
            }

            // every peak FlashLFQ used to trace the feature is flagged, and nothing else is
            Assert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.IndexedPeak).Distinct().Count(),
                window.Scans.Sum(s => s.IsPeakfindingPeak.Count(f => f)));

            // and the window holds far more than those, which is the whole point of the output
            Assert.IsTrue(window.PeakCount > window.Scans.Sum(s => s.IsPeakfindingPeak.Count(f => f)));

            // a wider window is a superset of a narrower one
            PeakWindow narrowWindow = PeakWindow.Create(peak, dataFile, map, mzExpansion: 0.5);
            Assert.IsTrue(narrowWindow.PeakCount <= window.PeakCount);
            Assert.AreEqual(window.Scans.Length, narrowWindow.Scans.Length);
            var wideMzs = window.Scans.SelectMany(s => s.Mz).ToHashSet();
            Assert.IsTrue(narrowWindow.Scans.SelectMany(s => s.Mz).All(mz => wideMzs.Contains(mz)));

            // peaks with no isotopic envelopes have no window, and neither do missing arguments
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var emptyPeak = new ChromatographicPeak(new Identification(mzml, "PEPTIDE", "PEPTIDE", 799.36, 94.1, 2,
                new List<ProteinGroup> { pg }), mzml);
            Assert.IsNull(PeakWindow.Create(emptyPeak, dataFile, map));
            Assert.IsNull(PeakWindow.Create(null, dataFile, map));
            Assert.IsNull(PeakWindow.Create(peak, null, map));
            Assert.IsNull(PeakWindow.Create(peak, dataFile, null));
        }

        /// <summary>
        /// Reading the scans directly has to produce the same peaks the indexing engine would have. The comparison
        /// stays clear of the window edges: the index stores m/z narrowed to float, so a peak sitting within float
        /// rounding of a boundary can legitimately fall on either side of it depending on which path is taken.
        /// </summary>
        [Test]
        public static void TestPeakWindowMatchesIndexingEngine()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            ChromatographicPeak peak = QuantifySlicedMzml(mzml, out _);

            MsDataFile dataFile = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            dataFile.LoadAllStaticData();
            int[] map = PeakWindow.BuildMs1ScanNumberMap(dataFile);
            PeakWindow window = PeakWindow.Create(peak, dataFile, map, mzExpansion: 0.5);

            PeakIndexingEngine engine = PeakIndexingEngine.InitializeIndexingEngine(mzml);
            var indexedPeaks = engine.GetPeaksInRange(window.MinMz, window.MaxMz,
                peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex),
                peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex));

            const double edgeMargin = 1e-3;
            var fromFile = window.Scans
                .SelectMany(s => s.Mz.Select(mz => (s.OneBasedScanNumber, Mz: (float)mz)))
                .Where(p => p.Mz > window.MinMz + edgeMargin && p.Mz < window.MaxMz - edgeMargin)
                .OrderBy(p => p.OneBasedScanNumber).ThenBy(p => p.Mz)
                .ToList();
            var fromIndex = indexedPeaks
                .Select(p => (OneBasedScanNumber: map[p.ZeroBasedScanIndex], Mz: p.M))
                .Where(p => p.Mz > window.MinMz + edgeMargin && p.Mz < window.MaxMz - edgeMargin)
                .OrderBy(p => p.OneBasedScanNumber).ThenBy(p => p.Mz)
                .ToList();

            Assert.IsTrue(fromFile.Count > 0);
            CollectionAssert.AreEqual(fromIndex, fromFile);

            engine.ClearIndex();
        }

        #endregion

        #region Output file

        [Test]
        public static void TestWritePeakWindows()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            ChromatographicPeak peak = QuantifySlicedMzml(mzml, out FlashLfqResults results);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindows.tsv");
            results.WritePeakWindows(outputPath, silent: true);

            Assert.IsTrue(File.Exists(outputPath));
            string[] lines = File.ReadAllLines(outputPath);
            Assert.AreEqual(PeakWindow.TabSeparatedHeader, lines[0]);
            Assert.IsTrue(lines.Length > 1);

            string[] headerColumns = lines[0].Split('\t');
            int mzColumn = Array.IndexOf(headerColumns, "MZ");
            int intensityColumn = Array.IndexOf(headerColumns, "Intensity");
            int scanRtColumn = Array.IndexOf(headerColumns, "Scan Retention Time");
            int scanNumberColumn = Array.IndexOf(headerColumns, "Scan Number");
            int minMzColumn = Array.IndexOf(headerColumns, "Window Min MZ");
            int maxMzColumn = Array.IndexOf(headerColumns, "Window Max MZ");
            int rtStartColumn = Array.IndexOf(headerColumns, "Peak RT Start");
            int rtEndColumn = Array.IndexOf(headerColumns, "Peak RT End");
            int peakfindingColumn = Array.IndexOf(headerColumns, "Is Peakfinding Peak");
            int fullSequenceColumn = Array.IndexOf(headerColumns, "Full Sequence");

            MsDataFile dataFile = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            dataFile.LoadAllStaticData();

            int peakfindingRows = 0;
            foreach (string line in lines.Skip(1))
            {
                string[] cells = line.Split('\t');
                Assert.AreEqual(headerColumns.Length, cells.Length);
                Assert.AreEqual("EGFQVADGPLYR", cells[fullSequenceColumn]);

                double mz = double.Parse(cells[mzColumn], CultureInfo.InvariantCulture);
                double rt = double.Parse(cells[scanRtColumn], CultureInfo.InvariantCulture);
                Assert.IsTrue(mz >= double.Parse(cells[minMzColumn], CultureInfo.InvariantCulture));
                Assert.IsTrue(mz <= double.Parse(cells[maxMzColumn], CultureInfo.InvariantCulture));
                Assert.IsTrue(rt >= double.Parse(cells[rtStartColumn], CultureInfo.InvariantCulture));
                Assert.IsTrue(rt <= double.Parse(cells[rtEndColumn], CultureInfo.InvariantCulture));
                Assert.IsTrue(double.Parse(cells[intensityColumn], CultureInfo.InvariantCulture) > 0);

                // the scan number addresses the raw file directly, and names an MS1 scan at the row's retention time
                int scanNumber = int.Parse(cells[scanNumberColumn], CultureInfo.InvariantCulture);
                MsDataScan sourceScan = dataFile.GetOneBasedScan(scanNumber);
                Assert.AreEqual(1, sourceScan.MsnOrder);
                Assert.AreEqual(sourceScan.RetentionTime, rt);

                if (bool.Parse(cells[peakfindingColumn]))
                {
                    peakfindingRows++;
                }
            }

            // the rows FlashLFQ actually quantified are flagged, and they are a strict subset of what was written
            Assert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.IndexedPeak).Distinct().Count(), peakfindingRows);
            Assert.IsTrue(lines.Length - 1 > peakfindingRows);

            File.Delete(outputPath);
        }

        [Test]
        public static void TestWritePeakWindowsNullPathIsANoOp()
        {
            var results = new FlashLfqResults(new List<SpectraFileInfo>(), new List<Identification>());
            Assert.DoesNotThrow(() => results.WritePeakWindows(null, silent: true));
        }

        #endregion
    }
}
