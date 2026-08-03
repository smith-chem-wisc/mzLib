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
using Test.FileReadingTests;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using StringAssert = NUnit.Framework.Legacy.StringAssert;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;

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

        /// <summary>
        /// Writes a synthetic mzML holding one MS1 scan per supplied peak list, and returns a SpectraFileInfo for it.
        /// </summary>
        private static SpectraFileInfo WriteSyntheticMzml(string fileName, IReadOnlyList<double[]> mzPerScan)
        {
            var scans = new MsDataScan[mzPerScan.Count];
            for (int s = 0; s < mzPerScan.Count; s++)
            {
                double[] intensities = mzPerScan[s].Select((_, i) => 1e5 * (i + 1)).ToArray();
                scans[s] = new MsDataScan(new MzSpectrum(mzPerScan[s].ToArray(), intensities, false), s + 1, 1, true,
                    Polarity.Positive, 1.0 + s * 0.1, new MzRange(50, 1600), "f", MZAnalyzerType.Orbitrap,
                    intensities.Sum(), 1.0, null, "scan=" + (s + 1));
            }

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, fileName);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans), path, false);
            return new SpectraFileInfo(path, "a", 0, 0, 0);
        }

        /// <summary>
        /// Builds a ChromatographicPeak directly, without running the engine, so that a test can put its isotopic
        /// envelopes in exactly the scans it wants - including scans that do not exist.
        /// </summary>
        private static ChromatographicPeak BuildPeak(SpectraFileInfo file, string sequence,
            params (double Mz, int ZeroBasedScanIndex, double RetentionTime)[] envelopePeaks)
        {
            return BuildMultiChargePeak(file, sequence,
                envelopePeaks.Select(p => (p.Mz, 2, p.ZeroBasedScanIndex, p.RetentionTime)).ToArray());
        }

        /// <summary>
        /// As BuildPeak, but each envelope carries its own charge state, so a test can build the peak FlashLFQ
        /// produces for a peptide traced at more than one charge.
        /// </summary>
        private static ChromatographicPeak BuildMultiChargePeak(SpectraFileInfo file, string sequence,
            params (double Mz, int ChargeState, int ZeroBasedScanIndex, double RetentionTime)[] envelopePeaks)
        {
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var id = new Identification(file, sequence, sequence, 1350.65681, 1.0, 2, new List<ProteinGroup> { pg });
            var peak = new ChromatographicPeak(id, file);

            foreach (var envelopePeak in envelopePeaks)
            {
                peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(
                    new IndexedMassSpectralPeak(envelopePeak.Mz, 1e6, envelopePeak.ZeroBasedScanIndex,
                        envelopePeak.RetentionTime), envelopePeak.ChargeState, 1e6, 1.0));
            }
            peak.CalculateIntensityForThisFeature(false);
            return peak;
        }

        private static MsDataFile LoadDataFile(SpectraFileInfo file)
        {
            MsDataFile dataFile = MsDataFileReader.GetDataFile(file.FullFilePathWithExtension);
            dataFile.LoadAllStaticData();
            return dataFile;
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

            int chargeState = peak.Apex.ChargeState;
            var envelopes = peak.IsotopicEnvelopes.Where(e => e.ChargeState == chargeState).ToList();
            PeakWindow window = PeakWindow.Create(peak, dataFile, map, chargeState, mzExpansion: 2.0);
            Assert.AreEqual(chargeState, window.ChargeState);

            // the window is bounded by the charge state's observed m/z values, expanded by mzExpansion on either side
            Assert.AreEqual(envelopes.Min(e => e.IndexedPeak.M) - 2.0, window.MinMz, 1e-6);
            Assert.AreEqual(envelopes.Max(e => e.IndexedPeak.M) + 2.0, window.MaxMz, 1e-6);

            // and it covers exactly the scans that charge state was observed in
            int firstMs1Index = envelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex);
            int lastMs1Index = envelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex);
            Assert.AreEqual(lastMs1Index - firstMs1Index + 1, window.Scans.Length);
            Assert.AreEqual(map[firstMs1Index], window.Scans.First().OneBasedScanNumber);
            Assert.AreEqual(map[lastMs1Index], window.Scans.Last().OneBasedScanNumber);

            // retention times are the scans' own, taken from the file
            Assert.AreEqual(dataFile.GetOneBasedScan(map[firstMs1Index]).RetentionTime, window.MinRetentionTime);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[lastMs1Index]).RetentionTime, window.MaxRetentionTime);
            Assert.AreEqual(dataFile.GetOneBasedScan(
                    map[envelopes.MaxBy(e => e.Intensity).IndexedPeak.ZeroBasedScanIndex]).RetentionTime,
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

            // every peak FlashLFQ used to trace this charge state is flagged, and nothing else is
            Assert.AreEqual(envelopes.Select(e => e.IndexedPeak).Distinct().Count(),
                window.Scans.Sum(s => s.IsPeakfindingPeak.Count(f => f)));

            // and the window holds far more than those, which is the whole point of the output
            Assert.IsTrue(window.PeakCount > window.Scans.Sum(s => s.IsPeakfindingPeak.Count(f => f)));

            // a wider window is a superset of a narrower one
            PeakWindow narrowWindow = PeakWindow.Create(peak, dataFile, map, chargeState, mzExpansion: 0.5);
            Assert.IsTrue(narrowWindow.PeakCount <= window.PeakCount);
            Assert.AreEqual(window.Scans.Length, narrowWindow.Scans.Length);
            var wideMzs = window.Scans.SelectMany(s => s.Mz).ToHashSet();
            Assert.IsTrue(narrowWindow.Scans.SelectMany(s => s.Mz).All(mz => wideMzs.Contains(mz)));

            // peaks with no isotopic envelopes have no window, and neither do missing arguments or a charge state
            // the peak was never traced at
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var emptyPeak = new ChromatographicPeak(new Identification(mzml, "PEPTIDE", "PEPTIDE", 799.36, 94.1, 2,
                new List<ProteinGroup> { pg }), mzml);
            Assert.IsNull(PeakWindow.Create(emptyPeak, dataFile, map, chargeState));
            Assert.IsNull(PeakWindow.Create(null, dataFile, map, chargeState));
            Assert.IsNull(PeakWindow.Create(peak, null, map, chargeState));
            Assert.IsNull(PeakWindow.Create(peak, dataFile, null, chargeState));
            Assert.IsNull(PeakWindow.Create(peak, dataFile, map, chargeState: 97));

            // and the same holds for the per-charge-state enumeration
            Assert.IsFalse(PeakWindow.CreateForEachChargeState(emptyPeak, dataFile, map).Any());
            Assert.IsFalse(PeakWindow.CreateForEachChargeState(null, dataFile, map).Any());
            Assert.IsFalse(PeakWindow.CreateForEachChargeState(peak, null, map).Any());
            Assert.IsFalse(PeakWindow.CreateForEachChargeState(peak, dataFile, null).Any());
        }

        /// <summary>
        /// A ChromatographicPeak holds one isotopic envelope per scan per charge state, and a peptide's m/z at one
        /// charge is nowhere near its m/z at another. Bounding every envelope at once would produce a single window
        /// spanning the gap between the charge states - hundreds of daltons of unrelated spectrum - so each charge
        /// state gets its own window instead.
        /// </summary>
        [Test]
        public static void TestPeakWindowIsBoundedPerChargeState()
        {
            // one peak every 0.5 m/z from 400 to 1400, so a window's width is directly readable from its peak count
            var mzs = Enumerable.Range(0, 2001).Select(i => 400.0 + i * 0.5).ToArray();
            SpectraFileInfo file = WriteSyntheticMzml("peakWindowMultiCharge.mzML", new[] { mzs, mzs, mzs });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = PeakWindow.BuildMs1ScanNumberMap(dataFile);

            // EGFQVADGPLYR, 1350.66 Da, traced at z=2 (676.3 m/z) and z=3 (451.2 m/z) - both entirely ordinary
            ChromatographicPeak peak = BuildMultiChargePeak(file, "EGFQVADGPLYR",
                (676.3, 2, 0, 1.0), (451.2, 3, 0, 1.0),
                (676.3, 2, 1, 1.1), (451.2, 3, 1, 1.1),
                (676.3, 2, 2, 1.2), (451.2, 3, 2, 1.2));
            Assert.AreEqual(2, peak.NumChargeStatesObserved);

            var windows = PeakWindow.CreateForEachChargeState(peak, dataFile, map).ToList();

            // one window per charge state, in ascending charge order
            Assert.AreEqual(2, windows.Count);
            CollectionAssert.AreEqual(new[] { 2, 3 }, windows.Select(w => w.ChargeState).ToArray());

            // each window is mzExpansion wide on either side of its own charge state's m/z, and nowhere near the other
            // indexed peaks narrow m/z to float, so the bounds carry float precision - hence the loose deltas
            const double floatDelta = 1e-3;
            foreach (PeakWindow chargeWindow in windows)
            {
                Assert.AreEqual(2 * PeakWindow.DefaultMzExpansion, chargeWindow.MaxMz - chargeWindow.MinMz, floatDelta);
                Assert.AreEqual(3, chargeWindow.Scans.Length);

                // 1 m/z of a spectrum with a peak every 0.5 m/z: 2 peaks per scan, over 3 scans
                Assert.AreEqual(6, chargeWindow.PeakCount);
                Assert.AreEqual(6, chargeWindow.ToTsvRows(1).Count());
            }

            Assert.AreEqual(675.8, windows[0].MinMz, floatDelta);
            Assert.AreEqual(676.8, windows[0].MaxMz, floatDelta);
            Assert.AreEqual(450.7, windows[1].MinMz, floatDelta);
            Assert.AreEqual(451.7, windows[1].MaxMz, floatDelta);

            // the regression this guards: a single window bounding both charge states would have run from 450.7 to
            // 676.8 m/z - 226 Da - and swept in the whole spectrum between them
            Assert.IsTrue(windows.Sum(w => w.PeakCount) < 50);
        }

        /// <summary>
        /// Each charge state's window is centred on that charge state's own apex scan, not on the peak's overall
        /// apex, which can belong to a different charge state and lie outside the window entirely.
        /// </summary>
        [Test]
        public static void TestPeakWindowApexIsPerChargeState()
        {
            var mzs = new[] { 499.9, 500.0, 500.1, 599.9, 600.0, 600.1 };
            SpectraFileInfo file = WriteSyntheticMzml("peakWindowPerChargeApex.mzML", new[] { mzs, mzs, mzs });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = PeakWindow.BuildMs1ScanNumberMap(dataFile);

            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var id = new Identification(file, "PEPTIDE", "PEPTIDE", 1350.65681, 1.0, 2, new List<ProteinGroup> { pg });
            var peak = new ChromatographicPeak(id, file);

            // z=2 is most intense in the first scan; z=3 is most intense in the last one
            peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(new IndexedMassSpectralPeak(500.0, 9e6, 0, 1.0), 2, 9e6, 1.0));
            peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(new IndexedMassSpectralPeak(500.0, 1e6, 2, 1.2), 2, 1e6, 1.0));
            peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(new IndexedMassSpectralPeak(600.0, 1e5, 0, 1.0), 3, 1e5, 1.0));
            peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(new IndexedMassSpectralPeak(600.0, 5e6, 2, 1.2), 3, 5e6, 1.0));
            peak.CalculateIntensityForThisFeature(false);

            var windows = PeakWindow.CreateForEachChargeState(peak, dataFile, map).ToList();
            Assert.AreEqual(2, windows.Count);

            // the peak's own apex is the z=2 envelope in scan 1, but the z=3 window apexes in scan 3
            Assert.AreEqual(2, peak.Apex.ChargeState);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[0]).RetentionTime, windows[0].ApexRetentionTime, 1e-6);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[2]).RetentionTime, windows[1].ApexRetentionTime, 1e-6);
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
            int chargeState = peak.Apex.ChargeState;
            var envelopes = peak.IsotopicEnvelopes.Where(e => e.ChargeState == chargeState).ToList();
            PeakWindow window = PeakWindow.Create(peak, dataFile, map, chargeState, mzExpansion: 0.5);

            PeakIndexingEngine engine = PeakIndexingEngine.InitializeIndexingEngine(mzml);
            var indexedPeaks = engine.GetPeaksInRange(window.MinMz, window.MaxMz,
                envelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex),
                envelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex));

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

        /// <summary>
        /// A scan can hold nothing inside the window even though the peak elutes across it. Both ways that happens
        /// have to survive the binary search: the scan's peaks can all sit below the window, in which case the
        /// search clamps to the last peak and the start index has to be walked past the end of the array, or they
        /// can all sit above it, in which case it clamps to the first.
        /// </summary>
        [Test]
        public static void TestPeakWindowHandlesScansWithNothingInRange()
        {
            SpectraFileInfo file = WriteSyntheticMzml("peakWindowEmptyScans.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },  // scan 1: inside the window
                new[] { 100.0, 101.0 },         // scan 2: entirely below it
                new[] { 900.0, 901.0 },         // scan 3: entirely above it
                new[] { 499.9, 500.0, 500.1 },  // scan 4: inside the window
            });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = PeakWindow.BuildMs1ScanNumberMap(dataFile);

            // the peak is traced through the first and last scan, so the window spans all four
            ChromatographicPeak peak = BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0), (500.0, 3, 1.3));
            PeakWindow window = PeakWindow.CreateForEachChargeState(peak, dataFile, map).Single();

            Assert.AreEqual(499.5, window.MinMz, 1e-6);
            Assert.AreEqual(500.5, window.MaxMz, 1e-6);
            Assert.AreEqual(4, window.Scans.Length);

            // the scans with nothing in range are kept, holding no peaks
            Assert.AreEqual(3, window.Scans[0].Mz.Length);
            Assert.AreEqual(0, window.Scans[1].Mz.Length);
            Assert.AreEqual(0, window.Scans[2].Mz.Length);
            Assert.AreEqual(3, window.Scans[3].Mz.Length);
            Assert.AreEqual(6, window.PeakCount);

            foreach (PeakWindowScan scan in window.Scans)
            {
                Assert.AreEqual(scan.Mz.Length, scan.Intensity.Length);
                Assert.AreEqual(scan.Mz.Length, scan.IsPeakfindingPeak.Length);
            }

            // an empty scan contributes no rows, and the peakfinding peaks are still flagged
            Assert.AreEqual(6, window.ToTsvRows(1).Count());
            Assert.AreEqual(2, window.Scans.Sum(s => s.IsPeakfindingPeak.Count(f => f)));
        }

        [Test]
        public static void TestPeakWindowReturnsNullWhenScanIndexIsNotInTheMap()
        {
            SpectraFileInfo file = WriteSyntheticMzml("peakWindowShortFile.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },
                new[] { 499.9, 500.0, 500.1 },
            });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = PeakWindow.BuildMs1ScanNumberMap(dataFile);
            Assert.AreEqual(2, map.Length);

            // an envelope in a scan the file does not have cannot be translated to a scan number
            ChromatographicPeak beyondEnd = BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0), (500.0, 99, 1.3));
            Assert.IsNull(PeakWindow.Create(beyondEnd, dataFile, map, chargeState: 2));
            Assert.IsFalse(PeakWindow.CreateForEachChargeState(beyondEnd, dataFile, map).Any());
        }

        [Test]
        public static void TestPeakWindowToString()
        {
            SpectraFileInfo file = WriteSyntheticMzml("peakWindowToString.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },
                new[] { 499.9, 500.0, 500.1 },
            });

            MsDataFile dataFile = LoadDataFile(file);
            PeakWindow window = PeakWindow.CreateForEachChargeState(
                BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0), (500.0, 1, 1.1)),
                dataFile, PeakWindow.BuildMs1ScanNumberMap(dataFile)).Single();

            string description = window.ToString();
            StringAssert.Contains("PEPTIDE", description);
            StringAssert.Contains("z=2", description);
            StringAssert.Contains("6 peaks in 2 scans", description);
            StringAssert.Contains("m/z", description);
            StringAssert.Contains("min", description);
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
            int rtStartColumn = Array.IndexOf(headerColumns, "Window RT Start");
            int rtEndColumn = Array.IndexOf(headerColumns, "Window RT End");
            int peakfindingColumn = Array.IndexOf(headerColumns, "Is Peakfinding Peak");
            int fullSequenceColumn = Array.IndexOf(headerColumns, "Full Sequence");
            int chargeColumn = Array.IndexOf(headerColumns, "Peak Charge");

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

                // every row names the charge state of the window it came from, not the peak's apex charge
                int rowCharge = int.Parse(cells[chargeColumn], CultureInfo.InvariantCulture);
                Assert.IsTrue(peak.IsotopicEnvelopes.Any(e => e.ChargeState == rowCharge));

                if (bool.Parse(cells[peakfindingColumn]))
                {
                    peakfindingRows++;
                }
            }

            // the rows FlashLFQ actually quantified are flagged, and they are a strict subset of what was written
            Assert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.IndexedPeak).Distinct().Count(), peakfindingRows);
            Assert.IsTrue(lines.Length - 1 > peakfindingRows);

            // every charge state the peak was traced at contributed rows
            CollectionAssert.AreEquivalent(
                peak.IsotopicEnvelopes.Select(e => e.ChargeState).Distinct().OrderBy(z => z).ToArray(),
                lines.Skip(1).Select(l => int.Parse(l.Split('\t')[chargeColumn], CultureInfo.InvariantCulture))
                    .Distinct().OrderBy(z => z).ToArray());

            File.Delete(outputPath);
        }

        [Test]
        public static void TestWritePeakWindowsNullPathIsANoOp()
        {
            var results = new FlashLfqResults(new List<SpectraFileInfo>(), new List<Identification>());
            Assert.DoesNotThrow(() => results.WritePeakWindows(null, silent: true));
        }

        /// <summary>
        /// Files with no quantified peaks and peaks with no isotopic envelopes are both skipped. Runs with console
        /// output on, since that is how a caller would normally invoke it.
        /// </summary>
        [Test]
        public static void TestWritePeakWindowsSkipsWhatItCannotWrite()
        {
            SpectraFileInfo fileWithPeaks = WriteSyntheticMzml("peakWindowSkipping.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },
                new[] { 499.9, 500.0, 500.1 },
            });

            // this file is never read, because it has no quantified peaks - so it need not exist at all
            var fileWithoutPeaks = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory,
                "peakWindowNoSuchFile.mzML"), "a", 1, 0, 0);
            Assert.IsFalse(File.Exists(fileWithoutPeaks.FullFilePathWithExtension));

            var results = new FlashLfqResults(new List<SpectraFileInfo> { fileWithPeaks, fileWithoutPeaks },
                new List<Identification>());
            results.Peaks[fileWithPeaks].Add(BuildPeak(fileWithPeaks, "PEPTIDE", (500.0, 0, 1.0), (500.0, 1, 1.1)));
            results.Peaks[fileWithPeaks].Add(BuildPeak(fileWithPeaks, "NOENVELOPES"));

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindowsSkipping.tsv");
            Assert.DoesNotThrow(() => results.WritePeakWindows(outputPath, silent: false));

            string[] lines = File.ReadAllLines(outputPath);
            Assert.AreEqual(PeakWindow.TabSeparatedHeader, lines[0]);

            // only the peak that had envelopes produced rows
            Assert.AreEqual(6, lines.Length - 1);
            Assert.IsTrue(lines.Skip(1).All(l => l.Split('\t')[2] == "PEPTIDE"));

            File.Delete(outputPath);
        }

        /// <summary>
        /// A file whose peaks cannot be placed in any MS1 scan is reported and skipped rather than throwing.
        /// </summary>
        [Test]
        public static void TestWritePeakWindowsWithNoMs1Scans()
        {
            var ms2Scans = new MsDataScan[2];
            for (int s = 0; s < ms2Scans.Length; s++)
            {
                double[] mz = { 200.0, 300.0 };
                double[] intensity = { 1e5, 2e5 };
                ms2Scans[s] = new MsDataScan(new MzSpectrum(mz, intensity, false), s + 1, 2, true, Polarity.Positive,
                    1.0 + s * 0.1, new MzRange(50, 1600), "f", MZAnalyzerType.Orbitrap, intensity.Sum(), 1.0, null,
                    "scan=" + (s + 1), selectedIonMz: 500.0, dissociationType: DissociationType.HCD);
            }

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindowMs2Only.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(ms2Scans), path, false);
            var file = new SpectraFileInfo(path, "a", 0, 0, 0);

            var results = new FlashLfqResults(new List<SpectraFileInfo> { file }, new List<Identification>());
            results.Peaks[file].Add(BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0)));

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindowsMs2Only.tsv");
            Assert.DoesNotThrow(() => results.WritePeakWindows(outputPath, silent: false));

            // the header is written, but the file contributes no rows
            string[] lines = File.ReadAllLines(outputPath);
            Assert.AreEqual(1, lines.Length);
            Assert.AreEqual(PeakWindow.TabSeparatedHeader, lines[0]);

            File.Delete(outputPath);
        }

        #endregion
    }
}
