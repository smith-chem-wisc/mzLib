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
using System.Text.Json;
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
            return BuildPeak(file, sequence, true, envelopePeaks);
        }

        /// <param name="resolveApex"> when false the peak keeps its envelopes but never gets an apex, as happens
        /// before CalculateIntensityForThisFeature has run </param>
        private static ChromatographicPeak BuildPeak(SpectraFileInfo file, string sequence, bool resolveApex,
            params (double Mz, int ZeroBasedScanIndex, double RetentionTime)[] envelopePeaks)
        {
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var id = new Identification(file, sequence, sequence, 1350.65681, 1.0, 2, new List<ProteinGroup> { pg });
            var peak = new ChromatographicPeak(id, file);

            foreach (var envelopePeak in envelopePeaks)
            {
                peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(
                    new IndexedMassSpectralPeak(envelopePeak.Mz, 1e6, envelopePeak.ZeroBasedScanIndex,
                        envelopePeak.RetentionTime), 2, 1e6, 1.0));
            }
            if (resolveApex)
            {
                peak.CalculateIntensityForThisFeature(false);
            }
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

            PeakWindow window = PeakWindow.Create(peak, dataFile, map, mzExpansion: 2.0);

            // the window is bounded by the peak's observed m/z values, expanded by mzExpansion on either side
            Assert.AreEqual(peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.M) - 2.0, window.MinMz, 1e-6);
            Assert.AreEqual(peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.M) + 2.0, window.MaxMz, 1e-6);

            // and it covers exactly the scans the peak was observed in
            int firstMs1Index = peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex);
            int lastMs1Index = peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex);
            Assert.AreEqual(lastMs1Index - firstMs1Index + 1, window.ScanNumbers.Length);
            Assert.AreEqual(map[firstMs1Index], window.ScanNumbers.First());
            Assert.AreEqual(map[lastMs1Index], window.ScanNumbers.Last());

            // retention times are the scans' own, taken from the file
            Assert.AreEqual(dataFile.GetOneBasedScan(map[firstMs1Index]).RetentionTime, window.MinRetentionTime);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[lastMs1Index]).RetentionTime, window.MaxRetentionTime);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[peak.Apex.IndexedPeak.ZeroBasedScanIndex]).RetentionTime,
                window.ApexRetentionTime);

            // each scan's slice is the m/z window of that scan's spectrum, held as parallel arrays
            Assert.AreEqual(window.Mz.Sum(a => a.Length), window.PeakCount);
            Assert.AreEqual(window.ScanNumbers.Length, window.RetentionTimes.Length);
            Assert.AreEqual(window.ScanNumbers.Length, window.Mz.Length);
            Assert.AreEqual(window.ScanNumbers.Length, window.Intensity.Length);
            Assert.AreEqual(window.ScanNumbers.Length, window.PeakfindingIndices.Length);

            for (int i = 0; i < window.ScanNumbers.Length; i++)
            {
                MsDataScan sourceScan = dataFile.GetOneBasedScan(window.ScanNumbers[i]);
                Assert.AreEqual(sourceScan.RetentionTime, window.RetentionTimes[i]);
                Assert.AreEqual(window.Mz[i].Length, window.Intensity[i].Length);
                Assert.IsTrue(window.Mz[i].All(mz => mz >= window.MinMz && mz <= window.MaxMz));

                // the slice matches a brute force filter of the whole spectrum, so the binary search found it all
                var expected = sourceScan.MassSpectrum.XArray
                    .Where(mz => mz >= window.MinMz && mz <= window.MaxMz).Select(mz => (float)mz).ToArray();
                CollectionAssert.AreEqual(expected, window.Mz[i]);

                // peakfinding entries are indices into this scan's arrays, ascending and in range
                Assert.IsTrue(window.PeakfindingIndices[i].All(index => index >= 0 && index < window.Mz[i].Length));
                CollectionAssert.AreEqual(window.PeakfindingIndices[i].OrderBy(index => index).ToArray(),
                    window.PeakfindingIndices[i]);
            }

            // every peak FlashLFQ used to trace the feature is marked, and nothing else is
            int markedPeaks = window.PeakfindingIndices.Sum(a => a.Length);
            Assert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.IndexedPeak).Distinct().Count(), markedPeaks);

            // and the window holds far more than those, which is the whole point of the output
            Assert.IsTrue(window.PeakCount > markedPeaks);

            // a wider window is a superset of a narrower one
            PeakWindow narrowWindow = PeakWindow.Create(peak, dataFile, map, mzExpansion: 0.5);
            Assert.IsTrue(narrowWindow.PeakCount <= window.PeakCount);
            Assert.AreEqual(window.ScanNumbers.Length, narrowWindow.ScanNumbers.Length);
            var wideMzs = window.Mz.SelectMany(a => a).ToHashSet();
            Assert.IsTrue(narrowWindow.Mz.SelectMany(a => a).All(mz => wideMzs.Contains(mz)));

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
            var fromFile = window.ScanNumbers
                .SelectMany((scanNumber, i) => window.Mz[i].Select(mz => (OneBasedScanNumber: scanNumber, Mz: mz)))
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
            PeakWindow window = PeakWindow.Create(peak, dataFile, map);

            Assert.AreEqual(499.5, window.MinMz, 1e-6);
            Assert.AreEqual(500.5, window.MaxMz, 1e-6);
            Assert.AreEqual(4, window.ScanNumbers.Length);

            // the scans with nothing in range are kept, holding empty arrays
            Assert.AreEqual(3, window.Mz[0].Length);
            Assert.AreEqual(0, window.Mz[1].Length);
            Assert.AreEqual(0, window.Mz[2].Length);
            Assert.AreEqual(3, window.Mz[3].Length);
            Assert.AreEqual(6, window.PeakCount);

            for (int i = 0; i < window.ScanNumbers.Length; i++)
            {
                Assert.AreEqual(window.Mz[i].Length, window.Intensity[i].Length);
            }

            // the peakfinding peaks are still marked, in the scans that hold them
            Assert.AreEqual(2, window.PeakfindingIndices.Sum(a => a.Length));
            Assert.AreEqual(1, window.PeakfindingIndices[0].Single());
            Assert.AreEqual(0, window.PeakfindingIndices[1].Length);
            Assert.AreEqual(0, window.PeakfindingIndices[2].Length);
            Assert.AreEqual(1, window.PeakfindingIndices[3].Single());
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
            Assert.IsNull(PeakWindow.Create(beyondEnd, dataFile, map));
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
            PeakWindow window = PeakWindow.Create(BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0), (500.0, 1, 1.1)),
                dataFile, PeakWindow.BuildMs1ScanNumberMap(dataFile));

            string description = window.ToString();
            StringAssert.Contains("PEPTIDE", description);
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

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindows.jsonl");
            results.WritePeakWindows(outputPath, silent: true);

            Assert.IsTrue(File.Exists(outputPath));
            string[] lines = File.ReadAllLines(outputPath);

            // one line per quantified peak, each a complete JSON object
            Assert.AreEqual(1, lines.Length);

            MsDataFile dataFile = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            dataFile.LoadAllStaticData();

            using JsonDocument document = JsonDocument.Parse(lines[0]);
            JsonElement root = document.RootElement;

            // the metadata appears once for the whole peak, not once per MS1 peak
            Assert.AreEqual("sliced-mzml", root.GetProperty("fileName").GetString());
            Assert.AreEqual(1, root.GetProperty("peakId").GetInt32());
            Assert.AreEqual("EGFQVADGPLYR", root.GetProperty("baseSequence").GetString());
            Assert.AreEqual("EGFQVADGPLYR", root.GetProperty("fullSequence").GetString());
            Assert.AreEqual("MSMS", root.GetProperty("detectionType").GetString());
            Assert.AreEqual(peak.Apex.ChargeState, root.GetProperty("charge").GetInt32());

            double rtStart = root.GetProperty("rtStart").GetDouble();
            double rtEnd = root.GetProperty("rtEnd").GetDouble();
            double rtApex = root.GetProperty("rtApex").GetDouble();
            Assert.IsTrue(rtStart <= rtApex && rtApex <= rtEnd);

            double minMz = root.GetProperty("minMz").GetDouble();
            double maxMz = root.GetProperty("maxMz").GetDouble();

            int[] scanNumbers = root.GetProperty("scanNumbers").EnumerateArray().Select(e => e.GetInt32()).ToArray();
            double[] retentionTimes = root.GetProperty("retentionTimes").EnumerateArray().Select(e => e.GetDouble()).ToArray();
            var mzArrays = root.GetProperty("mz").EnumerateArray()
                .Select(a => a.EnumerateArray().Select(e => e.GetSingle()).ToArray()).ToArray();
            var intensityArrays = root.GetProperty("intensity").EnumerateArray()
                .Select(a => a.EnumerateArray().Select(e => e.GetSingle()).ToArray()).ToArray();
            var peakfindingIndices = root.GetProperty("peakfindingIndices").EnumerateArray()
                .Select(a => a.EnumerateArray().Select(e => e.GetInt32()).ToArray()).ToArray();

            // the arrays are parallel, one entry per scan
            Assert.AreEqual(scanNumbers.Length, retentionTimes.Length);
            Assert.AreEqual(scanNumbers.Length, mzArrays.Length);
            Assert.AreEqual(scanNumbers.Length, intensityArrays.Length);
            Assert.AreEqual(scanNumbers.Length, peakfindingIndices.Length);
            Assert.AreEqual(mzArrays.Sum(a => a.Length), root.GetProperty("peakCount").GetInt32());

            for (int i = 0; i < scanNumbers.Length; i++)
            {
                // the scan number addresses the raw file directly, and names an MS1 scan at the stated time
                MsDataScan sourceScan = dataFile.GetOneBasedScan(scanNumbers[i]);
                Assert.AreEqual(1, sourceScan.MsnOrder);
                Assert.AreEqual(sourceScan.RetentionTime, retentionTimes[i]);
                Assert.IsTrue(retentionTimes[i] >= rtStart && retentionTimes[i] <= rtEnd);

                Assert.AreEqual(mzArrays[i].Length, intensityArrays[i].Length);
                Assert.IsTrue(mzArrays[i].All(mz => mz >= minMz && mz <= maxMz));
                Assert.IsTrue(intensityArrays[i].All(intensity => intensity > 0));
                Assert.IsTrue(peakfindingIndices[i].All(index => index >= 0 && index < mzArrays[i].Length));
            }

            // the peaks FlashLFQ actually quantified are marked, and are a strict subset of what was written
            int markedPeaks = peakfindingIndices.Sum(a => a.Length);
            Assert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.IndexedPeak).Distinct().Count(), markedPeaks);
            Assert.IsTrue(mzArrays.Sum(a => a.Length) > markedPeaks);

            File.Delete(outputPath);
        }

        /// <summary>
        /// The reason for the format: the metadata is written once per peak rather than once per MS1 peak, so it
        /// stays a rounding error in the output no matter how much raw data the window covers.
        /// </summary>
        [Test]
        public static void TestPeakWindowJsonDoesNotRepeatMetadata()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            QuantifySlicedMzml(mzml, out FlashLfqResults results);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindowsMetadata.jsonl");
            results.WritePeakWindows(outputPath, mzExpansion: 2.0, silent: true);

            string json = File.ReadAllText(outputPath);
            using JsonDocument document = JsonDocument.Parse(json);
            int peakCount = document.RootElement.GetProperty("peakCount").GetInt32();

            // the sequence names each appear once for the whole peak, however many MS1 peaks it covers
            Assert.IsTrue(peakCount > 100);
            Assert.AreEqual(1, CountOccurrences(json, "\"fullSequence\""));
            Assert.AreEqual(1, CountOccurrences(json, "EGFQVADGPLYR\",\"detectionType\""));
            Assert.AreEqual(1, CountOccurrences(json, "\"peakId\""));

            File.Delete(outputPath);
        }

        /// <summary>
        /// A peak with no apex has no apex retention time, and NaN is not a JSON number - Utf8JsonWriter throws on
        /// it - so it has to be written as null or the whole file is lost to one unresolved peak.
        /// </summary>
        [Test]
        public static void TestPeakWithNoApexWritesNullRetentionTime()
        {
            SpectraFileInfo file = WriteSyntheticMzml("peakWindowNoApex.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },
                new[] { 499.9, 500.0, 500.1 },
            });

            var results = new FlashLfqResults(new List<SpectraFileInfo> { file }, new List<Identification>());
            results.Peaks[file].Add(BuildPeak(file, "PEPTIDE", false, (500.0, 0, 1.0), (500.0, 1, 1.1)));

            MsDataFile dataFile = LoadDataFile(file);
            PeakWindow window = PeakWindow.Create(results.Peaks[file].Single(), dataFile,
                PeakWindow.BuildMs1ScanNumberMap(dataFile));
            Assert.IsNaN(window.ApexRetentionTime);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindowsNoApex.jsonl");
            Assert.DoesNotThrow(() => results.WritePeakWindows(outputPath, silent: true));

            // the object still parses, with a null apex and a zero charge rather than a broken number
            string line = File.ReadAllLines(outputPath).Single();
            using JsonDocument document = JsonDocument.Parse(line);
            Assert.AreEqual(JsonValueKind.Null, document.RootElement.GetProperty("rtApex").ValueKind);
            Assert.AreEqual(0, document.RootElement.GetProperty("charge").GetInt32());
            Assert.AreEqual(6, document.RootElement.GetProperty("peakCount").GetInt32());

            File.Delete(outputPath);
        }

        private static int CountOccurrences(string haystack, string needle)
        {
            int count = 0;
            for (int i = haystack.IndexOf(needle, StringComparison.Ordinal); i >= 0;
                 i = haystack.IndexOf(needle, i + needle.Length, StringComparison.Ordinal))
            {
                count++;
            }
            return count;
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

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindowsSkipping.jsonl");
            Assert.DoesNotThrow(() => results.WritePeakWindows(outputPath, silent: false));

            // only the peak that had envelopes produced an object
            string[] lines = File.ReadAllLines(outputPath);
            Assert.AreEqual(1, lines.Length);

            using JsonDocument document = JsonDocument.Parse(lines[0]);
            Assert.AreEqual("PEPTIDE", document.RootElement.GetProperty("fullSequence").GetString());
            Assert.AreEqual(6, document.RootElement.GetProperty("peakCount").GetInt32());

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

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "peakWindowsMs2Only.jsonl");
            Assert.DoesNotThrow(() => results.WritePeakWindows(outputPath, silent: false));

            // the file contributes nothing, leaving an empty output rather than a malformed one
            Assert.AreEqual(0, new FileInfo(outputPath).Length);

            File.Delete(outputPath);
        }

        #endregion
    }
}
