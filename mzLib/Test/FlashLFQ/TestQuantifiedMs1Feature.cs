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
    internal class TestQuantifiedMs1Feature
    {
        private static string SlicedMzmlPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "FlashLFQ", "TestData", @"sliced-mzml.mzml");

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
        /// <param name="subdirectory"> written under the test directory when given, so that two files can share a
        /// base name </param>
        private static SpectraFileInfo WriteSyntheticMzml(string fileName, IReadOnlyList<double[]> mzPerScan,
            string subdirectory = null)
        {
            var scans = new MsDataScan[mzPerScan.Count];
            for (int s = 0; s < mzPerScan.Count; s++)
            {
                double[] intensities = mzPerScan[s].Select((_, i) => 1e5 * (i + 1)).ToArray();
                scans[s] = new MsDataScan(new MzSpectrum(mzPerScan[s].ToArray(), intensities, false), s + 1, 1, true,
                    Polarity.Positive, 1.0 + s * 0.1, new MzRange(50, 1600), "f", MZAnalyzerType.Orbitrap,
                    intensities.Sum(), 1.0, null, "scan=" + (s + 1));
            }

            string directory = subdirectory == null
                ? TestContext.CurrentContext.TestDirectory
                : Path.Combine(TestContext.CurrentContext.TestDirectory, subdirectory);
            Directory.CreateDirectory(directory);

            string path = Path.Combine(directory, fileName);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans), path, false);
            return new SpectraFileInfo(path, "a", 0, 0, 0);
        }

        /// <summary>
        /// Builds a ChromatographicPeak directly, without running the engine, so that a test can put its isotopic
        /// envelopes in exactly the scans it wants - including scans that do not exist. Every envelope is traced at
        /// charge 2.
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
            return BuildPeak(file, sequence, resolveApex,
                envelopePeaks.Select(e => (e.Mz, 2, e.ZeroBasedScanIndex, e.RetentionTime)).ToArray());
        }

        /// <summary>
        /// Builds a ChromatographicPeak whose envelopes can be spread over several charge states, so that a test
        /// can control how many windows the feature ends up with.
        /// </summary>
        private static ChromatographicPeak BuildPeak(SpectraFileInfo file, string sequence, bool resolveApex,
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

            int[] map = QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile);
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

        /// <summary>
        /// The map and the indexing engine's own scan metadata have to agree, since the zero based indices being
        /// translated were handed out by the engine. Both come from PeakIndexingEngine.GetScansToIndex so that they
        /// cannot drift: were they two copies of the same selection, every scan number in this output would go
        /// silently wrong the first time one of them changed.
        /// </summary>
        [Test]
        public static void TestMs1ScanNumberMapAgreesWithTheIndexingEngine()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            QuantifySlicedMzml(mzml, out FlashLfqResults results);

            MsDataFile dataFile = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            dataFile.LoadAllStaticData();

            // the engine leaves its scan metadata on the results, and it survives the index being cleared
            Assert.IsNotNull(results.Ms1ScanInfo);
            ScanInfo[] scanInfo = results.Ms1ScanInfo[mzml];

            CollectionAssert.AreEqual(QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile),
                scanInfo.Select(s => s.OneBasedScanNumber).ToArray());

            // the zero based indices are the positions in that same array, which is what makes it a map
            for (int i = 0; i < scanInfo.Length; i++)
            {
                Assert.AreEqual(i, scanInfo[i].ZeroBasedScanIndex);
                Assert.AreEqual(1, scanInfo[i].MsnOrder);
            }
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
            int[] map = QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile);

            QuantifiedMs1Feature feature = QuantifiedMs1Feature.Create(peak, dataFile, map, mzExpansion: 2.0);

            // one window per charge state the peak was observed at, ascending
            CollectionAssert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.ChargeState).Distinct().OrderBy(z => z).ToArray(),
                feature.ChargeStates.ToArray());
            Assert.AreEqual(peak.Apex.ChargeState, feature.ApexChargeState);

            foreach (PeakWindow window in feature.Windows)
            {
                var envelopes = peak.IsotopicEnvelopes.Where(e => e.ChargeState == window.ChargeState).ToList();
                Assert.AreSame(window, feature.GetWindow(window.ChargeState));

                // the window is bounded by the m/z values observed at its own charge state, expanded by
                // mzExpansion on either side - not by the m/z values of the peak's other charge states
                Assert.AreEqual(envelopes.Min(e => e.IndexedPeak.M) - 2.0, window.MinMz, 1e-6);
                Assert.AreEqual(envelopes.Max(e => e.IndexedPeak.M) + 2.0, window.MaxMz, 1e-6);

                // and it covers exactly the scans that charge state was observed in
                int firstMs1Index = envelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex);
                int lastMs1Index = envelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex);
                Assert.AreEqual(lastMs1Index - firstMs1Index + 1, window.ScanNumbers.Length);
                Assert.AreEqual(map[firstMs1Index], window.ScanNumbers.First());
                Assert.AreEqual(map[lastMs1Index], window.ScanNumbers.Last());

                // retention times are the scans' own, taken from the file
                Assert.AreEqual(dataFile.GetOneBasedScan(map[firstMs1Index]).RetentionTime, window.MinRetentionTime);
                Assert.AreEqual(dataFile.GetOneBasedScan(map[lastMs1Index]).RetentionTime, window.MaxRetentionTime);

                // the window's apex is the apex of its own charge state
                IsotopicEnvelope chargeApex = envelopes.MaxBy(e => e.Intensity);
                Assert.AreEqual(dataFile.GetOneBasedScan(map[chargeApex.IndexedPeak.ZeroBasedScanIndex]).RetentionTime,
                    window.ApexRetentionTime);
                Assert.AreEqual(chargeApex.Intensity, window.ApexIntensity, 1e-6);

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

                // every peak FlashLFQ used to trace this charge state is marked, and nothing else is
                Assert.AreEqual(envelopes.Select(e => e.IndexedPeak).Distinct().Count(),
                    window.PeakfindingIndices.Sum(a => a.Length));
            }

            // the feature's bounds are those of its windows taken together
            Assert.AreEqual(feature.Windows.Min(w => w.MinMz), feature.MinMz);
            Assert.AreEqual(feature.Windows.Max(w => w.MaxMz), feature.MaxMz);
            Assert.AreEqual(feature.Windows.Min(w => w.MinRetentionTime), feature.MinRetentionTime);
            Assert.AreEqual(feature.Windows.Max(w => w.MaxRetentionTime), feature.MaxRetentionTime);
            Assert.AreEqual(feature.Windows.Sum(w => w.PeakCount), feature.PeakCount);
            Assert.AreEqual(dataFile.GetOneBasedScan(map[peak.Apex.IndexedPeak.ZeroBasedScanIndex]).RetentionTime,
                feature.ApexRetentionTime);

            // and the feature holds far more peaks than the ones FlashLFQ traced, which is the whole point
            int markedPeaks = feature.Windows.Sum(w => w.PeakfindingIndices.Sum(a => a.Length));
            Assert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.IndexedPeak).Distinct().Count(), markedPeaks);
            Assert.IsTrue(feature.PeakCount > markedPeaks);

            // a wider window is a superset of a narrower one
            QuantifiedMs1Feature narrow = QuantifiedMs1Feature.Create(peak, dataFile, map, mzExpansion: 0.5);
            Assert.AreEqual(feature.Windows.Count, narrow.Windows.Count);
            Assert.IsTrue(narrow.PeakCount <= feature.PeakCount);
            for (int i = 0; i < narrow.Windows.Count; i++)
            {
                Assert.AreEqual(feature.Windows[i].ScanNumbers.Length, narrow.Windows[i].ScanNumbers.Length);
                var wideMzs = feature.Windows[i].Mz.SelectMany(a => a).ToHashSet();
                Assert.IsTrue(narrow.Windows[i].Mz.SelectMany(a => a).All(mz => wideMzs.Contains(mz)));
            }

            // peaks with no isotopic envelopes have no feature, and neither do missing arguments
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var emptyPeak = new ChromatographicPeak(new Identification(mzml, "PEPTIDE", "PEPTIDE", 799.36, 94.1, 2,
                new List<ProteinGroup> { pg }), mzml);
            Assert.IsNull(QuantifiedMs1Feature.Create(emptyPeak, dataFile, map));
            Assert.IsNull(QuantifiedMs1Feature.Create(null, dataFile, map));
            Assert.IsNull(QuantifiedMs1Feature.Create(peak, null, map));
            Assert.IsNull(QuantifiedMs1Feature.Create(peak, dataFile, null));
        }

        /// <summary>
        /// The reason a window covers one charge state: the charge states of a peptide sit at unrelated m/z values,
        /// and a single rectangle enclosing all of them would drag in every species in between.
        /// </summary>
        [Test]
        public static void TestChargeStatesGetSeparateWindows()
        {
            // the 3+ cluster near 334, the 2+ cluster near 500, and an unrelated species at 400 sitting between them
            SpectraFileInfo file = WriteSyntheticMzml("ms1FeatureMultiCharge.mzML", new[]
            {
                new[] { 333.9, 334.0, 334.1, 400.0, 499.9, 500.0, 500.1 },
                new[] { 333.9, 334.0, 334.1, 400.0, 499.9, 500.0, 500.1 },
            });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile);

            // 2+ is traced across both scans, 3+ only in the first
            ChromatographicPeak peak = BuildPeak(file, "PEPTIDE", true,
                (500.0, 2, 0, 1.0), (500.0, 2, 1, 1.1), (334.0, 3, 0, 1.0));
            QuantifiedMs1Feature feature = QuantifiedMs1Feature.Create(peak, dataFile, map);

            Assert.AreEqual(2, feature.Windows.Count);
            CollectionAssert.AreEqual(new[] { 2, 3 }, feature.ChargeStates.ToArray());

            PeakWindow doublyCharged = feature.GetWindow(2);
            Assert.AreEqual(499.5, doublyCharged.MinMz, 1e-6);
            Assert.AreEqual(500.5, doublyCharged.MaxMz, 1e-6);
            Assert.AreEqual(2, doublyCharged.ScanNumbers.Length);
            Assert.AreEqual(6, doublyCharged.PeakCount);

            PeakWindow triplyCharged = feature.GetWindow(3);
            Assert.AreEqual(333.5, triplyCharged.MinMz, 1e-6);
            Assert.AreEqual(334.5, triplyCharged.MaxMz, 1e-6);
            Assert.AreEqual(1, triplyCharged.ScanNumbers.Length);
            Assert.AreEqual(3, triplyCharged.PeakCount);

            // neither window contains the other's charge state, nor the species sitting between them
            Assert.IsFalse(doublyCharged.Mz.SelectMany(a => a).Any(mz => mz < 499.5f));
            Assert.IsFalse(triplyCharged.Mz.SelectMany(a => a).Any(mz => mz > 334.5f));
            Assert.IsFalse(feature.Windows.SelectMany(w => w.Mz).SelectMany(a => a).Any(mz => mz == 400.0f));

            // the feature spans both windows, and holds nine peaks rather than the fourteen a single rectangle
            // from 333.5 to 500.5 would have swept up
            Assert.AreEqual(333.5, feature.MinMz, 1e-6);
            Assert.AreEqual(500.5, feature.MaxMz, 1e-6);
            Assert.AreEqual(9, feature.PeakCount);

            // a charge state with no window is simply absent
            Assert.IsNull(feature.GetWindow(4));

            // the peak's most intense envelope is 2+, since IsotopicEnvelope divides intensity by charge
            Assert.AreEqual(2, feature.ApexChargeState);
        }

        /// <summary>
        /// A window is a single charge state by construction, so envelopes from more than one are a caller bug
        /// rather than something to silently split.
        /// </summary>
        [Test]
        public static void TestPeakWindowRejectsMixedChargeStates()
        {
            SpectraFileInfo file = WriteSyntheticMzml("ms1FeatureMixedCharge.mzML", new[]
            {
                new[] { 333.9, 334.0, 334.1, 499.9, 500.0, 500.1 },
            });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile);
            ChromatographicPeak peak = BuildPeak(file, "PEPTIDE", true, (500.0, 2, 0, 1.0), (334.0, 3, 0, 1.0));

            var exception = Assert.Throws<ArgumentException>(() =>
                PeakWindow.Create(peak.IsotopicEnvelopes, dataFile, map));
            StringAssert.Contains("single charge state", exception.Message);

            // and the degenerate inputs come back null rather than throwing
            Assert.IsNull(PeakWindow.Create(null, dataFile, map));
            Assert.IsNull(PeakWindow.Create(new List<IsotopicEnvelope>(), dataFile, map));
            Assert.IsNull(PeakWindow.Create(peak.IsotopicEnvelopes.Take(1).ToList(), null, map));
            Assert.IsNull(PeakWindow.Create(peak.IsotopicEnvelopes.Take(1).ToList(), dataFile, null));
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
            SpectraFileInfo file = WriteSyntheticMzml("ms1FeatureEmptyScans.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },  // scan 1: inside the window
                new[] { 100.0, 101.0 },         // scan 2: entirely below it
                new[] { 900.0, 901.0 },         // scan 3: entirely above it
                new[] { 499.9, 500.0, 500.1 },  // scan 4: inside the window
            });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile);

            // the peak is traced through the first and last scan, so the window spans all four
            ChromatographicPeak peak = BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0), (500.0, 3, 1.3));
            PeakWindow window = QuantifiedMs1Feature.Create(peak, dataFile, map).Windows.Single();

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
        public static void TestPeakWindowIsDroppedWhenScanIndexIsNotInTheMap()
        {
            SpectraFileInfo file = WriteSyntheticMzml("ms1FeatureShortFile.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },
                new[] { 499.9, 500.0, 500.1 },
            });

            MsDataFile dataFile = LoadDataFile(file);
            int[] map = QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile);
            Assert.AreEqual(2, map.Length);

            // an envelope in a scan the file does not have cannot be translated to a scan number, so its charge
            // state produces no window - and with only the one charge state, no feature either
            ChromatographicPeak beyondEnd = BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0), (500.0, 99, 1.3));
            Assert.IsNull(QuantifiedMs1Feature.Create(beyondEnd, dataFile, map));

            // a peak whose other charge states are fine keeps them, and loses only the unmappable one
            ChromatographicPeak partiallyMappable = BuildPeak(file, "PEPTIDE", true,
                (500.0, 2, 0, 1.0), (500.0, 2, 1, 1.1), (500.0, 3, 99, 1.3));
            QuantifiedMs1Feature feature = QuantifiedMs1Feature.Create(partiallyMappable, dataFile, map);
            CollectionAssert.AreEqual(new[] { 2 }, feature.ChargeStates.ToArray());
        }

        [Test]
        public static void TestToString()
        {
            SpectraFileInfo file = WriteSyntheticMzml("ms1FeatureToString.mzML", new[]
            {
                new[] { 333.9, 334.0, 334.1, 499.9, 500.0, 500.1 },
                new[] { 333.9, 334.0, 334.1, 499.9, 500.0, 500.1 },
            });

            MsDataFile dataFile = LoadDataFile(file);
            ChromatographicPeak peak = BuildPeak(file, "PEPTIDE", true,
                (500.0, 2, 0, 1.0), (500.0, 2, 1, 1.1), (334.0, 3, 0, 1.0));
            QuantifiedMs1Feature feature = QuantifiedMs1Feature.Create(peak, dataFile,
                QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile));

            string featureDescription = feature.ToString();
            StringAssert.Contains("PEPTIDE", featureDescription);
            StringAssert.Contains("2 charge states (+2, +3)", featureDescription);
            StringAssert.Contains("9 peaks", featureDescription);
            StringAssert.Contains("m/z", featureDescription);
            StringAssert.Contains("min", featureDescription);

            string windowDescription = feature.GetWindow(2).ToString();
            StringAssert.Contains("+2", windowDescription);
            StringAssert.Contains("6 peaks in 2 scans", windowDescription);
            StringAssert.Contains("m/z", windowDescription);
            StringAssert.Contains("min", windowDescription);
        }

        #endregion

        #region Output file

        [Test]
        public static void TestWriteMs1Features()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            ChromatographicPeak peak = QuantifySlicedMzml(mzml, out FlashLfqResults results);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1Features.jsonl");
            results.WriteMs1Features(outputPath, silent: true);

            Assert.IsTrue(File.Exists(outputPath));
            string[] lines = File.ReadAllLines(outputPath);

            // one line per quantified peak, each a complete JSON object
            Assert.AreEqual(1, lines.Length);

            MsDataFile dataFile = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            dataFile.LoadAllStaticData();

            using JsonDocument document = JsonDocument.Parse(lines[0]);
            JsonElement root = document.RootElement;

            // the metadata appears once for the whole peak, not once per charge state and not once per MS1 peak
            Assert.AreEqual("sliced-mzml", root.GetProperty("fileName").GetString());
            Assert.AreEqual(1, root.GetProperty("featureId").GetInt32());
            Assert.AreEqual("EGFQVADGPLYR", root.GetProperty("baseSequence").GetString());
            Assert.AreEqual("EGFQVADGPLYR", root.GetProperty("fullSequence").GetString());
            Assert.AreEqual("MSMS", root.GetProperty("detectionType").GetString());
            Assert.AreEqual(peak.Apex.ChargeState, root.GetProperty("apexCharge").GetInt32());

            double rtStart = root.GetProperty("rtStart").GetDouble();
            double rtEnd = root.GetProperty("rtEnd").GetDouble();
            double rtApex = root.GetProperty("rtApex").GetDouble();
            Assert.IsTrue(rtStart <= rtApex && rtApex <= rtEnd);

            double featureMinMz = root.GetProperty("minMz").GetDouble();
            double featureMaxMz = root.GetProperty("maxMz").GetDouble();

            int[] charges = root.GetProperty("charges").EnumerateArray().Select(e => e.GetInt32()).ToArray();
            var windows = root.GetProperty("windows").EnumerateArray().ToArray();
            Assert.AreEqual(charges.Length, windows.Length);
            CollectionAssert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.ChargeState).Distinct().OrderBy(z => z).ToArray(),
                charges);

            int totalPeaks = 0;
            int markedPeaks = 0;

            for (int w = 0; w < windows.Length; w++)
            {
                JsonElement window = windows[w];
                Assert.AreEqual(charges[w], window.GetProperty("charge").GetInt32());

                double minMz = window.GetProperty("minMz").GetDouble();
                double maxMz = window.GetProperty("maxMz").GetDouble();
                Assert.IsTrue(minMz >= featureMinMz && maxMz <= featureMaxMz);

                double windowRtStart = window.GetProperty("rtStart").GetDouble();
                double windowRtEnd = window.GetProperty("rtEnd").GetDouble();
                double windowRtApex = window.GetProperty("rtApex").GetDouble();
                Assert.IsTrue(windowRtStart <= windowRtApex && windowRtApex <= windowRtEnd);
                Assert.IsTrue(windowRtStart >= rtStart && windowRtEnd <= rtEnd);
                Assert.IsTrue(window.GetProperty("apexIntensity").GetDouble() > 0);

                int[] scanNumbers = window.GetProperty("scanNumbers").EnumerateArray().Select(e => e.GetInt32()).ToArray();
                double[] retentionTimes = window.GetProperty("retentionTimes").EnumerateArray().Select(e => e.GetDouble()).ToArray();
                var mzArrays = window.GetProperty("mz").EnumerateArray()
                    .Select(a => a.EnumerateArray().Select(e => e.GetSingle()).ToArray()).ToArray();
                var intensityArrays = window.GetProperty("intensity").EnumerateArray()
                    .Select(a => a.EnumerateArray().Select(e => e.GetSingle()).ToArray()).ToArray();
                var peakfindingIndices = window.GetProperty("peakfindingIndices").EnumerateArray()
                    .Select(a => a.EnumerateArray().Select(e => e.GetInt32()).ToArray()).ToArray();

                // the arrays are parallel, one entry per scan
                Assert.AreEqual(scanNumbers.Length, retentionTimes.Length);
                Assert.AreEqual(scanNumbers.Length, mzArrays.Length);
                Assert.AreEqual(scanNumbers.Length, intensityArrays.Length);
                Assert.AreEqual(scanNumbers.Length, peakfindingIndices.Length);
                Assert.AreEqual(mzArrays.Sum(a => a.Length), window.GetProperty("peakCount").GetInt32());

                for (int i = 0; i < scanNumbers.Length; i++)
                {
                    // the scan number addresses the raw file directly, and names an MS1 scan at the stated time
                    MsDataScan sourceScan = dataFile.GetOneBasedScan(scanNumbers[i]);
                    Assert.AreEqual(1, sourceScan.MsnOrder);
                    Assert.AreEqual(sourceScan.RetentionTime, retentionTimes[i]);
                    Assert.IsTrue(retentionTimes[i] >= windowRtStart && retentionTimes[i] <= windowRtEnd);

                    Assert.AreEqual(mzArrays[i].Length, intensityArrays[i].Length);
                    Assert.IsTrue(mzArrays[i].All(mz => mz >= minMz && mz <= maxMz));
                    Assert.IsTrue(intensityArrays[i].All(intensity => intensity > 0));
                    Assert.IsTrue(peakfindingIndices[i].All(index => index >= 0 && index < mzArrays[i].Length));
                }

                totalPeaks += mzArrays.Sum(a => a.Length);
                markedPeaks += peakfindingIndices.Sum(a => a.Length);
            }

            // the peaks FlashLFQ actually quantified are marked, and are a strict subset of what was written
            Assert.AreEqual(totalPeaks, root.GetProperty("peakCount").GetInt32());
            Assert.AreEqual(peak.IsotopicEnvelopes.Select(e => e.IndexedPeak).Distinct().Count(), markedPeaks);
            Assert.IsTrue(totalPeaks > markedPeaks);

            File.Delete(outputPath);
        }

        /// <summary>
        /// The reason for the format: the metadata is written once per peak rather than once per charge state or
        /// per MS1 peak, so it stays a rounding error in the output no matter how much raw data the windows cover.
        /// </summary>
        [Test]
        public static void TestMs1FeatureJsonDoesNotRepeatMetadata()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            QuantifySlicedMzml(mzml, out FlashLfqResults results);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeaturesMetadata.jsonl");
            results.WriteMs1Features(outputPath, mzExpansion: 2.0, silent: true);

            string json = File.ReadAllText(outputPath);
            using JsonDocument document = JsonDocument.Parse(json);
            int peakCount = document.RootElement.GetProperty("peakCount").GetInt32();

            // the sequence names each appear once for the whole peak, however many MS1 peaks it covers
            Assert.IsTrue(peakCount > 100);
            Assert.AreEqual(1, CountOccurrences(json, "\"fullSequence\""));
            Assert.AreEqual(1, CountOccurrences(json, "EGFQVADGPLYR\",\"detectionType\""));
            Assert.AreEqual(1, CountOccurrences(json, "\"featureId\""));

            File.Delete(outputPath);
        }

        /// <summary>
        /// A peak with no apex has no apex retention time, and NaN is not a JSON number - Utf8JsonWriter throws on
        /// it - so it has to be written as null or the whole file is lost to one unresolved peak.
        /// </summary>
        [Test]
        public static void TestPeakWithNoApexWritesNullRetentionTime()
        {
            SpectraFileInfo file = WriteSyntheticMzml("ms1FeatureNoApex.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },
                new[] { 499.9, 500.0, 500.1 },
            });

            var results = new FlashLfqResults(new List<SpectraFileInfo> { file }, new List<Identification>());
            results.Peaks[file].Add(BuildPeak(file, "PEPTIDE", false, (500.0, 0, 1.0), (500.0, 1, 1.1)));

            MsDataFile dataFile = LoadDataFile(file);
            QuantifiedMs1Feature feature = QuantifiedMs1Feature.Create(results.Peaks[file].Single(), dataFile,
                QuantifiedMs1Feature.BuildMs1ScanNumberMap(dataFile));
            Assert.IsNaN(feature.ApexRetentionTime);
            Assert.AreEqual(0, feature.ApexChargeState);

            // the windows still have an apex of their own: a charge state's apex is its most intense envelope,
            // which does not depend on the peak's intensity having been calculated
            Assert.IsFalse(double.IsNaN(feature.Windows.Single().ApexRetentionTime));

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeaturesNoApex.jsonl");
            Assert.DoesNotThrow(() => results.WriteMs1Features(outputPath, silent: true));

            // the object still parses, with a null apex and a zero charge rather than a broken number
            string line = File.ReadAllLines(outputPath).Single();
            using JsonDocument document = JsonDocument.Parse(line);
            Assert.AreEqual(JsonValueKind.Null, document.RootElement.GetProperty("rtApex").ValueKind);
            Assert.AreEqual(0, document.RootElement.GetProperty("apexCharge").GetInt32());
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

        /// <summary>
        /// featureId joins this output to QuantifiedPeaks.tsv, whose File Name column is the base file name. Two
        /// SpectraFileInfo can share one - the same run registered as two fractions, or two directories holding a
        /// file of the same name - and are then indistinguishable in both outputs, so the ids have to be unique
        /// across all of them rather than restarting per SpectraFileInfo.
        /// </summary>
        [Test]
        public static void TestFeatureIdsAreUniqueWithinAFileName()
        {
            var mzPerScan = new[] { new[] { 499.9, 500.0, 500.1 }, new[] { 499.9, 500.0, 500.1 } };
            SpectraFileInfo firstFile = WriteSyntheticMzml("sharedName.mzML", mzPerScan, "ms1FeatureIdsA");
            SpectraFileInfo secondFile = WriteSyntheticMzml("sharedName.mzML", mzPerScan, "ms1FeatureIdsB");

            Assert.AreEqual(firstFile.FilenameWithoutExtension, secondFile.FilenameWithoutExtension);
            Assert.AreNotEqual(firstFile.FullFilePathWithExtension, secondFile.FullFilePathWithExtension);

            var results = new FlashLfqResults(new List<SpectraFileInfo> { firstFile, secondFile },
                new List<Identification>());
            foreach (SpectraFileInfo file in new[] { firstFile, secondFile })
            {
                results.Peaks[file].Add(BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0), (500.0, 1, 1.1)));
                results.Peaks[file].Add(BuildPeak(file, "PEPTIDEK", (500.0, 0, 1.0), (500.0, 1, 1.1)));
            }

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeaturesSharedName.jsonl");
            results.WriteMs1Features(outputPath, silent: true);

            var written = File.ReadAllLines(outputPath)
                .Select(line => JsonDocument.Parse(line).RootElement)
                .Select(root => (FileName: root.GetProperty("fileName").GetString(),
                    FeatureId: root.GetProperty("featureId").GetInt32()))
                .ToList();

            // all four peaks are written under the one base name, and no two share an id
            Assert.AreEqual(4, written.Count);
            Assert.IsTrue(written.All(w => w.FileName == "sharedName"));
            CollectionAssert.AreEquivalent(new[] { 1, 2, 3, 4 }, written.Select(w => w.FeatureId).ToArray());

            File.Delete(outputPath);
        }

        /// <summary>
        /// A file contributing only a few peaks is read scan by scan through a dynamic connection rather than being
        /// materialized in full, which it can only do when a run left its scan metadata behind. The two paths have
        /// to produce the same bytes, or the output would depend on how the results were assembled.
        /// </summary>
        [Test]
        public static void TestDynamicAndStaticReadsProduceTheSameOutput()
        {
            SpectraFileInfo mzml = new SpectraFileInfo(SlicedMzmlPath, "a", 0, 0, 0);
            ChromatographicPeak peak = QuantifySlicedMzml(mzml, out FlashLfqResults engineResults);

            // the engine's results carry the scan metadata, and this peak's windows cover far fewer scans than the
            // file holds, so the writer seeks to the scans it needs instead of loading the file
            Assert.IsNotNull(engineResults.Ms1ScanInfo);
            Assert.IsTrue(QuantifiedMs1Feature.CountScansToRead(peak) < engineResults.Ms1ScanInfo[mzml].Length);

            // the writer falls back to a static load for files that cannot be read a scan at a time, so pin down
            // that this file is not one of them - otherwise both halves of this test would take the same path and
            // the comparison below would be vacuous
            MsDataFile dynamicallyRead = MsDataFileReader.GetDataFile(SlicedMzmlPath);
            Assert.DoesNotThrow(() => dynamicallyRead.InitiateDynamicConnection());
            MsDataScan firstMs1Scan = dynamicallyRead.GetOneBasedScanFromDynamicConnection(
                engineResults.Ms1ScanInfo[mzml].First().OneBasedScanNumber);
            Assert.IsNotNull(firstMs1Scan);
            Assert.AreEqual(1, firstMs1Scan.MsnOrder);
            dynamicallyRead.CloseDynamicConnection();

            string dynamicPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeaturesDynamic.jsonl");
            engineResults.WriteMs1Features(dynamicPath, silent: true);

            // the same peak, in results with no run behind them: nothing to take the map from, so it is rebuilt
            // from the file and the scans are read out of a static load
            var handBuiltResults = new FlashLfqResults(new List<SpectraFileInfo> { mzml }, new List<Identification>());
            handBuiltResults.Peaks[mzml].Add(peak);
            Assert.IsNull(handBuiltResults.Ms1ScanInfo);
            string staticPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeaturesStatic.jsonl");
            handBuiltResults.WriteMs1Features(staticPath, silent: true);

            Assert.IsTrue(new FileInfo(dynamicPath).Length > 0);
            Assert.AreEqual(File.ReadAllText(staticPath), File.ReadAllText(dynamicPath));

            File.Delete(dynamicPath);
            File.Delete(staticPath);
        }

        [Test]
        public static void TestWriteMs1FeaturesNullPathIsANoOp()
        {
            var results = new FlashLfqResults(new List<SpectraFileInfo>(), new List<Identification>());
            Assert.DoesNotThrow(() => results.WriteMs1Features(null, silent: true));
        }

        /// <summary>
        /// Files with no quantified peaks and peaks with no isotopic envelopes are both skipped. Runs with console
        /// output on, since that is how a caller would normally invoke it.
        /// </summary>
        [Test]
        public static void TestWriteMs1FeaturesSkipsWhatItCannotWrite()
        {
            SpectraFileInfo fileWithPeaks = WriteSyntheticMzml("ms1FeatureSkipping.mzML", new[]
            {
                new[] { 499.9, 500.0, 500.1 },
                new[] { 499.9, 500.0, 500.1 },
            });

            // this file is never read, because it has no quantified peaks - so it need not exist at all
            var fileWithoutPeaks = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory,
                "ms1FeatureNoSuchFile.mzML"), "a", 1, 0, 0);
            Assert.IsFalse(File.Exists(fileWithoutPeaks.FullFilePathWithExtension));

            var results = new FlashLfqResults(new List<SpectraFileInfo> { fileWithPeaks, fileWithoutPeaks },
                new List<Identification>());
            results.Peaks[fileWithPeaks].Add(BuildPeak(fileWithPeaks, "PEPTIDE", (500.0, 0, 1.0), (500.0, 1, 1.1)));
            results.Peaks[fileWithPeaks].Add(BuildPeak(fileWithPeaks, "NOENVELOPES"));

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeaturesSkipping.jsonl");
            Assert.DoesNotThrow(() => results.WriteMs1Features(outputPath, silent: false));

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
        public static void TestWriteMs1FeaturesWithNoMs1Scans()
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

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeatureMs2Only.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(ms2Scans), path, false);
            var file = new SpectraFileInfo(path, "a", 0, 0, 0);

            var results = new FlashLfqResults(new List<SpectraFileInfo> { file }, new List<Identification>());
            results.Peaks[file].Add(BuildPeak(file, "PEPTIDE", (500.0, 0, 1.0)));

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ms1FeaturesMs2Only.jsonl");
            Assert.DoesNotThrow(() => results.WriteMs1Features(outputPath, silent: false));

            // the file contributes nothing, leaving an empty output rather than a malformed one
            Assert.AreEqual(0, new FileInfo(outputPath).Length);

            File.Delete(outputPath);
        }

        #endregion
    }
}
