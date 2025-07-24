using MassSpectrometry;
using MathNet.Numerics;
using MzLibUtil;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestTimsTofFileReader
    {
        // The timsTOF_snippet.d contains DDA_PASEF data
        public string _testDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "timsTOF_snippet.d");
        public TimsTofFileReader _testReader;
        public TimsDataScan _testMs2Scan;
        public TimsDataScan _testMs1Scan;
        public FilteringParams _filteringParams = new FilteringParams(numberOfPeaksToKeepPerWindow:200, minimumAllowedIntensityRatioToBasePeak: 0.01);

        [OneTimeSetUp]
        public void SetUp()
        {
            _testReader = new TimsTofFileReader(_testDataPath);
            _testReader.LoadAllStaticData(filteringParams: _filteringParams, maxThreads: 10);
            _testMs2Scan = (TimsDataScan)_testReader.Scans.Skip(1000).First(scan => scan.MsnOrder > 1);
            _testMs1Scan = (TimsDataScan)_testReader.Scans.Skip(500).First(scan => scan.MsnOrder == 1);
        }

        [Test]
        public void TestReadForMrmFile()
        {
            // This mrm file has been modified! There are actually 6437 scans per frame in this file. However, that takes a while to read in,
            // merge, and centroid. So I edited the SQL database (contained in the .tdf) file to show that each frame only contains 100 scans)
            // This resulted in a ~20x speed up
            // This probably implies that the reader should handle the data fundamentally differently (e.g., chunking up each frame into 100 scan segments)
            string mrmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "timsTOF_MRM.d");
            var reader = new TimsTofFileReader(mrmFilePath);
            reader.LoadAllStaticData(filteringParams: _filteringParams, maxThreads: 10);
            Assert.That(reader.NumSpectra, Is.EqualTo(909), "Number of spectra in the MRM file is not as expected.");

            reader.MrmScanArray[10] = null;
            reader.AssignScanNumbersToMrmScans();
            Assert.That(reader.Scans[10], Is.Not.Null); // Null scan should be removed in AssignScanNumbersToMrmScans
            Assert.That(reader.NumSpectra, Is.EqualTo(908));
            Assert.That(reader.Scans[^1].OneBasedScanNumber, Is.EqualTo(908));
            Assert.That(reader.Scans[0].IsolationWidth, Is.EqualTo(5).Within(0.01), "Isolation width of first scan is not as expected.");
            Assert.That(reader.Scans[0].IsolationMz, Is.EqualTo(627.52).Within(0.01), "Selected ion m/z of first scan is not as expected.");
            Assert.That(reader.Scans[0].HcdEnergy, Is.EqualTo("28"));

            reader.MrmScanArray = new TimsDataScan[reader.NumSpectra];
            reader.AssignScanNumbersToMrmScans();
            Assert.Pass(); // Make sure empty scan array doesn't crash
        }

        [Test]
        public static void TestDiaFileThrowsAppropriateException()
        {
            string diaFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "timsTOF_DIA.d");
            var reader = new TimsTofFileReader(diaFilePath);
            Assert.That(
                () => reader.LoadAllStaticData(filteringParams: new FilteringParams()),
                Throws.TypeOf<MzLibException>()
                    .With.Message.EqualTo("The timsTOF file contains unsupported scan mode: DIA. Only DDA-PASEF and MRM data is supported at this time.")
            );
        }

        [Test]
        public void TestDynamicConnection()
        {
            _testReader.InitiateDynamicConnection();
            _testReader.CloseDynamicConnection();
            _testReader.InitiateDynamicConnection();
            Assert.Pass();
        }

        [Test, NonParallelizable]
        public void TestCorruptedDataDoesntThrowException()
        {
            var tempReader = new TimsTofFileReader(_testDataPath);
            try
            {
                tempReader.InitiateDynamicConnection();
                FrameProxy frame = tempReader.FrameProxyFactory.GetFrameProxy(100);
                var scanOffsetsProperty = typeof(FrameProxy).GetField("_scanOffsets", BindingFlags.NonPublic | BindingFlags.Instance);
                int[] corruptedOffsets = new int[frame.NumberOfScans];
                int bigNumber = int.MaxValue / 3; // Use a large number to simulate corruption
                for (int i = 0; i < corruptedOffsets.Length; i++)
                {
                    corruptedOffsets[i] = bigNumber; // Set all offsets to a large number to simulate corruption
                }
                scanOffsetsProperty.SetValue(frame, corruptedOffsets); // Simulate corrupted data by setting scan offsets to null
                Assert.That(!frame.IsFrameValid());
                Assert.Throws<ArgumentException>(() => frame.GetScanIntensities(100));
            }
            finally
            {
                tempReader.CloseDynamicConnection();
            }
        }

        [Test, NonParallelizable]
        public void TestCorruptedDataProducesCorrectWarningMessage()
        {
            var tempReader = new TimsTofFileReader(_testDataPath);
            try
            {
                tempReader.InitiateDynamicConnection();
                var mockFactory = new MockFrameProxyFactory(tempReader.FrameProxyFactory);
                var x = typeof(TimsTofFileReader).GetProperty("FrameProxyFactory", BindingFlags.NonPublic | BindingFlags.Instance);
                Assert.That(x, Is.Not.Null, "FrameProxyFactory property is null. This should not happen.");
                x.SetValue(tempReader, mockFactory);
                //typeof(TimsTofFileReader).GetProperty("FrameProxyFactory", BindingFlags.NonPublic | BindingFlags.Instance).SetValue(tempReader, mockFactory);

                tempReader.BuildDDAScans(frameId: 1, null);
                var scan = new TimsDataScan(massSpectrum: null,
                    oneBasedScanNumber: -1, // This will be adjusted once all scans have been read
                    msnOrder: 2,
                    isCentroid: true,
                    polarity: Polarity.Positive,
                    retentionTime: 1,
                    scanWindowRange: new MzLibUtil.MzRange(200, 2000),
                    scanFilter: "xyz",
                    mzAnalyzer: MZAnalyzerType.TOF,
                    totalIonCurrent: -1, // Will be set later
                    injectionTime: 1,
                    noiseData: null,
                    nativeId: "frames=",
                    frameId: 2,
                    scanNumberStart: 1,
                    scanNumberEnd:900,
                    medianOneOverK0: 1, // Needs to be set later
                    precursorId: 1,
                    selectedIonMz:1,
                    selectedIonChargeStateGuess: 1,
                    selectedIonIntensity: 1,
                    isolationMZ: 1,
                    isolationWidth: 1,
                    dissociationType: DissociationType.CID,
                    oneBasedPrecursorScanNumber: -1, // This will be set later
                    selectedIonMonoisotopicGuessMz: 1,
                    hcdEnergy: "1",
                    frames:new List<long> { 2 });
                tempReader.PopulateSpectraForPasefScans(new List<TimsDataScan>() { scan }, new List<long> { 2 }, null);

                Assert.That(tempReader.FaultyFrameIds.Contains(1));
                Assert.That(tempReader.FaultyFrameIds.Contains(2));
                Assert.That(tempReader.Warnings, Is.EqualTo("The following frames were not read correctly: 1, 2"));
            }
            finally
            {
                tempReader.CloseDynamicConnection();
            }
        }

        [Test, NonParallelizable]
        public void TestGetPasefScanFromDynamicConnectionUsingFrameId()
        {
            var dynamicReader = new TimsTofFileReader(_testDataPath);
            dynamicReader.InitiateDynamicConnection();
            var dynamicScan = dynamicReader
                .GetScanFromPrecursorAndFrameIdFromDynamicConnection((int)_testMs2Scan.PrecursorId, (int)_testMs2Scan.FrameId, _filteringParams);
            Assert.IsNotNull(dynamicScan);

            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs2Scan.PrecursorId), "PrecursorId values are not equal.");
            Assert.That(dynamicScan.ScanNumberStart, Is.EqualTo(_testMs2Scan.ScanNumberStart), "ScanStart values are not equal.");
            Assert.That(dynamicScan.ScanNumberEnd, Is.EqualTo(_testMs2Scan.ScanNumberEnd), "ScanEnd values are not equal.");
            Assert.That(dynamicScan.OneOverK0, Is.EqualTo(_testMs2Scan.OneOverK0), "ScanMedian values are not equal.");
            Assert.That(dynamicScan.IsolationMz, Is.EqualTo(_testMs2Scan.IsolationMz), "IsolationMz values are not equal.");
            Assert.That(dynamicScan.IsolationWidth, Is.EqualTo(_testMs2Scan.IsolationWidth), "IsolationWidth values are not equal.");
            Assert.That(dynamicScan.HcdEnergy, Is.EqualTo(_testMs2Scan.HcdEnergy), "CollisionEnergy values are not equal.");
            Assert.That(dynamicScan.SelectedIonMZ, Is.EqualTo(_testMs2Scan.SelectedIonMZ), "MostAbundantPrecursorMz values are not equal.");
            Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz, Is.EqualTo(_testMs2Scan.SelectedIonMonoisotopicGuessMz), "PrecursorMonoisotopicMz values are not equal.");
            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs2Scan.PrecursorId), "PrecursorID values are not equal.");
            Assert.That(dynamicScan.MassSpectrum, Is.EqualTo(_testMs2Scan.MassSpectrum), "Mass spectra are not equal");
        }

        [Test, NonParallelizable]
        public void TestGetMs1ScanFromDynamicConnectionUsingFrameId()
        {
            var dynamicReader = new TimsTofFileReader(_testDataPath);
            dynamicReader.InitiateDynamicConnection();
            var dynamicScan = dynamicReader
                .GetScanFromPrecursorAndFrameIdFromDynamicConnection((int)_testMs1Scan.PrecursorId, (int)_testMs1Scan.FrameId, _filteringParams);
            Assert.IsNotNull(dynamicScan);

            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs1Scan.PrecursorId), "PrecursorId values are not equal.");
            Assert.That(dynamicScan.ScanNumberStart, Is.EqualTo(_testMs1Scan.ScanNumberStart), "ScanStart values are not equal.");
            Assert.That(dynamicScan.ScanNumberEnd, Is.EqualTo(_testMs1Scan.ScanNumberEnd), "ScanEnd values are not equal.");
            Assert.That(dynamicScan.OneOverK0, Is.EqualTo(_testMs1Scan.OneOverK0), "ScanMedian values are not equal.");
            Assert.That(dynamicScan.IsolationMz, Is.EqualTo(_testMs1Scan.IsolationMz), "IsolationMz values are not equal.");
            Assert.That(dynamicScan.IsolationWidth, Is.EqualTo(_testMs1Scan.IsolationWidth), "IsolationWidth values are not equal.");
            Assert.That(dynamicScan.HcdEnergy, Is.EqualTo(_testMs1Scan.HcdEnergy), "CollisionEnergy values are not equal.");
            Assert.That(dynamicScan.SelectedIonMZ, Is.EqualTo(_testMs1Scan.SelectedIonMZ), "MostAbundantPrecursorMz values are not equal.");
            Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz, Is.EqualTo(_testMs1Scan.SelectedIonMonoisotopicGuessMz), "PrecursorMonoisotopicMz values are not equal.");
            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs1Scan.PrecursorId), "PrecursorID values are not equal.");
            Assert.That(dynamicScan.MassSpectrum, Is.EqualTo(_testMs1Scan.MassSpectrum), "Mass spectra are not equal");
        }

        [Test, NonParallelizable]
        public void TestGetScanFromDynamicConnectionUsingOneBasedScanNumber()
        {
            var dynamicReader = new TimsTofFileReader(_testDataPath);
            dynamicReader.InitiateDynamicConnection();
            var scanBeforeCast = dynamicReader.GetOneBasedScanFromDynamicConnection(_testMs1Scan.OneBasedScanNumber, _filteringParams);
            var dynamicScan = scanBeforeCast as TimsDataScan;
            Assert.IsNotNull(dynamicScan);

            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs1Scan.PrecursorId), "PrecursorId values are not equal.");
            Assert.That(dynamicScan.ScanNumberStart, Is.EqualTo(_testMs1Scan.ScanNumberStart), "ScanStart values are not equal.");
            Assert.That(dynamicScan.ScanNumberEnd, Is.EqualTo(_testMs1Scan.ScanNumberEnd), "ScanEnd values are not equal.");
            Assert.That(dynamicScan.OneOverK0, Is.EqualTo(_testMs1Scan.OneOverK0), "ScanMedian values are not equal.");
            Assert.That(dynamicScan.IsolationMz, Is.EqualTo(_testMs1Scan.IsolationMz), "IsolationMz values are not equal.");
            Assert.That(dynamicScan.IsolationWidth, Is.EqualTo(_testMs1Scan.IsolationWidth), "IsolationWidth values are not equal.");
            Assert.That(dynamicScan.HcdEnergy, Is.EqualTo(_testMs1Scan.HcdEnergy), "CollisionEnergy values are not equal.");
            Assert.That(dynamicScan.SelectedIonMZ, Is.EqualTo(_testMs1Scan.SelectedIonMZ), "MostAbundantPrecursorMz values are not equal.");
            Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz, Is.EqualTo(_testMs1Scan.SelectedIonMonoisotopicGuessMz), "PrecursorMonoisotopicMz values are not equal.");
            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs1Scan.PrecursorId), "PrecursorID values are not equal.");
            Assert.That(dynamicScan.MassSpectrum, Is.EqualTo(_testMs1Scan.MassSpectrum), "Mass spectra are not equal");
            Assert.That(dynamicScan.OneBasedScanNumber, Is.EqualTo(_testMs1Scan.OneBasedScanNumber));


            scanBeforeCast = dynamicReader.GetOneBasedScanFromDynamicConnection(_testMs2Scan.OneBasedScanNumber, _filteringParams);
            dynamicScan = scanBeforeCast as TimsDataScan;
            Assert.IsNotNull(dynamicScan);

            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs2Scan.PrecursorId), "PrecursorId values are not equal.");
            Assert.That(dynamicScan.ScanNumberStart, Is.EqualTo(_testMs2Scan.ScanNumberStart), "ScanStart values are not equal.");
            Assert.That(dynamicScan.ScanNumberEnd, Is.EqualTo(_testMs2Scan.ScanNumberEnd), "ScanEnd values are not equal.");
            Assert.That(dynamicScan.OneOverK0, Is.EqualTo(_testMs2Scan.OneOverK0), "ScanMedian values are not equal.");
            Assert.That(dynamicScan.IsolationMz, Is.EqualTo(_testMs2Scan.IsolationMz), "IsolationMz values are not equal.");
            Assert.That(dynamicScan.IsolationWidth, Is.EqualTo(_testMs2Scan.IsolationWidth), "IsolationWidth values are not equal.");
            Assert.That(dynamicScan.HcdEnergy, Is.EqualTo(_testMs2Scan.HcdEnergy), "CollisionEnergy values are not equal.");
            Assert.That(dynamicScan.SelectedIonMZ, Is.EqualTo(_testMs2Scan.SelectedIonMZ), "MostAbundantPrecursorMz values are not equal.");
            Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz, Is.EqualTo(_testMs2Scan.SelectedIonMonoisotopicGuessMz), "PrecursorMonoisotopicMz values are not equal.");
            Assert.That(dynamicScan.PrecursorId, Is.EqualTo(_testMs2Scan.PrecursorId), "PrecursorID values are not equal.");
            Assert.That(dynamicScan.MassSpectrum, Is.EqualTo(_testMs2Scan.MassSpectrum), "Mass spectra are not equal");
            Assert.That(dynamicScan.OneBasedScanNumber, Is.EqualTo(_testMs2Scan.OneBasedScanNumber));
        }

        [Test]
        public void TestTwoPointerMerge()
        {
            uint[] indices1 = new uint[] { 1, 3, 5, 7, 9, 11 };
            uint[] indices2 = new uint[] { 0, 2, 4, 6, 8, 10 };

            int[] intensities1 = new int[] { 1, 3, 5, 7, 9, 11 };
            int[] intensities2 = new int[] { 0, 2, 4, 6, 8, 10 };

            int[] intendedOutput = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

            var mergerOutput = TofSpectraMerger.TwoPointerMerge(indices1, indices2, intensities1, intensities2);

            Assert.That(mergerOutput.Intensities, Is.EqualTo(intendedOutput));
            Assert.That(mergerOutput.Indices.Select(i => (int)i).ToArray(), Is.EqualTo(intendedOutput));

            indices2 = new uint[] { 0, 2, 4, 6, 8, 10, 12, 13, 14, 15, 16 };

            intensities2 = new int[] { 0, 2, 4, 6, 8, 10, 12, 13, 14, 15, 16 };

            intendedOutput = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };

            mergerOutput = TofSpectraMerger.TwoPointerMerge(indices1, indices2, intensities1, intensities2);

            Assert.That(mergerOutput.Intensities, Is.EqualTo(intendedOutput));
            Assert.That(mergerOutput.Indices.Select(i => (int)i).ToArray(), Is.EqualTo(intendedOutput));
        }

        [Test]
        public void TestCollapse()
        {
            uint[] indices = new uint[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
            int[] intensities = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

            List<int> intendedIdx = new List<int> { 1, 4, 7, 10 };
            List<int> intendedIntensities = new List<int> { 3, 12, 21, 30 };

            var collapsedOutput = TofSpectraMerger.CollapseArrays(indices, intensities);

            Assert.That(collapsedOutput.Indices, Is.EqualTo(intendedIdx));
            Assert.That(collapsedOutput.Intensities, Is.EqualTo(intendedIntensities));


            indices = new uint[] { 0, 1, 2, 3, 4, 5, 6, 7,  9, 11 };
            intensities = new int[] { 0, 1, 2, 3, 4, 5, 6, 7,  9, 11 };

            intendedIdx = new List<int> { 1, 4, 6, 9 };
            intendedIntensities = new List<int> { 3, 12, 13, 20 };

            collapsedOutput = TofSpectraMerger.CollapseArrays(indices, intensities);

            Assert.That(collapsedOutput.Indices, Is.EqualTo(intendedIdx));
            Assert.That(collapsedOutput.Intensities, Is.EqualTo(intendedIntensities));

            indices = new uint[] { 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 11, 18, 523, 1000, 1000, 1000 };
            intensities = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 11, 18, 523, 1000, 1000, 1000 };

            intendedIdx = new List<int> { 1, 4, 6, 11, 18, 523, 1000 };
            intendedIntensities = new List<int> { 3, 12, 13, 31, 18, 523, 3000 };

            collapsedOutput = TofSpectraMerger.CollapseArrays(indices, intensities);

            Assert.That(collapsedOutput.Indices, Is.EqualTo(intendedIdx));
            Assert.That(collapsedOutput.Intensities, Is.EqualTo(intendedIntensities));
        }

        [Test]
        public void TestConstructor()
        {
            var reader = MsDataFileReader.GetDataFile(_testDataPath);
            Assert.That(reader, !Is.Null);
        }

        [Test]
        public void TestFileDoesntExist()
        {
            string fakePath = "fakePath.d";
            Assert.Throws<FileNotFoundException>(() =>
                MsDataFileReader.GetDataFile(fakePath));

            TimsTofFileReader reader = new TimsTofFileReader(fakePath);

            Assert.Throws<FileNotFoundException>(() =>
              reader.LoadAllStaticData());
        }


        [Test]
        public void TestLoadAllStaticData()
        {
            Assert.That(_testReader.NumSpectra, Is.EqualTo(4096));

            Assert.That(_testMs2Scan.Polarity == Polarity.Positive);
            Assert.That(_testMs2Scan.DissociationType == DissociationType.CID);
            Assert.That(_testMs2Scan.TotalIonCurrent == 25130);
            Assert.That(_testMs2Scan.NativeId == "frames=64-64;scans=410-435");
            Assert.That(_testMs2Scan.SelectedIonMZ, Is.EqualTo(739.3668).Within(0.001));
            Assert.That(_testMs2Scan.MsnOrder == 2);
            Assert.That(_testMs2Scan.IsCentroid);
            Assert.That(_testMs2Scan.ScanNumberStart == 410);
            Assert.That(_testMs2Scan.OneOverK0, Is.EqualTo(1.0424).Within(0.0001));
            Assert.That(_testReader.Scans.All(s => s != null));
            Assert.That(_testMs1Scan.RetentionTime, Is.EqualTo(11.28).Within(0.01));
        }

        [Test]
        public void TestOneBasedPrecursor()
        {
            TimsDataScan ms1Scan = (TimsDataScan)_testReader.GetOneBasedScan((int)_testMs2Scan.OneBasedPrecursorScanNumber);

            Assert.AreEqual(_testMs2Scan.PrecursorId, ms1Scan.PrecursorId);
            // Check that the child and parent scan are both looking at the same timsScans (i.e., the same region in the ion-mobility dimension)
            Assert.AreEqual(_testMs2Scan.ScanNumberStart, ms1Scan.ScanNumberStart);
            Assert.AreEqual(_testMs2Scan.ScanNumberEnd, ms1Scan.ScanNumberEnd);
            Assert.AreEqual(_testMs2Scan.OneOverK0, ms1Scan.OneOverK0);

        }

        [Test]
        public void TestSpectraMerger()
        {
            double[] mz1 = new double[] { 1, 3, 5, 7, 9 };
            double[] mz2 = new double[] { 2, 4, 6, 8, 10 };

            int[] intensity1 = new int[] { 1, 3, 5, 7, 9 };
            int[] intensity2 = new int[] { 2, 4, 6, 8, 10 };

            MzSpectrum outSpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(
                new List<double[]> { mz1, mz2 },
                new List<int[]> { intensity1, intensity2 });

            Assert.AreEqual(outSpectrum.Size, 10);
            CollectionAssert.AreEqual(outSpectrum.XArray, new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
        }

        [Test]
        public void TestSpectraMerger2()
        {
            double[] mz1 = new double[] { 1, 3, 5, 7, 9, 10 };
            double[] mz2 = new double[] { 2, 4, 6, 8, 10 };

            int[] intensity1 = new int[] { 1, 3, 5, 7, 9, 10 };
            int[] intensity2 = new int[] { 2, 4, 6, 8, 10 };

            MzSpectrum outSpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(
                new List<double[]> { mz1, mz2 },
                new List<int[]> { intensity1, intensity2 });

            Assert.AreEqual(outSpectrum.Size, 10);
            CollectionAssert.AreEqual(outSpectrum.XArray, new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
            CollectionAssert.AreEqual(outSpectrum.YArray, new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 20 });
        }

        [Test]
        public void TestSpectraMerger3()
        {
            double[] mz1 = new double[] { 1, 4, 7, 10 };
            double[] mz2 = new double[] { 2, 5, 8 };
            double[] mz3 = new double[] { 3, 6, 9 };

            int[] intensity1 = new int[] { 1, 4, 7, 10 };
            int[] intensity2 = new int[] { 2, 5, 8 };
            int[] intensity3 = new int[] { 3, 6, 9 };

            MzSpectrum outSpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(
                new List<double[]> { mz1, mz2, mz3 },
                new List<int[]> { intensity1, intensity2, intensity3 });

            Assert.AreEqual(outSpectrum.Size, 10);
            CollectionAssert.AreEqual(outSpectrum.XArray, new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
            CollectionAssert.AreEqual(outSpectrum.YArray, new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
        }

        // Test that weighted averaging works when two peaks are close together
        [Test]
        public void TestSpectraMerger4()
        {
            double[] mz1 = new double[] { 1, 3, 5, 7, 9 };
            double[] mz2 = new double[] { 2, 4, 6, 8, 10 };
            double[] mz3 = new double[] { 1 + 1e-6, 2 + 1e-6, 11 + 1e-6 };
 
            int[] intensity1 = new int[] { 1, 3, 5, 7, 9 };
            int[] intensity2 = new int[] { 2, 4, 6, 8, 10 };
            int[] intensity3 = new int[] { 10, 10, 11 };

            MzSpectrum outSpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(
                new List<double[]> { mz1, mz2, mz3 },
                new List<int[]> { intensity1, intensity2, intensity3 });

            Assert.AreEqual(outSpectrum.Size, 11);
            // Peaks (mz = 1, intensity = 1) and (mz = 1+1e-6, intensity = 10) are close together, so they should be averaged
            // Same thing for (mz = 2, intensity = 2) and (mz = 2+1e-6, intensity = 10) 
            CollectionAssert.AreEqual(outSpectrum.XArray.Select(mz => mz.Round(7)).ToArray(),
                new double[] { 1 + 9e-7, 2 + 8e-7, 3, 4, 5, 6, 7, 8, 9, 10, 11 + 1e-6 });
            CollectionAssert.AreEqual(outSpectrum.YArray, new double[] { 11, 12, 3, 4, 5, 6, 7, 8, 9, 10, 11 });
        }
    }

    internal class MockFrameProxyFactory : FrameProxyFactory
    {
        public MockFrameProxyFactory(FrameProxyFactory realFactory) : base(realFactory.FramesTable, realFactory.FileHandle, realFactory.FileLock, realFactory.MaxIndex)
        {
            // Mock constructor does not need to do anything special
        }
        internal override FrameProxy GetFrameProxy(long frameId)
        {
            return null;
        }
    }
}
