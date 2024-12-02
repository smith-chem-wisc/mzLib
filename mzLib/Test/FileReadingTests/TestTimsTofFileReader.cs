using IO.ThermoRawFileReader;
using MassSpectrometry;
using MathNet.Numerics;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Readers;
using Readers.Bruker;
using Readers.Bruker;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestTimsTofFileReader
    {

        public string _testDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "timsTOF_snippet.d");
        public TimsTofFileReader _testReader;

        [Test]
        public void SetUp()
        {
            DateTime time = DateTime.Now;

            _testReader = new TimsTofFileReader(_testDataPath);
            _testReader.LoadAllStaticData();

            DateTime time2 = DateTime.Now;
            var timeDiff = time2 - time;
            int placeholder = 0;
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
        }

        [Test]
        public static void TestThreeMinuteReader()
        {
            Stopwatch watch = Stopwatch.StartNew();
            var memoryBeforeReading = GC.GetTotalMemory(true) / 1000;

            FilteringParams filter = new FilteringParams(
                numberOfPeaksToKeepPerWindow: 200,
                minimumAllowedIntensityRatioToBasePeak: 0.005,
                windowWidthThomsons: null,
                numberOfWindows: null,
                normalizePeaksAcrossAllWindows: false,
                applyTrimmingToMs1: false,
                applyTrimmingToMsMs: true);

            string filePath = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\50ng_K562_extreme_3min.d";
            var test = new TimsTofFileReader(filePath).LoadAllStaticData(filteringParams: filter, maxThreads: 10);
            //TimsTofFileReader testAsTimmyReader = (TimsTofFileReader)test;
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\new_version_release_no_MS1centroiding.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                output.WriteLine(test.Scans.Length + " scans read from file.");
                output.WriteLine(memoryBeforeReading.ToString() + " kB used before reading.");
                GC.Collect();
                output.WriteLine((GC.GetTotalMemory(true) / 1000).ToString() + " kB used after reading.");
            }

            //test.ExportAsMzML(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\50ng_K562_extreme_3min_10ppm_Centroid.mzML", writeIndexed: true);
            //Assert.Pass();
        }

        [Test]
        public void TestLoadAllStaticData()
        {
            Assert.That(_testReader.NumSpectra, Is.EqualTo(4096));
            Assert.That(_testReader.Scans[50].Polarity == Polarity.Positive);
            Assert.That(_testReader.Scans[50].DissociationType == DissociationType.CID);
            Assert.That(_testReader.Scans[50].TotalIonCurrent == 24722);
            Assert.That(_testReader.Scans[50].NativeId == "frames=2-2;scans=354-379");
            Assert.That(_testReader.Scans[50].SelectedIonMZ, Is.EqualTo(876.915).Within(0.001));
            Assert.That(_testReader.Scans[50].MsnOrder == 2);
            Assert.That(_testReader.Scans[50].IsCentroid);
        }

        [Test]
        public void TestOneBasedPrecursor()
        {
            TimsDataScan ms2Scan = (TimsDataScan)_testReader.Scans[50];
            TimsDataScan ms1Scan = (TimsDataScan)_testReader.GetOneBasedScan((int)ms2Scan.OneBasedPrecursorScanNumber);

            Assert.AreEqual(ms2Scan.PrecursorId, ms1Scan.PrecursorId);
            // Check that the child and parent scan are both looking at the same timsScans (i.e., the same region in the ion-mobility dimension)
            Assert.AreEqual(ms2Scan.ScanNumberStart, ms1Scan.ScanNumberStart);
            Assert.AreEqual(ms2Scan.ScanNumberEnd, ms1Scan.ScanNumberEnd);

        }

        // TODO: Ask Nic about complicated test case implementations
        //[TestCase(
        //    new List<double[]> { new double[] { 1, 3, 5, 7, 9 }, new double[] { 2, 4, 6, 8, 10 } }
        //    )]
        //public void TestSpectraMerger(
        //    List<double[]> mzArrays, 
        //    List<int[]> intensityArrays, 
        //    int expectedSize, 
        //    double[] expectedXArray,
        //    double[] expectedYArray)

        [Test]
        public void TestSpectraMerger()
        {
            double[] mz1 = new double[] { 1, 3, 5, 7, 9 };
            double[] mz2 = new double[] { 2, 4, 6, 8, 10 };

            int[] intensity1 = new int[] { 1, 3, 5, 7, 9 };
            int[] intensity2 = new int[] { 2, 4, 6, 8, 10 };

            MzSpectrum outSpectrum = TofSpectraMerger.MergesMs1Spectra(
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

            MzSpectrum outSpectrum = TofSpectraMerger.MergesMs1Spectra(
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

            MzSpectrum outSpectrum = TofSpectraMerger.MergesMs1Spectra(
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

            MzSpectrum outSpectrum = TofSpectraMerger.MergesMs1Spectra(
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
}
