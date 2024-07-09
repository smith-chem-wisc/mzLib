using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Readers.Bruker;
using Readers.Bruker.TimsTofReader;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestTimsTofFileReader
    {

        [Test]
        public static void TestReadingRealLocalData()
        {
            Stopwatch watch = Stopwatch.StartNew();
            var memoryBeforeReading = GC.GetTotalMemory(true) / 1000;

            FilteringParams filter = new FilteringParams(
                numberOfPeaksToKeepPerWindow: 200,
                minimumAllowedIntensityRatioToBasePeak: 0.5,
                windowWidthThomsons: null,
                numberOfWindows: null,
                normalizePeaksAcrossAllWindows: false,
                applyTrimmingToMs1: false,
                applyTrimmingToMsMs: true);

            //MsDataFile brukerData = MsDataFileReader.GetDataFile(@"C:\Users\Alex\Downloads\transfer_292991_files_907ddd5f\data_files\T03797_AurEl3_trap1_CMB-1380_1_GC1_1_4093.mzML").LoadAllStaticData();
            //string filePath = @"C:\Users\Alex\Documents\timsTOF Data\timsTOF_User_Example_file\data_files\T03797_AurEl3_trap1_CMB-1380_1_GC1_1_4093.d";

            string filePath = @"D:\timsTof_DDA_PXD016870\c10_1_GF3_01_1619.d";
            var test = new TimsTofFileReader(filePath).LoadAllStaticData(filteringParams: filter, maxThreads: 10);
            //TimsTofFileReader testAsTimmyReader = (TimsTofFileReader)test;
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"D:\timsTof_DDA_PXD016870\c10_1_GF3_01_1619_reader_results.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                output.WriteLine(test.Scans.Length + " scans read from file.");
                output.WriteLine(memoryBeforeReading.ToString() + " kB used before reading.");
                GC.Collect();
                output.WriteLine((GC.GetTotalMemory(true) / 1000).ToString() + " kB used after reading.");
            }

            // This is failing during write, giving the following error:
            // System.InvalidOperationException : Nullable object must have a value.
            // MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(MsDataFile myMsDataFile, String outputFile, Boolean writeIndexed) line 496
            test.ExportAsMzML(@"D:\timsTof_DDA_PXD016870\c10_1_GF3_01_1619_from_mzLib_reader.mzML", writeIndexed: true);
            Assert.Pass();
            
        }

        [Test]
        public static void TestOneMinuteReader()
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

            string filePath = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min.d";
            var test = new TimsTofFileReader(filePath).LoadAllStaticData(filteringParams: filter, maxThreads: 10);
            //TimsTofFileReader testAsTimmyReader = (TimsTofFileReader)test;
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min_results_5_5ppm.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                output.WriteLine(test.Scans.Length + " scans read from file.");
                output.WriteLine(memoryBeforeReading.ToString() + " kB used before reading.");
                GC.Collect();
                output.WriteLine((GC.GetTotalMemory(true) / 1000).ToString() + " kB used after reading.");
            }

            test.ExportAsMzML(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min_5_5ppmCentroid.mzML", writeIndexed: true);
            Assert.Pass();

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

            StreamWriter output = new StreamWriter(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\50ng_K562_extreme_3min_Results_10ppm.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                output.WriteLine(test.Scans.Length + " scans read from file.");
                output.WriteLine(memoryBeforeReading.ToString() + " kB used before reading.");
                GC.Collect();
                output.WriteLine((GC.GetTotalMemory(true) / 1000).ToString() + " kB used after reading.");
            }

            test.ExportAsMzML(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\50ng_K562_extreme_3min_10ppm_Centroid.mzML", writeIndexed: true);
            Assert.Pass();
        }

        [Test]
        public static void TestNinetyMinuteReader()
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

            string filePath = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\20191021_K562_200ng_90min_Slot1-1_01_4786.d";
            var test = new TimsTofFileReader(filePath).LoadAllStaticData(filteringParams: filter, maxThreads: 10);
            //TimsTofFileReader testAsTimmyReader = (TimsTofFileReader)test;
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\20191021_K562_200ng_90min_Slot1-1_01_4786_results.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                output.WriteLine(test.Scans.Length + " scans read from file.");
                output.WriteLine(memoryBeforeReading.ToString() + " kB used before reading.");
                GC.Collect();
                output.WriteLine((GC.GetTotalMemory(true) / 1000).ToString() + " kB used after reading.");
            }

            // This is failing during write, giving the following error:
            // System.InvalidOperationException : Nullable object must have a value.
            // MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(MsDataFile myMsDataFile, String outputFile, Boolean writeIndexed) line 496
            test.ExportAsMzML(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\20191021_K562_200ng_90min_Slot1-1_01_4786_mzLlib.mzML", writeIndexed: true);
            Assert.Pass();

        }

        [Test]
        public static void TestReadingRealLocalDataTenThreads()
        {
            Stopwatch watch = Stopwatch.StartNew();
            var memoryBeforeReading = GC.GetTotalMemory(true) / 1000;

            //MsDataFile brukerData = MsDataFileReader.GetDataFile(@"C:\Users\Alex\Downloads\transfer_292991_files_907ddd5f\data_files\T03797_AurEl3_trap1_CMB-1380_1_GC1_1_4093.mzML").LoadAllStaticData();
            string filePath = @"C:\Users\Alex\Documents\timsTOF Data\timsTOF_User_Example_file\data_files\T03797_AurEl3_trap1_CMB-1380_1_GC1_1_4093.d";
            var test = new TimsTofFileReader(filePath).LoadAllStaticData(maxThreads: 10);
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"C:\Users\Alex\Documents\timsTOF Data\timsTOF_User_Example_file\data_files\runtimeTenThreads.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                output.WriteLine(test.Scans.Length + " scans read from file.");
                output.WriteLine(memoryBeforeReading.ToString() + " kB used before reading.");
                output.WriteLine((GC.GetTotalMemory(true) / 1000).ToString() + " kB used after reading.");
            }

            Assert.Pass();
        }

        [Test]
        public static void TestReadingRawFile()
        {
            Stopwatch watch = Stopwatch.StartNew();
            var memoryBeforeReading = GC.GetTotalMemory(true) / 1000;

            //MsDataFile brukerData = MsDataFileReader.GetDataFile(@"C:\Users\Alex\Downloads\transfer_292991_files_907ddd5f\data_files\T03797_AurEl3_trap1_CMB-1380_1_GC1_1_4093.mzML").LoadAllStaticData();
            string filePath = @"D:\SingleCellDataSets\Organoid\raw_files\HFL1SC_2_Healthy_CH1_I3.raw";
            //var test = new ThermoRawFileReader(filePath).LoadAllStaticData();
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"C:\Users\Alex\Documents\timsTOF Data\timsTOF_User_Example_file\data_files\rawRuntime.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                //output.WriteLine(test.Scans.Length + " scans read from file.");
                output.WriteLine(memoryBeforeReading.ToString() + " kB used before reading.");
                output.WriteLine((GC.GetTotalMemory(true) / 1000).ToString() + " kB used after reading.");
            }


            Assert.Pass();
        }

        // TODO: Ask Nic about complicated test case implementations
        //[TestCase(
        //    new List<double[]> { new double[] { 1, 3, 5, 7, 9 }, new double[] { 2, 4, 6, 8, 10 } }
        //    )]
        //public static void TestSpectraMerger(
        //    List<double[]> mzArrays, 
        //    List<int[]> intensityArrays, 
        //    int expectedSize, 
        //    double[] expectedXArray,
        //    double[] expectedYArray)

        [Test]
        public static void TestSpectraMerger()
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
        public static void TestSpectraMerger2()
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
        public static void TestSpectraMerger3()
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

        [Test]
        public static void TestSpectraMerger4()
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
            CollectionAssert.AreEqual(outSpectrum.XArray, new double[] { 1 + 5e-7, 2 + 5e-7, 3, 4, 5, 6, 7, 8, 9, 10, 11 + 1e-6 });
            CollectionAssert.AreEqual(outSpectrum.YArray, new double[] { 11, 12, 3, 4, 5, 6, 7, 8, 9, 10, 11 });
        }
    }
}
