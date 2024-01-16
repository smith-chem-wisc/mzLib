using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
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

            //MsDataFile brukerData = MsDataFileReader.GetDataFile(@"C:\Users\Alex\Downloads\transfer_292991_files_907ddd5f\data_files\T03797_AurEl3_trap1_CMB-1380_1_GC1_1_4093.mzML").LoadAllStaticData();
            string filePath = @"C:\Users\Alex\Documents\timsTOF Data\timsTOF_User_Example_file\data_files\T03797_AurEl3_trap1_CMB-1380_1_GC1_1_4093.d";
            var test = new TimsTofFileReader(filePath).LoadAllStaticData(maxThreads: 10);
            //TimsTofFileReader testAsTimmyReader = (TimsTofFileReader)test;
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"C:\Users\Alex\Documents\timsTOF Data\timsTOF_User_Example_file\data_files\MS1_Parallel10_500_plusMS2_Parallel_SpectraMergeChange.txt");
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
            var test = new ThermoRawFileReader(filePath).LoadAllStaticData();
            watch.Stop();
            var elapsedTimeSecond = watch.ElapsedMilliseconds / 1000;

            StreamWriter output = new StreamWriter(@"C:\Users\Alex\Documents\timsTOF Data\timsTOF_User_Example_file\data_files\rawRuntime.txt");
            using (output)
            {
                output.WriteLine(elapsedTimeSecond.ToString() + " seconds to read in the file.");
                output.WriteLine(test.Scans.Length + " scans read from file.");
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

            MzSpectrum outSpectrum = TofSpectraMerger.MergeSpectra(
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

            MzSpectrum outSpectrum = TofSpectraMerger.MergeSpectra(
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

            MzSpectrum outSpectrum = TofSpectraMerger.MergeSpectra(
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

            MzSpectrum outSpectrum = TofSpectraMerger.MergeSpectra(
                new List<double[]> { mz1, mz2, mz3 },
                new List<int[]> { intensity1, intensity2, intensity3 });

            Assert.AreEqual(outSpectrum.Size, 11);
            CollectionAssert.AreEqual(outSpectrum.XArray, new double[] { 1 + 5e-7, 2 + 5e-7, 3, 4, 5, 6, 7, 8, 9, 10, 11 + 1e-6 });
            CollectionAssert.AreEqual(outSpectrum.YArray, new double[] { 11, 12, 3, 4, 5, 6, 7, 8, 9, 10, 11 });
        }
    }
}
