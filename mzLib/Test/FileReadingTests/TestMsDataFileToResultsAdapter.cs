using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestMsDataFileToResultsAdapter
    {
        [Test]
        [TestCase("DataFiles/small.RAW", 48, "Thermo nativeID format")]
        [TestCase("DataFiles/sliced_ethcd.raw", 6, "Thermo nativeID format")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", 142, "Thermo nativeID format")]
        [TestCase("DataFiles/tester.mzML", 7, null)]
        [TestCase("DataFiles/tester.mgf", 5, "no nativeID format")]
        public static void TestAdapterIsEquivalent(string filePath, int expectedScanCount, string sourceFormat)
        {
            List<MsDataScan> scans;
            List<MsDataScan> factoryScans;
            MsDataFile dataFile;
            MsDataFileToResultFileAdapter adapter;

            switch (Path.GetExtension(filePath).ToLower())
            {
                case ".raw":
                    dataFile = new IO.ThermoRawFileReader.ThermoRawFileReader(filePath);
                    scans = dataFile.GetAllScansList();

                    adapter = FileReader.ReadFile<MsDataFileToResultFileAdapter>(filePath);
                    factoryScans = adapter.Results;
                    Assert.That(adapter.FileType, Is.EqualTo(SupportedFileType.ThermoRaw));
                    break;

                case ".mzml":
                    dataFile = new IO.MzML.Mzml(filePath);
                    scans = dataFile.GetAllScansList();

                    adapter = FileReader.ReadFile<MsDataFileToResultFileAdapter>(filePath);
                    factoryScans = adapter.Results;
                    Assert.That(adapter.FileType, Is.EqualTo(SupportedFileType.MzML));
                    break;

                case ".mgf":
                    dataFile = new IO.Mgf.Mgf(filePath);
                    scans = dataFile.GetAllScansList();

                    adapter = FileReader.ReadFile<MsDataFileToResultFileAdapter>(filePath);
                    factoryScans = adapter.Results;
                    Assert.That(adapter.FileType, Is.EqualTo(SupportedFileType.Mgf));
                    break;

                default:
                    throw new MzLibException("File type not needed to test for backwards compatibility");
            }

            var sourceFile = dataFile.GetSourceFile();
            Assert.That(scans.Count, Is.EqualTo(expectedScanCount));
            Assert.That(sourceFile.NativeIdFormat, Is.EqualTo(sourceFormat));

            sourceFile = adapter.GetSourceFile();
            Assert.That(factoryScans.Count, Is.EqualTo(expectedScanCount));
            Assert.That(sourceFile.NativeIdFormat, Is.EqualTo(sourceFormat));
            Assert.That(filePath, Is.EqualTo(adapter.FilePath));
            Assert.That(adapter.Software, Is.EqualTo(Software.MassSpecFile));
        }


        [Test]
        public static void TestDynamicConnection()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var path1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var factoryConstruction = FileReader.ReadFile<MsDataFileToResultFileAdapter>(path1);
            factoryConstruction.InitiateDynamicConnection();

            var path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "testFileWMS2.raw");
            var constructorConstruction = new MsDataFileToResultFileAdapter(path2);
            constructorConstruction.LoadResults();
            constructorConstruction.InitiateDynamicConnection();


            var a = factoryConstruction.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(a != null);

            var b = constructorConstruction.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(b != null);

            Assert.That(a.MassSpectrum.XArray.Length != b.MassSpectrum.XArray.Length);

            a = factoryConstruction.GetOneBasedScanFromDynamicConnection(10000);
            factoryConstruction.CloseDynamicConnection();
            constructorConstruction.CloseDynamicConnection();

            Console.WriteLine($"Analysis time for TestDynamicConnectionRawFileReader: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
            Assert.That(a == null);
        }

        [Test]
        [TestCase("small.raw", "a.mzML",  0)]
        [TestCase("small.raw", "a.mzML",  1)]
        public void TestWriting(string filePath, string outfile, int loop)
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", filePath);
            outfile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", outfile);

            MsDataFileToResultFileAdapter dataFile;
            if (loop == 0)
                dataFile = new MsDataFileToResultFileAdapter(path);
            else if (loop == 1)
                dataFile = FileReader.ReadFile<MsDataFileToResultFileAdapter>(path);
            else
                throw new ArgumentOutOfRangeException();
            dataFile.LoadResults();

            if (File.Exists(outfile))
                File.Delete(outfile);
            dataFile.WriteResults(outfile);
            Assert.That(File.Exists(outfile));

            MsDataFile readInFile = MsDataFileReader.GetDataFile(outfile).LoadAllStaticData();
            Assert.That(readInFile.Scans.Length, Is.EqualTo(dataFile.Results.Count));
            for (var i = 0; i < readInFile.Scans.Length; i++)
            {
                var readInScan = readInFile.Scans[i];
                var writtenScan = dataFile.Results[i];
                Assert.That(readInScan.OneBasedScanNumber.Equals(writtenScan.OneBasedScanNumber));
                Assert.That(readInScan.MsnOrder.Equals(writtenScan.MsnOrder));
                Assert.That(readInScan.IsCentroid.Equals(writtenScan.IsCentroid));
                Assert.That(readInScan.MassSpectrum.Equals(writtenScan.MassSpectrum));
            }

            File.Delete(outfile);
            Assert.That(!File.Exists(outfile));
        }

        [Test]
        public void TestWritingNoExtension()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var outfile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "a");
            var dataFile = FileReader.ReadFile<MsDataFileToResultFileAdapter>(path);
            dataFile.LoadAllStaticData();

            if (File.Exists(outfile))
                File.Delete(outfile);
            dataFile.WriteResults(outfile);
            Assert.That(File.Exists(outfile + ".mzML"));
            File.Delete(outfile + ".mzML");
        }
    }
}
