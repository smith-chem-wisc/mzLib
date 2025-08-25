using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using MzLibUtil;
using System.Linq;
using System.Collections.Concurrent;
using System;
using System.Threading.Tasks;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestReaders
    {
        [Test]
        [TestCase("DataFiles/small.RAW", 48, "Thermo nativeID format")]
        [TestCase("DataFiles/sliced_ethcd.raw", 6, "Thermo nativeID format")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", 142, "Thermo nativeID format")]
        [TestCase("DataFiles/tester.mzML", 7, null)]
        [TestCase("DataFiles/tester.mgf", 5, "no nativeID format")]
        public static void TestLoadingRawFilesAndSourceFiles(string filePath, int expectedScanCount, string sourceFormat)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            MsDataFile datafile = MsDataFileReader.GetDataFile(spectraPath);
            List<MsDataScan> scans = datafile.GetAllScansList();
            Assert.That(scans.Count == expectedScanCount);

            SourceFile file = datafile.GetSourceFile();
            Assert.That(file.NativeIdFormat == sourceFormat);
        }

        [Test]
        [TestCase("DataFiles/small.RAW", 48, "Thermo nativeID format")]
        [TestCase("DataFiles/sliced_ethcd.raw", 6, "Thermo nativeID format")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", 142, "Thermo nativeID format")]
        [TestCase("DataFiles/tester.mzML", 7, null)]
        [TestCase("DataFiles/tester.mgf", 5, "no nativeID format")]
        public static void EnsureBackwardsCompatibility(string filePath, int expectedScanCount, string sourceFormat)
        {
            List<MsDataScan> scans;
            MsDataFile dataFile;
            switch (Path.GetExtension(filePath).ToLower())
            {
                case ".raw":
                    dataFile = IO.ThermoRawFileReader.ThermoRawFileReader.LoadAllStaticData(filePath);
                    scans = dataFile.GetAllScansList();
                    break;

                case ".mzml":
                    dataFile = IO.MzML.Mzml.LoadAllStaticData(filePath);
                    scans = dataFile.GetAllScansList();
                    break;

                case ".mgf":
                    dataFile = IO.Mgf.Mgf.LoadAllStaticData(filePath);
                    scans = dataFile.GetAllScansList();
                    break;

                default:
                    throw new MzLibException("File type not needed to test for backwards compatibility");
            }

            var sourceFile = dataFile.GetSourceFile();
            Assert.That(scans.Count, Is.EqualTo(expectedScanCount));
            Assert.That(sourceFile.NativeIdFormat, Is.EqualTo(sourceFormat));
        }

        [Test]
        [TestCase("DataFiles/small.RAW", 48, "Thermo nativeID format")]
        [TestCase("DataFiles/sliced_ethcd.raw", 6, "Thermo nativeID format")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", 142, "Thermo nativeID format")]
        [TestCase("DataFiles/tester.mzML", 7, null)]
        [TestCase("DataFiles/tester.mgf", 5, "no nativeID format")]
        public static void EnsureBackwardsCompatibilityConstrucor(string filePath, int expectedScanCount, string sourceFormat)
        {
            List<MsDataScan> scans;
            MsDataFile dataFile;
            switch (Path.GetExtension(filePath).ToLower())
            {
                case ".raw":
                    dataFile = new IO.ThermoRawFileReader.ThermoRawFileReader(filePath);
                    scans = dataFile.GetAllScansList();
                    break;

                case ".mzml":
                    dataFile = new IO.MzML.Mzml(filePath);
                    scans = dataFile.GetAllScansList();
                    break;
                
                case ".mgf":
                    dataFile = new IO.Mgf.Mgf(filePath);
                    scans = dataFile.GetAllScansList();
                    break;

                default:
                    throw new MzLibException("File type not needed to test for backwards compatibility");
            }

            var sourceFile = dataFile.GetSourceFile();
            Assert.That(scans.Count, Is.EqualTo(expectedScanCount));
            Assert.That(sourceFile.NativeIdFormat, Is.EqualTo(sourceFormat));
        }

        [TestCase("DataFiles/small.RAW")]
        [TestCase("DataFiles/sliced_ethcd.raw")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml")]
        [TestCase("DataFiles/tester.mzML")]
        [TestCase("DataFiles/tester.mgf")]
        public static void TestDynamicConnection_IsThreadSafe(string filePath)
        {
            int numThreads = 8;
            var dataFile = MsDataFileReader.GetDataFile(filePath).LoadAllStaticData();
            var scanNumbers = dataFile.Scans.Select(scan => scan.OneBasedScanNumber).ToList();

            var exceptions = new ConcurrentBag<Exception>();
            var results = new ConcurrentDictionary<int, MsDataScan>();

            dataFile.InitiateDynamicConnection();

            Parallel.ForEach(scanNumbers, new ParallelOptions { MaxDegreeOfParallelism = numThreads }, scanNumber =>
            {
                try
                {
                    var scan = dataFile.GetOneBasedScanFromDynamicConnection(scanNumber);
                    results[scanNumber] = scan;
                }
                catch (Exception ex)
                {
                    exceptions.Add(ex);
                }
            });

            dataFile.CloseDynamicConnection();

            Assert.That(exceptions, Is.Empty, "Exceptions occurred during concurrent scan reading.");
            Assert.That(results.Count, Is.EqualTo(scanNumbers.Count), "Not all scans were read successfully.");
            Assert.That(results.Values.All(scan => scan != null), Is.True, "Null scan(s) returned.");
        }

    }
}
