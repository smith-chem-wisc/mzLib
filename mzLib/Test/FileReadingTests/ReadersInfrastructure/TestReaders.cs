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

namespace Test.FileReadingTests.ReadersInfrastructure
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
                    dataFile = ThermoRawFileReader.LoadAllStaticData(filePath);
                    scans = dataFile.GetAllScansList();
                    break;

                case ".mzml":
                    dataFile = Mzml.LoadAllStaticData(filePath);
                    scans = dataFile.GetAllScansList();
                    break;

                case ".mgf":
                    dataFile = Mgf.LoadAllStaticData(filePath);
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

        [Test]
        [TestCase("DataFiles/small.RAW", 48, "Thermo nativeID format")]
        [TestCase("DataFiles/sliced_ethcd.raw", 6, "Thermo nativeID format")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", 142, "Thermo nativeID format")]
        [TestCase("DataFiles/tester.mzML", 7, null)]
        public static void TestLoadingDynamicAfterStaticLoading(string filePath, int expectedScanCount, string sourceFormat)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            MsDataFile datafile = MsDataFileReader.GetDataFile(spectraPath);
            datafile.InitiateDynamicConnection();
            var dynScan = datafile.GetOneBasedScanFromDynamicConnection(1);
            datafile.CloseDynamicConnection();

            datafile.LoadAllStaticData();
            var staticScan = datafile.GetOneBasedScan(1);
            var dynScanAfterLoad = datafile.GetOneBasedScanFromDynamicConnection(1);

            Assert.That(dynScan.MassSpectrum, Is.EqualTo(dynScanAfterLoad.MassSpectrum));
            Assert.That(dynScan.MassSpectrum, Is.EqualTo(staticScan.MassSpectrum));
        }

        /// <summary>
        /// Exercises the Initiate -> Get -> Close -> Get lifecycle on the genuinely dynamic path
        /// (no prior LoadAllStaticData, so CheckIfScansLoaded() does not short-circuit the read).
        /// Covers: reading before the connection is created throws, an in-range scan is returned,
        /// an out-of-range scan number returns null (documented contract), and reading after the
        /// connection is closed throws cleanly instead of crashing on a disposed connection.
        /// </summary>
        [Test]
        [TestCase("DataFiles/small.RAW", 48)]
        [TestCase("DataFiles/sliced_ethcd.raw", 6)]
        public static void TestDynamicConnectionLifecycle(string filePath, int scanCount)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            var dataFile = MsDataFileReader.GetDataFile(spectraPath);

            // Reading before InitiateDynamicConnection must throw - the connection has not been created.
            Assert.Throws<MzLibException>(() => dataFile.GetOneBasedScanFromDynamicConnection(1));

            dataFile.InitiateDynamicConnection();

            // An in-range scan comes back through the dynamic connection.
            var scan = dataFile.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(scan, Is.Not.Null);
            Assert.That(scan.OneBasedScanNumber, Is.EqualTo(1));

            // An out-of-range scan number returns null rather than throwing.
            Assert.That(dataFile.GetOneBasedScanFromDynamicConnection(scanCount + 1), Is.Null);

            dataFile.CloseDynamicConnection();

            // After close the connection is disposed-but-non-null; the guard must catch it and throw
            // an MzLibException rather than crashing inside the vendor layer on RunHeaderEx.
            Assert.Throws<MzLibException>(() => dataFile.GetOneBasedScanFromDynamicConnection(1));

            // A redundant close is harmless.
            Assert.DoesNotThrow(() => dataFile.CloseDynamicConnection());
        }

        /// <summary>
        /// Concurrent reads over a shared dynamic connection, without a prior static load, so the reads
        /// actually go through the dynamic path guarded by DynamicReadingLock. Each returned scan must
        /// carry the scan number it was requested with - a race on the stateful, non-thread-safe
        /// IRawDataPlus handle would hand back a scan whose header came from a different spectrum.
        /// </summary>
        [Test]
        [TestCase("DataFiles/small.RAW")]
        [TestCase("DataFiles/sliced_ethcd.raw")]
        public static void TestDynamicConnection_ConcurrentReads_NoStaticLoad(string filePath)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            var dataFile = MsDataFileReader.GetDataFile(spectraPath);
            dataFile.InitiateDynamicConnection();

            // Derive the scan range from the connection itself, not from a static load.
            int scanCount = dataFile.GetMsOrderByScanInDynamicConnection().Length;
            var scanNumbers = Enumerable.Range(1, scanCount).ToList();

            var exceptions = new ConcurrentBag<Exception>();
            var results = new ConcurrentDictionary<int, MsDataScan>();

            Parallel.ForEach(scanNumbers, new ParallelOptions { MaxDegreeOfParallelism = 8 }, scanNumber =>
            {
                try
                {
                    results[scanNumber] = dataFile.GetOneBasedScanFromDynamicConnection(scanNumber);
                }
                catch (Exception ex)
                {
                    exceptions.Add(ex);
                }
            });

            dataFile.CloseDynamicConnection();

            Assert.That(exceptions, Is.Empty, "Exceptions occurred during concurrent dynamic reads.");
            Assert.That(results.Count, Is.EqualTo(scanNumbers.Count), "Not all scans were read successfully.");
            Assert.That(results.All(kvp => kvp.Value != null && kvp.Value.OneBasedScanNumber == kvp.Key), Is.True,
                "A concurrent read returned a scan whose number does not match the request (connection race).");
        }
    }
}
