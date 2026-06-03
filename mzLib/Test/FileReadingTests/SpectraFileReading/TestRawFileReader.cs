using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestRawFileReader
    {
        [Test]
        public void TestFileDoesntExist()
        {
            string fakePath = "fakePath.raw";
            var reader = MsDataFileReader.GetDataFile(fakePath);
            NUnit.Framework.Assert.Throws<FileNotFoundException>(() =>
            {
                reader.InitiateDynamicConnection();
            });
        }

        #region Testing Exceptions

        [Test]
        public void TestRawFileReaderFileNotFoundException()
        {
            var fakeRawFile = "asdasd.raw";

            var ex = NUnit.Framework.Assert.Throws<FileNotFoundException>(() => MsDataFileReader.GetDataFile(fakeRawFile).LoadAllStaticData());

            NUnit.Framework.Assert.That(ex.Message, Is.EqualTo(new FileNotFoundException().Message));
        }

        #endregion


        [Test]
        public void TestScanDescription()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "ScanDescriptionTestData.raw");
            var scans = MsDataFileReader.GetDataFile(filePath).GetAllScansList();
            var ms1Scans = scans.Where(x => x.MsnOrder == 1).ToList();
            var ms2Scans = scans.Where(x => x.MsnOrder == 2).ToList();

            ms1Scans.ForEach(x => NUnit.Framework.Assert.That(x.ScanDescription, Is.EqualTo(null)));
            ms2Scans.ForEach(x => NUnit.Framework.Assert.That(x.ScanDescription, Is.EqualTo("Testing2")));
        }

        /// <summary>
        /// Tests LoadAllStaticData for ThermoRawFileReader
        /// </summary>
        /// <param name="infile"></param>
        /// <param name="outfile1"></param>
        /// <param name="outfile2"></param>
        [Test]
        [TestCase("testFileWMS2.raw", "a.mzML", "aa.mzML")]
        [TestCase("small.raw", "a.mzML", "aa.mzML")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw", "a.mzML", "aa.mzML")]
        public static void TestLoadAllStaticDataRawFileReader(string infile, string outfile1, string outfile2)
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", infile);
            outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", outfile1);
            outfile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", outfile2);

            string dummyPath = "aaljienxmbelsiemxmbmeba.raw";
            NUnit.Framework.Assert.Throws<FileNotFoundException>(() =>
            {
                var dummyReader = MsDataFileReader.GetDataFile(dummyPath);
                dummyReader.LoadAllStaticData();
            });

            // testing load with multiple threads 
            var parallelReader = MsDataFileReader.GetDataFile(path);
            parallelReader.LoadAllStaticData(null, maxThreads: 4);

            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData(null, maxThreads: 1);
            reader.ExportAsMzML(outfile1, false);
            reader.LoadAllStaticData();
            reader.ExportAsMzML(outfile2, true);
            var readerMzml = MsDataFileReader.GetDataFile(outfile2);
            readerMzml.LoadAllStaticData();
            Console.WriteLine($"Analysis time for TestLoadAllStaticDataRawFileReader({infile}): {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public void TestThermoGetSourceFile()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var reader = MsDataFileReader.GetDataFile(path);
            SourceFile sf = reader.GetSourceFile();
            NUnit.Framework.Assert.That(sf.NativeIdFormat, Is.EqualTo(@"Thermo nativeID format"));
        }

        /// <summary>
        /// Tests the dynamic connection for thermorawfilereader
        /// </summary>
        [Test]
        public static void TestDynamicConnectionRawFileReader()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var path1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var thermoDynamic1 = MsDataFileReader.GetDataFile(path1);
            thermoDynamic1.InitiateDynamicConnection();

            var path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "testFileWMS2.raw");
            var thermoDynamic2 = MsDataFileReader.GetDataFile(path2);
            thermoDynamic2.InitiateDynamicConnection();

            var msOrders = thermoDynamic1.GetMsOrderByScanInDynamicConnection();
            NUnit.Framework.Assert.That(msOrders != null && msOrders.Length > 0);

            var a = thermoDynamic1.GetOneBasedScanFromDynamicConnection(1);
            NUnit.Framework.Assert.That(a != null);

            var b = thermoDynamic2.GetOneBasedScanFromDynamicConnection(1);
            NUnit.Framework.Assert.That(b != null);

            NUnit.Framework.Assert.That(a.MassSpectrum.XArray.Length != b.MassSpectrum.XArray.Length);

            a = thermoDynamic1.GetOneBasedScanFromDynamicConnection(10000);
            thermoDynamic1.CloseDynamicConnection();
            thermoDynamic2.CloseDynamicConnection();

            Console.WriteLine($"Analysis time for TestDynamicConnectionRawFileReader: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
            NUnit.Framework.Assert.That(a == null);
        }

        [Test]
        public static void TestDynamicConnectionRawFileReader_AfterStaticLoading()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var path1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var thermoDynamic1 = MsDataFileReader.GetDataFile(path1).LoadAllStaticData();
            thermoDynamic1.InitiateDynamicConnection();

            var path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "testFileWMS2.raw");
            var thermoDynamic2 = MsDataFileReader.GetDataFile(path2).LoadAllStaticData();
            thermoDynamic2.InitiateDynamicConnection();

            var msOrders = thermoDynamic1.GetMsOrderByScanInDynamicConnection();
            NUnit.Framework.Assert.That(msOrders != null && msOrders.Length > 0);

            var a = thermoDynamic1.GetOneBasedScanFromDynamicConnection(1);
            NUnit.Framework.Assert.That(a != null);

            var b = thermoDynamic2.GetOneBasedScanFromDynamicConnection(1);
            NUnit.Framework.Assert.That(b != null);

            NUnit.Framework.Assert.That(a.MassSpectrum.XArray.Length != b.MassSpectrum.XArray.Length);

            a = thermoDynamic1.GetOneBasedScanFromDynamicConnection(10000);
            thermoDynamic1.CloseDynamicConnection();
            thermoDynamic2.CloseDynamicConnection();

            Console.WriteLine($"Analysis time for TestDynamicConnectionRawFileReader: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
            NUnit.Framework.Assert.That(a == null);
        }

        /// <summary>
        /// Tests peak filtering for ThermoRawFileReader
        /// </summary>
        /// <param name="infile"></param>
        [Test]
        [TestCase("testFileWMS2.raw")]
        [TestCase("small.raw")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw")]
        public static void TestPeakFilteringRawFileReader(string infile)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            var filterParams = new FilteringParams(200, 0.01, 0, 1, false, true, true);

            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", infile);
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData(filterParams, maxThreads: 1);
            var rawScans = reader.GetAllScansList();
            foreach (var scan in rawScans)
            {
                NUnit.Framework.Assert.That(scan.MassSpectrum.XArray.Length <= 200);
            }

            string outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", Path.GetFileNameWithoutExtension(infile) + ".mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, outfile1, false);
            var mzml = MsDataFileReader.GetDataFile(outfile1);
            mzml.LoadAllStaticData(filterParams, maxThreads: 1);

            var mzmlScans = mzml.GetAllScansList();
            for (int i = 0; i < mzmlScans.Count; i++)
            {
                var mzmlScan = mzmlScans[i];
                var rawScan = rawScans[i];

                for (int j = 0; j < mzmlScan.MassSpectrum.XArray.Length; j++)
                {
                    double roundedRawMz = Math.Round(rawScan.MassSpectrum.XArray[j], 4);
                    double roundedMzmlMz = Math.Round(mzmlScan.MassSpectrum.XArray[j], 4);

                    // XArray is rounded to the 4th digit during CreateAndWrite
                    Assert.AreEqual(roundedMzmlMz, roundedRawMz);

                    double roundedMzmlIntensity = Math.Round(mzmlScan.MassSpectrum.XArray[j], 0);
                    double roundedRawIntensity = Math.Round(rawScan.MassSpectrum.XArray[j], 0);

                    Assert.AreEqual(roundedMzmlIntensity, roundedRawIntensity);
                }
            }

            Console.WriteLine($"Analysis time for TestPeakFilteringRawFileReader: {stopwatch.Elapsed.Hours}h " +
                $"{stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        /// <summary>
        /// Test Thermo License for ThermoRawFileReader
        /// </summary>
        [Test]
        public static void TestThermoLicence()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var licence = ThermoRawFileReaderLicence.ThermoLicenceText;
            NUnit.Framework.Assert.That(licence.Length > 100);

            Console.WriteLine($"Analysis time for TestThermoLicence: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        /// <summary>
        /// Test that raw files can be opened dynamically in ThermoRawFileReader
        /// </summary>
        /// <param name="fileName"></param>
        [Test]
        [TestCase("small.RAW")]
        [TestCase("testFileWMS2.raw")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw")]
        public static void TestDynamicRaw(string fileName)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

            var staticRaw = MsDataFileReader.GetDataFile(filePath);
            staticRaw.LoadAllStaticData();
            staticRaw.InitiateDynamicConnection();

            foreach (MsDataScan staticScan in staticRaw.GetAllScansList())
            {
                MsDataScan dynamicScan = staticRaw.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                Assert.IsFalse(staticScan.MassSpectrum.YArray.Contains(0));
                Assert.IsFalse(dynamicScan.MassSpectrum.YArray.Contains(0));
                NUnit.Framework.Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                NUnit.Framework.Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);
                NUnit.Framework.Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                NUnit.Framework.Assert.That(dynamicScan.Polarity == staticScan.Polarity);
                NUnit.Framework.Assert.That(dynamicScan.ScanWindowRange.Minimum == staticScan.ScanWindowRange.Minimum);
                NUnit.Framework.Assert.That(dynamicScan.ScanWindowRange.Maximum == staticScan.ScanWindowRange.Maximum);
                NUnit.Framework.Assert.That(dynamicScan.ScanFilter == staticScan.ScanFilter);
                NUnit.Framework.Assert.That(dynamicScan.NativeId == staticScan.NativeId);
                NUnit.Framework.Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                NUnit.Framework.Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                NUnit.Framework.Assert.That(dynamicScan.InjectionTime == staticScan.InjectionTime);
                NUnit.Framework.Assert.That(dynamicScan.NoiseData == staticScan.NoiseData);

                NUnit.Framework.Assert.That(dynamicScan.IsolationMz == staticScan.IsolationMz);
                NUnit.Framework.Assert.That(dynamicScan.SelectedIonChargeStateGuess == staticScan.SelectedIonChargeStateGuess);
                NUnit.Framework.Assert.That(dynamicScan.SelectedIonIntensity == staticScan.SelectedIonIntensity);
                NUnit.Framework.Assert.That(dynamicScan.SelectedIonMZ == staticScan.SelectedIonMZ);
                NUnit.Framework.Assert.That(dynamicScan.DissociationType == staticScan.DissociationType);
                NUnit.Framework.Assert.That(dynamicScan.IsolationWidth == staticScan.IsolationWidth);
                NUnit.Framework.Assert.That(dynamicScan.OneBasedPrecursorScanNumber == staticScan.OneBasedPrecursorScanNumber);
                NUnit.Framework.Assert.That(dynamicScan.SelectedIonMonoisotopicGuessIntensity == staticScan.SelectedIonMonoisotopicGuessIntensity);
                NUnit.Framework.Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz == staticScan.SelectedIonMonoisotopicGuessMz);

                if (dynamicScan.IsolationRange != null || staticScan.IsolationRange != null)
                {
                    NUnit.Framework.Assert.That(dynamicScan.IsolationRange.Minimum == staticScan.IsolationRange.Minimum);
                    NUnit.Framework.Assert.That(dynamicScan.IsolationRange.Maximum == staticScan.IsolationRange.Maximum);
                }

                NUnit.Framework.Assert.That(dynamicScan.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                NUnit.Framework.Assert.That(dynamicScan.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);

                for (int i = 0; i < staticScan.MassSpectrum.XArray.Length; i++)
                {
                    double staticMz = staticScan.MassSpectrum.XArray[i];
                    double staticIntensity = staticScan.MassSpectrum.YArray[i];

                    double dynamicMz = dynamicScan.MassSpectrum.XArray[i];
                    double dynamicIntensity = dynamicScan.MassSpectrum.YArray[i];

                    NUnit.Framework.Assert.That(dynamicMz == staticMz);
                    NUnit.Framework.Assert.That(dynamicIntensity == staticIntensity);
                }
            }
        }

        /// <summary>
        /// Tests that you can read EtHCD files in ThermoRawFileReader
        /// </summary>
        [Test]
        public static void TestEthcdReading()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "sliced_ethcd.raw");
            var spectra = MsDataFileReader.GetDataFile(filePath);
            spectra.LoadAllStaticData(null, 1);
            var hcdScan = spectra.GetOneBasedScan(5);
            NUnit.Framework.Assert.That(hcdScan.DissociationType == DissociationType.HCD);
            NUnit.Framework.Assert.That(hcdScan.HcdEnergy, Is.EqualTo("36.00"));
            var ethcdScan = spectra.GetOneBasedScan(6);
            NUnit.Framework.Assert.That(ethcdScan.DissociationType == DissociationType.EThcD);
            NUnit.Framework.Assert.That(ethcdScan.HcdEnergy, Is.EqualTo("25.00"));
        }

        [Test]
        public static void TestCompensationVoltageReading()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "TestCompensationVoltageReading.raw");
            var spectra = MsDataFileReader.GetDataFile(filePath);
            spectra.LoadAllStaticData();
            var availableCvValues = spectra
                .GetAllScansList()
                .Select(i => i.CompensationVoltage)
                .Distinct()
                .OrderByDescending(i => i)
                .ToArray();
            double?[] expected =  new double?[] {-45d, -60d}; 

        Assert.AreEqual(expected, availableCvValues); 

        }
    }
}