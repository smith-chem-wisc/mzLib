using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestRawFileReader
    {
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

            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            var a = ThermoRawFileReader.LoadAllStaticData(path, maxThreads: 1);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, outfile1, false);
            var aa = Mzml.LoadAllStaticData(outfile1);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(aa, outfile2, true);
            Mzml.LoadAllStaticData(outfile2);
            Console.WriteLine($"Analysis time for TestLoadAllStaticDataRawFileReader({infile}): {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
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
            var dynamicConnection1 = new ThermoDynamicData(path1);

            var path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "testFileWMS2.raw");
            var dynamicConnection2 = new ThermoDynamicData(path2);

            var msOrders = dynamicConnection1.MsOrdersByScan;
            Assert.That(msOrders != null && msOrders.Length > 0);

            var a = dynamicConnection1.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(a != null);

            var b = dynamicConnection2.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(b != null);

            Assert.That(a.MassSpectrum.XArray.Length != b.MassSpectrum.XArray.Length);

            a = dynamicConnection1.GetOneBasedScanFromDynamicConnection(10000);
            Assert.That(a == null);

            dynamicConnection1.CloseDynamicConnection();
            dynamicConnection2.CloseDynamicConnection();

            Console.WriteLine($"Analysis time for TestDynamicConnectionRawFileReader: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
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

            var a = ThermoRawFileReader.LoadAllStaticData(path, filterParams, maxThreads: 1);
            var rawScans = a.GetAllScansList();
            foreach (var scan in rawScans)
            {
                Assert.That(scan.MassSpectrum.XArray.Length <= 200);
            }

            string outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", Path.GetFileNameWithoutExtension(infile) + ".mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, outfile1, false);
            var mzml = Mzml.LoadAllStaticData(outfile1, filterParams, maxThreads: 1);

            var mzmlScans = mzml.GetAllScansList();
            for (int i = 0; i < mzmlScans.Count; i++)
            {
                var mzmlScan = mzmlScans[i];
                var rawScan = rawScans[i];

                for (int j = 0; j < mzmlScan.MassSpectrum.XArray.Length; j++)
                {
                    double roundedMzmlMz = Math.Round(mzmlScan.MassSpectrum.XArray[j], 2);
                    double roundedRawMz = Math.Round(rawScan.MassSpectrum.XArray[j], 2);

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
            Assert.That(licence.Length > 100);

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

            ThermoRawFileReader staticRaw = ThermoRawFileReader.LoadAllStaticData(filePath);
            ThermoDynamicData dynamicRaw = new ThermoDynamicData(filePath);

            foreach (MsDataScan staticScan in staticRaw.GetAllScansList())
            {
                MsDataScan dynamicScan = dynamicRaw.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                Assert.IsFalse(staticScan.MassSpectrum.YArray.Contains(0));
                Assert.IsFalse(dynamicScan.MassSpectrum.YArray.Contains(0));
                Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);
                Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                Assert.That(dynamicScan.Polarity == staticScan.Polarity);
                Assert.That(dynamicScan.ScanWindowRange.Minimum == staticScan.ScanWindowRange.Minimum);
                Assert.That(dynamicScan.ScanWindowRange.Maximum == staticScan.ScanWindowRange.Maximum);
                Assert.That(dynamicScan.ScanFilter == staticScan.ScanFilter);
                Assert.That(dynamicScan.NativeId == staticScan.NativeId);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.InjectionTime == staticScan.InjectionTime);
                Assert.That(dynamicScan.NoiseData == staticScan.NoiseData);

                Assert.That(dynamicScan.IsolationMz == staticScan.IsolationMz);
                Assert.That(dynamicScan.SelectedIonChargeStateGuess == staticScan.SelectedIonChargeStateGuess);
                Assert.That(dynamicScan.SelectedIonIntensity == staticScan.SelectedIonIntensity);
                Assert.That(dynamicScan.SelectedIonMZ == staticScan.SelectedIonMZ);
                Assert.That(dynamicScan.DissociationType == staticScan.DissociationType);
                Assert.That(dynamicScan.IsolationWidth == staticScan.IsolationWidth);
                Assert.That(dynamicScan.OneBasedPrecursorScanNumber == staticScan.OneBasedPrecursorScanNumber);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessIntensity == staticScan.SelectedIonMonoisotopicGuessIntensity);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz == staticScan.SelectedIonMonoisotopicGuessMz);

                if (dynamicScan.IsolationRange != null || staticScan.IsolationRange != null)
                {
                    Assert.That(dynamicScan.IsolationRange.Minimum == staticScan.IsolationRange.Minimum);
                    Assert.That(dynamicScan.IsolationRange.Maximum == staticScan.IsolationRange.Maximum);
                }

                Assert.That(dynamicScan.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                Assert.That(dynamicScan.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);

                for (int i = 0; i < staticScan.MassSpectrum.XArray.Length; i++)
                {
                    double staticMz = staticScan.MassSpectrum.XArray[i];
                    double staticIntensity = staticScan.MassSpectrum.YArray[i];

                    double dynamicMz = dynamicScan.MassSpectrum.XArray[i];
                    double dynamicIntensity = dynamicScan.MassSpectrum.YArray[i];

                    Assert.That(dynamicMz == staticMz);
                    Assert.That(dynamicIntensity == staticIntensity);
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
            var spectra = ThermoRawFileReader.LoadAllStaticData(filePath, null, 1);
            var hcdScan = spectra.GetOneBasedScan(5);
            Assert.That(hcdScan.DissociationType == DissociationType.HCD);
            var ethcdScan = spectra.GetOneBasedScan(6);
            Assert.That(ethcdScan.DissociationType == DissociationType.EThcD);
        }
    }
}