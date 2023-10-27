using System;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using MzLibUtil;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMgf
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestLoadMgf()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "ThereIsNothingHerePleaseDoNotGenerateThisFile.mgf");
            var reader = MsDataFileReader.GetDataFile(path);

            try
            {
                reader.LoadAllStaticData();
                Assert.IsTrue(false);
            }
            catch
            {
                //woohoo, there was an exception!
            }
            reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester.mgf"));
            reader.LoadAllStaticData();
            var ya = reader.GetOneBasedScan(14);

            Assert.AreEqual(192, ya.MassSpectrum.Size);
            Assert.AreEqual(2, ya.MsnOrder);
            Assert.AreEqual(14, ya.OneBasedScanNumber);
            Assert.AreEqual(Polarity.Positive, ya.Polarity);
            Assert.AreEqual(0.26666666666666666, ya.RetentionTime);
            Assert.AreEqual(571.806916, ya.IsolationMz);
            Assert.AreEqual(571.806916, ya.SelectedIonMZ);
            Assert.AreEqual(2, ya.SelectedIonChargeStateGuess);
            Assert.AreEqual(571.806916, ya.SelectedIonMonoisotopicGuessMz);
            Assert.AreEqual(1294963.5999999996, ya.TotalIonCurrent);
            Assert.AreEqual(110.0719, ya.ScanWindowRange.Minimum);
            Assert.AreEqual(1038.8018, ya.ScanWindowRange.Maximum);
            var ya2 = reader.GetOneBasedScan(20).MassSpectrum;
            Assert.AreEqual(165, ya2.Size);
            var ya3 = reader.GetOneBasedScan(2).MassSpectrum;
            Assert.AreEqual(551, ya3.Size);
        }

        [Test]
        public static void TestLoadMgfTabSeparated()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Tab_separated_peak_list.mgf");
            var dataReader = MsDataFileReader.GetDataFile(path);
            dataReader.LoadAllStaticData();

            var ya = dataReader.GetOneBasedScan(2);

            Assert.AreEqual(19, ya.MassSpectrum.Size);
            Assert.AreEqual(2, ya.MsnOrder);
            Assert.AreEqual(2, ya.OneBasedScanNumber);
            Assert.AreEqual(Polarity.Positive, ya.Polarity);
            Assert.That(ya.RetentionTime, Is.EqualTo(15.393).Within(0.1));
            Assert.That(ya.IsolationMz, Is.EqualTo(354.8).Within(0.1));
            Assert.That(ya.SelectedIonMZ, Is.EqualTo(354.8).Within(0.1));
            Assert.That(ya.SelectedIonChargeStateGuess, Is.EqualTo(2));
            Assert.That(ya.SelectedIonMonoisotopicGuessMz, Is.EqualTo(354.8).Within(0.1));
            Assert.That(ya.TotalIonCurrent, Is.EqualTo(1737).Within(0.1));
            Assert.That(ya.ScanWindowRange.Minimum, Is.EqualTo(227.787).Within(0.1));
            Assert.That(ya.ScanWindowRange.Maximum, Is.EqualTo(565.64).Within(0.1));
        }

        [Test]
        public void TestMgfDynamicConnection()
        {
            var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester.mgf"));
            var fileDoesntExistReader = MsDataFileReader.GetDataFile("fakeFile.mgf");

            Assert.Throws<FileNotFoundException>(() =>
            {
                fileDoesntExistReader.InitiateDynamicConnection();
            });
            IFilteringParams filter = new FilteringParams(1, 0.01);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester.mgf");
            var readerReal = MsDataFileReader.GetDataFile(path);
            readerReal.InitiateDynamicConnection();
            readerReal.GetOneBasedScanFromDynamicConnection(2, filter);
        }

        [Test]
        public void EliminateZeroIntensityPeaksFromMgfOnFileLoad()
        {
            //read the mgf file. zero intensity peaks should be eliminated during read
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "withZeros.mgf");

            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            //insure that read scans contain no zero intensity peaks
            Assert.IsFalse(reader.GetAllScansList()[0].MassSpectrum.YArray.Contains(0));
            Assert.IsFalse(reader.GetAllScansList()[1].MassSpectrum.YArray.Contains(0));

            reader.InitiateDynamicConnection();
            MsDataScan dynamicScan1 = reader.GetOneBasedScanFromDynamicConnection(1);
            MsDataScan dynamicScan2 = reader.GetOneBasedScanFromDynamicConnection(2);
            Assert.IsFalse(dynamicScan1.MassSpectrum.YArray.Contains(0));
            Assert.IsFalse(dynamicScan2.MassSpectrum.YArray.Contains(0));
            reader.CloseDynamicConnection();
        }


        [Test]
        public static void TestLoadCorruptMgf()
        {
            //tester_corrupt.mgf is extracted from tester.mgf except it contains empty lines or unknown words. You can compare the two files and find the differences.
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester_corrupt.mgf");
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            var ya = reader.GetOneBasedScan(14);

            Assert.AreEqual(192, ya.MassSpectrum.Size);
            Assert.AreEqual(2, ya.MsnOrder);
            Assert.AreEqual(14, ya.OneBasedScanNumber);
            Assert.AreEqual(Polarity.Positive, ya.Polarity);
            Assert.AreEqual(0.26666666666666666, ya.RetentionTime);
            Assert.AreEqual(571.806916, ya.IsolationMz);
            Assert.AreEqual(571.806916, ya.SelectedIonMZ);
            Assert.AreEqual(2, ya.SelectedIonChargeStateGuess);
            Assert.AreEqual(571.806916, ya.SelectedIonMonoisotopicGuessMz);
            Assert.AreEqual(1294963.5999999996, ya.TotalIonCurrent);
            Assert.AreEqual(110.0719, ya.ScanWindowRange.Minimum);
            Assert.AreEqual(1038.8018, ya.ScanWindowRange.Maximum);
        }


        [Test]
        [TestCase("tester.mgf")]
        [TestCase("SmallCalibratibleYeast.mgf")]
        public static void TestDynamicMgf(string fileName)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData();
            reader.InitiateDynamicConnection();

            foreach (MsDataScan staticScan in reader.GetAllScansList())
            {
                MsDataScan dynamicScan = reader.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);

                if (!double.IsNaN(dynamicScan.RetentionTime) || !double.IsNaN(staticScan.RetentionTime))
                {
                    Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                }

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

        [Test]
        public void TestGetByteOffsetAtCurrentPositionReaderNullBranch()
        {
            // create a stream reader that will generate a null
            StreamReader streamReader = null;
            Assert.Throws<MzLibException>(() => TextFileReading.GetByteOffsetAtCurrentPosition(streamReader));
        }
    }
}