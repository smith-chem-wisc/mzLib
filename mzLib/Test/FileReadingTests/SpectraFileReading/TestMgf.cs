using System;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Readers;
using MzLibUtil;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.FileReadingTests.SpectraFileReading
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
        public void SkipsEmptySpectra()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "withEmptySpectra.mgf");
            var file = MsDataFileReader.GetDataFile(path);

            file.LoadAllStaticData();

            var scans = file.GetAllScansList();

            // Skipped two empty spectra, so only one scan should be loaded
            NUnit.Framework.Assert.That(scans.Count, Is.EqualTo(1));
            NUnit.Framework.Assert.That(scans[0].OneBasedScanNumber, Is.EqualTo(25501));
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
            NUnit.Framework.Assert.That(ya.RetentionTime, Is.EqualTo(15.393).Within(0.1));
            NUnit.Framework.Assert.That(ya.IsolationMz, Is.EqualTo(354.8).Within(0.1));
            NUnit.Framework.Assert.That(ya.SelectedIonMZ, Is.EqualTo(354.8).Within(0.1));
            NUnit.Framework.Assert.That(ya.SelectedIonChargeStateGuess, Is.EqualTo(2));
            NUnit.Framework.Assert.That(ya.SelectedIonMonoisotopicGuessMz, Is.EqualTo(354.8).Within(0.1));
            NUnit.Framework.Assert.That(ya.TotalIonCurrent, Is.EqualTo(1737).Within(0.1));
            NUnit.Framework.Assert.That(ya.ScanWindowRange.Minimum, Is.EqualTo(227.787).Within(0.1));
            NUnit.Framework.Assert.That(ya.ScanWindowRange.Maximum, Is.EqualTo(565.64).Within(0.1));
        }

        [Test]
        public void TestMgfDynamicConnection()
        {
            var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester.mgf"));
            var fileDoesntExistReader = MsDataFileReader.GetDataFile("fakeFile.mgf");

            NUnit.Framework.Assert.Throws<FileNotFoundException>(() =>
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
        public void TestMgfDynamicConnection_AfterStaticLoading()
        {
            var filter = new FilteringParams(1, 0.01);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester.mgf");
            var readerReal = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filter);
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

            // A separate reader, with no static load: reusing the one above would return the cached
            // scans and check the static path twice rather than the dynamic parse.
            var dynamicReader = MsDataFileReader.GetDataFile(path);
            dynamicReader.InitiateDynamicConnection();
            MsDataScan dynamicScan1 = dynamicReader.GetOneBasedScanFromDynamicConnection(1);
            MsDataScan dynamicScan2 = dynamicReader.GetOneBasedScanFromDynamicConnection(2);
            Assert.IsFalse(dynamicScan1.MassSpectrum.YArray.Contains(0));
            Assert.IsFalse(dynamicScan2.MassSpectrum.YArray.Contains(0));
            dynamicReader.CloseDynamicConnection();
        }

        [Test]
        public void ReadsPrecursorIntensityWhenPresent()
        {
            //read the mgf file. zero intensity peaks should be eliminated during read
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "withIntensity.mgf");

            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();

            var scans = reader.GetAllScansList();
            NUnit.Framework.Assert.That(scans.Count, Is.EqualTo(1));
            NUnit.Framework.Assert.That(scans[0].SelectedIonIntensity, Is.Not.Null);
            NUnit.Framework.Assert.That(scans[0].SelectedIonIntensity, Is.EqualTo(47641904.0).Within(0.00001));
        }

        /// <summary>
        /// Pins the sign parsed off the CHARGE line. The polarity assertions are not independent
        /// evidence -- MsDataScan's constructor normalizes the charge sign to the polarity, and Mgf
        /// derives that polarity from the same charge -- so the scan count and retention times are
        /// asserted as well, which the normalization cannot manufacture.
        /// </summary>
        [Test]
        public void NegativeModeSetsCorrectCharge_FromMgfChargeLine()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "negativeModeCharge.mgf");

            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();

            var scans = reader.GetAllScansList();
            NUnit.Framework.Assert.That(scans.Count, Is.EqualTo(2), path + " did not parse as two scans.");
            NUnit.Framework.Assert.That(scans.Select(s => s.RetentionTime).Distinct().Count(), Is.EqualTo(2),
                "Both scans carry the same retention time, so the two BEGIN IONS blocks did not parse separately.");

            var negative = reader.GetOneBasedScan(1);
            NUnit.Framework.Assert.That(negative.SelectedIonChargeStateGuess, Is.EqualTo(-2),
                "CHARGE=2- did not parse as -2.");
            NUnit.Framework.Assert.That(negative.Polarity, Is.EqualTo(Polarity.Negative),
                "Scan polarity is not negative.");

            var positive = reader.GetOneBasedScan(2);
            NUnit.Framework.Assert.That(positive.SelectedIonChargeStateGuess, Is.EqualTo(3),
                "CHARGE=3+ did not parse as 3.");
            NUnit.Framework.Assert.That(positive.Polarity, Is.EqualTo(Polarity.Positive),
                "Scan polarity is not positive.");
        }

        /// <summary>
        /// The dangerous half of the optional sign suffix, kept in its own fixture so the assertion
        /// that fails is the VALUE. Stripping the last character unconditionally dropped a digit, so
        /// CHARGE=12 read as 1 -- and 1 is a valid charge, so nothing surfaced.
        /// </summary>
        [Test]
        public void UnsignedMultiDigitChargeKeepsEveryDigit()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "unsignedMultiDigitCharge.mgf");

            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();

            NUnit.Framework.Assert.That(reader.GetOneBasedScan(1).SelectedIonChargeStateGuess, Is.EqualTo(12),
                "CHARGE=12 did not parse as 12; an unconditional suffix strip reads it as 1.");
            NUnit.Framework.Assert.That(reader.GetOneBasedScan(2).SelectedIonChargeStateGuess, Is.EqualTo(13),
                "CHARGE=13 did not parse as 13; an unconditional suffix strip reads it as 1.");
            NUnit.Framework.Assert.That(reader.GetAllScansList().Select(s => s.Polarity),
                Is.All.EqualTo(Polarity.Positive), "An unsigned charge is positive.");
        }

        /// <summary>
        /// The loud half: a single-digit unsigned charge left nothing to parse, so the whole file
        /// failed to load. Separate fixture, so this failing cannot mask the silent case above.
        /// </summary>
        [Test]
        public void UnsignedSingleDigitChargeLoads()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "unsignedSingleDigitCharge.mgf");

            var reader = MsDataFileReader.GetDataFile(path);
            NUnit.Framework.Assert.DoesNotThrow(() => reader.LoadAllStaticData(),
                "CHARGE=2 left an empty string to parse, so the file did not load at all.");

            NUnit.Framework.Assert.That(reader.GetOneBasedScan(1).SelectedIonChargeStateGuess, Is.EqualTo(2),
                "CHARGE=2 did not parse as 2.");
            NUnit.Framework.Assert.That(reader.GetOneBasedScan(1).Polarity, Is.EqualTo(Polarity.Positive),
                "An unsigned charge is positive.");
        }


        /// <summary>
        /// Regression for MetaMorpheus issue #2366: a malformed line threw IndexOutOfRangeException out
        /// of the reader, so one bad line made the whole file unreadable rather than being skipped the
        /// way tester_corrupt.mgf's unknown words already are.
        ///
        /// malformedLines.mgf carries one such line per scan: a peak with an m/z but no intensity, a
        /// header with no "=", a header with "=" and nothing after it, and a line that starts with a
        /// digit but is not a peak. Each scan's well-formed peaks must survive.
        /// </summary>
        [Test]
        public static void TestLoadMgfWithMalformedLines()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "malformedLines.mgf");
            var reader = MsDataFileReader.GetDataFile(path);

            NUnit.Framework.Assert.DoesNotThrow(() => reader.LoadAllStaticData(),
                "one malformed line must not make the file unreadable");

            NUnit.Framework.Assert.That(reader.NumSpectra, Is.EqualTo(4));

            // scan 1: "202.456" has no intensity and is dropped; the two complete peaks remain
            var peakWithoutIntensity = reader.GetOneBasedScan(1);
            NUnit.Framework.Assert.That(peakWithoutIntensity.MassSpectrum.Size, Is.EqualTo(2));
            NUnit.Framework.Assert.That(peakWithoutIntensity.MassSpectrum.XArray,
                Is.EqualTo(new[] { 201.123, 203.789 }).Within(1e-9));

            // scan 2: bare "CHARGE" carries no value, so the default charge stands
            var chargeWithoutEquals = reader.GetOneBasedScan(2);
            NUnit.Framework.Assert.That(chargeWithoutEquals.MassSpectrum.Size, Is.EqualTo(2));
            NUnit.Framework.Assert.That(chargeWithoutEquals.SelectedIonChargeStateGuess, Is.EqualTo(2));

            // scan 3: "CHARGE=" is empty, and "9notanumber here" starts with a digit but is not a peak
            var emptyChargeAndNonNumericPeak = reader.GetOneBasedScan(3);
            NUnit.Framework.Assert.That(emptyChargeAndNonNumericPeak.MassSpectrum.Size, Is.EqualTo(2));
            NUnit.Framework.Assert.That(emptyChargeAndNonNumericPeak.MassSpectrum.XArray,
                Is.EqualTo(new[] { 401.5, 402.5 }).Within(1e-9));

            // scan 4: valueless RTINSECONDS and SCANS leave their defaults rather than throwing
            var headersWithoutValues = reader.GetOneBasedScan(4);
            NUnit.Framework.Assert.That(headersWithoutValues.MassSpectrum.Size, Is.EqualTo(2));
            NUnit.Framework.Assert.That(headersWithoutValues.SelectedIonChargeStateGuess, Is.EqualTo(3));
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
        [TestCase("negativeModeCharge.mgf")]
        [TestCase("unsignedMultiDigitCharge.mgf")]
        public static void TestDynamicMgf(string fileName)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

            // Two readers, deliberately. GetOneBasedScanFromDynamicConnection short-circuits to the
            // static cache when scans are already loaded, so a single reader that has called
            // LoadAllStaticData returns the very same MsDataScan instance and every comparison below
            // is a scan against itself. Keeping the dynamic reader free of a static load is what makes
            // this exercise the dynamic parse.
            var staticReader = MsDataFileReader.GetDataFile(filePath);
            staticReader.LoadAllStaticData();

            var dynamicReader = MsDataFileReader.GetDataFile(filePath);
            dynamicReader.InitiateDynamicConnection();

            foreach (MsDataScan staticScan in staticReader.GetAllScansList())
            {
                MsDataScan dynamicScan = dynamicReader.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                NUnit.Framework.Assert.That(dynamicScan, Is.Not.Null,
                    $"No scan {staticScan.OneBasedScanNumber} came back from the dynamic connection.");
                NUnit.Framework.Assert.That(ReferenceEquals(dynamicScan, staticScan), Is.False,
                    "The dynamic read returned the statically cached instance, so nothing below is being compared.");

                NUnit.Framework.Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                NUnit.Framework.Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);

                if (!double.IsNaN(dynamicScan.RetentionTime) || !double.IsNaN(staticScan.RetentionTime))
                {
                    NUnit.Framework.Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                }

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

            dynamicReader.CloseDynamicConnection();
        }

        [Test]
        public void TestGetByteOffsetAtCurrentPositionReaderNullBranch()
        {
            // create a stream reader that will generate a null
            StreamReader streamReader = null;
            NUnit.Framework.Assert.Throws<MzLibException>(() => TextFileReading.GetByteOffsetAtCurrentPosition(streamReader));
        }
    }
}