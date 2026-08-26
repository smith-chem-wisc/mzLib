using System;
using System.Globalization;
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
        /// Builds an MS1 + MS2 file in memory, writes it, and reads it back. MGF is a lossy container --
        /// analyzer type, dissociation type and scan filters have nowhere to go, and the reader drops peaks
        /// below 0.01 intensity and sorts by m/z -- so this asserts the subset the format can actually carry
        /// rather than object equality, which could only pass by being weakened until it proved nothing.
        /// </summary>
        [Test]
        public void MgfRoundTripPreservesScanContentThatTheFormatCanCarry()
        {
            var ms1 = new MsDataScan(
                new MzSpectrum(new[] { 300.1, 400.2, 500.3 }, new[] { 1000.0, 2000.0, 3000.0 }, false),
                oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.5, scanWindowRange: new MzRange(300, 501), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 6000.0, injectionTime: null,
                noiseData: null, nativeId: "scan=1");

            var ms2 = new MsDataScan(
                new MzSpectrum(new[] { 110.05, 220.11 }, new[] { 500.0, 750.0 }, false),
                oneBasedScanNumber: 2, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.6, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1250.0, injectionTime: null,
                noiseData: null, nativeId: "scan=2", selectedIonMz: 571.8069,
                selectedIonChargeStateGuess: 2, selectedIonIntensity: 999999.0);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfRoundTrip.mgf");
            new GenericMsDataFile(new[] { ms1, ms2 }, null).ExportAsMgf(path);

            try
            {
                var reread = MsDataFileReader.GetDataFile(path);
                reread.LoadAllStaticData();
                var scans = reread.GetAllScansList();

                NUnit.Framework.Assert.That(scans.Count, Is.EqualTo(2));

                var readMs1 = scans[0];
                NUnit.Framework.Assert.That(readMs1.MsnOrder, Is.EqualTo(1));
                NUnit.Framework.Assert.That(readMs1.OneBasedScanNumber, Is.EqualTo(1));
                NUnit.Framework.Assert.That(readMs1.SelectedIonMZ, Is.Null);
                NUnit.Framework.Assert.That(readMs1.RetentionTime, Is.EqualTo(1.5).Within(1e-9));
                NUnit.Framework.Assert.That(readMs1.MassSpectrum.XArray, Is.EqualTo(new[] { 300.1, 400.2, 500.3 }).Within(1e-9));
                NUnit.Framework.Assert.That(readMs1.MassSpectrum.YArray, Is.EqualTo(new[] { 1000.0, 2000.0, 3000.0 }).Within(1e-9));

                var readMs2 = scans[1];
                NUnit.Framework.Assert.That(readMs2.MsnOrder, Is.EqualTo(2));
                NUnit.Framework.Assert.That(readMs2.OneBasedScanNumber, Is.EqualTo(2));
                NUnit.Framework.Assert.That(readMs2.SelectedIonMZ, Is.EqualTo(571.8069).Within(1e-9));
                NUnit.Framework.Assert.That(readMs2.SelectedIonChargeStateGuess, Is.EqualTo(2));
                NUnit.Framework.Assert.That(readMs2.SelectedIonIntensity, Is.EqualTo(999999.0).Within(1e-6));
                NUnit.Framework.Assert.That(readMs2.RetentionTime, Is.EqualTo(1.6).Within(1e-9));
                NUnit.Framework.Assert.That(readMs2.MassSpectrum.XArray, Is.EqualTo(new[] { 110.05, 220.11 }).Within(1e-9));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The Matrix Science specification requires no whitespace around '=', all parameters ahead of the
        /// fragment peaks, and '.' as the decimal separator regardless of machine locale. The culture is
        /// switched to one that uses ',' for decimals so a missing InvariantCulture would show up here.
        /// </summary>
        [Test]
        public void MgfWriterEmitsSpecConformantText()
        {
            var scan = new MsDataScan(
                new MzSpectrum(new[] { 110.05, 220.11 }, new[] { 500.0, 750.0 }, false),
                oneBasedScanNumber: 7, msnOrder: 2, isCentroid: true, polarity: Polarity.Negative,
                retentionTime: 2.0, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1250.0, injectionTime: null,
                noiseData: null, nativeId: "scan=7", selectedIonMz: 350.5,
                selectedIonChargeStateGuess: 2, selectedIonIntensity: null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfSpecConformance.mgf");
            var originalCulture = System.Threading.Thread.CurrentThread.CurrentCulture;

            try
            {
                System.Threading.Thread.CurrentThread.CurrentCulture = new CultureInfo("de-DE");
                new GenericMsDataFile(new[] { scan }, null).ExportAsMgf(path);

                string[] lines = File.ReadAllLines(path);

                NUnit.Framework.Assert.That(lines[0], Is.EqualTo("BEGIN IONS"));
                NUnit.Framework.Assert.That(lines[^1], Is.EqualTo("END IONS"));

                // a negative-polarity scan must carry the '-' suffix, and decimals must not be localised
                NUnit.Framework.Assert.That(lines, Has.Member("CHARGE=2-"));
                NUnit.Framework.Assert.That(lines, Has.Member("MSLEVEL=2"));
                NUnit.Framework.Assert.That(lines, Has.Member("SCANS=7"));
                NUnit.Framework.Assert.That(lines, Has.Member("PEPMASS=350.5"));
                NUnit.Framework.Assert.That(lines, Has.Member("RTINSECONDS=120"));
                NUnit.Framework.Assert.That(lines, Has.Member("110.05 500"));

                // every parameter must precede the first peak line
                int lastParameter = Array.FindLastIndex(lines, l => l.Contains('='));
                int firstPeak = Array.FindIndex(lines, l => l.Length > 0 && char.IsDigit(l[0]));
                NUnit.Framework.Assert.That(lastParameter, Is.LessThan(firstPeak));
            }
            finally
            {
                System.Threading.Thread.CurrentThread.CurrentCulture = originalCulture;
                File.Delete(path);
            }
        }

        [Test]
        public void MgfWriterRejectsANullFile()
        {
            NUnit.Framework.Assert.Throws<ArgumentNullException>(
                () => MgfMethods.WriteMgf(null, Path.Combine(TestContext.CurrentContext.TestDirectory, "unused.mgf")));
        }

        /// <summary>
        /// A file straight from MsDataFileReader has not loaded its scans yet, so the writer has to load
        /// them itself. It also has a FilePath, which is where TITLE takes its name from in preference to
        /// the output path.
        /// </summary>
        [Test]
        public void MgfWriterLoadsAnUnloadedFileAndNamesTitleFromItsPath()
        {
            var source = MsDataFileReader.GetDataFile(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester.mgf"));
            NUnit.Framework.Assert.That(source.Scans, Is.Null, "precondition: the file must not be loaded yet");

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfFromUnloaded.mgf");

            try
            {
                source.ExportAsMgf(path);

                string[] lines = File.ReadAllLines(path);
                NUnit.Framework.Assert.That(lines.Count(l => l == "BEGIN IONS"), Is.EqualTo(source.GetAllScansList().Count));

                // named for the source file, not the destination
                NUnit.Framework.Assert.That(lines.First(l => l.StartsWith("TITLE=")), Does.Contain("tester"));
                NUnit.Framework.Assert.That(lines.First(l => l.StartsWith("TITLE=")), Does.Not.Contain("mgfFromUnloaded"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// Exercises the writer's three "nothing to say" paths: a null scan and an empty spectrum are both
        /// skipped rather than emitted as peakless blocks, an unknown retention time omits RTINSECONDS
        /// rather than writing NaN, and a null NativeId still produces a well-formed TITLE.
        /// </summary>
        [Test]
        public void MgfWriterSkipsEmptyScansAndOmitsUnknownRetentionTime()
        {
            var empty = new MsDataScan(
                new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false),
                oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.0, scanWindowRange: null, scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 0.0, injectionTime: null,
                noiseData: null, nativeId: "scan=1");

            var noRetentionTime = new MsDataScan(
                new MzSpectrum(new[] { 150.0 }, new[] { 42.0 }, false),
                oneBasedScanNumber: 2, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: double.NaN, scanWindowRange: null, scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 42.0, injectionTime: null,
                noiseData: null, nativeId: null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfSparse.mgf");

            try
            {
                new GenericMsDataFile(new[] { null, empty, noRetentionTime }, null).ExportAsMgf(path);

                string[] lines = File.ReadAllLines(path);

                // only the one scan that has peaks is written
                NUnit.Framework.Assert.That(lines.Count(l => l == "BEGIN IONS"), Is.EqualTo(1));
                NUnit.Framework.Assert.That(lines, Has.Member("SCANS=2"));
                NUnit.Framework.Assert.That(lines, Has.No.Member("SCANS=1"));

                // NaN retention time drops the line rather than writing "RTINSECONDS=NaN"
                NUnit.Framework.Assert.That(lines.Any(l => l.StartsWith("RTINSECONDS")), Is.False);

                // a null NativeId still yields a parseable TITLE
                NUnit.Framework.Assert.That(lines.First(l => l.StartsWith("TITLE=")), Does.Contain("NativeID:\"\""));
            }
            finally
            {
                File.Delete(path);
            }
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

            NUnit.Framework.Assert.That(reader.NumSpectra, Is.EqualTo(8));

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

            // Scans 5-8 cover the header values that were PRESENT but unparseable. Guarding on the
            // "=" alone only proved a value existed, so these still took the whole file down -- as
            // FormatException instead of IndexOutOfRangeException, and MSLEVEL with no "=" was not
            // guarded at all and still threw IndexOutOfRangeException, the exception this fixture
            // exists to rule out.

            // scan 5: "PEPMASS=" is empty
            var emptyPepmass = reader.GetOneBasedScan(5);
            NUnit.Framework.Assert.That(emptyPepmass.MassSpectrum.Size, Is.EqualTo(2));

            // scan 6: non-numeric PEPMASS, RTINSECONDS and SCANS, and "CHARGE=+" whose sign strip
            // leaves an empty string
            var nonNumericHeaders = reader.GetOneBasedScan(6);
            NUnit.Framework.Assert.That(nonNumericHeaders.MassSpectrum.Size, Is.EqualTo(2));
            NUnit.Framework.Assert.That(nonNumericHeaders.MassSpectrum.XArray,
                Is.EqualTo(new[] { 701.5, 702.5 }).Within(1e-9));

            // scan 7: bare "MSLEVEL" with no "=" -- the one site with no value guard at all
            var mslevelWithoutEquals = reader.GetOneBasedScan(7);
            NUnit.Framework.Assert.That(mslevelWithoutEquals.MassSpectrum.Size, Is.EqualTo(2));

            // scan 8: "MSLEVEL=" is empty, and the peak list carries "(901.5)" and "1400.0-".
            // Those two are REJECTED rather than read as negatives: NumberStyles.Any would have
            // accepted both as -901.5 and -1400.0, inventing peaks the old code threw on.
            var parenthesisedAndTrailingSign = reader.GetOneBasedScan(8);
            NUnit.Framework.Assert.That(parenthesisedAndTrailingSign.MassSpectrum.Size, Is.EqualTo(1),
                "a parenthesised m/z and a trailing-sign intensity must be skipped, not negated");
            NUnit.Framework.Assert.That(parenthesisedAndTrailingSign.MassSpectrum.XArray,
                Is.EqualTo(new[] { 903.5 }).Within(1e-9));
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

        /// <summary>
        /// Polarity is only carried by CHARGE, and CHARGE is written only where a charge guess exists, so
        /// scans without one come back positive. The point of the test is that both read paths say so:
        /// a static load and a random-access dynamic read of the same scan must agree, since the dynamic
        /// reader cannot see what preceded the block it seeks to. Covers MS1 and, equally, an MS2 whose
        /// charge guess is null.
        /// </summary>
        [Test]
        public void MgfStaticAndDynamicReadsAgreeOnPolarity()
        {
            var ms1 = new MsDataScan(
                new MzSpectrum(new[] { 300.1, 400.2 }, new[] { 1000.0, 2000.0 }, false),
                oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Negative,
                retentionTime: 1.0, scanWindowRange: new MzRange(300, 401), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 3000.0, injectionTime: null,
                noiseData: null, nativeId: "scan=1");

            // Negative MS2 with a charge guess: CHARGE carries the sign, so this one survives.
            var ms2WithCharge = new MsDataScan(
                new MzSpectrum(new[] { 110.05, 220.11 }, new[] { 500.0, 750.0 }, false),
                oneBasedScanNumber: 2, msnOrder: 2, isCentroid: true, polarity: Polarity.Negative,
                retentionTime: 1.1, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1250.0, injectionTime: null,
                noiseData: null, nativeId: "scan=2", selectedIonMz: 571.8069,
                selectedIonChargeStateGuess: 2, selectedIonIntensity: 999999.0);

            // Negative MS2 with no charge guess: no CHARGE line is written, so the sign is unrecoverable.
            var ms2NoCharge = new MsDataScan(
                new MzSpectrum(new[] { 130.07, 240.13 }, new[] { 400.0, 650.0 }, false),
                oneBasedScanNumber: 3, msnOrder: 2, isCentroid: true, polarity: Polarity.Negative,
                retentionTime: 1.2, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1050.0, injectionTime: null,
                noiseData: null, nativeId: "scan=3", selectedIonMz: 481.2,
                selectedIonChargeStateGuess: null, selectedIonIntensity: 5000.0);

            // A second MS1 after the negative MS2, which is where an inherited polarity would have shown up.
            var trailingMs1 = new MsDataScan(
                new MzSpectrum(new[] { 310.1, 410.2 }, new[] { 1100.0, 2100.0 }, false),
                oneBasedScanNumber: 4, msnOrder: 1, isCentroid: true, polarity: Polarity.Negative,
                retentionTime: 1.3, scanWindowRange: new MzRange(300, 411), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 3200.0, injectionTime: null,
                noiseData: null, nativeId: "scan=4");

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfPolarityAgreement.mgf");
            new GenericMsDataFile(new[] { ms1, ms2WithCharge, ms2NoCharge, trailingMs1 }, null).ExportAsMgf(path);

            try
            {
                var staticFile = MsDataFileReader.GetDataFile(path);
                staticFile.LoadAllStaticData();
                var staticScans = staticFile.GetAllScansList();

                var dynamicReader = MsDataFileReader.GetDataFile(path);
                dynamicReader.InitiateDynamicConnection();
                try
                {
                    for (int scanNumber = 1; scanNumber <= 4; scanNumber++)
                    {
                        MsDataScan dynamicScan = dynamicReader.GetOneBasedScanFromDynamicConnection(scanNumber);
                        MsDataScan staticScan = staticScans[scanNumber - 1];

                        NUnit.Framework.Assert.That(dynamicScan.Polarity, Is.EqualTo(staticScan.Polarity),
                            $"scan {scanNumber}: the static and dynamic readers must not disagree on polarity");
                    }
                }
                finally
                {
                    dynamicReader.CloseDynamicConnection();
                }

                // Only the scan that wrote a CHARGE line keeps its sign; the format cannot carry the rest.
                NUnit.Framework.Assert.That(staticScans[1].Polarity, Is.EqualTo(Polarity.Negative),
                    "the MS2 with a charge guess wrote CHARGE=2- and must read back negative");
                NUnit.Framework.Assert.That(staticScans[0].Polarity, Is.EqualTo(Polarity.Positive),
                    "an MS1 has no CHARGE line, so its polarity is not recoverable from the format");
                NUnit.Framework.Assert.That(staticScans[2].Polarity, Is.EqualTo(Polarity.Positive),
                    "an MS2 with a null charge guess writes no CHARGE either, so it behaves like the MS1");
                NUnit.Framework.Assert.That(staticScans[3].Polarity, Is.EqualTo(Polarity.Positive),
                    "position in the file must not change the answer; polarity is not inherited");
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// A file whose every scan is empty would produce an mgf with no blocks, which the reader cannot
        /// load -- it indexes by the last scan number and throws on an empty array. Better to refuse to
        /// write it than to hand back a file that fails on the way in.
        /// </summary>
        [Test]
        public void MgfWriterRefusesToWriteAFileWithNoBlocks()
        {
            var empty = new MsDataScan(
                new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false),
                oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.0, scanWindowRange: new MzRange(300, 400), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 0.0, injectionTime: null,
                noiseData: null, nativeId: "scan=1");

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfNoBlocks.mgf");
            try
            {
                var exception = NUnit.Framework.Assert.Throws<MzLibException>(
                    () => new GenericMsDataFile(new[] { empty }, null).ExportAsMgf(path));
                NUnit.Framework.Assert.That(exception.Message, Does.Contain("no spectra"));

                // Refusing before opening the writer means no half-written file is left behind.
                NUnit.Framework.Assert.That(File.Exists(path), Is.False);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

    }
}