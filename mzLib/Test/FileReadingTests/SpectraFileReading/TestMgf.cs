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

        private static MsDataScan Ms1(int scanNumber, double[] mz, double[] intensity, double tic) =>
            new MsDataScan(new MzSpectrum(mz, intensity, false),
                oneBasedScanNumber: scanNumber, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.5, scanWindowRange: new MzRange(mz[0] - 1, mz[^1] + 1), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: tic, injectionTime: null,
                noiseData: null, nativeId: "scan=" + scanNumber);

        private static MsDataScan Ms2WithFullPrecursorMetadata(int scanNumber, int precursorScanNumber) =>
            new MsDataScan(new MzSpectrum(new[] { 110.05, 220.11 }, new[] { 500.0, 750.0 }, false),
                oneBasedScanNumber: scanNumber, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.6, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1250.0, injectionTime: null,
                noiseData: null, nativeId: "scan=" + scanNumber, selectedIonMz: 571.8069,
                selectedIonChargeStateGuess: 2, selectedIonIntensity: 999999.0,
                isolationMZ: 572.0, isolationWidth: 3.0, dissociationType: DissociationType.HCD,
                oneBasedPrecursorScanNumber: precursorScanNumber, selectedIonMonoisotopicGuessMz: 571.8069);

        /// <summary>
        /// The five fields the format has no standard home for. Without them a search reading the file back
        /// cannot deconvolute a precursor, because IsolationRange stays null and there is no MS1 to point at.
        /// </summary>
        [Test]
        public void MgfRoundTripPreservesTheExtensionHeaders()
        {
            var file = new GenericMsDataFile(
                new[] { Ms1(1, new[] { 300.1, 400.2 }, new[] { 1000.0, 2000.0 }, 6000.0),
                        Ms2WithFullPrecursorMetadata(2, 1) }, null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfExtensionHeaders.mgf");
            file.ExportAsMgf(path);

            try
            {
                var reread = MsDataFileReader.GetDataFile(path);
                reread.LoadAllStaticData();
                var readMs2 = reread.GetAllScansList()[1];

                NUnit.Framework.Assert.That(readMs2.DissociationType, Is.EqualTo(DissociationType.HCD));
                NUnit.Framework.Assert.That(readMs2.OneBasedPrecursorScanNumber, Is.EqualTo(1));
                NUnit.Framework.Assert.That(readMs2.IsolationWidth, Is.EqualTo(3.0).Within(1e-9));
                NUnit.Framework.Assert.That(readMs2.IsolationMz, Is.EqualTo(572.0).Within(1e-9));
                NUnit.Framework.Assert.That(readMs2.TotalIonCurrent, Is.EqualTo(1250.0).Within(1e-9));

                // the point of all four precursor headers: without them this is null and deconvolution
                // returns nothing
                NUnit.Framework.Assert.That(readMs2.IsolationRange, Is.Not.Null);
                NUnit.Framework.Assert.That(readMs2.IsolationRange.Minimum, Is.EqualTo(570.5).Within(1e-9));
                NUnit.Framework.Assert.That(readMs2.IsolationRange.Maximum, Is.EqualTo(573.5).Within(1e-9));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// An empty header is worse than an absent one, because a reader cannot tell "unknown" from "zero".
        /// TIC is the exception: every scan has one, so it is written for MS1 blocks too.
        /// </summary>
        [Test]
        public void MgfWriterOmitsPrecursorHeadersItHasNoValueFor()
        {
            var bare = new MsDataScan(new MzSpectrum(new[] { 110.05 }, new[] { 500.0 }, false),
                oneBasedScanNumber: 2, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.6, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 500.0, injectionTime: null,
                noiseData: null, nativeId: "scan=2", selectedIonMz: 571.8069,
                selectedIonChargeStateGuess: 2, selectedIonIntensity: null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfAbsentHeaders.mgf");
            new GenericMsDataFile(new[] { Ms1(1, new[] { 300.1 }, new[] { 1000.0 }, 1000.0), bare }, null)
                .ExportAsMgf(path);

            try
            {
                string[] lines = File.ReadAllLines(path);

                foreach (string key in new[] { "ACTIVATIONMETHOD", "PRECURSORSCAN", "ISOLATIONWIDTH", "ISOLATIONMZ" })
                {
                    NUnit.Framework.Assert.That(lines.Any(l => l.StartsWith(key)), Is.False, key + " should be absent, not empty");
                }

                // the precursor headers describe a precursor, so they never appear on an MS1 block --
                // but TIC does
                int ms1Start = Array.IndexOf(lines, "BEGIN IONS");
                int ms1End = Array.IndexOf(lines, "END IONS");
                string[] ms1Block = lines[ms1Start..ms1End];
                NUnit.Framework.Assert.That(ms1Block, Has.Member("TIC=1000"));
                NUnit.Framework.Assert.That(ms1Block.Any(l => l.StartsWith("PEPMASS")), Is.False);
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The recorded TIC is not the sum of the written peaks -- it predates centroiding and thresholding.
        /// A file without the header must still produce a usable value, which is what every mgf written
        /// before this header existed relies on.
        /// </summary>
        [Test]
        public void MgfReaderPrefersTheTicHeaderAndFallsBackToThePeakSum()
        {
            // TIC deliberately unequal to 500 + 750, so the two sources are distinguishable
            var scan = new MsDataScan(new MzSpectrum(new[] { 110.05, 220.11 }, new[] { 500.0, 750.0 }, false),
                oneBasedScanNumber: 1, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.6, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 99999.0, injectionTime: null,
                noiseData: null, nativeId: "scan=1", selectedIonMz: 571.8, selectedIonChargeStateGuess: 2,
                selectedIonIntensity: null);

            string withTic = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfWithTic.mgf");
            string withoutTic = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfWithoutTic.mgf");
            new GenericMsDataFile(new[] { scan }, null).ExportAsMgf(withTic);
            File.WriteAllLines(withoutTic, File.ReadAllLines(withTic).Where(l => !l.StartsWith("TIC=")));

            try
            {
                var read = MsDataFileReader.GetDataFile(withTic);
                read.LoadAllStaticData();
                NUnit.Framework.Assert.That(read.GetAllScansList()[0].TotalIonCurrent, Is.EqualTo(99999.0).Within(1e-9));

                var readLegacy = MsDataFileReader.GetDataFile(withoutTic);
                readLegacy.LoadAllStaticData();
                NUnit.Framework.Assert.That(readLegacy.GetAllScansList()[0].TotalIonCurrent, Is.EqualTo(1250.0).Within(1e-9));
            }
            finally
            {
                File.Delete(withTic);
                File.Delete(withoutTic);
            }
        }

        /// <summary>
        /// These come from files mzLib did not necessarily write, so one unparseable value must leave the
        /// field unset rather than take the whole file down.
        /// </summary>
        [Test]
        public void MgfReaderIgnoresUnparseableExtensionHeaders()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfGarbageHeaders.mgf");
            File.WriteAllLines(path, new[]
            {
                "BEGIN IONS", "TITLE=junk", "MSLEVEL=2", "PEPMASS=571.8", "CHARGE=2+", "SCANS=1",
                "TIC=not-a-number", "ACTIVATIONMETHOD=Telekinesis", "PRECURSORSCAN=three",
                "ISOLATIONWIDTH=wide", "ISOLATIONMZ=", "110.05 500", "END IONS"
            });

            try
            {
                var read = MsDataFileReader.GetDataFile(path);
                NUnit.Framework.Assert.DoesNotThrow(() => read.LoadAllStaticData());
                var scan = read.GetAllScansList()[0];

                NUnit.Framework.Assert.That(scan.DissociationType, Is.EqualTo(DissociationType.Unknown));
                NUnit.Framework.Assert.That(scan.OneBasedPrecursorScanNumber, Is.Null);
                NUnit.Framework.Assert.That(scan.IsolationWidth, Is.Null);
                NUnit.Framework.Assert.That(scan.TotalIonCurrent, Is.EqualTo(500.0).Within(1e-9));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// Regression for the reader gating peak trimming on ApplyTrimmingToMsMs alone and applying it to
        /// every block. Trimming an MS1 precursor scan strips its isotope envelopes, and a search reading
        /// the file back then finds no precursors at all.
        /// </summary>
        [Test]
        public void MgfReaderTrimsByMsLevelRatherThanTrimmingEveryBlock()
        {
            double[] mz = Enumerable.Range(0, 100).Select(i => 300.0 + i).ToArray();
            double[] intensity = Enumerable.Range(0, 100).Select(i => 1000.0 + i).ToArray();

            var ms2 = new MsDataScan(new MzSpectrum(mz, intensity, false),
                oneBasedScanNumber: 2, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.6, scanWindowRange: new MzRange(299, 400), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensity.Sum(), injectionTime: null,
                noiseData: null, nativeId: "scan=2", selectedIonMz: 571.8,
                selectedIonChargeStateGuess: 2, selectedIonIntensity: null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfTrimByMsLevel.mgf");
            new GenericMsDataFile(new[] { Ms1(1, mz, intensity, intensity.Sum()), ms2 }, null).ExportAsMgf(path);

            try
            {
                var trimMsMsOnly = new FilteringParams(numberOfPeaksToKeepPerWindow: 5, numberOfWindows: 1,
                    applyTrimmingToMs1: false, applyTrimmingToMsMs: true);
                var read = MsDataFileReader.GetDataFile(path);
                read.LoadAllStaticData(trimMsMsOnly);
                var scans = read.GetAllScansList();

                NUnit.Framework.Assert.That(scans[0].MassSpectrum.Size, Is.EqualTo(100), "MS1 must be untouched when ApplyTrimmingToMs1 is false");
                NUnit.Framework.Assert.That(scans[1].MassSpectrum.Size, Is.EqualTo(5), "MS2 must still be trimmed");

                var trimBoth = new FilteringParams(numberOfPeaksToKeepPerWindow: 5, numberOfWindows: 1,
                    applyTrimmingToMs1: true, applyTrimmingToMsMs: true);
                var readBoth = MsDataFileReader.GetDataFile(path);
                readBoth.LoadAllStaticData(trimBoth);

                NUnit.Framework.Assert.That(readBoth.GetAllScansList()[0].MassSpectrum.Size, Is.EqualTo(5), "MS1 must be trimmed when asked");
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// "Telekinesis" is the case Enum.TryParse already rejects, so it cannot fail against the bug that
        /// matters. TryParse also accepts the NUMERIC form of an enum, and numeric activation codes are
        /// exactly what turns up in files mzLib did not write.
        ///
        /// The two numeric cases fail differently and both need pinning. "999" parses to an undefined
        /// DissociationType of 999, which no switch matches and every default arm swallows. "3" is the
        /// dangerous one: it parses to SID, a real member -- not a skipped field but a confidently wrong
        /// one, and Enum.IsDefined does NOT catch it, because SID is defined. Only rejecting the numeric
        /// form does.
        ///
        /// "HCD,CID" is the third shape, and it defeats IsDefined for the same reason: TryParse accepts a
        /// comma-separated list for any enum and ORs the members, and since CID is 0 the result is plain
        /// HCD. A file naming two activations would be read as having named one.
        /// </summary>
        [Test]
        [TestCase("999")]
        [TestCase("3")]
        [TestCase("-1")]
        [TestCase("Telekinesis")]
        [TestCase("HCD,CID")]
        public void MgfReaderTreatsAnUnrecognisedActivationMethodAsUnknown(string activationMethod)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "mgfActivation" + activationMethod.Replace("-", "neg") + ".mgf");
            File.WriteAllLines(path, new[]
            {
                "BEGIN IONS", "TITLE=junk", "MSLEVEL=2", "PEPMASS=571.8", "CHARGE=2+", "SCANS=1",
                "ACTIVATIONMETHOD=" + activationMethod, "110.05 500", "END IONS"
            });

            try
            {
                var read = MsDataFileReader.GetDataFile(path);
                read.LoadAllStaticData();

                NUnit.Framework.Assert.That(read.GetAllScansList()[0].DissociationType,
                    Is.EqualTo(DissociationType.Unknown),
                    "an activation method we cannot read must be Unknown, never a value we half-recognised");
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// A written PRECURSORSCAN may only name a block this file actually contains. WriteMgf filters
        /// scansToWrite twice -- empty spectra, and includeMs1Scans -- and the header was emitted from the
        /// scan alone, so it could point at a block that was dropped.
        ///
        /// Read back, that populates OneBasedPrecursorScanNumber where master left null, and
        /// MsDataFile.GetOneBasedScan is a raw Scans[n - 1], so the missing slot returns null rather than
        /// throwing and the failure surfaces far from the mgf that caused it.
        ///
        /// The assertion is on the file and on the round trip rather than on a downstream export, because
        /// mgf-to-mzML export of an MS1-less file throws for a separate, pre-existing reason that has
        /// nothing to do with this header.
        /// </summary>
        [Test]
        public void MgfWriterOmitsPrecursorScanWhenThatScanWasNotWritten()
        {
            var file = new GenericMsDataFile(
                new[] { Ms1(1, new[] { 300.1, 400.2 }, new[] { 1000.0, 2000.0 }, 3000.0),
                        Ms2WithFullPrecursorMetadata(2, 1) }, null);

            string ms2Only = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfDanglingPrecursor.mgf");
            file.ExportAsMgf(ms2Only, includeMs1Scans: false);

            try
            {
                NUnit.Framework.Assert.That(File.ReadAllLines(ms2Only).Any(l => l.StartsWith("PRECURSORSCAN=")), Is.False,
                    "MS1 scan 1 was not written, so nothing may point at it");

                var read = MsDataFileReader.GetDataFile(ms2Only);
                read.LoadAllStaticData();
                NUnit.Framework.Assert.That(read.GetAllScansList()[0].OneBasedPrecursorScanNumber, Is.Null,
                    "a precursor pointer that resolves to nothing is worse than an absent one");

                // The control: when the MS1 IS written, the pointer is written too and still round-trips.
                string withMs1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfLivePrecursor.mgf");
                file.ExportAsMgf(withMs1);
                try
                {
                    var readBoth = MsDataFileReader.GetDataFile(withMs1);
                    readBoth.LoadAllStaticData();
                    NUnit.Framework.Assert.That(readBoth.GetAllScansList()[1].OneBasedPrecursorScanNumber, Is.EqualTo(1),
                        "dropping the header entirely would be the wrong fix");
                }
                finally
                {
                    File.Delete(withMs1);
                }
            }
            finally
            {
                File.Delete(ms2Only);
            }
        }

        /// <summary>
        /// The same dangling pointer without the new flag, which is why it had to be fixed at the source
        /// rather than inside the includeMs1Scans branch: WriteMgf has always dropped scans with no peaks,
        /// so a peakless MS1 leaves the MS2 pointing at a block that was never written, at default settings.
        /// </summary>
        [Test]
        public void MgfWriterOmitsPrecursorScanWhenThePrecursorHadNoPeaks()
        {
            var peaklessMs1 = new MsDataScan(new MzSpectrum(new double[0], new double[0], false),
                oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.5, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 0.0, injectionTime: null,
                noiseData: null, nativeId: "scan=1");

            var file = new GenericMsDataFile(
                new[] { peaklessMs1, Ms2WithFullPrecursorMetadata(2, 1) }, null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfPeaklessPrecursor.mgf");
            file.ExportAsMgf(path);

            try
            {
                NUnit.Framework.Assert.That(File.ReadAllLines(path).Count(l => l == "BEGIN IONS"), Is.EqualTo(1),
                    "the peakless MS1 is dropped by the pre-existing empty-spectrum filter");
                NUnit.Framework.Assert.That(File.ReadAllLines(path).Any(l => l.StartsWith("PRECURSORSCAN=")), Is.False,
                    "so nothing may point at it, flag or no flag");

                var read = MsDataFileReader.GetDataFile(path);
                read.LoadAllStaticData();
                NUnit.Framework.Assert.That(read.GetAllScansList()[0].OneBasedPrecursorScanNumber, Is.Null);
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The specification requires a PEPMASS in every block and an MS1 has none, so a file containing
        /// MS1 blocks is out of spec -- MSToolkit, and therefore Comet, exits on the first one.
        /// </summary>
        [Test]
        public void MgfWriterCanOmitMs1ScansToStayWithinTheSpecification()
        {
            var file = new GenericMsDataFile(
                new[] { Ms1(1, new[] { 300.1, 400.2 }, new[] { 1000.0, 2000.0 }, 3000.0),
                        Ms2WithFullPrecursorMetadata(2, 1) }, null);

            string withMs1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfWithMs1.mgf");
            string ms2Only = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfMs2Only.mgf");
            file.ExportAsMgf(withMs1);
            file.ExportAsMgf(ms2Only, includeMs1Scans: false);

            try
            {
                // default is unchanged: MS1 present, and it is the block that has no PEPMASS
                var readDefault = MsDataFileReader.GetDataFile(withMs1);
                readDefault.LoadAllStaticData();
                NUnit.Framework.Assert.That(readDefault.GetAllScansList().Count, Is.EqualTo(2));
                NUnit.Framework.Assert.That(File.ReadAllLines(withMs1).Count(l => l == "BEGIN IONS"), Is.EqualTo(2));

                var readMs2Only = MsDataFileReader.GetDataFile(ms2Only);
                readMs2Only.LoadAllStaticData();
                var scans = readMs2Only.GetAllScansList();
                NUnit.Framework.Assert.That(scans.Count, Is.EqualTo(1));
                NUnit.Framework.Assert.That(scans[0].MsnOrder, Is.EqualTo(2));

                // the whole point: every block carries the mandatory PEPMASS
                string[] lines = File.ReadAllLines(ms2Only);
                NUnit.Framework.Assert.That(lines.Count(l => l == "BEGIN IONS"),
                    Is.EqualTo(lines.Count(l => l.StartsWith("PEPMASS="))));
            }
            finally
            {
                File.Delete(withMs1);
                File.Delete(ms2Only);
            }
        }

        /// <summary>
        /// Dropping the MS1 scans can empty the file, and a blockless mgf is one the reader cannot load.
        /// </summary>
        [Test]
        public void MgfWriterRefusesToWriteWhenOmittingMs1LeavesNothing()
        {
            var ms1Only = new GenericMsDataFile(
                new[] { Ms1(1, new[] { 300.1, 400.2 }, new[] { 1000.0, 2000.0 }, 3000.0) }, null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "mgfNothingLeft.mgf");

            try
            {
                var ex = NUnit.Framework.Assert.Throws<MzLibException>(
                    () => ms1Only.ExportAsMgf(path, includeMs1Scans: false));
                NUnit.Framework.Assert.That(ex.Message, Does.Contain("no spectra"));
                NUnit.Framework.Assert.That(File.Exists(path), Is.False);
            }
            finally
            {
                if (File.Exists(path)) { File.Delete(path); }
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