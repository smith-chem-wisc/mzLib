using System.Diagnostics;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestBruker
    {
        static string _centroidPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
            "centroid_1x_MS1_4x_autoMS2.d");
        private string _profilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
            "profile_1x_MS1_4x_autoMS2.d");
        private string _profileAndCentroid = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
            "profile_and_centroid_1x_MS1_4x_autoMS2.d"); 

        [Test]
        public void TestConstructors()
        {
            var reader = MsDataFileReader.GetDataFile(_centroidPath); 
            Assert.That(reader, !Is.Null);
        }

        [Test]
        public void TestFileDoesntExist()
        {
            string fakePath = "fakePath.d";
            Assert.Throws<FileNotFoundException>(() =>
                MsDataFileReader.GetDataFile(fakePath));
        }

        [Test]
        public void TestLoadAllStaticDataCentroid()
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(_centroidPath).LoadAllStaticData();
            Assert.That(brukerData.NumSpectra, Is.EqualTo(5));
            Assert.That(brukerData.Scans[1].Polarity == Polarity.Positive);
            Assert.That(brukerData.Scans[1].DissociationType == DissociationType.CID);
            Assert.That(brukerData.Scans[1].TotalIonCurrent == 346d);
            Assert.That(brukerData.Scans[1].NativeId == "scan=2");
            Assert.That(brukerData.Scans[1].SelectedIonMZ, Is.EqualTo(721.86865).Within(0.001));
            Assert.That(brukerData.Scans[1].MsnOrder == 2);
            Assert.That(brukerData.Scans[1].IsCentroid);
            // Isolation width is read from the (Per)SpectrumVariables table (Quadrupole_IsolationResolution_Act).
            // Without it, IsolationRange is null and precursor deconvolution finds nothing (empty search).
            Assert.That(brukerData.Scans[1].IsolationWidth, Is.EqualTo(10).Within(1e-6));
            Assert.That(brukerData.Scans[1].IsolationRange, Is.Not.Null);
            Assert.That(brukerData.Scans[1].IsolationRange.Minimum, Is.EqualTo(721.86865 - 5).Within(0.001));
            Assert.That(brukerData.Scans[1].IsolationRange.Maximum, Is.EqualTo(721.86865 + 5).Within(0.001));
            // Retention time is converted from Bruker seconds to mzLib minutes (4.232 s -> 0.07053 min).
            Assert.That(brukerData.Scans[1].RetentionTime, Is.EqualTo(4.232 / 60.0).Within(1e-4));
            // This fixture's precursor charge state is recorded as 0 (undetermined) -> exposed as null, not 0.
            Assert.That(brukerData.Scans[1].SelectedIonChargeStateGuess, Is.Null);
        }

        [Test]

        public void TestLoadAllStaticDataProfile()
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(_profilePath).LoadAllStaticData();
            Assert.That(brukerData.NumSpectra, Is.EqualTo(5));
            Assert.That(brukerData.Scans[1].Polarity == Polarity.Positive);
            Assert.That(brukerData.Scans[1].DissociationType == DissociationType.CID);
            Assert.That(!brukerData.Scans[1].IsCentroid);
            Assert.That(brukerData.Scans[1].TotalIonCurrent == -1d);
            Assert.That(brukerData.Scans[1].NativeId == "scan=2");
            Assert.That(brukerData.Scans[1].SelectedIonMZ, Is.EqualTo(716.58715).Within(0.001));
            Assert.That(brukerData.Scans[1].MsnOrder == 2);
        }

        [Test]
        public void TestLoadAllStaticDataProfileAndCentroid()
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(_profileAndCentroid).LoadAllStaticData();
            Assert.That(brukerData.NumSpectra, Is.EqualTo(5));
            Assert.That(brukerData.Scans[1].Polarity == Polarity.Positive);
            Assert.That(brukerData.Scans[1].DissociationType == DissociationType.CID);
            // If centroided and profile are both present, we default to centroid data. 
            Assert.That(brukerData.Scans[1].IsCentroid);
            Assert.That(brukerData.Scans[1].TotalIonCurrent == 210d);
            Assert.That(brukerData.Scans[1].NativeId == "scan=2");
            Assert.That(brukerData.Scans[1].SelectedIonMZ, Is.EqualTo(1280.748901).Within(0.001));
            Assert.That(brukerData.Scans[1].MsnOrder == 2);
        }

        [Test]
        public void TestGetSourceFile()
        {
            var sourceFile = MsDataFileReader.GetDataFile(_centroidPath).GetSourceFile();
            Assert.That(_centroidPath == sourceFile.Uri.OriginalString);
            Assert.That(sourceFile.FileName == "analysis.baf");
            Assert.That(sourceFile.MassSpectrometerFileFormat == "mzML format");
            Assert.That(sourceFile.NativeIdFormat == "scan number only nativeID format");
        }

        [Test]
        public void TestDynamicConnection()
        {
            MsDataFile brukerReader = MsDataFileReader.GetDataFile(_centroidPath);
            brukerReader.InitiateDynamicConnection();
            var scan = brukerReader.GetOneBasedScanFromDynamicConnection(2);
            
            Assert.That(scan.Polarity == Polarity.Positive);
            Assert.That(scan.DissociationType == DissociationType.CID);
            Assert.That(scan.TotalIonCurrent == 346d);
            Assert.That(scan.NativeId == "scan=2");
            Assert.That(scan.SelectedIonMZ, Is.EqualTo(721.86865).Within(0.001));
            Assert.That(scan.MsnOrder == 2);
            Assert.That(scan.IsCentroid);
            // The dynamic path must populate isolation width too, otherwise IsolationRange is null.
            Assert.That(scan.IsolationWidth, Is.EqualTo(10).Within(1e-6));
            Assert.That(scan.IsolationRange, Is.Not.Null);
        }

        [Test]
        public void TestDynamicConnection_AfterStaticLoading()
        {
            MsDataFile brukerReader = MsDataFileReader.GetDataFile(_centroidPath);
            brukerReader.LoadAllStaticData();
            brukerReader.InitiateDynamicConnection();
            var scan = brukerReader.GetOneBasedScanFromDynamicConnection(2);

            Assert.That(scan.Polarity == Polarity.Positive);
            Assert.That(scan.DissociationType == DissociationType.CID);
            Assert.That(scan.TotalIonCurrent == 346d);
            Assert.That(scan.NativeId == "scan=2");
            Assert.That(scan.SelectedIonMZ, Is.EqualTo(721.86865).Within(0.001));
            Assert.That(scan.MsnOrder == 2);
            Assert.That(scan.IsCentroid);
        }

        [Test]
        public void TestDynamicConnectionToAllScans()
        {
            MsDataFile brukerReader = MsDataFileReader.GetDataFile(_centroidPath);
            brukerReader.InitiateDynamicConnection();
            int counter = 5;
            while (counter > 0)
            {
                Assert.DoesNotThrow(delegate
                {
                    brukerReader.GetOneBasedScanFromDynamicConnection(counter);
                });
                counter--; 
            }
        }

        [Test]
        public void TestOpenAndCloseConnection()
        {
            MsDataFile brukerReader = MsDataFileReader.GetDataFile(_centroidPath);
            Assert.DoesNotThrow(delegate
            {
                brukerReader.InitiateDynamicConnection();
            });
            Assert.DoesNotThrow(delegate
            {
                brukerReader.CloseDynamicConnection();
            });
        }

        [Test]
        public void TestPeakFiltering()
        {
            FilteringParams filteringParams = new(null, 0.5);
            var scan = MsDataFileReader.GetDataFile(_centroidPath).LoadAllStaticData(filteringParams).Scans[0];
            Assert.That(scan.MassSpectrum.XArray.Length == 1);
        }

        // Loading with a FilteringParams that has ApplyTrimmingToMsMs = true exercises the window-trimming
        // branches of GetSpectraData for both centroid (line) and profile spectra. Profile data previously
        // threw IndexOutOfRangeException here (it indexed profileMzs[^0], i.e. one past the end); centroid
        // was unaffected. This mirrors how MetaMorpheus's MyFileManager loads files (it always supplies a
        // FilteringParams unless DissociationType is LowCID).
        //
        // Reaching the fixed line without throwing is not sufficient on its own: a future regression could
        // avoid the exception yet silently drop the boundary peak or hand back an empty spectrum. So beyond
        // "doesn't throw", assert the trimmed MS/MS spectra actually have content and that their retained
        // m/z bounds sit inside the untrimmed m/z range (the profileMzs[0] / profileMzs[^1] values the fix
        // uses as the window bounds).
        [Test]
        [TestCase("centroid_1x_MS1_4x_autoMS2.d")]
        [TestCase("profile_1x_MS1_4x_autoMS2.d")]
        [TestCase("profile_and_centroid_1x_MS1_4x_autoMS2.d")]
        public void TestLoadAllStaticDataWithMsMsTrimming(string fileName)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);
            // 7th argument (applyTrimmingToMsMs) = true; the others mirror a typical MetaMorpheus setup.
            var filteringParams = new FilteringParams(200, 0.01, null, 1, false, false, true);

            MsDataFile brukerData = null;
            Assert.DoesNotThrow(() =>
                brukerData = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams, 1));
            Assert.That(brukerData.NumSpectra, Is.EqualTo(5));

            // Baseline load with no trimming: captures each MS/MS scan's full m/z range, i.e. the exact
            // first/last values the trimming window is built from.
            var untrimmed = MsDataFileReader.GetDataFile(path).LoadAllStaticData();

            var ms2Scans = brukerData.Scans.Where(s => s.MsnOrder == 2).ToList();
            Assert.That(ms2Scans, Is.Not.Empty);
            foreach (var scan in ms2Scans)
            {
                double[] trimmed = scan.MassSpectrum.XArray;
                double[] full = untrimmed.GetOneBasedScan(scan.OneBasedScanNumber).MassSpectrum.XArray;

                Assert.That(full, Is.Not.Empty, $"scan {scan.OneBasedScanNumber} unexpectedly had no peaks to trim");
                // Content survived trimming (catches "returns an empty spectrum").
                Assert.That(trimmed, Is.Not.Empty, $"scan {scan.OneBasedScanNumber} was trimmed to an empty spectrum");
                // Peaks stay sorted ascending and positive: [0] and [^1] really are the min/max m/z bounds.
                Assert.That(trimmed[0], Is.GreaterThan(0));
                Assert.That(trimmed, Is.Ordered);
                // Retained first/last m/z fall within the untrimmed window (catches a bad boundary bound).
                Assert.That(trimmed[0], Is.GreaterThanOrEqualTo(full[0]));
                Assert.That(trimmed[^1], Is.LessThanOrEqualTo(full[^1]));
                // Trimming can only remove peaks, never add them, and honors the per-window cap of 200.
                Assert.That(trimmed.Length, Is.LessThanOrEqualTo(full.Length));
                Assert.That(trimmed.Length, Is.LessThanOrEqualTo(200));
            }
        }

        // Regression guard for the "search returns nothing" bug: the reader never populated IsolationWidth,
        // so MsDataScan.IsolationRange was null for every MS2 scan, GetIsolatedMassesAndCharges returned an
        // empty list, and every consumer (e.g. MetaMorpheus) dropped every MS2 scan -> zero identifications.
        // Every MS2 scan must expose a positive isolation width and a non-null isolation range that brackets
        // its precursor m/z. This is the reader's contract; with it, GetIsolatedMassesAndCharges (and hence
        // precursor deconvolution in search engines) can run at all. Previously IsolationWidth was never set,
        // so IsolationRange was null and the method short-circuited to an empty list for every scan.
        [Test]
        [TestCase("centroid_1x_MS1_4x_autoMS2.d")]
        [TestCase("profile_1x_MS1_4x_autoMS2.d")]
        [TestCase("profile_and_centroid_1x_MS1_4x_autoMS2.d")]
        public void TestMs2ScansHaveUsableIsolationRange(string fileName)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);
            MsDataFile brukerData = MsDataFileReader.GetDataFile(path).LoadAllStaticData();

            var ms2Scans = brukerData.Scans.Where(s => s.MsnOrder == 2).ToList();
            Assert.That(ms2Scans, Is.Not.Empty);

            foreach (var scan in ms2Scans)
            {
                Assert.That(scan.IsolationWidth, Is.Not.Null.And.GreaterThan(0),
                    $"scan {scan.OneBasedScanNumber} has no isolation width");
                Assert.That(scan.IsolationRange, Is.Not.Null,
                    $"scan {scan.OneBasedScanNumber} has a null isolation range");
                Assert.That(scan.IsolationRange.Contains(scan.SelectedIonMZ.Value), Is.True,
                    $"scan {scan.OneBasedScanNumber} isolation range does not contain its precursor m/z");
                // IsolationRange being non-null is exactly what GetIsolatedMassesAndCharges requires to run;
                // when it is null (the old behavior) the method returns an empty list and no precursor is found.
                Assert.That(scan.GetIsolatedMassesAndCharges(brukerData.GetOneBasedScan(
                    scan.OneBasedPrecursorScanNumber.Value).MassSpectrum,
                    new ClassicDeconvolutionParameters(1, 6, 20, 3)), Is.Not.Null);
            }
        }

        // Peak trimming must respect the per-MS-level FilteringParams flags. The reader previously trimmed
        // every scan whenever ApplyTrimmingToMsMs was set, so it trimmed MS1 precursor scans even when
        // ApplyTrimmingToMs1 was false (as MetaMorpheus sets it). That stripped precursor isotope envelopes
        // from the MS1 and precursor deconvolution then found nothing for most MS2 scans.
        [Test]
        public void TestMs1NotTrimmedWhenApplyTrimmingToMs1IsFalse()
        {
            // The profile fixture has a rich MS1 (many peaks), so trimming is clearly observable.
            var raw = MsDataFileReader.GetDataFile(_profilePath).LoadAllStaticData();
            int rawMs1Peaks = raw.GetOneBasedScan(1).MassSpectrum.XArray.Length;   // scan 1 is MS1
            int rawMs2Peaks = raw.GetOneBasedScan(2).MassSpectrum.XArray.Length;   // scan 2 is MS2
            Assert.That(rawMs1Peaks, Is.GreaterThan(1), "fixture MS1 needs >1 peak to make trimming observable");
            Assert.That(rawMs2Peaks, Is.GreaterThan(1));

            // Aggressive trim (keep 200 peaks/window, 1 window) with MS1 trimming DISABLED, MS/MS ENABLED --
            // this is the shape MetaMorpheus uses.
            var ms1Off = new FilteringParams(200, null, null, 1, false,
                applyTrimmingToMs1: false, applyTrimmingToMsMs: true, applyTrimmingToMsN: true);
            var withMs1Off = MsDataFileReader.GetDataFile(_profilePath).LoadAllStaticData(ms1Off);
            // MS1 must be untouched (this is the regression: it used to be trimmed anyway)...
            Assert.That(withMs1Off.GetOneBasedScan(1).MassSpectrum.XArray.Length, Is.EqualTo(rawMs1Peaks));
            // ...while MS2 is trimmed.
            Assert.That(withMs1Off.GetOneBasedScan(2).MassSpectrum.XArray.Length, Is.LessThan(rawMs2Peaks));

            // With MS1 trimming ENABLED, the MS1 scan is trimmed too.
            var ms1On = new FilteringParams(200, null, null, 1, false,
                applyTrimmingToMs1: true, applyTrimmingToMsMs: true, applyTrimmingToMsN: true);
            var withMs1On = MsDataFileReader.GetDataFile(_profilePath).LoadAllStaticData(ms1On);
            Assert.That(withMs1On.GetOneBasedScan(1).MassSpectrum.XArray.Length, Is.LessThan(rawMs1Peaks));
        }

        // Bruker stores retention time in seconds; mzLib's convention is minutes. The reader must convert.
        [Test]
        public void TestRetentionTimeConvertedToMinutes()
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(_centroidPath).LoadAllStaticData();
            // Raw Bruker Rt values (seconds) for this fixture are 2.165, 4.232, 6.287, 8.341, 10.396.
            double[] expectedSeconds = { 2.165, 4.232, 6.287, 8.341, 10.396 };
            for (int i = 0; i < expectedSeconds.Length; i++)
            {
                Assert.That(brukerData.Scans[i].RetentionTime, Is.EqualTo(expectedSeconds[i] / 60.0).Within(1e-4));
            }
        }
    }
}
