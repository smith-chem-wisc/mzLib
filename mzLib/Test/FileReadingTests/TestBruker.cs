using System.IO;
using MassSpectrometry;
using NUnit; 
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
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
            var reader = MsDataFileReader.GetDataFile(fakePath);
            Assert.Throws<FileNotFoundException>(() =>
                reader.InitiateDynamicConnection());
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

        [Test]
        public void TestFileNotFoundExceptionThrown()
        {
            MsDataFile brukerReader = MsDataFileReader.GetDataFile("notrealfile.d");
            Assert.Throws<FileNotFoundException>(delegate
            {
                brukerReader.LoadAllStaticData(); 
            });
        }
    }
}
