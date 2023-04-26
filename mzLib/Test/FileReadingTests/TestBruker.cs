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
        [Test]
        [TestCase(@"D:\BurkerFileSupport\SmallFiles\centroid_1x_MS1_4x_autoMS2.d")]
        public void TestConstructors(string path)
        {
            var reader = MsDataFileReader.GetDataFile(path); 
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
        [TestCase(@"D:\BurkerFileSupport\SmallFiles\centroid_1x_MS1_4x_autoMS2.d")]
        public void TestLoadAllStaticDataCentroid(string path)
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(path).LoadAllStaticData();
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
        [TestCase(@"D:\BurkerFileSupport\SmallFiles\profile_1x_MS1_4x_autoMS2.d")]
        public void TestLoadAllStaticDataProfile(string path)
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(path).LoadAllStaticData();
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
        [TestCase(@"D:\BurkerFileSupport\SmallFiles\profile_and_centroid_1x_MS1_4x_autoMS2.d")]
        public void TestLoadAllStaticDataProfileAndCentroid(string path)
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(path).LoadAllStaticData();
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
        [TestCase(@"D:\BurkerFileSupport\SmallFiles\centroid_1x_MS1_4x_autoMS2.d")]
        public void TestGetSourceFile(string path)
        {
            var sourceFile = MsDataFileReader.GetDataFile(path).GetSourceFile();
            Assert.That(path == sourceFile.Uri.OriginalString);
            Assert.That(sourceFile.FileName == "analysis.baf");
            Assert.That(sourceFile.MassSpectrometerFileFormat == "mzML format");
            Assert.That(sourceFile.NativeIdFormat == "scan number only nativeID format");
        }

        [Test]
        [TestCase(@"D:\BurkerFileSupport\SmallFiles\centroid_1x_MS1_4x_autoMS2.d")]
        public void TestDynamicConnection(string path)
        {
            MsDataFile brukerReader = MsDataFileReader.GetDataFile(path);
            brukerReader.InitiateDynamicConnection();
            var scan = brukerReader.GetOneBasedScan(2);
            
            Assert.That(scan.Polarity == Polarity.Positive);
            Assert.That(scan.DissociationType == DissociationType.CID);
            Assert.That(scan.TotalIonCurrent == 346d);
            Assert.That(scan.NativeId == "scan=2");
            Assert.That(scan.SelectedIonMZ, Is.EqualTo(721.86865).Within(0.001));
            Assert.That(scan.MsnOrder == 2);
            Assert.That(scan.IsCentroid);
        }

        [Test]
        [TestCase(@"D:\BurkerFileSupport\SmallFiles\centroid_1x_MS1_4x_autoMS2.d")]
        public void TestPeakFiltering(string path)
        {
            FilteringParams filteringParams = new(null, 0.5);
            var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams).Scans[0];
            Assert.That(scan.MassSpectrum.XArray.Length == 1);
        }
    }
}
