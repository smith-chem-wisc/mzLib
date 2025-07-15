using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using Path = System.IO.Path;

namespace Test.FileReadingTests
{
    [TestFixture]
    public sealed class TestReaderConnection
    {
        [Test]
        [TestCase(@"DataFiles/sliced_ethcd.mzML", "sliced_ethcd.mzML")]
        [TestCase(@"DataFiles/sliced_ethcd.raw", "sliced_ethcd.raw")]
        [TestCase(@"DataFiles/small.RAW", "small.RAW")]
        [TestCase(@"DataFiles/SmallCalibratibleYeast.mzml", "SmallCalibratibleYeast.mzml")]
        public static void TestReaderClosesConnection(string filePath, string fileName)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            MsDataFile datafile = MsDataFileReader.GetDataFile(spectraPath);

            List<MsDataScan> scans = datafile.GetAllScansList();

            string movingDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "FileReadingTests/MoveHere");

            Directory.CreateDirectory(movingDirectory);

            File.Copy(spectraPath, movingDirectory + '/' + fileName);

            Directory.Delete(movingDirectory, true);

            Assert.Pass();
        }
    }
}
