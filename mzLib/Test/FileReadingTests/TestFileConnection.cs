using System;
using System.Collections.Generic;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System.IO;
using System.Printing.IndexedProperties;
using System.Threading;
using System.Windows.Shapes;
using Path = System.IO.Path;

namespace Test.FileReadingTests
{
    [TestFixture]
    public sealed class TestReaderConnection
    {
        [Test]
        //[Parallelizable(ParallelScope.All)]
        [TestCase(@"FileReadingTests/TestConnectionFiles/sliced_ethcd.mzML", "sliced_ethcd.mzML")]
        [TestCase(@"FileReadingTests/TestConnectionFiles/sliced_ethcd.raw", "sliced_ethcd.raw")]
        [TestCase(@"FileReadingTests/TestConnectionFiles/small.RAW", "small.RAW")]
        [TestCase(@"FileReadingTests/TestConnectionFiles/SmallCalibratibleYeast.mzml", "SmallCalibratibleYeast.mzml")]
        public static void TestReaderClosesConnection(string filePath, string fileName)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            MsDataFile datafile = MsDataFileReader.GetDataFile(spectraPath);

            List<MsDataScan> scans = datafile.GetAllScansList();

            string movingDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "FileReadingTests/MoveHere");
            
            Directory.CreateDirectory(movingDirectory);

            File.Move(spectraPath, movingDirectory + '/' + fileName);

            Directory.Delete("FileReadingTests/MoveHere", true);

            Assert.Pass();
        }
    }
}
