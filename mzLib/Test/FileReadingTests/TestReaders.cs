using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestReaders
    {
        [Test]
        [TestCase("DataFiles/small.RAW", 48, "Thermo nativeID format")]
        [TestCase("DataFiles/sliced_ethcd.raw", 6, "Thermo nativeID format")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", 142, "Thermo nativeID format")]
        [TestCase("DataFiles/tester.mzML", 7, null)]
        [TestCase("DataFiles/tester.mgf", 5, "no nativeID format")]
        public static void TestLoadingRawFilesAndSourceFiles(string filePath, int expectedScanCount, string sourceFormat)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            MsDataFile datafile = MsDataFileReader.GetDataFile(spectraPath);
            List<MsDataScan> scans = datafile.GetAllScansList();
            Assert.That(scans.Count == expectedScanCount);

            SourceFile file = datafile.GetSourceFile();
            Assert.That(file.NativeIdFormat == sourceFormat);
        }
    }
}
