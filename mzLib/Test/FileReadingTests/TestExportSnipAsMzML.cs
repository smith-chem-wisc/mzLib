using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    public class TestExportSnipAsMzML
    {
        private string _dataPath;
        private MsDataFile _originalFile;

        [SetUp]
        public void Setup()
        {
            _dataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast.mzml");
            _originalFile = MsDataFileReader.GetDataFile(_dataPath);
        }


        /// <summary>
        /// The test data file goes 
        /// 1     -> Ms1
        /// 2-5   -> MS2
        /// 6     -> MS1
        /// 7-12  -> MS2
        /// 13    -> MS1
        /// 14-16 -> MS2
        /// 17    -> MS1
        /// 18-21 -> MS2
        /// ...
        /// 141   -> MS2
        /// 142   -> MS1
        /// </summary>
        [TestCase(1, 10, 10, TestName = "Snip_StartsAtMS1_EndsBeforeMS1")]
        [TestCase(2, 11, 6, TestName = "Snip_StartsAfterMS1_EndsAfterMS1")]
        [TestCase(1, 142, 142, TestName = "Snip_WholeFile")]
        [TestCase(1, 1, 1, TestName = "Snip_SingleMS1")]
        [TestCase(1, 2, 2, TestName = "Snip_MS1AndFirstMS2")]
        [TestCase(10, 20, 8, TestName = "Snip_MiddleSection")]
        [TestCase(141, 142, 1, TestName = "Snip_LastTwoScans")]
        public void ExportSnipAsMzML_ScanNumberingIsCorrect(int startScan, int endScan, int expectedCount)
        {
            // Export snip
            _originalFile.ExportSnipAsMzML(startScan, endScan);

            // Build expected output path
            string outPath = _dataPath.Replace(".mzml", $"_snip_{startScan}-{endScan}.mzML");

            Assert.That(File.Exists(outPath), $"Output file {outPath} was not created.");

            // Read snipped file
            var snippedFile = MsDataFileReader.GetDataFile(outPath);
            var snippedScans = snippedFile.GetAllScansList();

            // Check scan count
            Assert.That(snippedScans.Count, Is.EqualTo(expectedCount), "Unexpected number of scans in snipped file.");

            // Check scan numbers are 1-based and sequential
            var scanNumbers = snippedScans.Select(s => s.OneBasedScanNumber).ToArray();
            var expectedScanNumbers = Enumerable.Range(1, expectedCount).ToArray();
            Assert.That(scanNumbers, Is.EqualTo(expectedScanNumbers), "Scan numbers are not sequential and 1-based.");

            // Check precursor scan numbers (if any MS2s)
            foreach (var scan in snippedScans.Where(s => s.MsnOrder > 1 && s.OneBasedPrecursorScanNumber.HasValue))
            {
                Assert.That(scan.OneBasedPrecursorScanNumber.Value, Is.GreaterThanOrEqualTo(1), "Precursor scan number is less than 1.");
                Assert.That(scan.OneBasedPrecursorScanNumber.Value, Is.LessThanOrEqualTo(expectedCount), "Precursor scan number is out of range.");
            }

            File.Delete(outPath);
        }
    }
}
