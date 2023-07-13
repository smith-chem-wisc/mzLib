using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using CsvHelper;
using MzLibUtil;
using Newtonsoft.Json;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestMsFeature
    {

        private static string directoryPath;

        [OneTimeSetUp]
        public void SetUp()
        {
            directoryPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ReadingWritingTests");
            Directory.CreateDirectory(directoryPath);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(directoryPath, true);
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature", 1, 7)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms2Feature_FlashDeconvOpenMs3.0.0_ms2.feature", 2, 1)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature", 1, 4)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms2Feature_TopFDv1.6.2_ms2.feature", 2, 5)]
        public void TestFeaturesLoadAndCountIsCorrect(string path, int type, int featureCount)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);

            switch (type)
            {
                case 1:
                    {
                        var featureFile = new Ms1FeatureFile(filePath);
                        Assert.That(featureFile.Count() == featureCount);
                        Assert.That(featureFile.Results.Count() == featureCount);
                        Assert.That(featureFile.CanRead(path));

                        featureFile = FileReader.ReadFile<Ms1FeatureFile>(filePath);
                        Assert.That(featureFile.Count() == featureCount);
                        Assert.That(featureFile.Results.Count() == featureCount);
                        break;
                    }
                case 2:
                    {
                        var featureFile = new Ms2FeatureFile(filePath);
                        Assert.That(featureFile.Count() == featureCount);
                        Assert.That(featureFile.Results.Count() == featureCount);
                        Assert.That(featureFile.CanRead(path));

                        featureFile = FileReader.ReadFile<Ms2FeatureFile>(filePath);
                        Assert.That(featureFile.Count() == featureCount);
                        Assert.That(featureFile.Results.Count() == featureCount);
                        break;
                    }
            }
        }

        [Test]
        public static void TestFlashDeconvMs1FeatureFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature");
            var ms1Features = new Ms1FeatureFile(filePath);
            var first = ms1Features.First();
            var last = ms1Features.Last();

            Assert.That(first.SampleId, Is.EqualTo(0));
            Assert.That(first.Id, Is.EqualTo(1));
            Assert.That(first.Mass, Is.EqualTo(10835.9).Within(0.00001));
            Assert.That(first.Intensity, Is.EqualTo(1.21e+10).Within(0.00001));
            Assert.That(first.RetentionTimeBegin, Is.EqualTo(2375.98));
            Assert.That(first.RetentionTimeEnd, Is.EqualTo(2398.21));
            Assert.That(first.RetentionTimeApex, Is.EqualTo(2390.8));
            Assert.That(first.IntensityApex, Is.EqualTo(null));
            Assert.That(first.ChargeStateMin, Is.EqualTo(7));
            Assert.That(first.ChargeStateMax, Is.EqualTo(18));
            Assert.That(first.FractionIdMin, Is.EqualTo(0));
            Assert.That(first.FractionIdMax, Is.EqualTo(0));

            Assert.That(last.SampleId, Is.EqualTo(0));
            Assert.That(last.Id, Is.EqualTo(16533));
            Assert.That(last.Mass, Is.EqualTo(11299.4).Within(0.00001));
            Assert.That(last.Intensity, Is.EqualTo(813423));
            Assert.That(last.RetentionTimeBegin, Is.EqualTo(5391.7));
            Assert.That(last.RetentionTimeEnd, Is.EqualTo(5393.7));
            Assert.That(last.RetentionTimeApex, Is.EqualTo(5392.7));
            Assert.That(last.IntensityApex, Is.EqualTo(null));
            Assert.That(last.ChargeStateMin, Is.EqualTo(7));
            Assert.That(last.ChargeStateMax, Is.EqualTo(18));
            Assert.That(last.FractionIdMin, Is.EqualTo(0));
            Assert.That(last.FractionIdMax, Is.EqualTo(0));
        }

        [Test]
        public static void TestFlashDeconvMs2FeatureFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\Ms2Feature_FlashDeconvOpenMs3.0.0_ms2.feature");
            var ms1Features = FileReader.ReadFile<Ms2FeatureFile>(filePath);

            var first = ms1Features.First();
            var last = ms1Features.Last();

            Assert.That(first.ZeroBasedScanNumber, Is.EqualTo(-1));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(0));
            Assert.That(first.FractionId, Is.EqualTo(0));
            Assert.That(first.FilePath, Is.EqualTo(
                    "D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(first.Scans, Is.EqualTo(-1));
            Assert.That(first.Ms1Id, Is.EqualTo(-1));
            Assert.That(first.Ms1Scans, Is.EqualTo(-1));
            Assert.That(first.PrecursorMass, Is.EqualTo(11299.4).Within(0.00001));
            Assert.That(first.PrecursorIntensity, Is.EqualTo(813423).Within(0.00001));
            Assert.That(first.FractionFeatureId, Is.EqualTo(16533));
            Assert.That(first.FractionFeatureIntensity, Is.EqualTo(813423).Within(0.00001));
            Assert.That(first.FractionFeatureScore, Is.EqualTo(-1000).Within(0.00001));
            Assert.That(first.FractionFeatureApex, Is.EqualTo(5392.7).Within(0.00001));
            Assert.That(first.SampleFeatureId, Is.EqualTo(16533));
            Assert.That(first.SampleFeatureIntensity, Is.EqualTo(813423).Within(0.00001));

            Assert.That(last.ZeroBasedScanNumber, Is.EqualTo(-1));
            Assert.That(last.OneBasedScanNumber, Is.EqualTo(0));
            Assert.That(last.FractionId, Is.EqualTo(0));
            Assert.That(last.FilePath, Is.EqualTo(
                "D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(last.Scans, Is.EqualTo(-1));
            Assert.That(last.Ms1Id, Is.EqualTo(-1));
            Assert.That(last.Ms1Scans, Is.EqualTo(-1));
            Assert.That(last.PrecursorMass, Is.EqualTo(11299.4).Within(0.00001));
            Assert.That(last.PrecursorIntensity, Is.EqualTo(813423).Within(0.00001));
            Assert.That(last.FractionFeatureId, Is.EqualTo(16533));
            Assert.That(last.FractionFeatureIntensity, Is.EqualTo(813423).Within(0.00001));
            Assert.That(last.FractionFeatureScore, Is.EqualTo(-1000).Within(0.00001));
            Assert.That(last.FractionFeatureApex, Is.EqualTo(5392.7).Within(0.00001));
            Assert.That(last.SampleFeatureId, Is.EqualTo(16533));
            Assert.That(last.SampleFeatureIntensity, Is.EqualTo(813423).Within(0.00001));

        }

        [Test]
        public static void TestTopFDMs1FeatureFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature");
            var ms1Features = FileReader.ReadFile<Ms1FeatureFile>(filePath);

            var first = ms1Features.First();
            var last = ms1Features.Last();

            Assert.That(first.SampleId, Is.EqualTo(0));
            Assert.That(first.Id, Is.EqualTo(0));
            Assert.That(first.Mass, Is.EqualTo(10835.85272090354).Within(0.00001));
            Assert.That(first.Intensity, Is.EqualTo(10849947123.04).Within(0.00001));
            Assert.That(first.RetentionTimeBegin, Is.EqualTo(2372.27));
            Assert.That(first.RetentionTimeEnd, Is.EqualTo(2401.92));
            Assert.That(first.RetentionTimeApex, Is.EqualTo(2390.8));
            Assert.That(first.IntensityApex, Is.EqualTo(912795138.8));
            Assert.That(first.ChargeStateMin, Is.EqualTo(7));
            Assert.That(first.ChargeStateMax, Is.EqualTo(17));
            Assert.That(first.FractionIdMin, Is.EqualTo(0));
            Assert.That(first.FractionIdMax, Is.EqualTo(0));

            Assert.That(last.SampleId, Is.EqualTo(0));
            Assert.That(last.Id, Is.EqualTo(20553));
            Assert.That(last.Mass, Is.EqualTo(1841.06822).Within(0.00001));
            Assert.That(last.Intensity, Is.EqualTo(3279.65));
            Assert.That(last.RetentionTimeBegin, Is.EqualTo(3122.55));
            Assert.That(last.RetentionTimeEnd, Is.EqualTo(3122.55));
            Assert.That(last.RetentionTimeApex, Is.EqualTo(3122.55));
            Assert.That(last.IntensityApex, Is.EqualTo(3279.65));
            Assert.That(last.ChargeStateMin, Is.EqualTo(1));
            Assert.That(last.ChargeStateMax, Is.EqualTo(1));
            Assert.That(last.FractionIdMin, Is.EqualTo(0));
            Assert.That(last.FractionIdMax, Is.EqualTo(0));
        }

        [Test]
        public static void TestTopFDMs2FeatureFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Ms2Feature_TopFDv1.6.2_ms2.feature");
            var ms1Features = new Ms2FeatureFile(filePath);
            var first = ms1Features.First();
            var last = ms1Features.Last();

            Assert.That(first.ZeroBasedScanNumber, Is.EqualTo(0));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(1));
            Assert.That(first.FractionId, Is.EqualTo(0));
            Assert.That(first.FilePath, Is.EqualTo(
                    "D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(first.Scans, Is.EqualTo(614));
            Assert.That(first.Ms1Id, Is.EqualTo(612));
            Assert.That(first.Ms1Scans, Is.EqualTo(613));
            Assert.That(first.PrecursorMass, Is.EqualTo(6937.649520000001).Within(0.00001));
            Assert.That(first.PrecursorIntensity, Is.EqualTo(145784.16).Within(0.00001));
            Assert.That(first.FractionFeatureId, Is.EqualTo(20141));
            Assert.That(first.FractionFeatureIntensity, Is.EqualTo(141018972.3899998).Within(0.00001));
            Assert.That(first.FractionFeatureScore, Is.EqualTo(-1000).Within(0.00001));
            Assert.That(first.FractionFeatureApex, Is.EqualTo(1321.02).Within(0.00001));
            Assert.That(first.SampleFeatureId, Is.EqualTo(20141));
            Assert.That(first.SampleFeatureIntensity, Is.EqualTo(141018972.3899998).Within(0.00001));

            Assert.That(last.ZeroBasedScanNumber, Is.EqualTo(2808));
            Assert.That(last.OneBasedScanNumber, Is.EqualTo(2809));
            Assert.That(last.FractionId, Is.EqualTo(0));
            Assert.That(last.FilePath, Is.EqualTo(
                "D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(last.Scans, Is.EqualTo(5555));
            Assert.That(last.Ms1Id, Is.EqualTo(2745));
            Assert.That(last.Ms1Scans, Is.EqualTo(5554));
            Assert.That(last.PrecursorMass, Is.EqualTo(11299.3603).Within(0.00001));
            Assert.That(last.PrecursorIntensity, Is.EqualTo(117824.17).Within(0.00001));
            Assert.That(last.FractionFeatureId, Is.EqualTo(20402));
            Assert.That(last.FractionFeatureIntensity, Is.EqualTo(10451720.37999999).Within(0.00001));
            Assert.That(last.FractionFeatureScore, Is.EqualTo(-1000).Within(0.00001));
            Assert.That(last.FractionFeatureApex, Is.EqualTo(5267.55).Within(0.00001));
            Assert.That(last.SampleFeatureId, Is.EqualTo(20402));
            Assert.That(last.SampleFeatureIntensity, Is.EqualTo(10451720.37999999).Within(0.00001));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature")]
        public static void TestMs1FeatureReadWrite(string filePath)
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            var testOutputPath = Path.Combine(directoryPath, "featureOut_ms1.feature");

            // load in file
            Ms1FeatureFile deconFile = FileReader.ReadFile<Ms1FeatureFile>(testFilePath);

            // write and reread file
            deconFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<Ms1FeatureFile>(testOutputPath);
            Assert.That(File.Exists(testOutputPath));

            // check are equivalent
            for (int i = 0; i < deconFile.Results.Count; i++)
            {
                var original = JsonConvert.SerializeObject(deconFile.Results[i]);
                var written = JsonConvert.SerializeObject(writtenDeconFile.Results[i]);
                Assert.That(original, Is.EqualTo(written));
            }

            // test writer still works without specifying extensions
            File.Delete(testOutputPath);
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "featureOut");
            deconFile.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms2Feature_FlashDeconvOpenMs3.0.0_ms2.feature")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms2Feature_TopFDv1.6.2_ms2.feature")]
        public static void TestMs2FeatureReadWrite(string filePath)
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            var testOutputPath = Path.Combine(directoryPath, "featureOut_ms2.feature");

            // load in file
            Ms2FeatureFile deconFile = FileReader.ReadFile<Ms2FeatureFile>(testFilePath);

            // write and reread file
            deconFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<Ms2FeatureFile>(testOutputPath);
            Assert.That(File.Exists(testOutputPath));

            // check are equivalent
            for (int i = 0; i < deconFile.Results.Count; i++)
            {
                var original = JsonConvert.SerializeObject(deconFile.Results[i]);
                var written = JsonConvert.SerializeObject(writtenDeconFile.Results[i]);
                Assert.That(original, Is.EqualTo(written));
            }

            // test writer still works without specifying extensions
            File.Delete(testOutputPath);
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "featureOut");
            deconFile.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }
    }
}
