using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Newtonsoft.Json;
using NUnit.Framework;
using pepXML.Generated;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestFlashDeconFiles
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
        [TestCase(@"FileReadingTests\ExternalFileTypes\FlashDeconvMs1Tsv_OpenMs3.0.0_ms1.tsv", 8)]
        public void TestMs1TsvFileLoadsAndCountCorrect(string path, int count)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            FlashDeconvMs1TsvFile file = new FlashDeconvMs1TsvFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FlashDeconvTsv_OpenMs3.0.0.tsv", 5)]
        public void TestFlashDeconTsvFileLoadsAndCountCorrect(string path, int count)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            FlashDeconvTsvFile file = new FlashDeconvTsvFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        public void TestFlashDeconMs1TsvFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FlashDeconvMs1Tsv_OpenMs3.0.0_ms1.tsv");
            var ms1TsvFile = new FlashDeconvMs1TsvFile(filePath);
            var first = ms1TsvFile.First();
            var last = ms1TsvFile.Last();


            Assert.That(first.Index, Is.EqualTo(1));
            Assert.That(first.FileName, Is.EqualTo("D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(first.ZeroBasedScanNumber, Is.EqualTo(-1));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(0));
            Assert.That(first.Decoy, Is.EqualTo(0));
            Assert.That(first.RetentionTime, Is.EqualTo(1.236057));
            Assert.That(first.MassCountInSpec, Is.EqualTo(58));
            Assert.That(first.AverageMass, Is.EqualTo(609.530436));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(609.177191));
            Assert.That(first.SumIntensity, Is.EqualTo(2069.13));
            Assert.That(first.MinCharge, Is.EqualTo(1));
            Assert.That(first.MaxCharge, Is.EqualTo(1));
            Assert.That(first.PeakCount, Is.EqualTo(2));
            Assert.That(first.IsotopeCosine, Is.EqualTo(0.979821));
            Assert.That(first.ChargeScore, Is.EqualTo(1));
            Assert.That(first.MassSNR, Is.EqualTo(5.68748));
            Assert.That(first.ChargeSNR, Is.EqualTo(5.68748));
            Assert.That(first.RepresentativeCharge, Is.EqualTo(1));
            Assert.That(first.RepresentativeMzMin, Is.EqualTo(610.185324));
            Assert.That(first.RepresentativeMzMax, Is.EqualTo(611.185295));
            Assert.That(first.QScore, Is.EqualTo(0.622148));
            Assert.That(first.Qvalue, Is.EqualTo(0.0333624));
            Assert.That(first.QvalueWithIsotopeDecoyOnly, Is.EqualTo(0.0211648));
            Assert.That(first.QvalueWithNoiseDecoyOnly, Is.EqualTo(0.00020122));
            Assert.That(first.QvalueWithChargeDecoyOnly, Is.EqualTo(0.0119964));
            Assert.That(first.PerChargeIntensity, Is.EqualTo(new List<double> { 2067.13 }));
            Assert.That(first.PerIsotopeIntensity, Is.EqualTo(new List<double> { 1330.02, 737.107 }));

            Assert.That(last.Index, Is.EqualTo(58));
            Assert.That(last.FileName, Is.EqualTo("D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(last.ZeroBasedScanNumber, Is.EqualTo(-1));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(0));
            Assert.That(last.Decoy, Is.EqualTo(2));
            Assert.That(last.RetentionTime, Is.EqualTo(5392.704300));
            Assert.That(last.MassCountInSpec, Is.EqualTo(58));
            Assert.That(last.AverageMass, Is.EqualTo(6271.571495));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(6267.708258));
            Assert.That(last.SumIntensity, Is.EqualTo(750.03));
            Assert.That(last.MinCharge, Is.EqualTo(4));
            Assert.That(last.MaxCharge, Is.EqualTo(6));
            Assert.That(last.PeakCount, Is.EqualTo(0));
            Assert.That(last.IsotopeCosine, Is.EqualTo(0.861777));
            Assert.That(last.ChargeScore, Is.EqualTo(1));
            Assert.That(last.MassSNR, Is.EqualTo(0.140764));
            Assert.That(last.ChargeSNR, Is.EqualTo(0.0997334));
            Assert.That(last.RepresentativeCharge, Is.EqualTo(6));
            Assert.That(last.RepresentativeMzMin, Is.EqualTo(1046.205181));
            Assert.That(last.RepresentativeMzMax, Is.EqualTo(1046.486308));
            Assert.That(last.QScore, Is.EqualTo(0.0662357));
            Assert.That(last.Qvalue, Is.EqualTo(1));
            Assert.That(last.QvalueWithIsotopeDecoyOnly, Is.EqualTo(1));
            Assert.That(last.QvalueWithNoiseDecoyOnly, Is.EqualTo(1));
            Assert.That(last.QvalueWithChargeDecoyOnly, Is.EqualTo(1));
            Assert.That(last.PerChargeIntensity, Is.EqualTo(new List<double> { 126.762, 219.947, 399.321 }));
            Assert.That(last.PerIsotopeIntensity, Is.EqualTo(new List<double> { 0, 126.762, 196.26, 203.061, 219.947 }));
        }


        [Test]
        public void TestFlashDeconvTsvFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FlashDeconvTsv_OpenMs3.0.0.tsv");
            var ms1TsvFile = new FlashDeconvTsvFile(filePath);
            var first = ms1TsvFile.First();
            var last = ms1TsvFile.Last();


            Assert.That(first.FeatureIndex, Is.EqualTo(0));
            Assert.That(first.FilePath,
                Is.EqualTo(
                    "D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(10835.850545));
            Assert.That(first.AverageMass, Is.EqualTo(10842.596816));
            Assert.That(first.MassCount, Is.EqualTo(7));
            Assert.That(first.RetentionTimeBegin, Is.EqualTo(2375.98));
            Assert.That(first.RetentionTimeEnd, Is.EqualTo(2398.21));
            Assert.That(first.RetentionTimeDuration, Is.EqualTo(22.233));
            Assert.That(first.RetentionTimeApex, Is.EqualTo(2390.8));
            Assert.That(first.IntensityOfAllPeaks, Is.EqualTo(1.21145e+10));
            Assert.That(first.IntensityOfMostAbundantPeak, Is.EqualTo(3.31923e+09));
            Assert.That(first.FeatureQuantity, Is.EqualTo(4.30496e+10));
            Assert.That(first.ChargeStateMin, Is.EqualTo(7));
            Assert.That(first.ChargeStateMax, Is.EqualTo(18));
            Assert.That(first.ChargeCount, Is.EqualTo(12));
            Assert.That(first.IsotopeCosineSimilarity, Is.EqualTo(0.997606));
            Assert.That(first.MaxQScore, Is.EqualTo(0.730994));
            Assert.That(first.IntensityPerChargeState,
                Is.EqualTo(new List<double>
                {
                    1.29036e+07, 4.50673e+07, 4.91774e+07, 1.29783e+08, 4.30799e+08, 1.42702e+09, 2.86165e+09,
                    3.56044e+09, 2.62395e+09, 8.8017e+08, 8.57345e+07, 7.77948e+06
                }));
            Assert.That(first.IntensityPerIsotope,
                Is.EqualTo(new List<double>
                {
                    2.56518e+07, 1.07828e+08, 3.55397e+08, 7.20686e+08, 1.29869e+09, 1.81351e+09,
                    1.93354e+09, 1.77206e+09, 1.4457e+09, 1.05851e+09, 7.10388e+08, 4.02286e+08,
                    2.45191e+08, 1.313e+08, 5.85071e+07, 3.52401e+07
                }
                ));


            Assert.That(last.FeatureIndex, Is.EqualTo(16531));
            Assert.That(last.FilePath,
                Is.EqualTo(
                    "D:/Averaging/Rep1CalibCentroidAverageCentroid/id_02-17-20_jurkat_td_rep1_fract2-calib-centroided-averaged-centroided.mzML"));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(618.236581));
            Assert.That(last.AverageMass, Is.EqualTo(618.603738));
            Assert.That(last.MassCount, Is.EqualTo(4));
            Assert.That(last.RetentionTimeBegin, Is.EqualTo(504.821));
            Assert.That(last.RetentionTimeEnd, Is.EqualTo(530.084));
            Assert.That(last.RetentionTimeDuration, Is.EqualTo(25.2625));
            Assert.That(last.RetentionTimeApex, Is.EqualTo(528.399));
            Assert.That(last.IntensityOfAllPeaks, Is.EqualTo(885.019));
            Assert.That(last.IntensityOfMostAbundantPeak, Is.EqualTo(225.488));
            Assert.That(last.FeatureQuantity, Is.EqualTo(5521.45));
            Assert.That(last.ChargeStateMin, Is.EqualTo(1));
            Assert.That(last.ChargeStateMax, Is.EqualTo(1));
            Assert.That(last.ChargeCount, Is.EqualTo(1));
            Assert.That(last.IsotopeCosineSimilarity, Is.EqualTo(0.919741));
            Assert.That(last.MaxQScore, Is.EqualTo(0.0723919));
            Assert.That(last.IntensityPerChargeState,
                Is.EqualTo(new List<double>
                {
                    877.019
                }));
            Assert.That(last.IntensityPerIsotope,
                Is.EqualTo(new List<double>
                {
                    464.488, 412.531
                }));
        }

        [Test]
        public static void TestFlashDeconMs1TsvReadWrite()
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FlashDeconvMs1Tsv_OpenMs3.0.0_ms1.tsv");
            var testOutputPath = Path.Combine(directoryPath, "ms1TsvOut_ms1.tsv");

            // load in file
            FlashDeconvMs1TsvFile deconFile = FileReader.ReadFile<FlashDeconvMs1TsvFile>(testFilePath);

            // write and reread file
            deconFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<FlashDeconvMs1TsvFile>(testOutputPath);
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
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "ms1TsvOut");
            deconFile.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }


        [Test]
        public static void TestFlashDeconTsvReadWrite()
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FlashDeconvTsv_OpenMs3.0.0.tsv");
            var testOutputPath = Path.Combine(directoryPath, "tsvOut.tsv");

            // load in file
            FlashDeconvTsvFile deconFile = FileReader.ReadFile<FlashDeconvTsvFile>(testFilePath);

            // write and reread file
            deconFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<FlashDeconvTsvFile>(testOutputPath);
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
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "tsvOut");
            deconFile.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }
    }
}
