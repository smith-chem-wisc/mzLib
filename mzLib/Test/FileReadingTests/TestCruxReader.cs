using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Newtonsoft.Json;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestCruxReader
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
        [TestCase(@"FileReadingTests\ExternalFileTypes\crux.txt", 14)]
        public void TestCruxResultsLoadsAndCountCorrect(string path, int recordCount)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            CruxResultFile file = new CruxResultFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(recordCount));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\crux.txt", 14)]
        public static void TestCruxResultsFromGenericReader(string path, int recordCount)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            var constructedFile = new CruxResultFile(filePath);
            var genericFile = FileReader.ReadFile<CruxResultFile>(filePath);

            Assert.That(genericFile.Count(), Is.EqualTo(recordCount));
            Assert.That(genericFile.Count(), Is.EqualTo(constructedFile.Count()));
            Assert.That(genericFile.FilePath, Is.EqualTo(constructedFile.FilePath));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\crux.txt")]
        public void TestCruxResultsFirstAndLastAreCorrect(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            var file = new CruxResultFile(filePath);

            var first = file.First();
            var last = file.Last();

            Assert.That(first.FilePath, Is.EqualTo(@"/hdd/data/PXD005590/B02_21_161103_D4_HCD_OT_4ul.raw.mzXML"));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(14674));
            Assert.That(first.Charge, Is.EqualTo(3));
            Assert.That(first.RetentionTime, Is.EqualTo(2747.6599));
            Assert.That(first.PrecursorMz, Is.EqualTo(1075.1815));
            Assert.That(first.NeutralMass, Is.EqualTo(3222.5227));
            Assert.That(first.PeptideMass, Is.EqualTo(3222.5222));
            Assert.That(first.DeltaCn, Is.EqualTo(0.84335566));
            Assert.That(first.XCorrScore, Is.EqualTo(6.4364114));
            Assert.That(first.XCorrRank, Is.EqualTo(1));
            Assert.That(first.TailorScore, Is.EqualTo(1.9659604));
            Assert.That(first.TdcQValue, Is.EqualTo(0.0000038850189).Within(1E-6));
            Assert.That(first.BAndYIonsMatched, Is.EqualTo(51));
            Assert.That(first.BAndYIonsTotal, Is.EqualTo(116));
            Assert.That(first.BAndYIonsFraction, Is.EqualTo(0.43965518));
            Assert.That(first.BAndYIonRepeatMatch, Is.EqualTo(0));
            Assert.That(first.BaseSequence, Is.EqualTo("RPQYSNPPVQGEVMEGADNQGAGEQGRPVR"));
            Assert.That(first.FullSequence, Is.EqualTo("RPQYSNPPVQGEVMEGADNQGAGEQGRPVR"));
            Assert.That(first.ProteinId, Is.EqualTo("sp|P67809|YBOX1_HUMAN(205)"));
            Assert.That(first.FlankingAa, Is.EqualTo("RQ"));
            Assert.That(first.FileNameWithoutExtension, Is.EqualTo("B02_21_161103_D4_HCD_OT_4ul.raw"));
            Assert.That(first.Accession, Is.EqualTo("P67809"));

            Assert.That(last.FilePath, Is.EqualTo(@"/hdd/data/PXD005590/B02_20_161103_E4_HCD_OT_4ul.raw.mzXML"));
            Assert.That(last.OneBasedScanNumber, Is.EqualTo(19752));
            Assert.That(last.Charge, Is.EqualTo(3));
            Assert.That(last.RetentionTime, Is.EqualTo(3220.75));
            Assert.That(last.PrecursorMz, Is.EqualTo(827.3577));
            Assert.That(last.NeutralMass, Is.EqualTo(2479.0515));
            Assert.That(last.PeptideMass, Is.EqualTo(2479.0457));
            Assert.That(last.DeltaCn, Is.EqualTo(0.85362118));
            Assert.That(last.XCorrScore, Is.EqualTo(5.682076));
            Assert.That(last.XCorrRank, Is.EqualTo(1));
            Assert.That(last.TailorScore, Is.EqualTo(1.8917845));
            Assert.That(last.TdcQValue, Is.EqualTo(0.00000388501896).Within(1E-6));
            Assert.That(last.BAndYIonsMatched, Is.EqualTo(32));
            Assert.That(last.BAndYIonsTotal, Is.EqualTo(92));
            Assert.That(last.BAndYIonsFraction, Is.EqualTo(0.34782609));
            Assert.That(last.BAndYIonRepeatMatch, Is.EqualTo(0));
            Assert.That(last.BaseSequence, Is.EqualTo("QDHPSSMGVYGQESGGFSGPGENR"));
            Assert.That(last.FullSequence, Is.EqualTo("QDHPSSMGVYGQESGGFSGPGENR"));
            Assert.That(last.ProteinId, Is.EqualTo("sp|Q01844|EWS_HUMAN(269)"));
            Assert.That(last.FlankingAa, Is.EqualTo("RS"));
            Assert.That(last.FileNameWithoutExtension, Is.EqualTo("B02_20_161103_E4_HCD_OT_4ul.raw"));
            Assert.That(last.Accession, Is.EqualTo("Q01844"));
        }


        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\crux.txt")]
        public void TestCruxResultsWriteResults(string path)
        {
            // load in original
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            var original = new CruxResultFile(filePath);

            // write out original
            var outputPath = Path.Combine(directoryPath, "cruxResults.csv");
            original.WriteResults(outputPath);
            Assert.That(File.Exists(outputPath));

            // read in new original
            var written = new CruxResultFile(outputPath);
            Assert.That(written.Count(), Is.EqualTo(original.Count()));

            // check are equivalent
            for (int i = 0; i < original.Count(); i++)
            {
                var oldRecord = JsonConvert.SerializeObject(original.Results[i]);
                var newRecord = JsonConvert.SerializeObject(written.Results[i]);
                Assert.That(oldRecord, Is.EqualTo(newRecord));
            }

            // test writer still works without specifying extensions
            var outputPathWithoutExtension = Path.Combine(directoryPath, "cruxResults");
            original.WriteResults(outputPathWithoutExtension);
            Assert.That(File.Exists(outputPathWithoutExtension + ".csv"));

            var writtenWithoutExtension = new CruxResultFile(outputPathWithoutExtension + ".csv");
            Assert.That(writtenWithoutExtension.Count(), Is.EqualTo(original.Count()));
        }
    }
}
