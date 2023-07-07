using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestTopFDFiles
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
        [TestCase(@"FileReadingTests\ExternalFileTypes\TopFDmzrt_v1.6.2.mzrt.csv", 4)]
        public void TestTopFDMzRTLoadsAndCountCorrect(string path, int recordCount)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            TopFDMzrtFile file = new TopFDMzrtFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(recordCount));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\TopFDmzrt_v1.6.2.mzrt.csv")]
        public void TestTopFDMzRTFirstAndLastAreCorrect(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            var file = new TopFDMzrtFile(filePath);

            var first = file.First();
            var last = file.Last();

            Assert.That(first.Id, Is.EqualTo(0));
            Assert.That(first.FractionId, Is.EqualTo(0));
            Assert.That(first.EnvelopeNumber, Is.EqualTo(5));
            Assert.That(first.Mass, Is.EqualTo(10835.9));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(1548.99));
            Assert.That(first.Charge, Is.EqualTo(7));
            Assert.That(first.Intensity, Is.EqualTo(8.57517e+06));
            Assert.That(first.MzMin, Is.EqualTo(1549.03));
            Assert.That(first.MzMax, Is.EqualTo(1551.23));
            Assert.That(first.RetentionTimeBegin, Is.EqualTo(39.6613));
            Assert.That(first.RetentionTimeEnd, Is.EqualTo(39.9085));
            Assert.That(first.ColorHexcode, Is.EqualTo("#FF0000"));
            Assert.That(first.Opacity, Is.EqualTo(0.1));
            Assert.That(first.PromexScore, Is.EqualTo(40.1092));

            Assert.That(last.Id, Is.EqualTo(20553));
            Assert.That(last.FractionId, Is.EqualTo(0));
            Assert.That(last.EnvelopeNumber, Is.EqualTo(1));
            Assert.That(last.Mass, Is.EqualTo(1841.07));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(1842.08));
            Assert.That(last.Charge, Is.EqualTo(1));
            Assert.That(last.Intensity, Is.EqualTo(3279.65));
            Assert.That(last.MzMin, Is.EqualTo(1841.98));
            Assert.That(last.MzMax, Is.EqualTo(1846.18));
            Assert.That(last.RetentionTimeBegin, Is.EqualTo(52.0425));
            Assert.That(last.RetentionTimeEnd, Is.EqualTo(52.0425));
            Assert.That(last.ColorHexcode, Is.EqualTo("#FF0000"));
            Assert.That(last.Opacity, Is.EqualTo(0.1));
            Assert.That(last.PromexScore, Is.EqualTo(-1000));
        }

        [Test]
        public static void TestTopFdMzRteReadWrite()
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\TopFDmzrt_v1.6.2.mzrt.csv");
            var testOutputPath = Path.Combine(directoryPath, "topfdmzrt.mzrt.csv");

            // load in file
            TopFDMzrtFile deconFile = FileReader.ReadFile<TopFDMzrtFile>(testFilePath);

            // write and reread file
            deconFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<TopFDMzrtFile>(testOutputPath);
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
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "topfdmzrt");
            deconFile.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }
    }
}
