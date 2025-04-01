using Newtonsoft.Json;
using NUnit.Framework;
using Readers;
using Readers.SpectrumLibraries;
using System.Diagnostics.CodeAnalysis;
using System.IO;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSpectrumLibraryReader
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
        public static void TestMspStyleSpectrumReadWrite()
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\mspStyleSpectrumLibary.msp");
            var testOutputPath = Path.Combine(directoryPath, "mspStyleSpectrumLibary_out.msp");

            // load in file
            MspSpectrumFile mspSpectrumFile = FileReader.ReadFile<MspSpectrumFile>(testFilePath);
            Assert.That(mspSpectrumFile.Results.Count, Is.EqualTo(10));

            // write and reread file
            mspSpectrumFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<MspSpectrumFile>(testOutputPath);
            Assert.That(File.Exists(testOutputPath));

            // check are equivalent
            for (int i = 0; i < mspSpectrumFile.Results.Count; i++)
            {
                var original = JsonConvert.SerializeObject(mspSpectrumFile.Results[i]);
                var written = JsonConvert.SerializeObject(writtenDeconFile.Results[i]);
                Assert.That(original, Is.EqualTo(written));
            }

            // test writer still works without specifying extensions
            File.Delete(testOutputPath);
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "ms1TsvOut");
            mspSpectrumFile.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }

        [Test]
        public static void TestMsFraggerSpeclibReadWrite()
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\MsFraggerSpecLibExample.speclib");
            var testOutputPath = Path.Combine(directoryPath, "MsFraggerSpecLibExample_out.speclib");

            // load in file
            MsFraggerSpeclibFile msFraggerSpeclibFile = FileReader.ReadFile<MsFraggerSpeclibFile>(testFilePath);
            Assert.That(msFraggerSpeclibFile.OriginalRecords.Count, Is.EqualTo(344));
            Assert.That(msFraggerSpeclibFile.Results.Count, Is.EqualTo(29));

            // write and reread file
            msFraggerSpeclibFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<MsFraggerSpeclibFile>(testOutputPath);
            Assert.That(File.Exists(testOutputPath));

            // check are equivalent
            for (int i = 0; i < msFraggerSpeclibFile.Results.Count; i++)
            {
                var original = JsonConvert.SerializeObject(msFraggerSpeclibFile.Results[i]);
                var written = JsonConvert.SerializeObject(writtenDeconFile.Results[i]);
                Assert.That(original, Is.EqualTo(written));
            }

            // test writer still works without specifying extensions
            File.Delete(testOutputPath);
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "MsFraggerSpecLibExample_out");
            msFraggerSpeclibFile.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }


        [Test]
        public static void TestMspEntryWithNterminalMod()
        {
            //TOTO
            // [common bio mod]PEPTIDE read does not work
        }
    }
}
