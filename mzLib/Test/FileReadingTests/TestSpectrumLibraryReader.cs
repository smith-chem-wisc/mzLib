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
                @"FileReadingTests\ExternalFileTypes\mspStyleSpectrumLibary");
            var testOutputPath = Path.Combine(directoryPath, "ms1TsvOut_ms1.tsv");

            // load in file
            LibrarySpectrumFile mspSpectrumFile = FileReader.ReadFile<LibrarySpectrumFile>(testFilePath);

            // write and reread file
            mspSpectrumFile.WriteResults(testOutputPath);
            var writtenDeconFile = FileReader.ReadFile<FlashDeconvMs1TsvFile>(testOutputPath);
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
    }
}
