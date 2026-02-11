using NUnit.Framework;
using Readers.ExternalResults.ResultFiles;
using Readers;
using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestCasanovoMzTab
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
        public void TestCasanovoMzTabLoadsAndCountCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\Casanovo_5.0.0.mztab");
            var file = new CasanovoMzTabFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(5));
            Assert.That(file.CanRead(filePath));

            // Header captured and ms_run parsed
            Assert.That(file.HeaderLines.Any(l => l.StartsWith("MTD\tmzTab-version")));
            Assert.That(file.MsRunDictionary.ContainsKey("ms_run[1]"));
            Assert.That(Path.GetFileName(file.MsRunDictionary["ms_run[1]"]).EndsWith(".mgf", StringComparison.InvariantCultureIgnoreCase));
        }

        [Test]
        public void TestCasanovoMzTabFirstAndLastCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\Casanovo_5.0.0.mztab");
            var file = new CasanovoMzTabFile(filePath);
            var first = file.First();
            var last = file.Last();

            // First record
            Assert.That(first.Id, Is.EqualTo(1));
            Assert.That(first.ModifiedSequence, Is.EqualTo("AGAHLQGGAK"));
            Assert.That(first.BaseSequence, Is.EqualTo("AGAHLQGGAK"));
            Assert.That(first.FullSequence, Is.EqualTo("AGAHLQGGAK"));
            Assert.That(first.Charge, Is.EqualTo(2));
            Assert.That(first.FileNameWithoutExtension, Is.EqualTo("GODWIN_022123_TISSUE_DDA_SF12"));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(1));
            Assert.That(first.AminoAcidScores.Length, Is.EqualTo(10));
            Assert.That(first.SpectraRef, Is.EqualTo("ms_run[1]:index=0"));
            Assert.That(first.RetentionTime, Is.Null);
            Assert.That(first.Start, Is.Null);
            Assert.That(first.End, Is.Null);

            // Last record
            Assert.That(last.Id, Is.EqualTo(5));
            Assert.That(last.ModifiedSequence, Is.EqualTo("RQEFEM[Oxidation]K"));
            Assert.That(last.BaseSequence, Is.EqualTo("RQEFEMK"));
            Assert.That(last.FullSequence, Is.EqualTo("RQEFEM[Common Variable:Oxidation on M]K"));
            Assert.That(last.Charge, Is.EqualTo(2));
            Assert.That(last.FileNameWithoutExtension, Is.EqualTo("GODWIN_022123_TISSUE_DDA_SF12"));
            Assert.That(last.OneBasedScanNumber, Is.EqualTo(7));
            Assert.That(last.AminoAcidScores.Length, Is.EqualTo(7));
            Assert.That(last.SpectraRef, Is.EqualTo("ms_run[1]:index=6"));
            Assert.That(last.RetentionTime, Is.Null);
            Assert.That(last.Start, Is.Null);
            Assert.That(last.End, Is.Null);
        }

        [Test]
        public void TestCasanovoMzTabReadWriteRoundTrip()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\Casanovo_5.0.0.mztab");
            string outPath = Path.Combine(directoryPath, "testCasanovoMzTab");

            var file = new CasanovoMzTabFile(filePath);
            file.WriteResults(outPath);

            var outFile = new CasanovoMzTabFile(outPath + ".mztab");
            Assert.That(outFile.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file.Count(); i++)
            {
                var original = JsonSerializer.Serialize(
                    file.ElementAt(i),
                    new JsonSerializerOptions
                    {
                        NumberHandling = System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals
                    });
                var written = JsonSerializer.Serialize(
                    outFile.ElementAt(i),
                    new JsonSerializerOptions
                    {
                        NumberHandling = System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals
                    });
                Assert.That(original, Is.EqualTo(written));
            }

            // Ensure file written contains PSH/PSM prefixes
            var lines = File.ReadAllLines(outPath + ".mztab");
            Assert.That(lines.Any(l => l.StartsWith("PSH\t")));
            Assert.That(lines.Any(l => l.StartsWith("PSM\t")));
        }

        [Test]
        public void Test_CommaDelimitedToDoubleArrayTypeConverter_ReadWrite()
        {
            var conv = new CommaDelimitedToDoubleArrayTypeConverter();
            var input = "0.1,0.2,,0.3"; // includes an empty value
            var arr = (double[])conv.ConvertFromString(input, row: null, memberMapData: null);
            Assert.That(arr, Is.EqualTo(new[] { 0.1, 0.2, 0.3 }).Within(1e-9));

            var outStr = conv.ConvertToString(new[] { 1.0, 2.5, 3.75 }, row: null, memberMapData: null);
            Assert.That(outStr, Is.EqualTo("1,2.5,3.75"));
        }

        [Test]
        public void Test_ModificationConverter_GetClosestMod_And_ParseFullSequence()
        {
            // GetClosestMod picks Oxidation on M
            var oxM = ModificationConverter.GetClosestMod("Oxidation", 'M', null);
            Assert.That(oxM.IdWithMotif, Does.Contain("Oxidation"));
            Assert.That(oxM.IdWithMotif, Does.Contain(" on M"));
            Assert.That(oxM.MonoisotopicMass ?? 15.9949, Is.EqualTo(15.9949).Within(0.01));

            // Parse MetaMorpheus-style full sequence
            var fullSeq = "KVHGSLARAGKVRGQTPKVAKQEKKKKKTGRAKRRM[Common Variable:Oxidation on M]QYNRRFVNVVPTFGKKKGPNANS";
            var dict = ModificationConverter.GetModificationDictionaryFromFullSequence(fullSeq);
            Assert.That(dict.Count, Is.EqualTo(1));
            var only = dict.Single();
            Assert.That(only.Value.IdWithMotif, Does.Contain("Oxidation"));
            Assert.That(only.Value.MonoisotopicMass ?? 15.9949, Is.EqualTo(15.9949).Within(0.01));

            // Parse non-MM style short sequence (fallback path)
            var dict2 = ModificationConverter.GetModificationDictionaryFromFullSequence("AM[Oxidation]A");
            Assert.That(dict2.Count, Is.EqualTo(1));
            Assert.That(dict2.Keys.Single(), Is.EqualTo(3)); // M at position 3
            Assert.That(dict2.Values.Single().IdWithMotif, Does.Contain("Oxidation"));
        }
    }
}
