using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using Readers.ExternalResults.IndividualResultRecords;
using Readers.ExternalResults.ResultFiles;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// The test file is a Pytheas match output printed from the GUI_version_data_paper data
    /// (Lumos_Orbi.mgf searched against a 16S E. coli digest). It holds 4744 match lines spread
    /// over 822 PRECURSOR_ION groups, covering light and heavy isotopologues, ambiguous base
    /// calls (UXG), modified bases (CC[mG]CG), multi-location matches, and decoy rows.
    /// Distinct from the newer tab-delimited readers above, the reader here has to parse a
    /// space-delimited custom format line by line.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestPytheasResultFile
    {
        private const string PytheasResultPath = @"FileReadingTests\ExternalFileTypes\match_output_Lumos_Orbi.txt";

        private static string TestFilePath =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, PytheasResultPath);

        private string _outputDirectory;

        [OneTimeSetUp]
        public void SetUp()
        {
            _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\PytheasReadWriteTests");
            Directory.CreateDirectory(_outputDirectory);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(_outputDirectory, true);
        }

        [Test]
        public void TestPytheasResultLoadsAndCountCorrect()
        {
            PytheasResultFile file = new PytheasResultFile(TestFilePath);

            Assert.That(file.Count(), Is.EqualTo(4744));
            Assert.That(file.CanRead(TestFilePath));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.PytheasResult));
            Assert.That(file.Software, Is.EqualTo(Software.Pytheas));
        }

        [Test]
        public void TestPytheasResultFirstRecordCorrect()
        {
            PytheasResult first = new PytheasResultFile(TestFilePath).First();

            Assert.That(first.PrecursorIon, Is.EqualTo("485.563080"));
            Assert.That(first.MeasuredMz, Is.EqualTo(485.563080));
            Assert.That(first.RetentionTime, Is.EqualTo(14.673));
            Assert.That(first.TheoreticalMz, Is.EqualTo(485.56246141));
            Assert.That(first.OffsetPpm, Is.EqualTo(1.3));
            Assert.That(first.SpScore, Is.EqualTo(0.239));
            Assert.That(first.DspScore, Is.EqualTo(0.0));
            Assert.That(first.Rank, Is.EqualTo(1));
            Assert.That(first.Ms2MatchCount, Is.EqualTo(5));
            Assert.That(first.Isotope, Is.EqualTo("light"));
            Assert.That(first.Length, Is.EqualTo(3));
            Assert.That(first.Charge, Is.EqualTo(-2));
            Assert.That(first.Sequence, Is.EqualTo("CCG"));
            Assert.That(first.SequenceModification, Is.EqualTo("-"));
            Assert.That(first.FivePrimeEnd, Is.EqualTo("OH"));
            Assert.That(first.ThreePrimeEnd, Is.EqualTo("P"));
            Assert.That(first.MoleculeLocation,
                Is.EqualTo("16S.ecoli,400,402;16S.ecoli,525,527;16S.ecoli,1140,1142;16S.ecoli,896,898;16S.ecoli,810,812"));
            Assert.That(first.Ms2Matches, Does.StartWith("892.167114(0.8ppm)[680]:892.166417[M-P(-1)]"));
            Assert.That(first.Score, Is.EqualTo(0.239));
            Assert.That(first.DetailedScore, Does.StartWith("SCORE=0.239"));
        }

        [Test]
        public void TestPytheasResultLastRecordCorrect()
        {
            PytheasResult last = new PytheasResultFile(TestFilePath).Last();

            Assert.That(last.PrecursorIon, Is.EqualTo("1677.156250"));
            Assert.That(last.MeasuredMz, Is.EqualTo(1677.156250));
            Assert.That(last.RetentionTime, Is.EqualTo(17.395));
            Assert.That(last.TheoreticalMz, Is.EqualTo(1677.180824070));
            Assert.That(last.OffsetPpm, Is.EqualTo(-14.7));
            Assert.That(last.SpScore, Is.EqualTo(0.043));
            Assert.That(last.DspScore, Is.EqualTo(0.73));
            Assert.That(last.Rank, Is.EqualTo(6));
            Assert.That(last.Ms2MatchCount, Is.EqualTo(8));
            Assert.That(last.Isotope, Is.EqualTo("heavy"));
            Assert.That(last.Length, Is.EqualTo(5));
            Assert.That(last.Charge, Is.EqualTo(-1));
            Assert.That(last.Sequence, Is.EqualTo("AAAXG"));
            Assert.That(last.SequenceModification, Is.EqualTo("AAACG|AAAUG"));
            Assert.That(last.FivePrimeEnd, Is.EqualTo("OH"));
            Assert.That(last.ThreePrimeEnd, Is.EqualTo("P"));
            Assert.That(last.MoleculeLocation,
                Is.EqualTo("16S.ecoli,160,164;16S.ecoli,694,698;16S.ecoli,1080,1084"));
            Assert.That(last.Ms2Matches, Does.StartWith("1309.125854(-9.0ppm)[86]:1309.137650[c4(-1)]"));
            Assert.That(last.Score, Is.EqualTo(0.043));
            Assert.That(last.DetailedScore, Does.StartWith("SCORE=0.043"));
        }

        [Test]
        public void TestPytheasResultIsDetectedByHeaderRegardlessOfFileName()
        {
            string renamedPath = Path.Combine(_outputDirectory, "some_unrelated_name.tsv");
            File.Copy(TestFilePath, renamedPath, true);

            Assert.That(renamedPath.ParseFileType(), Is.EqualTo(SupportedFileType.PytheasResult));

            IResultFile read = FileReader.ReadResultFile(renamedPath);
            Assert.That(read, Is.TypeOf<PytheasResultFile>());
            Assert.That(read.FileType, Is.EqualTo(SupportedFileType.PytheasResult));
            Assert.That(((PytheasResultFile)read).Count(), Is.EqualTo(4744));
        }

        [Test]
        public void TestPytheasResultIsDetectedByHeaderInTxtExtension()
        {
            string renamedPath = Path.Combine(_outputDirectory, "renamed.txt");
            File.Copy(TestFilePath, renamedPath, true);

            Assert.That(renamedPath.ParseFileType(), Is.EqualTo(SupportedFileType.PytheasResult));
        }

        [Test]
        public void TestPytheasResultReadWriteRoundTrip()
        {
            PytheasResultFile original = new PytheasResultFile(TestFilePath);
            string outputPath = Path.Combine(_outputDirectory, "roundTrip.txt");

            original.WriteResults(outputPath);
            Assert.That(File.Exists(outputPath));

            PytheasResultFile rewritten = new PytheasResultFile(outputPath);
            Assert.That(rewritten.Count(), Is.EqualTo(original.Count()));

            foreach (var (before, after) in original.Zip(rewritten))
            {
                Assert.That(after.PrecursorIon, Is.EqualTo(before.PrecursorIon));
                Assert.That(after.RawLine, Is.EqualTo(before.RawLine));
                Assert.That(after.MeasuredMz, Is.EqualTo(before.MeasuredMz));
                Assert.That(after.SpScore, Is.EqualTo(before.SpScore));
                Assert.That(after.Sequence, Is.EqualTo(before.Sequence));
                Assert.That(after.DetailedScore, Is.EqualTo(before.DetailedScore));
            }

            Assert.That(outputPath.ParseFileType(), Is.EqualTo(SupportedFileType.PytheasResult));
        }

        [Test]
        public void TestPytheasResultHeaderLinesRoundTrip()
        {
            PytheasResultFile original = new PytheasResultFile(TestFilePath);
            _ = original.Count();
            Assert.That(original.HeaderLines.Count, Is.EqualTo(21));

            string outputPath = Path.Combine(_outputDirectory, "headerRoundTrip.txt");
            original.WriteResults(outputPath);

            PytheasResultFile rewritten = new PytheasResultFile(outputPath);
            _ = rewritten.Count();
            Assert.That(rewritten.HeaderLines, Is.EqualTo(original.HeaderLines));
        }

        /// <summary>
        /// The measured m/z sometimes carries a trailing '*' (phosphor/ambiguity flag in Pytheas'
        /// notation). The numeric part must parse, and the raw line must keep the flag so it
        /// round-trips.
        /// </summary>
        [Test]
        public void TestStarredMeasuredMzParses()
        {
            PytheasResult starred = new PytheasResultFile(TestFilePath)
                .First(r => r.RawLine.StartsWith("491.536102*"));

            Assert.That(starred.MeasuredMz, Is.EqualTo(491.536102));
            Assert.That(starred.RawLine, Does.StartWith("491.536102*"));
            Assert.That(starred.SpScore, Is.EqualTo(0.223));
        }

        [Test]
        public void TestMzLibExceptionOnTooFewColumns()
        {
            Assert.Throws<MzLibException>(() => PytheasResult.Parse("too few tokens", "100.0"));
        }
    }
}