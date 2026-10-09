using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;

namespace Test.FileReadingTests.InternalFileReading
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestTranscriptGroupFromTsv
    {
        private static string FixturePath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\ExternalFileTypes\MetaMorpheus_RNA_AllQuantifiedTranscriptGroups.tsv");

        [Test]
        public void ReadsEveryRowThroughBothEntryPoints()
        {
            var direct = new TranscriptGroupFromTsvFile(FixturePath);
            var factory = FileReader.ReadFile<TranscriptGroupFromTsvFile>(FixturePath);

            Assert.That(direct.Count(), Is.EqualTo(3));
            Assert.That(factory.Count(), Is.EqualTo(3));
            Assert.That(direct.CanRead(FixturePath));
            Assert.That(direct.FileType, Is.EqualTo(SupportedFileType.MetaMorpheusQuantifiedTranscriptGroups));
            Assert.That(direct.Software, Is.EqualTo(Software.MetaMorpheus));
        }

        [Test]
        public void ReadsMetaMorpheusColumnsAndVerbatimSampleLabels()
        {
            var row = new TranscriptGroupFromTsvFile(FixturePath).Single(r => r.ProteinGroupName == "FLuc");

            Assert.Multiple(() =>
            {
                Assert.That(row.Organism, Is.EqualTo("standard"));
                Assert.That(row.FullName, Is.EqualTo("FLuc"));
                Assert.That(row.NumberOfMembers, Is.EqualTo(1));
                Assert.That(row.NumberOfSequences, Is.EqualTo(231));
                Assert.That(row.QValue, Is.EqualTo(0));
                Assert.That(row.DecoyContaminantTarget, Is.EqualTo("T"));
                Assert.That(row.SampleGroups, Has.Count.EqualTo(37));
            });

            var sample = row.SampleGroups["1:1_1"];
            Assert.That(sample.SpectralCount, Is.EqualTo(355));
            Assert.That(sample.Intensity, Is.EqualTo(89077470.71004736));
            Assert.That(sample.CountOccupancyText, Is.Null);
            Assert.That(sample.CountOccupancy.Sites, Is.Empty);
        }

        [Test]
        public void ParsesRnaModificationOccupancyCells()
        {
            var group = new TranscriptGroupFromTsvFile(FixturePath).Single(r => r.ProteinGroupName == "20mer2");
            var sample = group.SampleGroups["1:1_1"];
            var countSites = sample.CountOccupancy.Sites.ToArray();

            Assert.That(countSites.Select(s => s.Position), Is.EqualTo(new[] { 2, 8, 13, 18 }));
            Assert.That(countSites[0].ModificationIdWithMotif, Is.EqualTo("2'-O-methyluridine on U"));
            Assert.That((countSites[0].Numerator, countSites[0].Denominator, countSites[0].Fraction),
                Is.EqualTo((22.0, 22.0, 1.0)));

            var intensitySite = sample.IntensityOccupancy.Sites.First();
            Assert.That(intensitySite.ModificationIdWithMotif, Is.EqualTo("2'-O-methyluridine on U"));
            Assert.That(intensitySite.Fraction, Is.EqualTo(1.0));
            Assert.That(intensitySite.Numerator, Is.EqualTo(4.279e6));
            Assert.That(intensitySite.Denominator, Is.EqualTo(4.279e6));
        }

        [Test]
        public void ReadsMultipleMembersDecoysContaminantsAndZeroPositionOccupancy()
        {
            string directory = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(TestTranscriptGroupFromTsv) + "Synthetic");
            Directory.CreateDirectory(directory);
            string path = Path.Combine(directory, "Synthetic_TranscriptGroups.tsv");
            File.WriteAllLines(path,
            [
                "Protein Accession\tProtein Decoy/Contaminant/Target\tProtein QValue\tCountOccupancy_1:1_1",
                "T1|T2\tT\t0\tpos0[5'-phosphate on X,info:fraction=1.00(1/1)]|pos2[2'-O-methyluridine on U,info:fraction=1.00(2/2)]",
                "DECOY_T3\tD\t0.01\t",
                "CONTAM_T4\tC\t0.02\t"
            ]);

            try
            {
                var rows = new TranscriptGroupFromTsvFile(path).ToArray();
                Assert.That(rows[0].Accessions, Is.EqualTo(new[] { "T1", "T2" }));
                Assert.That(rows[0].SampleGroups["1:1_1"].CountOccupancy.Entities, Has.Count.EqualTo(2));
                var nTerminus = rows[0].SampleGroups["1:1_1"].CountOccupancy.Sites.First();
                Assert.That(nTerminus.IsNTerminus);
                Assert.That(nTerminus.ModificationIdWithMotif, Is.EqualTo("5'-phosphate on X"));
                Assert.That(rows[1].IsDecoy && !rows[1].IsContaminant);
                Assert.That(rows[2].IsContaminant && !rows[2].IsDecoy);
            }
            finally
            {
                Directory.Delete(directory, true);
            }
        }

        [TestCase("AllTranscriptGroups.tsv")]
        [TestCase("Sample1_TranscriptGroups.tsv")]
        [TestCase("ALLQUANTIFIEDTRANSCRIPTGROUPS.TSV")]
        public void RecognizesTranscriptGroupFileNames(string fileName)
        {
            string directory = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(TestTranscriptGroupFromTsv));
            Directory.CreateDirectory(directory);
            string path = Path.Combine(directory, fileName);
            File.WriteAllText(path, "Protein Accession\tProtein Decoy/Contaminant/Target\tProtein QValue\nT1\tT\t0.004\n");

            try
            {
                Assert.That(path.ParseFileType(), Is.EqualTo(SupportedFileType.MetaMorpheusQuantifiedTranscriptGroups));
                var result = (TranscriptGroupFromTsvFile)FileReader.ReadResultFile(path);
                Assert.That(result.Single().QValue, Is.EqualTo(0.004));
                Assert.That(FileReader.ReadFile<TranscriptGroupFromTsvFile>(path).Single().ProteinGroupName, Is.EqualTo("T1"));
            }
            finally
            {
                Directory.Delete(directory, true);
            }
        }

        [Test]
        public void ReportsMalformedInputWithPathAndRefusesWriting()
        {
            string directory = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(TestTranscriptGroupFromTsv) + "Errors");
            Directory.CreateDirectory(directory);
            string path = Path.Combine(directory, "Broken_TranscriptGroups.tsv");
            File.WriteAllText(path, "Protein Accession\tProtein Decoy/Contaminant/Target\tProtein QValue\nT1\tT\tnot-a-number\n");

            try
            {
                var exception = Assert.Throws<MzLibException>(() => new TranscriptGroupFromTsvFile(path).LoadResults());
                Assert.That(exception!.Message, Does.Contain(path));
                Assert.Throws<NotSupportedException>(() => new TranscriptGroupFromTsvFile(FixturePath).WriteResults("x"));
            }
            finally
            {
                Directory.Delete(directory, true);
            }
        }
    }
}
