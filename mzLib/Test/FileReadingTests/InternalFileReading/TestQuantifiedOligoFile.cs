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
    internal class TestQuantifiedOligoFile
    {
        private static string FixturePath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\ExternalFileTypes\MetaMorpheus_RNA_AllQuantifiedOligos.tsv");

        [Test]
        public void ReadsFixedAndPerSampleColumnsThroughBothEntryPoints()
        {
            var direct = new QuantifiedOligoFile(FixturePath);
            var factory = FileReader.ReadFile<QuantifiedOligoFile>(FixturePath);
            var row = direct.First();

            Assert.Multiple(() =>
            {
                Assert.That(direct.Count(), Is.EqualTo(492));
                Assert.That(factory.Count(), Is.EqualTo(492));
                Assert.That(direct.FileType, Is.EqualTo(SupportedFileType.FlashLFQQuantifiedOligo));
                Assert.That(direct.CanRead(FixturePath));
                Assert.That(row.BaseSequence, Is.EqualTo("AAAAAAAAACUCG"));
                Assert.That(row.Sequence, Is.EqualTo("AAAAAAAAACUCG"));
                Assert.That(row.ProteinGroups, Is.EqualTo("MALAT1"));
                Assert.That(row.Organism, Is.EqualTo("plasmid"));
                Assert.That(row.Samples, Has.Count.EqualTo(37));
                Assert.That(row.PeakOrder, Is.Null);
            });

            var mbr = row.Samples.Values.Single(s => s.DetectionType == "MBR");
            Assert.That(mbr.Intensity, Is.GreaterThan(0));
            Assert.That(mbr.RetentionTime, Is.Null);
        }

        [Test]
        public void PreservesModifiedRnaSequenceAndNotDetectedZero()
        {
            var file = new QuantifiedOligoFile(FixturePath);
            Assert.That(file.Any(r => r.Sequence.Contains("[Digestion Termini:Cyclic Phosphate on X]")));

            var notDetected = file.SelectMany(r => r.Samples.Values).First(s => s.DetectionType == "NotDetected");
            Assert.That(notDetected.Intensity, Is.EqualTo(0.0));
        }

        [Test]
        public void BlankSampleValuesAreNull()
        {
            string directory = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(TestQuantifiedOligoFile));
            Directory.CreateDirectory(directory);
            string path = Path.Combine(directory, "QuantifiedOligos.tsv");
            File.WriteAllLines(path,
            [
                "Sequence\tBase Sequence\tIntensity_f1\tDetection Type_f1",
                "ACGU\tACGU\t\t"
            ]);

            try
            {
                var sample = new QuantifiedOligoFile(path).Single().Samples["f1"];
                Assert.That((sample.Intensity, sample.DetectionType), Is.EqualTo(((double?)null, (string?)null)));
            }
            finally
            {
                Directory.Delete(directory, true);
            }
        }

        [Test]
        public void ReportsMalformedInputWithPathAndRefusesWriting()
        {
            string directory = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(TestQuantifiedOligoFile) + "Errors");
            Directory.CreateDirectory(directory);
            string path = Path.Combine(directory, "Broken_QuantifiedOligos.tsv");
            File.WriteAllLines(path, ["Sequence\tBase Sequence\tIntensity_f1", "ACGU\tACGU\tnot-a-number"]);

            try
            {
                var exception = Assert.Throws<MzLibException>(() => new QuantifiedOligoFile(path).LoadResults());
                Assert.That(exception!.Message, Does.Contain(path));
                Assert.Throws<NotSupportedException>(() => new QuantifiedOligoFile(FixturePath).WriteResults("x"));
            }
            finally
            {
                Directory.Delete(directory, true);
            }
        }
    }
}
