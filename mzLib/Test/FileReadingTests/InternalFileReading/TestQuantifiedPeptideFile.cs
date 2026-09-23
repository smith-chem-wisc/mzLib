using NUnit.Framework;
using Readers;
using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;

namespace Test.FileReadingTests.InternalFileReading
{
    /// <summary>
    /// The fixture is the header and first three rows of the AllQuantifiedPeptides.tsv MetaMorpheus 1.1.11
    /// wrote for the public PRIDE dataset PXD036557 (18 files): FlashLFQ's peptide format, with one MBR
    /// detection and NotDetected cells written as intensity 0.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestQuantifiedPeptideFile
    {
        private static string FixturePath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\ExternalFileTypes\MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv");

        [Test]
        public void ReadsTheFixedAndPerSampleColumns()
        {
            var file = FileReader.ReadFile<QuantifiedPeptideFile>(FixturePath);
            Assert.That(file.Count(), Is.EqualTo(3));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.FlashLFQQuantifiedPeptide));

            var row = file.First();
            Assert.That(row.BaseSequence, Is.EqualTo("CACASHVAK"));
            Assert.That(row.Sequence, Does.StartWith("[Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]"));
            Assert.That(row.ProteinGroups, Is.EqualTo("Q16643"));
            Assert.That(row.Samples, Has.Count.EqualTo(18));
            Assert.That(row.PeakOrder, Is.Null);

            var mbr = row.Samples["QE-002118_GM6_a-calib"];
            Assert.That((mbr.DetectionType, mbr.Intensity), Is.EqualTo(("MBR", (double?)143944.34375)));
            Assert.That(mbr.RetentionTime, Is.Null);
        }

        /// <summary>FlashLFQ writes 0, not a blank, for a peptide it did not quantify. The reader keeps the 0.</summary>
        [Test]
        public void NotDetectedIsAWrittenZero()
        {
            var notDetected = new QuantifiedPeptideFile(FixturePath).First()
                .Samples.Values.First(s => s.DetectionType == "NotDetected");
            Assert.That(notDetected.Intensity, Is.EqualTo(0.0));
        }

        [Test]
        public void ReadsTheIsoTrackerLayout()
        {
            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestQuantifiedPeptideFile");
            Directory.CreateDirectory(dir);
            string path = Path.Combine(dir, "IsoTracker_QuantifiedPeptides.tsv");
            File.WriteAllLines(path,
            [
                "Sequence\tBase Sequence\tPeak Order\tProtein Groups\tGene Names\tOrganism\tIntensity_f1\tRetentionTime (min)_f1\tDetection Type_f1",
                "PEPTIDE\tPEPTIDE\t2\tP1\tG1\tHuman\t1.5E+06\t42.5\tMSMS"
            ]);
            try
            {
                var row = new QuantifiedPeptideFile(path).Single();
                Assert.That(row.PeakOrder, Is.EqualTo(2));
                var f1 = row.Samples["f1"];
                Assert.That((f1.Intensity, f1.RetentionTime, f1.DetectionType), Is.EqualTo(((double?)1.5e6, (double?)42.5, "MSMS")));
            }
            finally
            {
                Directory.Delete(dir, true);
            }
        }

        [Test]
        public void WritingIsRefused()
        {
            Assert.Throws<NotSupportedException>(() => new QuantifiedPeptideFile(FixturePath).WriteResults("x"));
        }
    }
}
