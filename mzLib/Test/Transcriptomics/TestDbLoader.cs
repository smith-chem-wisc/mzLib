using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics
{
    [TestFixture]
    internal class TestDbLoader
    {
        public static string ModomicsUnmodifedFastaPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.fasta");


        [Test]
        public static void TestModomicsUnmodifiedFasta()
        {
            var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.None, false,
                out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(5));
            Assert.That(oligos.First().BaseSequence,
                Is.EqualTo("GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA"));
            Assert.That(oligos.First().Name, Is.EqualTo("tdbR00000010"));
            Assert.That(oligos.First().Accession, Is.EqualTo("SO:0000254"));
            Assert.That(oligos.First().Organism, Is.EqualTo("Escherichia coli"));
            Assert.That(oligos.First().DatabaseFilePath, Is.EqualTo(ModomicsUnmodifedFastaPath));
            Assert.That(oligos.First().IsContaminant, Is.False);
            Assert.That(oligos.First().IsDecoy, Is.False);
            Assert.That(oligos.First().AdditionalDatabaseFields!.Count, Is.EqualTo(5));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Id"], Is.EqualTo("1"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Type"], Is.EqualTo("tRNA"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Subtype"], Is.EqualTo("Ala"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Feature"], Is.EqualTo("VGC"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Cellular Localization"], Is.EqualTo("prokaryotic cytosol"));
        }

        [Test]
        public static void TestContaminantFollowsThrough()
        {
            var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.None, true,
                               out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(5));
            Assert.That(oligos.First().BaseSequence,
                               Is.EqualTo("GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA"));
            Assert.That(oligos.All(p => p.IsContaminant));
            Assert.That(oligos.All(p => !p.IsDecoy));
        }

        [Test]
        public static void TestNotGeneratingTargetsOrDecoys()
        {
            var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, false, DecoyType.None, true,
                out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestModomicsModifiedFasta()
        {

        }

        [Test]
        public static void Test__Fna()
        {

        }
    }
}
