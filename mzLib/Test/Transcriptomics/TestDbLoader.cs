using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Omics.Modifications;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
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
        public static void TestXmlWriter()
        {
            var rna = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.None, false, out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));

            var modString = "ID   Methylation\r\nMT   Biological\r\nPP   Anywhere.\r\nTG   G\r\nCF   C1H2\r\n" + @"//";
            var methylG = PtmListLoader.ReadModsFromString(modString, out List<(Modification, string)> modsOut).First();

            Dictionary<string,HashSet<Tuple<int, Modification>>> mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            mods.Add("SO:0000254", new HashSet<Tuple<int, Modification>>()
            {
                new Tuple<int, Modification>(1, methylG),
                new Tuple<int, Modification>(3, methylG)
            });

            string outpath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.xml");

            var xml = ProteinDbWriter.WriteXmlDatabase(mods, rna, outpath);
            var temp = RnaDbLoader.LoadRnaXML(outpath, true, DecoyType.None, false,
                new List<Modification>() { methylG }, new List<string>(), out var unknownMods);

            Assert.That(unknownMods.Count, Is.EqualTo(0));
            Assert.That(temp.Count, Is.EqualTo(5));
            var first = temp.First();
            var loadedMods    = first.OneBasedPossibleLocalizedModifications;
            Assert.That(loadedMods.Count, Is.EqualTo(2));
            Assert.That(loadedMods[1].Count, Is.EqualTo(1));
            Assert.That(loadedMods[3].Count, Is.EqualTo(1));
            Assert.That(loadedMods[1].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
            Assert.That(loadedMods[3].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
        }
    }
}
