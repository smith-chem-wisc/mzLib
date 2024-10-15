using NUnit.Framework;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework.Interfaces;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases.Transcriptomics;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestDecoyGeneration
    {
        public static string ModomicsUnmodifiedFastaPath => TestDbLoader.ModomicsUnmodifedFastaPath;

        [Test]
        public static void TestReverseDecoy_Simple()
        {
            var oligos = new List<RNA>()
            {
                new RNA("GUUCUG"),
                new RNA("GUGCUA"),
            };
            var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Reverse, 1);
            Assert.That(decoys.Count, Is.EqualTo(2));
            Assert.That(decoys[0].BaseSequence, Is.EqualTo("UCUUGG"));
            Assert.That(decoys[1].BaseSequence, Is.EqualTo("UCGUGA"));

            var example = oligos.First();
            Assert.That(decoys.All(p => !p.IsContaminant));
            Assert.That(decoys.All(p => p.IsDecoy));
            Assert.That(decoys.All(p => p.DatabaseFilePath == example.DatabaseFilePath));
            Assert.That(decoys.All(p => p.Organism == example.Organism));
            Assert.That(decoys.All(p => p.AdditionalDatabaseFields == example.AdditionalDatabaseFields));
            Assert.That(decoys.All(p => p.Accession == example.Accession));
            Assert.That(decoys.All(p => p.Name == example.Name));
            Assert.That(decoys.All(p => p.Length == example.Length));
            Assert.That(decoys.All(p => Equals(p.FivePrimeTerminus, example.FivePrimeTerminus)));
            Assert.That(decoys.All(p => Equals(p.ThreePrimeTerminus, example.ThreePrimeTerminus)));
            Assert.That(decoys.All(p => p.OneBasedPossibleLocalizedModifications.Count == example.OneBasedPossibleLocalizedModifications.Count));
        }

        [Test]
        [TestCase("GUACUG", 1, "UCAUGG", 5)]
        [TestCase("GUACUA", 2, "UCAUGA", 4)]
        [TestCase("GUACUA", 3, "UCAUGA", 3)]
        [TestCase("GUACUA", 4, "UCAUGA", 2)]
        [TestCase("GUCCAA", 5, "ACCUGA", 1)]
        [TestCase("GUUCUA", 6, "UCUUGA", 6)]
        public static void TestReverseDecoy_SimpleWithMods(string rnaSequence, int modPosition, string expectedDecoySequence, int expectedDecoyModPosition)
        {
            var mod = new Modification();
            var oligos = new List<RNA>()
            {
                new RNA(rnaSequence, null, null,
                    new Dictionary<int, List<Modification>>()
                        { { modPosition, new List<Modification>() { mod } } }),
            };
            var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Reverse, 1);
            Assert.That(decoys.Count, Is.EqualTo(1));

            var decoy = decoys.First();
            var originalRna = oligos.First();
            Assert.That(decoy.BaseSequence, Is.EqualTo(expectedDecoySequence));
            Assert.That(decoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
            Assert.That(decoy.OneBasedPossibleLocalizedModifications.First().Key, Is.EqualTo(expectedDecoyModPosition));
            Assert.That(decoy.OneBasedPossibleLocalizedModifications.First().Value.Count, Is.EqualTo(1));
            Assert.That(decoy.OneBasedPossibleLocalizedModifications.First().Value.First(), Is.EqualTo(mod));
            Assert.That(decoy.Name, Is.EqualTo(originalRna.Name));
            Assert.That(decoy.Accession, Is.EqualTo(originalRna.Accession));
            Assert.That(decoy.Organism, Is.EqualTo(originalRna.Organism));
            Assert.That(decoy.DatabaseFilePath, Is.EqualTo(originalRna.DatabaseFilePath));
            Assert.That(decoy.IsContaminant, Is.EqualTo(originalRna.IsContaminant));
            Assert.That(decoy.IsDecoy, Is.True);
            Assert.That(decoy.AdditionalDatabaseFields, Is.EqualTo(originalRna.AdditionalDatabaseFields));
            Assert.That(decoy.FivePrimeTerminus, Is.EqualTo(originalRna.FivePrimeTerminus));
            Assert.That(decoy.ThreePrimeTerminus, Is.EqualTo(originalRna.ThreePrimeTerminus));
        }

        [Test]
        public void TestReverseDecoy_FromDatabase()
        {
            int numSequences = 5;
            Dictionary<string, string> expectedSequences = new Dictionary<string, string>()
            {
                { "tdbR00000010", "CCACCUCGAUACGCCCUAGCUUGGCGUCUGGAGGACGCACGUUUCGUCCGCGAGAGGGUCGACUCGAUAUCGGGGA"},
                { "tdbR00000008", "CCACCUCGAUUCGCCCUAGCUUGGCGACUGGAGAACGUACGGUACGUUCGCGAGAGGGUCGACUCGAUAUCGGGGA"},
                { "tdbR00000356", "CCACGUAGGCCCUCCUAAGCUUGGAGGCUGGCGAGCCAAGCAUCGGCUCAUGAGAUAGGUCGACUCGAUGCCUACGA"},
                { "tdbR00000359", "CCGCGCGGGCUGUCCUAAGCUUGGACUCUGGAGACGGAGGCCUCCCGUCGCGAGAUAGGUCGACUCGAUGCCCGCGA"},
                { "tdbR00000358", "CCGCGCGGGACGUCCUAAGCUUGGACGCCGGGUGCUGAAUCUUCCAGCAACGAGAUAGGUUGACUCGAUUCCCGCGA"},
            };

            var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifiedFastaPath, true, DecoyType.Reverse, false,
                out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(numSequences * 2));
            Assert.That(oligos.Count(p => p.IsDecoy), Is.EqualTo(numSequences));
            Assert.That(oligos.Count(p => !p.IsDecoy), Is.EqualTo(numSequences));

            foreach (var targetDecoyGroup in oligos.GroupBy(p => p.Name))
            {
                Assert.That(targetDecoyGroup.Count(), Is.EqualTo(2));
                var target = targetDecoyGroup.First(p => !p.IsDecoy);
                var decoy = targetDecoyGroup.First(p => p.IsDecoy);
                var expectedSequence = expectedSequences[target.Name];

                Assert.That(target.FivePrimeTerminus, Is.EqualTo(decoy.FivePrimeTerminus));
                Assert.That(target.ThreePrimeTerminus, Is.EqualTo(decoy.ThreePrimeTerminus));
                Assert.That(target.AdditionalDatabaseFields, Is.EqualTo(decoy.AdditionalDatabaseFields));
                Assert.That(target.IsContaminant, Is.EqualTo(decoy.IsContaminant));
                Assert.That(target.DatabaseFilePath, Is.EqualTo(decoy.DatabaseFilePath));
                Assert.That(target.DatabaseFilePath, Is.EqualTo(ModomicsUnmodifiedFastaPath));
                Assert.That(target.Organism, Is.EqualTo(decoy.Organism));
                Assert.That(target.Accession, Is.EqualTo(decoy.Accession));
                Assert.That(target.Name, Is.EqualTo(decoy.Name));
                Assert.That(target.Length, Is.EqualTo(decoy.Length));
                Assert.That(target.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(decoy.OneBasedPossibleLocalizedModifications.Count));

                Assert.That(decoy.BaseSequence, Is.EqualTo(expectedSequence));
            }
        }


        // TODO: Implement these test once other decoy generation methods are availiable

        [Test]
        public void TestShuffledDecoy_Simple()
        {
            var oligos = new List<RNA>()
            {
                new RNA("GUACUG"),
                new RNA("GUACUA"),
            };
            Assert.Throws<NotImplementedException>(() =>
            {
                var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Shuffle);
            });


            //var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Shuffle);
            //Assert.That(decoys.Count, Is.EqualTo(2));
        }

        [Test]
        public void TestShuffledDecoy_SimpleWithMods()
        {
            var oligos = new List<RNA>()
            {
                new RNA("GUACUG"),
                new RNA("GUACUA"),
            };
            Assert.Throws<NotImplementedException>(() =>
            {
                var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Shuffle);
            });
            //var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Shuffle);
            //Assert.That(decoys.Count, Is.EqualTo(2));
        }

        [Test]
        public void TestShuffledDecoy_FromDatabase()
        {
            Assert.Throws<NotImplementedException>(() =>
            {
                var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifiedFastaPath, true, DecoyType.Shuffle, false, out var errors);
            });

            //var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifiedFastaPath, true, DecoyType.Shuffle, false, out var errors);
            //Assert.That(errors.Count, Is.EqualTo(0));
            //Assert.That(oligos.Count, Is.EqualTo(10));
        }

        [Test]
        public void TestSlideDecoy_Simple()
        {
            var oligos = new List<RNA>()
            {
                new RNA("GUACUG"),
                new RNA("GUACUA"),
            };
            Assert.Throws<NotImplementedException>(() =>
            {
                var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Slide);
            });

            //var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Slide);
            //Assert.That(decoys.Count, Is.EqualTo(2));
        }

        [Test]
        public void TestSlideDecoy_SimpleWithMods()
        {
            var oligos = new List<RNA>()
            {
                new RNA("GUACUG"),
                new RNA("GUACUA"),
            };

            Assert.Throws<NotImplementedException>(() =>
            {
                var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Slide);
            });

            //var decoys = RnaDecoyGenerator.GenerateDecoys(oligos, DecoyType.Slide);
            //Assert.That(decoys.Count, Is.EqualTo(2));
        }

        [Test]
        public void TestSlideDecoy_FromDatabase()
        {
            Assert.Throws<NotImplementedException>(() =>
            {
                var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifiedFastaPath, true, DecoyType.Shuffle, false, out var errors);
            });

            //var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifiedFastaPath, true, DecoyType.Slide, false, out var errors);
            //Assert.That(errors.Count, Is.EqualTo(0));
            //Assert.That(oligos.Count, Is.EqualTo(10));
        }


        [Test]
        public void TestCreateNew()
        {
            var mods = PtmListLoader.ReadModsFromString(
                "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
                out List<(Modification, string)> modsOut).ToList();
            var modDict = mods.ToDictionary(p => p.IdWithMotif, p => p);
            var oneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>()
            {
                { 1, new List<Modification>() { modDict["Sodium on A"] } },
                { 3, new List<Modification>() { modDict["Sodium on A"] } },
            };

            var rna = new RNA("GAACUG", "name", "accession", "organism", "databaseFilePath",
                null, null, oneBasedPossibleLocalizedModifications, false, false, new List<Tuple<string, string>>(),
                new Dictionary<string, string>());
            var oligos = rna
                .Digest(new RnaDigestionParams(maxMods: 1), new List<Modification>(), mods)
                .ToList();

            var clonedRna = rna.CreateNew(null, null, true);
            var clonedOligo =  oligos.First().CreateNew(null, null, true);

            // ensure they are identical except for the isDecoy field
            Assert.That(rna.BaseSequence, Is.EqualTo(clonedRna.BaseSequence));
            Assert.That(rna.OneBasedPossibleLocalizedModifications, Is.EqualTo(clonedRna.OneBasedPossibleLocalizedModifications));
            Assert.That(rna.IsDecoy, Is.Not.EqualTo(clonedRna.IsDecoy));

            Assert.That(oligos.First().BaseSequence, Is.EqualTo(clonedOligo.BaseSequence));
            Assert.That(oligos.First().OneBasedPossibleLocalizedModifications, Is.EqualTo(clonedOligo.OneBasedPossibleLocalizedModifications));
            Assert.That(oligos.First().Parent.IsDecoy, Is.Not.EqualTo(clonedOligo.Parent.IsDecoy));


            var newMods = new Dictionary<int, List<Modification>>()
            {
                { 1, new List<Modification>() { modDict["Sodium on A"] } },
                { 2, new List<Modification>() { modDict["Sodium on A"] } },
                { 3, new List<Modification>() { modDict["Sodium on A"] } },
            };
            clonedRna = rna.CreateNew("AAAAAA", newMods, null);
            clonedOligo = oligos.First().CreateNew("AAAAAA", newMods, null);

            Assert.That(rna.BaseSequence, Is.Not.EqualTo(clonedRna.BaseSequence));
            Assert.That(rna.OneBasedPossibleLocalizedModifications, Is.Not.EqualTo(clonedRna.OneBasedPossibleLocalizedModifications));
            Assert.That(rna.IsDecoy, Is.EqualTo(clonedRna.IsDecoy));

            Assert.That(oligos.First().BaseSequence, Is.Not.EqualTo(clonedOligo.BaseSequence));
            Assert.That(oligos.First().OneBasedPossibleLocalizedModifications, Is.Not.EqualTo(clonedOligo.OneBasedPossibleLocalizedModifications));
            Assert.That(oligos.First().Parent.IsDecoy, Is.EqualTo(clonedOligo.Parent.IsDecoy));
        }
    }
}
