using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class DecoyDigestionTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        [TestCase("cRAP_databaseGPTMD.xml")]
        [TestCase("uniprot_aifm1.fasta")]
        public static void TestDecoyScramblingIsReproducible(string fileName)
        {
            // Load in proteins
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", fileName);
            DecoyType decoyType = DecoyType.Reverse;
            List<Protein> proteins1 = null;
            List<Protein> proteins2 = null;
            if (fileName.Contains(".xml"))
            {
                proteins1 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out unknownModifications);
            }
            else if (fileName.Contains(".fasta"))
            {
                proteins1 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out unknownModifications);
            }
            else
            {
                NUnit.Framework.Assert.Fail("Unknown file type");
            }

            DigestionParams d = new DigestionParams(
                        maxMissedCleavages: 1,
                        minPeptideLength: 5,
                        initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            // Digest target proteins
            var pepsToReplace = proteins1.Where(p => !p.IsDecoy)
                .SelectMany(p => p.Digest(d, new List<Modification>(), new List<Modification>()).ToList())
                .Select(pep => pep.BaseSequence)
                .ToHashSet();

            // Ensure at least one decoy peptide from each protein is problematic and must be replaced
            var singleDecoyPeptides = proteins1
                .Where(p => p.IsDecoy)
                .Select(p => p.Digest(d, new List<Modification>(), new List<Modification>()).Skip(2).Take(1))
                .Select(pwsm => pwsm.First().BaseSequence)
                .ToHashSet();

            //modify targetpeptides in place
            pepsToReplace.UnionWith(singleDecoyPeptides);

            // Scramble every decoy from db1
            List<Protein> decoys1 = new();
            foreach (var protein in proteins1.Where(p => p.IsDecoy))
            {
                decoys1.Add(DecoySequenceValidator.ScrambleDecoyBioPolymer(protein, d, pepsToReplace));
            }
            // Scramble every decoy from db2
            List<Protein> decoys2 = new();
            foreach (var protein in proteins2.Where(p => p.IsDecoy))
            {
                decoys2.Add(DecoySequenceValidator.ScrambleDecoyBioPolymer(protein, d, pepsToReplace));
            }

            // check are equivalent lists of proteins
            Assert.AreEqual(decoys1.Count, decoys2.Count);
            foreach (var decoyPair in decoys1.Concat(decoys2).GroupBy(p => p.Accession))
            {
                Assert.AreEqual(2, decoyPair.Count());
                Assert.AreEqual(decoyPair.First().BaseSequence, decoyPair.Last().BaseSequence);
            }
        }

        [Test]
        public static void TestDecoyScramblerReplacesPeptides()
        {
            DigestionParams d = new DigestionParams(
                        maxMissedCleavages: 1,
                        minPeptideLength: 5,
                        initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            Protein target = new Protein("MEDEEKFVGYKYGVFK", "target");
            Protein decoy = new Protein("EEDEMKYGVFKFVGYK", "decoy");

            var targetPep = target.Digest(d, new List<Modification>(), new List<Modification>());
            var decoyPep = decoy.Digest(d, new List<Modification>(), new List<Modification>());

            HashSet<string> targetPepSeqs = targetPep.Select(p => p.FullSequence).ToHashSet();
            var offendingDecoys = decoyPep.Where(p => targetPepSeqs.Contains(p.FullSequence)).Select(d => d.FullSequence).ToList();

            Assert.AreEqual(2, offendingDecoys.Count);

            Protein scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs, offendingDecoys);
            var scrambledPep = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>());

            Assert.AreEqual(decoyPep.Count(), scrambledPep.Count());
            Assert.IsFalse(scrambledPep.Any(p => offendingDecoys.Contains(p.FullSequence)));

            // Check to make sure that decoy generation also works in no offending sequences are passed in
            scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs);
            scrambledPep = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>());

            Assert.AreEqual(decoyPep.Count(), scrambledPep.Count());
            Assert.IsFalse(scrambledPep.Any(p => offendingDecoys.Contains(p.FullSequence)));
        }
    }
}
