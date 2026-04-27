using Chemistry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SequenceConversion;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;
using Transcriptomics;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class PeptideWithSetModsEqualityTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        /// <summary>
        /// CRITICAL: Tests that peptides from different proteases are NOT equal even with identical sequences.
        /// This is essential for multiprotease parsimony in MetaMorpheus - without this distinction,
        /// peptides from different enzyme digests would be incorrectly collapsed during protein inference.
        /// </summary>
        [Test]
        public static void TestDifferentProteaseEquals()
        {
            Protein myProtein = new Protein("SEQUENCEK", "accession");

            DigestionParams digest1 = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            DigestionParams digest2 = new DigestionParams(protease: "Lys-C|P", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            PeptideWithSetModifications pep1 = myProtein.Digest(digest1, new List<Modification>(), new List<Modification>()).First();
            PeptideWithSetModifications pep2 = myProtein.Digest(digest2, new List<Modification>(), new List<Modification>()).First();

            Assert.That(pep1.FullSequence.Equals(pep2.FullSequence));
            Assert.That(pep1.Parent.Equals(pep2.Parent));
            Assert.That(!pep1.DigestionParams.DigestionAgent.Equals(pep2.DigestionParams.DigestionAgent));
            Assert.That(!pep1.Equals(pep2));
            Assert.That(!pep1.Equals((object)pep2));
            Assert.That(!pep1.GetHashCode().Equals(pep2.GetHashCode()));
        }

        /// <summary>
        /// CRITICAL: Tests type safety between PeptideWithSetModifications and OligoWithSetMods.
        /// Ensures these different biopolymer types are never incorrectly compared as equal,
        /// which would corrupt search results in multi-omics analyses.
        /// </summary>
        [Test]
        public static void TestPeptideOligoEquality()
        {
            var oligo = new OligoWithSetMods("GUACUG", []);
            var peptide = new PeptideWithSetModifications("PEPTIDE", []);

            Assert.That(!oligo.Equals(peptide));
            Assert.That(!peptide.Equals(oligo));
            Assert.That(!((IBioPolymerWithSetMods)oligo).Equals(peptide));
            Assert.That(!((IBioPolymerWithSetMods)peptide).Equals(oligo));
            Assert.That(!((object)oligo).Equals(peptide));
            Assert.That(!((object)peptide).Equals(oligo));
        }

        /// <summary>
        /// CRITICAL: Tests hash code generation for peptides with and without digestion params.
        /// Hash codes are essential for Dictionary/HashSet operations. A null return or
        /// inconsistent hashing would break peptide deduplication and grouping operations.
        /// </summary>
        [Test]
        public static void TestPeptideWithSetMod_GetHashCode()
        {
            PeptideWithSetModifications pep1 = new PeptideWithSetModifications("SEQUENCEK", new Dictionary<string, Modification>());
            int oneHashCode = pep1.GetHashCode();

            //if digestion params is not defined, the peptidewithsetmods should still return a hashcode.
            Assert.IsNotNull(oneHashCode);

            Protein myProtein = new Protein("SEQUENCEK", "accession");

            DigestionParams digest1 = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            PeptideWithSetModifications pep2 = myProtein.Digest(digest1, new List<Modification>(), new List<Modification>()).First();

            int twoHashCode = pep2.GetHashCode();

            //if digestion params IS defined, the peptidewithsetmods should  return a hashcode.
            Assert.IsNotNull(twoHashCode);
        }

        /// <summary>
        /// CRITICAL: Tests equality comparison for PeptideWithSetModifications.
        /// Correct equality is essential for HashSet/Dictionary operations, peptide
        /// deduplication, and protein inference. Tests same peptide, different proteins,
        /// different positions, and null comparisons.
        /// </summary>
        [Test]
        public static void TestPeptideWithSetModsEquals()
        {
            // Create two proteins
            Protein protein1 = new Protein("SEQUENCEK", "accession1");
            Protein protein2 = new Protein("SEQUENCEK", "accession2");

            // Create digestion parameters
            DigestionParams digestionParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            // Digest the proteins to get peptides
            PeptideWithSetModifications peptide1 = protein1.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            PeptideWithSetModifications peptide2 = protein2.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            // Test equality - same peptide
            Assert.IsTrue(peptide1.Equals(peptide1));

            // different peptide
            Assert.IsTrue(!peptide1.Equals(peptide2));
            Assert.IsTrue(!peptide1.Equals((object)peptide2));
            Assert.IsTrue(!peptide1.Equals((IBioPolymerWithSetMods)peptide2));
            Assert.AreNotEqual(peptide1.GetHashCode(), peptide2.GetHashCode());

            // Test inequality with different start residue
            PeptideWithSetModifications peptide3 = new PeptideWithSetModifications(protein1, digestionParams, 2, 9, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            Assert.IsFalse(peptide1.Equals(peptide3));

            // Test inequality with different parent accession
            PeptideWithSetModifications peptide4 = new PeptideWithSetModifications(protein2, digestionParams, 1, 9, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            Assert.IsFalse(peptide1.Equals(peptide4));

            // all fail on null
            Assert.That(!peptide1.Equals(null));
            Assert.That(!peptide1.Equals((object)null));
            Assert.That(!peptide1.Equals((PeptideWithSetModifications)null));
        }
    }
}
