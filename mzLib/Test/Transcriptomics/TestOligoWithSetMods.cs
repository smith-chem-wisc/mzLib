using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Omics.Modifications;
using Transcriptomics.Digestion;
using Transcriptomics;
using Omics;
using MassSpectrometry;
using Omics.Fragmentation;

namespace Test.Transcriptomics
{
    [ExcludeFromCodeCoverage]
    public static class TestOligoWithSetMods
    {
        [Test]
        [TestCase( 0, 1, 20.45)]
        [TestCase(1, 1, 20.45)]
        [TestCase( 0, 2, 20.45)]
        [TestCase(1, 2, 20.45)]
        [TestCase( 0, 5, 28.37)]
        [TestCase(1, 5, 28.37)]
        [TestCase( 0, 6, 28.37)]
        [TestCase(1, 6, 28.37)]
        public static void TestLocalize(int modsOnOligo, int indexOfMass, double massToLocalize)
        {
            var oligoWithSetMods = new RNA("GUACUG", 
                    oneBasedPossibleLocalizedModifications: new Dictionary<int, List<Modification>> { { 4, [TestDigestion.PotassiumAdducts[1]] } })
                .Digest(new RnaDigestionParams(), [], [])
                .ElementAt(modsOnOligo);

            Assert.That(oligoWithSetMods.AllModsOneIsNterminus.Count, Is.EqualTo(modsOnOligo));

            // Act
            var localizedOligo = oligoWithSetMods.Localize(indexOfMass - 2, massToLocalize);

            // Assert
            int expectedModificationCount;
            double expectedMass;
            if (modsOnOligo == 1) // if the oligo started with a mod
            {
                int indexOfOriginalMod = oligoWithSetMods.AllModsOneIsNterminus.Keys.First();

                // ensure original modification exist
                Assert.That(localizedOligo.AllModsOneIsNterminus.ContainsKey(indexOfOriginalMod));

                if (indexOfOriginalMod != indexOfMass) // Additional mass was added to a different location
                {
                    expectedModificationCount = modsOnOligo + 1;
                    expectedMass = massToLocalize;

                    // ensure original modification is still intact
                    Assert.That(oligoWithSetMods.OneBasedPossibleLocalizedModifications[indexOfOriginalMod][0].MonoisotopicMass, 
                        Is.EqualTo(localizedOligo.AllModsOneIsNterminus[indexOfOriginalMod].MonoisotopicMass));
                }
                else // Additional mass was added to the location of an existing modification
                {
                    expectedModificationCount = modsOnOligo;
                    expectedMass = massToLocalize + TestDigestion.PotassiumAdducts[1].MonoisotopicMass!.Value;

                    // ensure original modification has been altered
                    Assert.That(oligoWithSetMods.OneBasedPossibleLocalizedModifications[indexOfOriginalMod][0].MonoisotopicMass,
                        Is.Not.EqualTo(localizedOligo.AllModsOneIsNterminus[indexOfOriginalMod].MonoisotopicMass));
                }
            }
            else // oligo started with no modifications
            {
                expectedModificationCount = modsOnOligo + 1;
                expectedMass = massToLocalize;
            }


            Assert.That(expectedModificationCount, Is.EqualTo(localizedOligo.AllModsOneIsNterminus.Count));
            Assert.That(localizedOligo.AllModsOneIsNterminus.ContainsKey(indexOfMass));
            Assert.That(expectedMass, Is.EqualTo(localizedOligo.AllModsOneIsNterminus[indexOfMass].MonoisotopicMass));
        }

        [Test]
        public static void TestEquality()
        {
            var modDict = new Dictionary<int, List<Modification>> { { 4, [TestDigestion.PotassiumAdducts[1]] } };
            var oligoWithSetMods = new RNA("GUACUG",
                    oneBasedPossibleLocalizedModifications: modDict)
                .Digest(new RnaDigestionParams(), [], [])
                .ElementAt(1);

            IBioPolymerWithSetMods oligoWithSetMods2 = new RNA("GUACUG",
                    oneBasedPossibleLocalizedModifications: modDict)
                .Digest(new RnaDigestionParams(), [], [])
                .ElementAt(1);

            // same oligos
            Assert.That(oligoWithSetMods.Equals(oligoWithSetMods2));
            Assert.That(oligoWithSetMods.Equals((object)oligoWithSetMods2));
            Assert.That(oligoWithSetMods.Equals((OligoWithSetMods)oligoWithSetMods2));
            Assert.That(oligoWithSetMods.Equals(oligoWithSetMods));
            Assert.That(oligoWithSetMods.Equals((object)oligoWithSetMods));
            Assert.That(oligoWithSetMods.Equals((OligoWithSetMods)oligoWithSetMods));
            Assert.That(oligoWithSetMods.GetHashCode(), Is.EqualTo(oligoWithSetMods2.GetHashCode()));

            // all fail on null
            Assert.That(!oligoWithSetMods2.Equals(null));
            Assert.That(!oligoWithSetMods2.Equals((object)null));
            Assert.That(!oligoWithSetMods2.Equals((OligoWithSetMods)null));

            // Null parent checks
            oligoWithSetMods = new(oligoWithSetMods.FullSequence, modDict.ToDictionary(p => p.Value.First().IdWithMotif, p => p.Value.First()));
            oligoWithSetMods2 = new OligoWithSetMods(oligoWithSetMods.FullSequence, modDict.ToDictionary(p => p.Value.First().IdWithMotif, p => p.Value.First()));
            var oligoWithSetMods3 = new OligoWithSetMods(oligoWithSetMods.FullSequence + "AGAUA", modDict.ToDictionary(p => p.Value.First().IdWithMotif, p => p.Value.First()));

            // same oligo null parent
            Assert.That(oligoWithSetMods.Equals(oligoWithSetMods2));
            Assert.That(oligoWithSetMods.Equals((object)oligoWithSetMods2));
            Assert.That(oligoWithSetMods.Equals((OligoWithSetMods)oligoWithSetMods2));

            // different oligo null parent
            Assert.That(!oligoWithSetMods.Equals(oligoWithSetMods3));
            Assert.That(!oligoWithSetMods.Equals((object)oligoWithSetMods3));
            Assert.That(!oligoWithSetMods.Equals((IBioPolymerWithSetMods)oligoWithSetMods3));
        }

        [Test]
        [TestCase("GUACUG", "GUACUGGUACUG", "RNase A")]
        [TestCase("GUAGGAG", "GUAGCAG", "RNase A")]
        public static void TestInequality_DifferentParentSameDigestionProduct(string sequence1, string sequence2, string enzyme)
        {
            var digestionParams = new RnaDigestionParams(rnase: enzyme, minLength: 1, maxMissedCleavages: 0);

             var oligo1 = new RNA(sequence1,"rna1")
                .Digest(digestionParams, [], [])
                .First();

            var oligo2 = new RNA(sequence2, "rna3")
                .Digest(digestionParams, [], [])
                .First();

            Assert.That(oligo1, Is.Not.EqualTo(oligo2));
            Assert.That(oligo1.Equals(oligo1));
            Assert.That(oligo1, Is.Not.EqualTo((object)oligo2));
            Assert.That(oligo1.GetHashCode(), Is.Not.EqualTo(oligo2.GetHashCode()));
        }

        /// <summary>
        /// The purpose of this test is to ensure that two oligos digested from two different rnases are not equal even if their sequences are equal
        /// This is important for multiprotease parsimony in MetaMorpheus
        /// </summary>
        [Test]
        [TestCase("AUAGUCUGG", "RNase T1", "colicin_E5")]
        [TestCase("AUAGUCUGGGAUCUG",  "RNase T1", "colicin_E5")]
        public static void TestInequality_SameParentAndDigestionProduct_DifferentRnases(string sequence, string enzyme1, string enzyme2)
        {
            var digestionParams1 = new RnaDigestionParams(rnase: enzyme1, minLength: 1, maxMissedCleavages: 0);
            var digestionParams2 = new RnaDigestionParams(rnase: enzyme2, minLength: 1, maxMissedCleavages: 0);

            var oligo1 = new RNA(sequence)
                .Digest(digestionParams1, [], [])
                .ToArray();

            var oligo2 = new RNA(sequence)
                .Digest(digestionParams2, [], [])
                .ToArray();

            Assert.That(oligo1.Length, Is.Not.EqualTo(oligo2.Length));

            Assert.That(oligo1.First().BaseSequence, Is.EqualTo("AUAG"));
            Assert.That(oligo2.First().BaseSequence, Is.EqualTo("AUAG"));

            Assert.That(oligo1, Is.Not.EqualTo(oligo2));
            Assert.That(oligo1, Is.Not.EqualTo((object)oligo2));
            Assert.That(oligo1.GetHashCode(), Is.Not.EqualTo(oligo2.GetHashCode()));
        }

        [Test]
        public static void Fragment_WorksWhenParentIsNull()
        {
            // Arrange: create an OligoWithSetMods with no parent (using the string constructor)
            var baseSequence = "GUACUG";
            var mod = TestDigestion.PotassiumAdducts[1];
            var modDict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
            var oligo = new OligoWithSetMods(baseSequence, modDict);

            // Confirm parent is null
            Assert.That(oligo.Parent == null);

            // Act: call Fragment
            var products = new List<Product>();
            oligo.Fragment(DissociationType.CID, FragmentationTerminus.Both, products);

            // Assert: products are generated and not empty
            Assert.That(products, Is.Not.Null);
            Assert.That(products.Count, Is.GreaterThan(0), "Fragment should generate at least one product when parent is null.");

            // Optionally, check that all products are for the correct sequence length
            foreach (var product in products)
            {
                Assert.That(product.NeutralMass, Is.GreaterThan(0), "Product neutral mass should be positive.");
            }
        }
    }
}
