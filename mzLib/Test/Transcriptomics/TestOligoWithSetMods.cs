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
    }
}
