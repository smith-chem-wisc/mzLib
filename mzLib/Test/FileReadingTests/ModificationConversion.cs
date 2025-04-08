using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    public static class ModificationConversions
    {
        [Test]
        public static void ModsLoadWithoutCrash()
        {
            var allMods = ModificationConverter.AllKnownMods;
            Assert.That(allMods.Count, Is.GreaterThanOrEqualTo(3110)); // number at time of creation
        }

        [Test]
        [TestCase("Carbamidomethyl", "Carbamidomethyl", 'C')]
        [TestCase("O-(ADP-ribosyl)-L-serine", "ADP-ribosylation", 'S')]
        [TestCase("Acetyl", "Acetylation", 'K')]
        [TestCase("N6,N6,N6-trimethyl-L-lysine", "N6,N6,N6-trimethyllysine", 'K')]
        [TestCase("N,N,N-trimethyl-L-alanine", "N,N,N-trimethylalanine", 'A')]
        [TestCase("O4'-(phospho-5'-adenosine)-L-tyrosine", "Phosphoadenosine", 'Y')]
        [TestCase("O4'-(phopho-ChickenPotPie-adenosine)-L-tyrosine", "Phosphoadenosine", 'Y')]
        [TestCase("Acetyldeoxyhypusine", "Acetyldeoxyhypusine", 'K')]
        [TestCase("15N-oxobutanoic", "15N-oxobutanoic", 'T')]
        [TestCase("Amidine", "Amidine", 'T')]
        public static void ModsConvertToModification(string name, string expectedId, char residue)
        {
            var modification = ModificationConverter.GetClosestMod(name, residue);
            Assert.That(modification.OriginalId, Is.EqualTo(expectedId));

            if (name == "Amidine") // Edge Case where we want to pass in the wrong motif and check if it got the right one
                residue = 'X';

            var withMotif = $"{expectedId} on {residue}";
            Assert.That(modification.IdWithMotif, Does.Contain(withMotif));
        }
    }
}
