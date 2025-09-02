using System;
using System.IO;
using System.Linq;
using System.Reflection;
using Chemistry;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics
{
    [TestFixture]
    public class ModomicsLoaderTests
    {
        [Test]
        public void LoadModomics_ParsesEntriesCorrectly()
        {
            // Adjust resource name to match your project/namespace
            var resourceName = "UsefulProteomicsDatabases.Transcriptomics.modomics.json";

            var mods = ModomicsLoader.LoadModomics();

            Assert.That(mods, Is.Not.Empty);
            Assert.That(mods.All(m => !string.IsNullOrWhiteSpace(m.IdWithMotif)), Is.True);
            Assert.That(mods.All(m => m.IdWithMotif.Contains(" on ")), Is.True);
        }

        [Test]
        public void LoadModomics_SkipsUnknownBases()
        {
            var mods = ModomicsLoader.LoadModomics();

            // Should skip modifications with unknown bases (e.g. "X")
            Assert.That(mods.Any(m => m.IdWithMotif.Contains("on X")), Is.False);
        }

        [Test]
        public void LoadModomics_CachesResults()
        {
            var firstLoad = ModomicsLoader.LoadModomics();
            var secondLoad = ModomicsLoader.LoadModomics();
            Assert.That(ReferenceEquals(firstLoad, secondLoad), Is.True);
        }

        [Test]  
        public static void MethylsHaveCorrectFormula()
        {
            var mods = ModomicsLoader.LoadModomics();
            var singleMethylMods = mods.Where(m =>
                    System.Text.RegularExpressions.Regex.IsMatch(
                        m.IdWithMotif,
                        @"^m\d+[ACGU] on [ACGU]$"))
                    .ToList();

            // Exception to normal CH2 formula rule
            var m3C = singleMethylMods.FirstOrDefault(m => m.IdWithMotif.StartsWith("m3C on C"));
            Assert.That(m3C, Is.Not.Null, "Expected to find m3C modification");
            Assert.That(m3C.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C1H3")), Is.True, "m3C should have formula C1H3");
            singleMethylMods.Remove(m3C);

            var expectedFormula = ChemicalFormula.ParseFormula("C1H2");
            CollectionAssert.AreEqual(
                Enumerable.Repeat(expectedFormula, singleMethylMods.Count),
                singleMethylMods.Select(p => p.ChemicalFormula)
            );

            var diMethylMods = mods.Where(m =>
                    m.DatabaseReference["Modomics"].Any(p => p.Contains("dimethyl") && !p.Contains("inosine")) && !m.IdWithMotif.Contains("AA"))
                    .ToList();

            var phosphoDiMethyl = diMethylMods.Where(p => p.IdWithMotif.StartsWith("pm")).ToList();
            expectedFormula = ChemicalFormula.ParseFormula("C2H3O3P");
            CollectionAssert.AreEqual(
                Enumerable.Repeat(expectedFormula, phosphoDiMethyl.Count),
                phosphoDiMethyl.Select(p => p.ChemicalFormula));
            foreach (var pdm in phosphoDiMethyl)
                diMethylMods.Remove(pdm);

            expectedFormula = ChemicalFormula.ParseFormula("C2H4");
            CollectionAssert.AreEqual(
                Enumerable.Repeat(expectedFormula, diMethylMods.Count).ToList(),
                diMethylMods.Select(p => p.ChemicalFormula).ToList()
            );
        }
    }
}
