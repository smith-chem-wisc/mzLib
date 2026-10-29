using NUnit.Framework;
using Omics.Modifications;
using Omics.Modifications.IO;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Chemistry;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestTerminalModFixes
    {
        private static Modification Mod(string restriction)
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            return new Modification("Phospho", "", "test", "", motif, restriction, ChemicalFormula.ParseFormula("HO3P"));
        }

        [Test]
        [TestCase("5'-terminal.")]
        [TestCase("Oligo 5'-terminal.")]
        public static void TestReverseDecoy_FivePrimeModStaysOnFivePrimeEnd(string restriction)
        {
            var mod = Mod(restriction);
            var rna = new RNA("GUACUG", "rna1",
                new Dictionary<int, List<Modification>> { { 1, new List<Modification> { mod } } }, null, null, "", null, null);

            var decoy = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1).Single();

            Assert.That(decoy.OneBasedPossibleLocalizedModifications.Keys, Is.EquivalentTo(new[] { 1 }));
            Assert.That(decoy.OneBasedPossibleLocalizedModifications[1].Single(), Is.EqualTo(mod));
        }

        [Test]
        public static void TestReverseDecoy_ThreePrimeModStaysOnThreePrimeEnd()
        {
            var mod = Mod("3'-terminal.");
            var rna = new RNA("GUACUG", "rna1",
                new Dictionary<int, List<Modification>> { { 6, new List<Modification> { mod } } }, null, null, "", null, null);

            var decoy = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1).Single();

            Assert.That(decoy.OneBasedPossibleLocalizedModifications.Keys, Is.EquivalentTo(new[] { 6 }));
        }

        [Test]
        public static void TestReverseDecoy_MixedModsAtTerminusAreSplit()
        {
            ModificationMotif.TryGetMotif("G", out var g);
            var anywhere = new Modification("Methyl", "", "test", "", g, "Anywhere.", ChemicalFormula.ParseFormula("CH2"));
            var fivePrime = Mod("5'-terminal.");
            var rna = new RNA("GUACUA", "rna1",
                new Dictionary<int, List<Modification>> { { 1, new List<Modification> { fivePrime, anywhere } } }, null, null, "", null, null);

            var decoy = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1).Single();

            Assert.That(decoy.OneBasedPossibleLocalizedModifications[1].Single(), Is.EqualTo(fivePrime));
            Assert.That(decoy.OneBasedPossibleLocalizedModifications[6].Single(), Is.EqualTo(anywhere));
        }

        [Test]
        public static void TestModFits_FivePrimeTerminalOnlyAtPositionOne()
        {
            var mod = Mod("5'-terminal.");
            Assert.That(ModificationLocalization.ModFits(mod, "GUACUG", 1, 6, 1), Is.True);
            Assert.That(ModificationLocalization.ModFits(mod, "GUACUG", 1, 5, 2), Is.False);

            // proteins keep the initiator-Met rule
            ModificationMotif.TryGetMotif("X", out var x);
            var nTerm = new Modification("Acetyl", "", "test", "", x, "N-terminal.", ChemicalFormula.ParseFormula("C2H2O"));
            Assert.That(ModificationLocalization.ModFits(nTerm, "MPEPTIDE", 1, 7, 2), Is.True);
        }

        [Test]
        public static void TestReadMod_KeywordStartingWithHash()
        {
            // MODOMICS one-letter code of 2'-O-methylguanosine is '#', which G-PTM-D writes into stored mods
            string stored = string.Join(Environment.NewLine,
                "ID   2'-O-methylguanosine on G",
                "MT   Modomics",
                "TG   G",
                "PP   Anywhere.",
                "CF   CH2",
                "KW   # or 2'-O-methylguanosine or Gm",
                "//");

            var mods = ModificationLoader.ReadModsFromString(stored, out var warnings).ToList();

            Assert.That(mods.Count, Is.EqualTo(1));
            Assert.That(mods[0].Keywords.Single(), Is.EqualTo("# or 2'-O-methylguanosine or Gm"));
        }
    }
}
