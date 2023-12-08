using NUnit.Framework;
using Proteomics;
using System;
using Omics.Modifications;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestModFits
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
        public static void TestModFitss()
        {
            Protein protein = new Protein("M", null);
            int peptideOneBasedIndex = 1;
            int peptideLength = 1;
            int proteinOneBasedIndex = 1;

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            Modification attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);
            Assert.IsTrue(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("N", out motif);
            attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);
            Assert.IsFalse(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Mx", out motif);
            attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);
            Assert.IsFalse(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Mr", out motif);
            attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);
            Assert.IsFalse(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("xM", out motif);
            attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);
            Assert.IsFalse(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);

            protein = new Protein("MNRS", null);
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            Assert.IsFalse(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);

            protein = new Protein("MNRS", null);
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            Assert.IsFalse(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
            peptideOneBasedIndex = 2;
            peptideLength = 1;
            proteinOneBasedIndex = 2;
            Assert.IsTrue(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
            protein = new Protein("MNRN", null);
            Assert.IsFalse(ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
        }

        [Test]
        [TestCase("M", "X", true)]
        [TestCase("M", "J", false)]
        [TestCase("I", "J", true)]
        [TestCase("L", "X", true)]
        [TestCase("M", "B", false)]
        [TestCase("D", "B", true)]
        [TestCase("N", "B", true)]
        [TestCase("M", "Z", false)]
        [TestCase("E", "Z", true)]
        [TestCase("Q", "Z", true)]
        public static void TestAmbiguousModFits(string proteinSequence, string motifString, bool result)
        {
            Protein protein = new Protein(proteinSequence, null);
            int peptideOneBasedIndex = 1;
            int peptideLength = 1;
            int proteinOneBasedIndex = 1;
            ModificationMotif.TryGetMotif(motifString, out ModificationMotif motif);
            Modification attemptToLocalize = new Modification(_target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: double.NaN);
            Assert.AreEqual(result, ModificationLocalization.ModFits(attemptToLocalize, protein.BaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
        }
    }
}