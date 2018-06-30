using Chemistry;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using static Chemistry.PeriodicTable;

namespace Test
{
    [TestFixture]
    public static class MyPeptideTest
    {
        [Test]
        public static void TestGoodPeptide()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            var protease = new Protease("CustomizedProtease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();

            Assert.AreEqual(2, ye.Count);

            var pep1 = ye[0];
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
            {
                Assert.IsTrue(huh > 0);
            }
            var pep2 = ye[1];
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
            {
                Assert.IsTrue(huh > 0);
            }
        }

        [Test]
        public static void TestNoCleavage()
        {
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();
            var prot = new Protein("MNNNKQQQQ", null, null, null, new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(5, 6, "lala") });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            var ye = prot.Digest(digestionParams, fixedModifications, new List<ModificationWithMass>()).ToList();
            Assert.AreEqual(3, ye.Count);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null);
            var protease = new Protease("Custom Protease7", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(2, ye.Count);
            var pep1 = ye[0];
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
            {
                Assert.IsTrue(huh > 0);
            }

            var pep2 = ye[1];
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsTrue(double.IsNaN(cool[1]));
            Assert.IsTrue(double.IsNaN(cool[2]));
            Assert.IsTrue(cool.Length == 3);
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null);
            DigestionParams digestionParams = new DigestionParams(protease: "Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            variableModifications.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, TerminusLocalization.NProt, ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("pepNmod", null, motif, TerminusLocalization.NPep, ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("resMod", null, motif, TerminusLocalization.Any, ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("PepCmod", null, motif, TerminusLocalization.PepC, ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, TerminusLocalization.ProtC, ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));

            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();
            Assert.AreEqual(3 * 2 * 3, ye.Count);

            Assert.AreEqual("[H]M[H][H]", ye.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ye.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null);
            List<ModificationWithMass> fixedMods = new List<ModificationWithMass>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            fixedMods.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, TerminusLocalization.NProt, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("PepNmod", null, motif, TerminusLocalization.NPep, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("resMod", null, motif, TerminusLocalization.Any, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("PepCmod", null, motif, TerminusLocalization.PepC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, TerminusLocalization.ProtC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            var ok = prot.Digest(digestionParams, fixedMods, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(1, ok.Count);

            Assert.AreEqual("[:PepNmod]M[:resMod][:ProtCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestDigestIndices()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            Modification mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> {mod } },
                {8, new List<Modification> {mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict);
            DigestionParams digestionParams = new DigestionParams(
                "Custom Protease7",
                maxMissedCleavages: 0,
                minPeptideLength: 5,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestList = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            var ok1 = digestList[1];
            var ok2 = digestList[3];

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.AllModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.ContainsKey(3));
        }

        [Test]
        public static void TestDigestDecoy()
        {
            ModificationMotif.TryGetMotif("Abcdefg", out ModificationMotif motif);
            Modification mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> { mod } },
                {8, new List<Modification> { mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict, isDecoy: true);
            DigestionParams digestionParams = new DigestionParams(
                "Custom Protease7",
                maxMissedCleavages: 0,
                minPeptideLength: 5,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var digestedList = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            var ok1 = digestedList[1];
            var ok2 = digestedList[3];

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.AllModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.ContainsKey(3));

            prot = new Protein("MNNNNKRRRRR", null, null, null, modDict);
            ok1 = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            ok2 = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).Last();

            Assert.AreEqual(0, ok1.NumMods);
            Assert.IsFalse(ok1.AllModsOneIsNterminus.Any());
            Assert.AreEqual(0, ok2.NumMods);
            Assert.IsFalse(ok2.AllModsOneIsNterminus.Any());
        }

        [Test]
        public static void TestGoodPeptideWithLength()
        {
            var prot = new Protein("MNNNKQQQQMNNNKQQQQ", null);

            DigestionParams digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye1 = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 1, maxPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye2 = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 5, maxPeptideLength: 8, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye3 = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            Assert.AreEqual(3, ye.Count);
            Assert.AreEqual(2, ye1.Count);
            Assert.AreEqual(2, ye2.Count);
            Assert.AreEqual(1, ye3.Count);
        }
    }
}