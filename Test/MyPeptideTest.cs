using Chemistry;
using MassSpectrometry;
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
            List<ModificationGeneral> variableModifications = new List<ModificationGeneral>();
            var ye = prot.Digest(digestionParams, new List<ModificationGeneral>(), variableModifications).ToList();

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
            List<ModificationGeneral> fixedModifications = new List<ModificationGeneral>();
            var prot = new Protein("MNNNKQQQQ", null, null, null, new Dictionary<int, List<ModificationGeneral>>(), new List<ProteolysisProduct> { new ProteolysisProduct(5, 6, "lala") });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            var ye = prot.Digest(digestionParams, fixedModifications, new List<ModificationGeneral>()).ToList();
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
            var ye = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).ToList();

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
            var protease = new Protease("Custom Protease7", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(protease: "Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3);
            List<ModificationGeneral> variableModifications = new List<ModificationGeneral>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);

            variableModifications.Add(new ModificationGeneral(_Id: "ProtNmod", _Target: motif, _Position: "N-terminal.", _ChemicalFormula: ChemicalFormula.ParseFormula("H"), _MonoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationGeneral(_Id: "pepNmod", _Target: motif, _Position: "Peptide N-terminal.", _ChemicalFormula: ChemicalFormula.ParseFormula("H"), _MonoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationGeneral(_Id: "resMod", _Target: motif, _Position: "Anywhere.", _ChemicalFormula: ChemicalFormula.ParseFormula("H"), _MonoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationGeneral(_Id: "PepCmod", _Target: motif, _Position: "Peptide C-terminal.", _ChemicalFormula: ChemicalFormula.ParseFormula("H"), _MonoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationGeneral(_Id: "ProtCmod", _Target: motif, _Position: "C-terminal.", _ChemicalFormula: ChemicalFormula.ParseFormula("H"), _MonoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            var ye = prot.Digest(digestionParams, new List<ModificationGeneral>(), variableModifications).ToList();
            Assert.AreEqual(3 * 2 * 3, ye.Count);
            Assert.AreEqual("[H]M[H][H]", ye.Last().SequenceWithChemicalFormulas);

            double m1 = 5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass;
            double m2 = ye.Last().MonoisotopicMass;
            double m3 = m1 - m2;

            Assert.IsTrue(m3 < 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null);
            List<ModificationGeneral> fixedMods = new List<ModificationGeneral>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            fixedMods.Add(new ModificationGeneral("ProtNmod", null, null, null, motif, "Anywhere.", ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, null, null, null, null, null));
            fixedMods.Add(new ModificationGeneral("pepNmod", null, null, null, motif, "Peptide N-terminal.", ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, null, null, null, null, null));
            fixedMods.Add(new ModificationGeneral("resMod", null, null, null, motif, "Anywhere.", ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, null, null, null, null, null));
            fixedMods.Add(new ModificationGeneral("PepCmod", null, null, null, motif, "Peptide C-terminal.", ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, null, null, null, null, null));
            fixedMods.Add(new ModificationGeneral("ProtCmod", null, null, null, motif, "C-terminal.", ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, null, null, null, null, null));

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            var ok = prot.Digest(digestionParams, fixedMods, new List<ModificationGeneral>()).ToList();

            Assert.AreEqual(1, ok.Count);

            Assert.AreEqual("[:pepNmod]M[:resMod][:ProtCmod]", ok.First().Sequence);

            Assert.AreEqual("[H]M[H][H]", ok.First().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestDigestIndices()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            ModificationGeneral mod = new ModificationGeneral(null, null, null, null, motif, "Anywhere.", null, double.NaN, null, null, null, null, null, null);
            IDictionary<int, List<ModificationGeneral>> modDict = new Dictionary<int, List<ModificationGeneral>>
            {
                {2, new List<ModificationGeneral> {mod } },
                {8, new List<ModificationGeneral> {mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict);
            DigestionParams digestionParams = new DigestionParams(
                "Custom Protease7",
                maxMissedCleavages: 0,
                minPeptideLength: 5,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestList = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).ToList();
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
            ModificationGeneral mod = new ModificationGeneral(null, null, null, null, motif, "Anywhere.", null, double.NaN, null, null, null, null, null, null);
            IDictionary<int, List<ModificationGeneral>> modDict = new Dictionary<int, List<ModificationGeneral>>
            {
                {2, new List<ModificationGeneral> { mod } },
                {8, new List<ModificationGeneral> { mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict, isDecoy: true);
            DigestionParams digestionParams = new DigestionParams(
                "Custom Protease7",
                maxMissedCleavages: 0,
                minPeptideLength: 5,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var digestedList = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).ToList();
            var ok1 = digestedList[1];
            var ok2 = digestedList[3];

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.AllModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.ContainsKey(3));

            prot = new Protein("MNNNNKRRRRR", null, null, null, modDict);
            ok1 = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).First();
            ok2 = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).Last();

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
            var ye = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye1 = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 1, maxPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye2 = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 5, maxPeptideLength: 8, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye3 = prot.Digest(digestionParams, new List<ModificationGeneral>(), new List<ModificationGeneral>()).ToList();
            Assert.AreEqual(3, ye.Count);
            Assert.AreEqual(2, ye1.Count);
            Assert.AreEqual(2, ye2.Count);
            Assert.AreEqual(1, ye3.Count);
        }
    }
}