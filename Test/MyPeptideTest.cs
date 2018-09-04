using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
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
            var protease = new Protease("CustomizedProtease", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications = new List<Modification>();
            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();

            Assert.AreEqual(2, ye.Count);

            var pep1 = ye[0];
            Assert.IsTrue(pep1.MonoisotopicMass > 0);

            var test = pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            foreach (var huh in pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both))
            {
                Assert.IsTrue(huh.NeutralMass > 0);
            }
            var pep2 = ye[1];
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.Fragment(DissociationType.HCD, FragmentationTerminus.Both))
            {
                Assert.IsTrue(huh.NeutralMass > 0);
            }
        }

        [Test]
        public static void TestNoCleavage()
        {
            List<Modification> fixedModifications = new List<Modification>();
            var prot = new Protein("MNNNKQQQQ", null, null, null, new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(5, 6, "lala") });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            var ye = prot.Digest(digestionParams, fixedModifications, new List<Modification>()).ToList();
            Assert.AreEqual(3, ye.Count);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null);
            var protease = new Protease("Custom Protease7", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            Assert.AreEqual(2, ye.Count);
            var pep1 = ye[0];
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both))
            {
                Assert.IsTrue(huh.NeutralMass > 0);
            }

            var pep2 = ye[1];
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToArray();
            Assert.IsTrue(cool[0].NeutralMass > 0);
            Assert.IsTrue(cool[1].NeutralMass > 0);
            Assert.IsTrue(cool[3].NeutralMass > 0);
            Assert.IsTrue(double.IsNaN(cool[2].NeutralMass));
            Assert.IsTrue(double.IsNaN(cool[4].NeutralMass));
            Assert.IsTrue(double.IsNaN(cool[5].NeutralMass));
            Assert.IsTrue(cool.Length == 6);
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3); // if you pass Custom Protease7 this test gets really flakey.
            List<Modification> variableModifications = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);

            variableModifications.Add(new Modification(_originalId: "ProtNmod", _target: motif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "pepNmod", _target: motif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "resMod", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "PepCmod", _target: motif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "ProtCmod", _target: motif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(3 * 2 * 3, ye.Count);
            Assert.AreEqual("[H]M[H][H]", ye.Last().SequenceWithChemicalFormulas);

            double m1 = 5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass;

            m1 = Math.Round((double)m1, 9, MidpointRounding.AwayFromZero);

            double m2 = ye.Last().MonoisotopicMass;
            double m3 = m1 - m2;

            Assert.IsTrue(m3 < 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3); // if you pass Custom Protease7 this test gets really flakey.
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);

            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: motif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: motif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: motif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: motif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            var ok = prot.Digest(digestionParams, fixedMods, new List<Modification>()).ToList();

            Assert.AreEqual(1, ok.Count);

            Assert.AreEqual("[:pepNmod on M]M[:resMod on M][:ProtCmod on M]", ok.First().FullSequence);

            Assert.AreEqual("[H]M[H][H]", ok.First().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestDigestIndices()
        {
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("R", out ModificationMotif motifR);
            Modification modN = new Modification("myMod", null, "myModType", null, motifN, "Anywhere.", null, 10, null, null, null, null, null, null);
            Modification modR = new Modification("myMod", null, "myModType", null, motifR, "Anywhere.", null, 10, null, null, null, null, null, null);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> { modN } },
                {8, new List<Modification> { modR } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict, isDecoy: true);

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var digestedList = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            var ok1 = digestedList[1];
            var ok2 = digestedList[3];

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.AllModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.ContainsKey(3));
        }

        [Test]
        public static void TestDigestDecoy()
        {
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("R", out ModificationMotif motifR);
            Modification modN = new Modification("myMod", null, "myModType", null, motifN, "Anywhere.", null, 10, null, null, null, null, null, null);
            Modification modR = new Modification("myMod", null, "myModType", null, motifR, "Anywhere.", null, 10, null, null, null, null, null, null);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> { modN } },
                {8, new List<Modification> { modR } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict, isDecoy: true);

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var digestedList = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            var ok1 = digestedList[1];
            var ok2 = digestedList[3];

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.AllModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.ContainsKey(3));

            prot = new Protein("MNNNNKRRRRR", null, null, null, modDict);
            ok1 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            ok2 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).Last();

            Assert.AreEqual(0, ok1.NumMods);
            Assert.IsFalse(ok1.AllModsOneIsNterminus.Any());
            Assert.AreEqual(2, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.Any());
        }

        [Test]
        public static void TestGoodPeptideWithLength()
        {
            var prot = new Protein("MNNNKQQQQMNNNKQQQQ", null);

            DigestionParams digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye1 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 1, maxPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye2 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("Custom Protease7", maxMissedCleavages: 0, minPeptideLength: 5, maxPeptideLength: 8, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye3 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            Assert.AreEqual(3, ye.Count);
            Assert.AreEqual(2, ye1.Count);
            Assert.AreEqual(2, ye2.Count);
            Assert.AreEqual(1, ye3.Count);
        }

        [Test]
        public static void Test_ProteinDigest()
        {
            DigestionParams d = new DigestionParams(
                        maxMissedCleavages: 0,
                        minPeptideLength: 5,
                        initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);

            IDictionary<int, List<Modification>> oneBasedModification = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification>{ mod } }
            };

            Protein prot1 = new Protein("MEDEEK", "prot1", oneBasedModifications: oneBasedModification);

            var pep1 = prot1.Digest(d, new List<Modification>(), new List<Modification>()).First();
            var pep2 = prot1.Digest(d, new List<Modification>(), new List<Modification>()).Last();

            Assert.AreEqual("MEDEEK", pep1.FullSequence);
            Assert.AreEqual("MED[mt:mod1 on D]EEK", pep2.FullSequence);
        }

        [Test]
        /// <summary>
        /// Tests that a PeptideWithSetModifications object can be parsed correctly from a string, with mod info
        /// </summary>
        public static void TestReadPeptideFromString()
        {
            // set up the test

            ModificationMotif.TryGetMotif("T", out ModificationMotif target);

            Modification carbamidomethylOnC = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: target, _chemicalFormula: ChemicalFormula.ParseFormula("C2H3NO"));
            string sequence = "HQVC[Common Fixed:Carbamidomethyl on C]TPGGTTIAGLC[Common Fixed:Carbamidomethyl on C]VMEEK";

            // parse the peptide from the string
            PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification> { { carbamidomethylOnC.IdWithMotif, carbamidomethylOnC } });

            // test base sequence and full sequence
            Assert.That(peptide.BaseSequence == "HQVCTPGGTTIAGLCVMEEK");
            Assert.That(peptide.FullSequence == sequence);

            // test peptide mass
            Assert.That(Math.Round(peptide.MonoisotopicMass, 5) == 2187.01225);

            // test mods (correct id, correct number of mods, correct location of mods)
            Assert.That(peptide.AllModsOneIsNterminus.First().Value.IdWithMotif == "Carbamidomethyl on C");
            Assert.That(peptide.AllModsOneIsNterminus.Count == 2);
            Assert.That(new HashSet<int>(peptide.AllModsOneIsNterminus.Keys).SetEquals(new HashSet<int>() { 5, 16 }));

            // calculate fragments. just check that they exist and it doesn't crash
            List<Product> theoreticalFragments = peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            Assert.That(theoreticalFragments.Count > 0);
        }
    }
}