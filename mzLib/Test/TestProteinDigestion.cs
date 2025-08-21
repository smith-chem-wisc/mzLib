using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using UsefulProteomicsDatabases;
using static Chemistry.PeriodicTable;
using Stopwatch = System.Diagnostics.Stopwatch;
using MzLibUtil;
using Omics.BioPolymer;
using Omics;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestProteinDigestion
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
        public static void ProteaseLoader()
        {
            string path1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteaseFilesForLoadingTests", "TestProteases_badMod.tsv");
            string path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteaseFilesForLoadingTests", "TestProteases_badMod_dupName.tsv");
            string path3 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteaseFilesForLoadingTests", "TestProteases_dupName.tsv");
            string path4 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteaseFilesForLoadingTests", "TestProteases_Mod_dupName.tsv");
            var proteaseMods = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "ProteaseMods.txt"), out var errors).ToList();
            
            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path1, proteaseMods)); 
            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path2, proteaseMods));
            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path3, proteaseMods));
            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path4, proteaseMods));
        }

        [Test]
        public static void CNBrProteinDigestion()
        {
            var proteaseMods = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "ProteaseMods.txt"), out var errors).ToList();
            var prot = new Protein("PEPTIDEMPEPTIDEM", null);
            var prot2 = new Protein("MPEPTIDEMPEPTIDE", null);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DoubleProtease.tsv");
            NUnit.Framework.Assert.That(File.Exists(path));

            var proteaseDict = ProteaseDictionary.LoadProteaseDictionary(path, proteaseMods);
            ProteaseDictionary.Dictionary = ProteaseDictionary.LoadProteaseDictionary(Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "ProteolyticDigestion", "proteases.tsv"), proteaseMods);
            var protease1 = proteaseDict["CNBr"];
            DigestionParams digestionParams1 = new DigestionParams(
                protease: protease1.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);            
            List<Modification> variableModifications1 = new List<Modification>();

            var protease2 = proteaseDict["CNBr_old"];
            DigestionParams digestionParams2 = new DigestionParams(
                protease: protease2.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications2 = new List<Modification>();

            var protease3 = proteaseDict["CNBr_N"];
            DigestionParams digestionParams3 = new DigestionParams(
                protease: protease3.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications3 = new List<Modification>();

            var peps1 = prot.Digest(digestionParams1, new List<Modification>(), variableModifications1).ToList();
            var peps2 = prot.Digest(digestionParams2, new List<Modification>(), variableModifications2).ToList(); 
            var peps3 = prot2.Digest(digestionParams3, new List<Modification>(), variableModifications1).ToList();

            Assert.AreNotEqual(null, protease3.CleavageMod);
            Assert.AreEqual("M", protease3.CleavageMod.Target.ToString());


            Assert.AreNotEqual(peps3[0].MonoisotopicMass, peps3[1].MonoisotopicMass);

            Assert.AreEqual(882.39707781799996, peps3[1].MonoisotopicMass);
            Assert.AreEqual(930.400449121, peps3[0].MonoisotopicMass);


            Assert.AreEqual(null, protease2.CleavageMod);
            Assert.AreNotEqual(null, protease1.CleavageMod);
            Assert.AreEqual("M", protease1.CleavageMod.Target.ToString());

            Assert.AreEqual(peps1[1].MonoisotopicMass, peps2[1].MonoisotopicMass);
            Assert.AreEqual(peps1[1].MonoisotopicMass, peps2[0].MonoisotopicMass);
            Assert.AreEqual(peps2[0].MonoisotopicMass, peps2[1].MonoisotopicMass);
            Assert.AreNotEqual(peps1[0].MonoisotopicMass, peps1[1].MonoisotopicMass);
            Assert.AreNotEqual(peps1[0].MonoisotopicMass, peps2[0].MonoisotopicMass);
            Assert.AreNotEqual(peps1[0].MonoisotopicMass, peps2[1].MonoisotopicMass);

            Assert.AreEqual(882.39707781799996, peps1[0].MonoisotopicMass);
            Assert.AreEqual(930.400449121, peps1[1].MonoisotopicMass);

            
        }
       
        [Test]
        public static void TestGoodPeptide()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            var motifList = DigestionMotif.ParseDigestionMotifsFromString("K|");
            var protease = new Protease("CustomizedProtease", CleavageSpecificity.Full, null, null, motifList);
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

            var test = new List<Product>();
            pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both, test);

            foreach (var huh in test)
            {
                Assert.IsTrue(huh.NeutralMass > 0);
            }
            var pep2 = ye[1];
            pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both, test);
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in test)
            {
                Assert.IsTrue(huh.NeutralMass > 0);
            }
        }

        [Test]
        public static void TestNoCleavage()
        {
            List<Modification> fixedModifications = new List<Modification>();
            var prot = new Protein("MNNNKQQQQ", null, null, null, new Dictionary<int, List<Modification>>(), new List<TruncationProduct> { new TruncationProduct(5, 6, "lala") });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            var ye = prot.Digest(digestionParams, fixedModifications, new List<Modification>()).ToList();
            Assert.AreEqual(3, ye.Count);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null);
            var motifList = DigestionMotif.ParseDigestionMotifsFromString("K|");            
            var protease = new Protease("Custom Protease7", CleavageSpecificity.Full, null, null, motifList);
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

            var fragments = new List<Product>();
            pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);
            foreach (var huh in fragments)
            {
                Assert.IsTrue(huh.NeutralMass > 0);
            }

            var pep2 = ye[1];
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = new List<Product>();
            pep2.Fragment(DissociationType.HCD, FragmentationTerminus.Both, cool);
            Assert.IsTrue(cool[0].NeutralMass > 0);
            Assert.IsTrue(cool[1].NeutralMass > 0);
            Assert.IsTrue(cool[2].NeutralMass > 0);
            Assert.IsTrue(double.IsNaN(cool[3].NeutralMass));
            Assert.IsTrue(double.IsNaN(cool[4].NeutralMass));
            Assert.IsTrue(double.IsNaN(cool[5].NeutralMass));
            Assert.IsTrue(cool.Count == 6);
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
            Assert.AreEqual("[H]M[H]-[H]", ye.Last().SequenceWithChemicalFormulas);

            double m1 = 5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass;
            m1 = Math.Round(m1, 9, MidpointRounding.AwayFromZero);
            double m2 = ye.Last().MonoisotopicMass;
            double m3 = m1 - m2;
            Assert.IsTrue(m3 < 1e-9);
        }

        [Test]
        public static void TestPeptideDigestion_FixedModifications_ProtModsOverwritePepMods()
        {
            string baseSequence = "M";
            var prot = new Protein(baseSequence, null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3); // if you pass Custom Protease7 this test gets really flakey.
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: motif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: motif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: motif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: motif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            var ok = prot.Digest(digestionParams, fixedMods, new List<Modification>()).ToList();

            Assert.AreEqual(1, ok.Count);
            
            string expectedFullSequence = "[:ProtNmod on M]M[:resMod on M]-[:ProtCmod on M]";
            Assert.AreEqual(expectedFullSequence, ok.First().FullSequence);
            var mods = ok.First().AllModsOneIsNterminus;

            NUnit.Framework.Assert.That(IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(expectedFullSequence), Is.EqualTo(baseSequence));
            NUnit.Framework.Assert.That(IBioPolymerWithSetMods.DetermineFullSequence(baseSequence, mods), Is.EqualTo(expectedFullSequence));
            NUnit.Framework.Assert.That(ok.First().DetermineFullSequence(), Is.EqualTo(expectedFullSequence));

            Assert.AreEqual("[H]M[H]-[H]", ok.First().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideDigestion_FixedModifications_ProtModsOverwritePepMods_RandomizedModOrder()
        {
            var rand = new Random(42);
            string baseSequence = "M";
            var prot = new Protein(baseSequence, null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3); // if you pass Custom Protease7 this test gets really flakey.
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: motif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: motif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: motif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: motif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            // set expected values
            int expectedDigestionProducts = 1;
            string expectedFullSequence = "[:ProtNmod on M]M[:resMod on M]-[:ProtCmod on M]";
            string expectedSequenceWithChemicalFormulas = "[H]M[H]-[H]";
            double expectedMonoisotopicMass = 5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass;

            

            // randomly scramble all mods, digest, and ensure the answer is correct. 
            for (int i = 0; i < 10; i++)
            {
                var shuffledFixedMods = fixedMods.OrderBy(a => rand.Next()).ToList();
                var ok = prot.Digest(digestionParams, shuffledFixedMods, new List<Modification>()).ToList();
                var mods = ok.First().AllModsOneIsNterminus;

                NUnit.Framework.Assert.That(IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(expectedFullSequence), Is.EqualTo(baseSequence));
                NUnit.Framework.Assert.That(IBioPolymerWithSetMods.DetermineFullSequence(baseSequence, mods), Is.EqualTo(expectedFullSequence));
                NUnit.Framework.Assert.That(ok.First().DetermineFullSequence(), Is.EqualTo(expectedFullSequence));

                Assert.AreEqual(expectedDigestionProducts, ok.Count);
                Assert.AreEqual(expectedFullSequence, ok.First().FullSequence);
                Assert.AreEqual(expectedSequenceWithChemicalFormulas, ok.First().SequenceWithChemicalFormulas);
                Assert.AreEqual(expectedMonoisotopicMass, ok.Last().MonoisotopicMass, 1e-9);
            }
        }

        [Test]
        public static void TestPeptideDigestion_FixedModifications_ProtModsOverwritePepMods_TwoProducts()
        {
            var prot = new Protein("MKM", null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            ModificationMotif.TryGetMotif("K", out ModificationMotif kMotif);

            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: mMotif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: kMotif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: mMotif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: kMotif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod", _target: mMotif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: mMotif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: kMotif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: mMotif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: kMotif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            var ok = prot.Digest(digestionParams, fixedMods, new List<Modification>()).ToList();

            Assert.AreEqual(2, ok.Count);

            Assert.AreEqual("[:ProtNmod on M]M[:resMod on M]K-[:PepCmod on K]", ok.First().FullSequence);
            Assert.AreEqual("[:pepNmod on M]M[:resMod on M]-[:ProtCmod on M]", ok.Skip(1).First().FullSequence);

            Assert.AreEqual("[H]M[H]K-[H]", ok.First().SequenceWithChemicalFormulas);
            Assert.AreEqual("[H]M[H]-[H]", ok.Skip(1).First().SequenceWithChemicalFormulas);
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

            DigestionParams digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 0, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye1 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 0, minPeptideLength: 1, maxPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye2 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 0, minPeptideLength: 5, maxPeptideLength: 8, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
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

        /// <summary>
        /// We want to have protein digestion yield the same set of peptides regardless of the order their modifications are encoded in the XML.
        /// While all of the positions of the modifications are the same, the order of the modifications in the XML is different.
        /// </summary>
        [Test]
        public static void TestDigestionOfSameProteinFromDifferentXmls()
        {
            DigestionParams digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 2, minPeptideLength: 7, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif);
            Modification carbamidomethylOnC = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("C2H3NO"));
            var fixedModifications = new List<Modification> { carbamidomethylOnC };
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            Modification oxidationOnM = new Modification(_originalId: "Oxidation on M", _modificationType: "Common Variable", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("O"));
            var variableModifications = new List<Modification> { oxidationOnM };

            // Load in proteins
            var dbFive = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SingleEntry_ModOrder1.xml");
            var dbSix = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SingleEntry_ModOrder2.xml");

            var proteins5 = ProteinDbLoader.LoadProteinXML(dbFive, true, DecoyType.None, null, false, null, out var unknownModificationsFive);
            var proteins6 = ProteinDbLoader.LoadProteinXML(dbSix, true, DecoyType.None, null, false, null, out var unknownModificationsSix);

            var fiveMods = ProteinDbLoader.GetPtmListFromProteinXml(dbFive);
            var sixMods = ProteinDbLoader.GetPtmListFromProteinXml(dbSix);

            Assert.AreEqual(fiveMods.Count, sixMods.Count);
            CollectionAssert.AreEquivalent(fiveMods, sixMods);

            var peptides5 = proteins5.First().Digest(digestionParams, fixedModifications, variableModifications).ToList();
            var peptides6 = proteins6.First().Digest(digestionParams, fixedModifications, variableModifications).ToList();
            Assert.AreEqual(peptides5.Count, peptides6.Count);
            CollectionAssert.AreEqual(peptides5, peptides6);
        }

        [Test]
        [TestCase("cRAP_databaseGPTMD.xml")]
        [TestCase("uniprot_aifm1.fasta")]
        public static void TestDecoyScramblingIsReproducible(string fileName)
        {
            // Load in proteins
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", fileName);
            DecoyType decoyType = DecoyType.Reverse;
            List<Protein> proteins1 = null;
            List<Protein> proteins2 = null;
            if (fileName.Contains(".xml"))
            {
                proteins1 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out unknownModifications);
            }
            else if (fileName.Contains(".fasta"))
            {
                proteins1 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out unknownModifications);
            }
            else
            {
                NUnit.Framework.Assert.Fail("Unknown file type");
            }

            DigestionParams d = new DigestionParams(
                        maxMissedCleavages: 1,
                        minPeptideLength: 5,
                        initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            // Digest target proteins
            var pepsToReplace = proteins1.Where(p => !p.IsDecoy)
                .SelectMany(p => p.Digest(d, new List<Modification>(), new List<Modification>()).ToList())
                .Select(pep => pep.BaseSequence)
                .ToHashSet();

            // Ensure at least one decoy peptide from each protein is problematic and must be replaced
            var singleDecoyPeptides = proteins1
                .Where(p => p.IsDecoy)
                .Select(p => p.Digest(d, new List<Modification>(), new List<Modification>()).Skip(2).Take(1))
                .Select(pwsm => pwsm.First().BaseSequence)
                .ToHashSet();

            //modify targetpeptides in place
            pepsToReplace.UnionWith(singleDecoyPeptides);

            // Scramble every decoy from db1
            List<Protein> decoys1 = new();
            foreach (var protein in proteins1.Where(p => p.IsDecoy))
            {
                decoys1.Add(DecoySequenceValidator.ScrambleDecoyBioPolymer(protein, d, pepsToReplace));
            }
            // Scramble every decoy from db2
            List<Protein> decoys2 = new();
            foreach (var protein in proteins2.Where(p => p.IsDecoy))
            {
                decoys2.Add(DecoySequenceValidator.ScrambleDecoyBioPolymer(protein, d, pepsToReplace));
            }

            // check are equivalent lists of proteins
            Assert.AreEqual(decoys1.Count, decoys2.Count);
            foreach (var decoyPair in decoys1.Concat(decoys2).GroupBy(p => p.Accession))
            {
                Assert.AreEqual(2, decoyPair.Count());
                Assert.AreEqual(decoyPair.First().BaseSequence, decoyPair.Last().BaseSequence);
            }
        }

        [Test]
        public static void TestDecoyScramblerReplacesPeptides()
        {
            DigestionParams d = new DigestionParams(
                        maxMissedCleavages: 1,
                        minPeptideLength: 5,
                        initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            Protein target = new Protein("MEDEEKFVGYKYGVFK", "target");
            Protein decoy = new Protein("EEDEMKYGVFKFVGYK", "decoy");

            var targetPep = target.Digest(d, new List<Modification>(), new List<Modification>());
            var decoyPep = decoy.Digest(d, new List<Modification>(), new List<Modification>());

            HashSet<string> targetPepSeqs = targetPep.Select(p => p.FullSequence).ToHashSet();
            var offendingDecoys = decoyPep.Where(p => targetPepSeqs.Contains(p.FullSequence)).Select(d => d.FullSequence).ToList();

            Assert.AreEqual(2, offendingDecoys.Count);

            Protein scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy,  d, targetPepSeqs, offendingDecoys);
            var scrambledPep = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>());

            Assert.AreEqual(decoyPep.Count(), scrambledPep.Count());
            Assert.IsFalse(scrambledPep.Any(p => offendingDecoys.Contains(p.FullSequence)));

            // Check to make sure that decoy generation also works in no offending sequences are passed in
            scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs);
            scrambledPep = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>());

            Assert.AreEqual(decoyPep.Count(), scrambledPep.Count());
            Assert.IsFalse(scrambledPep.Any(p => offendingDecoys.Contains(p.FullSequence)));
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
            NUnit.Framework.Assert.That(peptide.BaseSequence == "HQVCTPGGTTIAGLCVMEEK");
            NUnit.Framework.Assert.That(peptide.FullSequence == sequence);

            // test peptide mass
            NUnit.Framework.Assert.That(Math.Round(peptide.MonoisotopicMass, 5) == 2187.01225);

            // test mods (correct id, correct number of mods, correct location of mods)
            NUnit.Framework.Assert.That(peptide.AllModsOneIsNterminus.First().Value.IdWithMotif == "Carbamidomethyl on C");
            NUnit.Framework.Assert.That(peptide.AllModsOneIsNterminus.Count == 2);
            NUnit.Framework.Assert.That(new HashSet<int>(peptide.AllModsOneIsNterminus.Keys).SetEquals(new HashSet<int>() { 5, 16 }));

            // calculate fragments. just check that they exist and it doesn't crash
            List<Product> theoreticalFragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theoreticalFragments);
            NUnit.Framework.Assert.That(theoreticalFragments.Count > 0);
        }

        [Test]
        public static void TestGlycoPeptide()
        {
            var prot = new Protein("MNNNYTKQQQQKS", null);
            var motifList = DigestionMotif.ParseDigestionMotifsFromString("K|");
            var protease = new Protease("CustomizedProtease_diffname", CleavageSpecificity.Full, null, null, motifList);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true);
            List<Modification> variableModifications = new List<Modification>();
            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();

            Assert.AreEqual(2, ye.Count);
        }

        [Test]
        public static void TestDigestionParamsClone()
        {
            DigestionParams digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true);

            DigestionParams digestionParamsClone = (DigestionParams)digestionParams.Clone();
            Assert.AreEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(digestionParams.FragmentationTerminus, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true,
                maxModificationIsoforms: 5,
                maxModsForPeptides: 6,
                maxPeptideLength: 7,
                searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.C,
                generateUnlabeledProteinsForSilac: false);

            digestionParamsClone = (DigestionParams)digestionParams.Clone();
            Assert.AreEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(digestionParams.FragmentationTerminus, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));
        }



        [Test]
        public static void TestDigestionParamsCloneWithNewTerminus()
        {
            DigestionParams digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true);

            DigestionParams digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.N);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.N, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true,
                maxModificationIsoforms: 5,
                maxModsForPeptides: 6,
                maxPeptideLength: 7,
                searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.None,
                generateUnlabeledProteinsForSilac: false);

            digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.N);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(ProteaseDictionary.Dictionary["singleN"], digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.N, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.C);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(ProteaseDictionary.Dictionary["singleC"], digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.C, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));
        }

        [Test]
        public static void TestWhenFixedModIsSamePositionAsUniProtModWithDigestion()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3); // if you pass Custom Protease7 this test gets really flakey.
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("S", out ModificationMotif serineMotif);
            ChemicalFormula ohFormula = ChemicalFormula.ParseFormula("OH");
            double ohMass = GetElement("O").PrincipalIsotope.AtomicMass + GetElement("H").PrincipalIsotope.AtomicMass;

            fixedMods.Add(new Modification(_originalId: "serineOhMod", _target: serineMotif, _locationRestriction: "Anywhere.", _chemicalFormula: ohFormula, _monoisotopicMass: ohMass));


            List<Protein> dbProteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml"), true, DecoyType.Reverse, UniProtPtms.Concat(fixedMods), false,
                new List<string>(), out Dictionary<string, Modification> un);

            Protein prot = dbProteins.First();

            var digestionProducts = prot.Digest(digestionParams, fixedMods, new List<Modification>()).ToList();
            var firstPeptideModifiedForms = digestionProducts.Where(p => p.BaseSequence == "MSGR").ToList();
            List<string> fullSequences = firstPeptideModifiedForms.Select(p => p.FullSequence).ToList();
            List<string> expectedFullSequences = new List<string>
            {
                "MS[:serineOhMod on S]GR",
                "MS[:serineOhMod on S]GR[UniProt:Asymmetric dimethylarginine on R]",
                "MS[:serineOhMod on S]GR[UniProt:Citrulline on R]",
                "MS[:serineOhMod on S]GR[UniProt:Omega-N-methylarginine on R]",
                "MS[:serineOhMod on S]GR[UniProt:Symmetric dimethylarginine on R]",
                "MS[UniProt:Phosphoserine on S]GR",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Asymmetric dimethylarginine on R]",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Citrulline on R]",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Omega-N-methylarginine on R]",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Symmetric dimethylarginine on R]"
            };

            CollectionAssert.AreEquivalent(expectedFullSequences, fullSequences);
        }
    }
}