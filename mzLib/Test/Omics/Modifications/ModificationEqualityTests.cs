// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ModificationEqualityTests.cs) is part of Proteomics.
//
// Proteomics is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Proteomics is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Proteomics. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.Modifications.IO;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class ModificationEqualityTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
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
        public static void Test_modificationsHashCode()
        {
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var mod1 = new Modification("mod", null, "type", null, motif, "Anywhere.", null, 1, null, null, null, null, null, null);
            var mod2 = new Modification("mod2", null, "type", null, motif, "Anywhere.", null, 10, null, null, null, null, null, null);

            Assert.AreNotEqual(mod1.GetHashCode(), mod2.GetHashCode());
            Assert.AreNotEqual(mod1, mod2);
            HashSet<Modification> myHashSet = new HashSet<Modification>
            {
                mod1,
                mod2
            };
            Assert.AreEqual(2, myHashSet.Count);
        }

        [Test]
        public static void Test_ModificationWithNoMassWritten()
        {
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var mod1 = new Modification(_originalId: "mod of M", _modificationType: "type", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: ChemicalFormula.ParseFormula("H").MonoisotopicMass);
            var mod1string = mod1.ToString();
            Assert.IsTrue(mod1string.Contains("MM"));
            var modAfterWriteRead = ModificationLoader.ReadModsFromString(mod1string + Environment.NewLine + "//", out var errors).First() as Modification;

            Assert.IsTrue(modAfterWriteRead.Equals(mod1));
        }

        [Test]
        public static void NameAndSites()
        {
            // Empty modification, has no name and by default has Sites = ModificationSites.Any
            OldSchoolModification a = new OldSchoolModification();
            OldSchoolModification b = new OldSchoolModification(a);
            Assert.AreEqual(" (Any)", b.NameAndSites);
        }

        [Test]
        public static void ModificationEquality()
        {
            // Empty modification, has no name and by default has Sites = ModificationSites.Any
            OldSchoolModification a = new OldSchoolModification();
            OldSchoolModification b = new OldSchoolModification();
            OldSchoolModification c = new OldSchoolModification(0, "c");
            OldSchoolModification d = new OldSchoolModification(0, "", ModificationSites.E);
            Assert.IsTrue(a.Equals(b));
            Assert.IsFalse(a.Equals(c));
            Assert.IsFalse(a.Equals(d));
        }

        [Test]
        public static void ModificationSitesTest()
        {
            // Empty modification, has no name and by default has Sites = ModificationSites.Any
            var a = ModificationSites.A | ModificationSites.E;
            Assert.AreEqual(ModificationSites.A | ModificationSites.E, a);
        }

        [Test]
        public static void Sites()
        {
            // Empty modification, has no name and by default has Sites = ModificationSites.Any
            var a = ModificationSites.A | ModificationSites.C | ModificationSites.E;
            Assert.IsTrue(a.ContainsSites(ModificationSites.E));

            Assert.IsTrue(a.ContainsSites(ModificationSites.A | ModificationSites.C));
            Assert.IsFalse(a.ContainsSites(ModificationSites.N));
            Assert.IsFalse(a.ContainsSites(ModificationSites.N | ModificationSites.C));
            var b = a.EnumerateActiveSites();
            Assert.IsTrue(b.Count() == 3);
        }

        [Test]
        public static void SilacLabelTest()
        {
            Protein originalProtein = new Protein("ACDEFGHIKAKAK", "TEST");
            List<PeptideWithSetModifications> originalDigest = originalProtein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).ToList();

            //Multiple SILAC labels
            Residue lysine = Residue.GetResidue('K');
            Residue arginine = Residue.GetResidue('R');
            Residue heavyLabel = new Residue("heavy", 'a', "aaa", ChemicalFormula.ParseFormula("C{13}6H12N2O"), ModificationSites.All); //Lysine +6
            Residue heavierLabel = new Residue("heavier", 'b', "bbb", ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //Lysine +8
            Residue heavyArginine = new Residue("heavyR", 'c', "ccc", ChemicalFormula.ParseFormula("C{13}6H5N{15}4O"), ModificationSites.All); //Arginine +10
            List<Residue> residuesToAdd = new List<Residue> { heavyLabel, heavierLabel, heavyArginine };
            Residue.AddNewResiduesToDictionary(residuesToAdd);
            List<SilacLabel> silacLabels = new List<SilacLabel>
            {
                new SilacLabel('K','a', heavyLabel.ThisChemicalFormula.Formula, heavyLabel.MonoisotopicMass - lysine.MonoisotopicMass),
                new SilacLabel('K','b', heavierLabel.ThisChemicalFormula.Formula, heavierLabel.MonoisotopicMass - lysine.MonoisotopicMass)
            };
            List<PeptideWithSetModifications> silacDigest = originalProtein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>(), silacLabels).ToList();
            Assert.IsTrue(originalDigest.Count * 3 == silacDigest.Count); //check that each peptide now has a light, heavy, and heavier compliment

            double silacPeptideLightMonoisotopicMass = silacDigest.Where(x => x.BaseSequence.Contains("K")).First().MonoisotopicMass;
            double silacPeptideHeavyMonoisotopicMass = silacDigest.Where(x => x.BaseSequence.Contains("a")).First().MonoisotopicMass;
            double silacPeptideHeavierMonoisotopicMass = silacDigest.Where(x => x.BaseSequence.Contains("b")).First().MonoisotopicMass;
            Assert.IsTrue(silacPeptideLightMonoisotopicMass != double.NaN); //if both NaN, then the mass comparison statements will always be true, because NaN + double = NaN
            Assert.IsTrue(silacPeptideHeavyMonoisotopicMass != double.NaN); //if both NaN, then the mass comparison statements will always be true, because NaN + double = NaN
            Assert.IsTrue(silacPeptideHeavierMonoisotopicMass != double.NaN); //if both NaN, then the mass comparison statements will always be true, because NaN + double = NaN

            Assert.IsTrue(Math.Round(silacPeptideLightMonoisotopicMass + heavyLabel.MonoisotopicMass, 5).Equals(Math.Round(silacPeptideHeavyMonoisotopicMass + lysine.MonoisotopicMass, 5))); //check that the residue masses were succesfully added
            Assert.IsTrue(Math.Round(silacPeptideHeavyMonoisotopicMass + heavierLabel.MonoisotopicMass, 5).Equals(Math.Round(silacPeptideHeavierMonoisotopicMass + heavyLabel.MonoisotopicMass, 5))); //check that the residue masses were succesfully added

            //code coverage
            SilacLabel testParameterlessConstructorForTomlsWithoutAnyRealTestAndMoreJustForCodeCoverage = new SilacLabel();
            Assert.IsTrue(testParameterlessConstructorForTomlsWithoutAnyRealTestAndMoreJustForCodeCoverage != null);

            Assert.IsTrue(silacLabels[0].AdditionalLabels == null);
            //The zero index label is K to a with a mass diff of 6
            //the one index label is K to b with a mass diff of 8
            silacLabels[0].AddAdditionalSilacLabel(new SilacLabel('D', 'c', heavyLabel.ThisChemicalFormula.Formula, heavyLabel.MonoisotopicMass - lysine.MonoisotopicMass - 1));
            Assert.IsTrue(silacLabels[0].AdditionalLabels.Count == 1);
            silacLabels[0].AddAdditionalSilacLabel(new SilacLabel('E', 'd', heavyLabel.ThisChemicalFormula.Formula, heavyLabel.MonoisotopicMass - lysine.MonoisotopicMass + 1));
            Assert.IsTrue(silacLabels[0].AdditionalLabels.Count == 2);
            silacDigest = originalProtein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>(), silacLabels).ToList();
            Assert.IsTrue(silacDigest.Count == 9);

            //Turnover
            silacLabels = new List<SilacLabel>
            {
                new SilacLabel('K','b', heavierLabel.ThisChemicalFormula.Formula, heavierLabel.MonoisotopicMass - lysine.MonoisotopicMass)
            };
            silacDigest = originalProtein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>(), silacLabels, (null, silacLabels[0])).ToList();
            Assert.IsTrue(silacDigest.Count == 14);
            silacDigest = originalProtein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>(), silacLabels, (silacLabels[0], null)).ToList();
            Assert.IsTrue(silacDigest.Count == 14);

            silacLabels = new List<SilacLabel>
            {
                new SilacLabel('K','a', heavyLabel.ThisChemicalFormula.Formula, heavyLabel.MonoisotopicMass - lysine.MonoisotopicMass),
                new SilacLabel('K','b', heavierLabel.ThisChemicalFormula.Formula, heavierLabel.MonoisotopicMass - lysine.MonoisotopicMass)
            };
            silacDigest = originalProtein.Digest(new DigestionParams(generateUnlabeledProteinsForSilac: false), new List<Modification>(), new List<Modification>(), silacLabels, (silacLabels[1], silacLabels[0])).ToList();
            Assert.IsTrue(silacDigest.Count == 14);

            originalProtein = new Protein("ACDEFGHIKARAK", "TEST");
            silacLabels[1].AddAdditionalSilacLabel(new SilacLabel('R', 'c', heavyArginine.ThisChemicalFormula.Formula, heavyArginine.MonoisotopicMass - arginine.MonoisotopicMass));
            silacDigest = originalProtein.Digest(new DigestionParams(generateUnlabeledProteinsForSilac: false), new List<Modification>(), new List<Modification>(), silacLabels, (silacLabels[1], silacLabels[0])).ToList();
            Assert.IsTrue(silacDigest.Count == 14);

            silacLabels = new List<SilacLabel>
            {
                new SilacLabel('K','a', heavyLabel.ThisChemicalFormula.Formula, heavyLabel.MonoisotopicMass - lysine.MonoisotopicMass),
                new SilacLabel('K','b', heavierLabel.ThisChemicalFormula.Formula, heavierLabel.MonoisotopicMass - lysine.MonoisotopicMass)
            };
            silacLabels[0].AddAdditionalSilacLabel(new SilacLabel('R', 'c', heavyArginine.ThisChemicalFormula.Formula, heavyArginine.MonoisotopicMass - arginine.MonoisotopicMass));
            silacDigest = originalProtein.Digest(new DigestionParams(generateUnlabeledProteinsForSilac: false), new List<Modification>(), new List<Modification>(), silacLabels, (silacLabels[1], silacLabels[0])).ToList();
            Assert.IsTrue(silacDigest.Count == 14);
        }

        [Test]
        public static void ModificationCollectionTest()
        {
            ModificationCollection a = new ModificationCollection(new OldSchoolModification(1, "Mod1"), new OldSchoolModification(2, "Mod2"));

            double lala = 0;
            IEnumerable aasdf = a;
            foreach (var jadfk in aasdf)
                lala += (jadfk as IHasMass).MonoisotopicMass;
            Assert.AreEqual(3, lala);

            Assert.AreEqual("Mod1 | Mod2", a.ToString());
            a.Add(new OldSchoolModification(3, "Mod3"));
            Assert.AreEqual("Mod1 | Mod2 | Mod3", a.ToString());
            Assert.IsTrue(a.Contains(new OldSchoolModification(2, "Mod2")));
            IHasMass[] myArray = new IHasMass[4];
            a.CopyTo(myArray, 1);
            Assert.AreEqual(3, myArray.Sum(b => b == null ? 0 : 1));
            Assert.AreEqual(3, a.Count());
            Assert.IsFalse(a.IsReadOnly);
            a.Remove(new OldSchoolModification(2, "Mod2"));
            Assert.AreEqual("Mod1 | Mod3", a.ToString());
            double ok = 0;
            foreach (var b in a)
                ok += b.MonoisotopicMass;
            Assert.AreEqual(4, ok);

            a.Clear();
            Assert.AreEqual("", a.ToString());
        }

        [Test]
        public static void ModificationCollectionTest2()
        {
            ModificationCollection a = new ModificationCollection(new OldSchoolModification(1, "Mod1"), new OldSchoolModification(2, "Mod2"));
            Assert.IsFalse(a.Remove(new OldSchoolModification(3, "Mod3")));
        }

        [Test]
        public static void ModificationWithMultiplePossibilitiesTest()
        {
            var m = new ModificationWithMultiplePossibilitiesCollection("My Iso Mod", ModificationSites.E);
            m.AddModification(new OldSchoolModification(1, "My Mod1a", ModificationSites.E));
            m.AddModification(new OldSchoolModification(2, "My Mod2b", ModificationSites.E));
            Assert.AreEqual(2, m.Count);
            Assert.AreEqual("My Mod2b", m[1].Name);
            Assert.Throws<MzLibException>(() => { m.AddModification(new OldSchoolModification(1, "gg", ModificationSites.R)); }, "Unable to add a modification with sites other than ModificationSites.E");
            Assert.IsTrue(m.Contains(new OldSchoolModification(2, "My Mod2b", ModificationSites.E)));
            double kk = 0;
            IEnumerable a = m;
            foreach (var b in a)
                kk += (b as OldSchoolModification).MonoisotopicMass;
            Assert.AreEqual(3, kk);
        }

        [Test]
        public static void ModificationCollectionTestTest()
        {
            OldSchoolModification mod1 = new OldSchoolModification(10, "mass 10 modification");
            OldSchoolModification mod2 = new OldSchoolModification(100, "mass 100 modification");
            OldSchoolModification mod3 = new OldSchoolModification(1000, "mass 1000 modification");
            ModificationCollection a = new ModificationCollection(mod1, mod2, mod3, mod1);
            ModificationCollection b = new ModificationCollection(mod1, mod3, mod1, mod2);
            Assert.IsTrue(a.Equals(b));
            ModificationCollection c = new ModificationCollection(mod1);
            Assert.IsFalse(c.Equals(b));
        }

        [Test]
        public static void ModificationSitesTest55()
        {
            Assert.IsTrue(ModificationSites.E.ContainsSites(ModificationSites.Any));
            Assert.IsFalse(ModificationSites.E.ContainsSites(ModificationSites.None));
            Assert.IsTrue(ModificationSites.None.ContainsSites(ModificationSites.None));
        }

        [Test]
        public static void ChemicalFormulaModification()
        {
            OldSchoolChemicalFormulaModification a = new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("OH"));
            OldSchoolChemicalFormulaModification b = new OldSchoolChemicalFormulaModification(a);
            Assert.AreEqual(a, b);
        }

        [Test]
        public static void ModificationCollectionScrambledEquals()
        {
            ModificationCollection a = new ModificationCollection(new OldSchoolModification(1, "Mod1"), new OldSchoolModification(2, "Mod2"));
            ModificationCollection b = new ModificationCollection(new OldSchoolModification(1, "Mod1"), new OldSchoolModification(3, "Mod3"));

            Assert.IsFalse(a.Equals(b));
        }

        [Test]
        public static void Test_modification_hash_set()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m1 = new Modification("23", null, "unknown", null, motif, null, null, 42.01, null, null, null, null, null, null);
            Modification m2 = new Modification("23", null, "unknown", null, motif, null, null, 42.01, null, null, null, null, null, null);
            HashSet<Modification> mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.AreEqual(1, mods.Count);
        }

        [Test]
        public static void Test_modification2_hash_set()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m1 = new Modification("id1", null, "modificationType", null, motif, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            Modification m2 = new Modification("id1", null, "modificationType", null, motif, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            m1.DatabaseReference.Add("key", new List<string> { "value" });
            m2.DatabaseReference.Add("key", new List<string> { "value" });
            HashSet<Modification> mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.True(m1.Equals(m2));
            Assert.AreEqual(1, mods.Count);
        }

        [Test]
        public static void Test_modification3_hash_set() // numerical tolerance is 1e-9 so these two mods need to evaluate as identical
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m1 = new Modification(_originalId: "id1", _modificationType: "modificationType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 1.11111d, _databaseReference: new Dictionary<string, IList<string>>(), _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.AnyActivationType, new List<double> { 2.222222 } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { DissociationType.AnyActivationType, new List<double> { 1.2233 } } });
            Modification m2 = new Modification(_originalId: "id1", _modificationType: "modificationType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 1.11111d - 1e-10, _databaseReference: new Dictionary<string, IList<string>>(), _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.AnyActivationType, new List<double> { 2.222222 + 1e-10 } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { DissociationType.AnyActivationType, new List<double> { 1.2233 } } });
            m1.DatabaseReference.Add("key", new List<string> { "value" });
            m2.DatabaseReference.Add("key", new List<string> { "value" });

            HashSet<Modification> mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.AreEqual(1, mods.Count);
            Assert.IsTrue(m1.Equals(m2));
        }

        [Test]
        public static void TestInvalidModificationHash()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m1 = new Modification(_originalId: "id1", _modificationType: "modificationType", _target: motif, _locationRestriction: "Anywhere.");
            Modification m2 = new Modification(_originalId: "id1", _modificationType: "modificationType", _target: motif, _locationRestriction: "Anywhere.");
            HashSet<Modification> mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.IsFalse(m1.ValidModification);
            Assert.IsFalse(m2.ValidModification);
            Assert.True(m1.Equals(m2));
            Assert.AreEqual(1, mods.Count);

            // test comparing invalid mods with null vs not-null MMs
            m1 = new Modification(_originalId: "id1", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 1);
            m2 = new Modification(_originalId: "id1", _target: motif, _locationRestriction: "Anywhere.");
            mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.IsFalse(m1.ValidModification);
            Assert.IsFalse(m2.ValidModification);
            Assert.False(m1.Equals(m2));
            Assert.AreEqual(2, mods.Count);

            // test comparing invalid mods with null vs not-null IDs
            m1 = new Modification(_originalId: "id1", _target: motif, _locationRestriction: "Anywhere.");
            m2 = new Modification(_target: motif, _locationRestriction: "Anywhere.");
            mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.IsFalse(m1.ValidModification);
            Assert.IsFalse(m2.ValidModification);
            Assert.False(m1.Equals(m2));
            Assert.AreEqual(2, mods.Count);
        }

        [Test]
        public static void TestUniprotNTerminalMod()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif1);
            Modification variableMod = new Modification(_originalId: "acetylation", _modificationType: "Variable", _target: motif1, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "N-terminal.");

            ModificationMotif.TryGetMotif("P", out ModificationMotif motif2);
            Modification uniprotMod = new Modification(_originalId: "acetylation", _modificationType: "UniProt", _target: motif2, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "N-terminal.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 1, new List<Modification> { uniprotMod } } };

            Protein protein = new Protein("PEPTIDE", "", oneBasedModifications: mods);
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>() { variableMod }).Where(p => p.AllModsOneIsNterminus.Count == 1).First();

            Assert.That(peptide.FullSequence == "[UniProt:acetylation on P]PEPTIDE");
        }

        [Test]
        public static void TestUniprotCTerminalMod()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif1);
            Modification variableMod = new Modification(_originalId: "acetylation", _modificationType: "Variable", _target: motif1, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "C-terminal.");

            ModificationMotif.TryGetMotif("E", out ModificationMotif motif2);
            Modification uniprotMod = new Modification(_originalId: "acetylation", _modificationType: "UniProt", _target: motif2, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "C-terminal.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 7, new List<Modification> { uniprotMod } } };

            Protein protein = new Protein("PEPTIDE", "", oneBasedModifications: mods);
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>() { variableMod }).Where(p => p.AllModsOneIsNterminus.Count == 1).First();

            Assert.That(peptide.FullSequence == "PEPTIDE-[UniProt:acetylation on E]");
        }

        [Test]
        public static void TestUniprotResidualMod()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif1);
            Modification variableMod = new Modification(_originalId: "acetylation", _modificationType: "Variable", _target: motif1, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "Anywhere.");

            ModificationMotif.TryGetMotif("T", out ModificationMotif motif2);
            Modification uniprotMod = new Modification(_originalId: "acetylation", _modificationType: "UniProt", _target: motif2, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "Anywhere.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { uniprotMod } } };

            Protein protein = new Protein("PEPTIDE", "", oneBasedModifications: mods);
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>() { variableMod }).Where(p => p.AllModsOneIsNterminus.Count == 1).First();

            Assert.That(peptide.FullSequence == "PEPT[UniProt:acetylation on T]IDE");
        }

        [Test]
        public static void CompareTo_SameModification_ReturnsZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.EqualTo(0));
        }

        [Test]
        public static void CompareTo_DifferentIdWithMotif_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod2", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
        }

        [Test]
        public static void CompareTo_DifferentModificationType_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type2", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
        }

        [Test]
        public static void CompareTo_DifferentTarget_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif1);
            ModificationMotif.TryGetMotif("B", out var motif2);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif1, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type1", "feature1", motif2, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
        }

        [Test]
        public static void CompareTo_DifferentLocationRestriction_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "C-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
        }

        [Test]
        public static void CompareTo_DifferentMonoisotopicMass_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 101.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(null), Is.EqualTo(1));
        }
    }
}
