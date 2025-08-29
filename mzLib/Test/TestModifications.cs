// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (TestFragments.cs) is part of Proteomics.
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
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Omics.Fragmentation;
using Omics.Modifications;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestModifications
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
            var modAfterWriteRead = PtmListLoader.ReadModsFromString(mod1string + Environment.NewLine + "//", out var errors).First() as Modification;

            Assert.IsTrue(modAfterWriteRead.Equals(mod1));
        }

        [Test]
        public void NameAndSites()
        {
            // Empty modification, has no name and by default has Sites = ModificationSites.Any
            OldSchoolModification a = new OldSchoolModification();
            OldSchoolModification b = new OldSchoolModification(a);
            Assert.AreEqual(" (Any)", b.NameAndSites);
        }

        [Test]
        public void ModificationEquality()
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
        public void ModificationSitesTest()
        {
            // Empty modification, has no name and by default has Sites = ModificationSites.Any
            var a = ModificationSites.A | ModificationSites.E;
            Assert.AreEqual(ModificationSites.A | ModificationSites.E, a);
        }

        [Test]
        public void Sites()
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
        public void SilacLabelTest()
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
        public void ModificationCollectionTest()
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
        public void ModificationCollectionTest2()
        {
            ModificationCollection a = new ModificationCollection(new OldSchoolModification(1, "Mod1"), new OldSchoolModification(2, "Mod2"));
            Assert.IsFalse(a.Remove(new OldSchoolModification(3, "Mod3")));
        }

        [Test]
        public void ModificationWithMultiplePossibilitiesTest()
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
        public void ModificationCollectionTestTest()
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
        public void ModificationSitesTest55()
        {
            Assert.IsTrue(ModificationSites.E.ContainsSites(ModificationSites.Any));
            Assert.IsFalse(ModificationSites.E.ContainsSites(ModificationSites.None));
            Assert.IsTrue(ModificationSites.None.ContainsSites(ModificationSites.None));
        }

        [Test]
        public void ChemicalFormulaModification()
        {
            OldSchoolChemicalFormulaModification a = new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("OH"));
            OldSchoolChemicalFormulaModification b = new OldSchoolChemicalFormulaModification(a);
            Assert.AreEqual(a, b);
        }

        [Test]
        public void ModificationCollectionScrambledEquals()
        {
            ModificationCollection a = new ModificationCollection(new OldSchoolModification(1, "Mod1"), new OldSchoolModification(2, "Mod2"));
            ModificationCollection b = new ModificationCollection(new OldSchoolModification(1, "Mod1"), new OldSchoolModification(3, "Mod3"));

            Assert.IsFalse(a.Equals(b));
        }

        [Test]
        public void Test_modification_hash_set()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m1 = new Modification("23", null, "unknown", null, motif, null, null, 42.01, null, null, null, null, null, null);
            Modification m2 = new Modification("23", null, "unknown", null, motif, null, null, 42.01, null, null, null, null, null, null);
            HashSet<Modification> mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.AreEqual(1, mods.Count);
        }

        [Test]
        public void Test_modification2_hash_set()
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
        public void Test_modification3_hash_set() // numerical tolerance is 1e-9 so these two mods need to evaluate as identical
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
        public void TestInvalidModificationHash()
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
        public void TestFragmentationNoMod()
        {
            // First we're checking to see if the fragment masses of an unmodified peptide a calculated correctly
            var prot = new Protein("PEPTIDE", null);
            DigestionParams digestionParams = new DigestionParams(

                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications = new List<Modification>();
            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();

            // check unmodified
            var unmodPeptide = ye.Where(p => p.AllModsOneIsNterminus.Count == 0).First();
            var fragments = new List<Product>();
            unmodPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);
            var myUnmodFragmentMasses = fragments.Select(v => (int)Math.Round(v.NeutralMass.ToMz(1), 1)).ToList();
            HashSet<int> expectedMzs = new HashSet<int> { 98, 227, 324, 425, 538, 653, 703, 574, 477, 376, 263, 148 };

            Assert.That(expectedMzs.SetEquals(myUnmodFragmentMasses));
        }

        [Test]
        public void TestFragmentationModNoNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "oxidation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("O1"), _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            // check unmodified
            var unmodPeptide = ye.Where(p => p.AllModsOneIsNterminus.Count == 0).First();
            var myUnmodFragments = new List<Product>();
            unmodPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myUnmodFragments);
            var neutralMasses = new List<double>();
            neutralMasses.AddRange(myUnmodFragments.Select(m => m.NeutralMass).ToList());
            var expectedMasses = new List<double> { 97, 226, 323, 424, 537, 652, 147, 262, 375, 476, 573, 702 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            var firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            var secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set without oxidation
            Assert.AreEqual(12, myUnmodFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);

            // with oxidation, no neutral loss
            var modPeptide = ye.Where(p => p.AllModsOneIsNterminus.Count == 1).First();

            var myModFragments = new List<Product>();
            modPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragments);
            neutralMasses = new List<double>();
            neutralMasses.AddRange(myModFragments.Select(m => m.NeutralMass).ToList());
            expectedMasses = new List<double> { 97, 226, 323, 440, 553, 668, 147, 262, 375, 492, 589, 718 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set with oxidation
            Assert.AreEqual(12, myUnmodFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);
        }

        [Test]
        public void Test_FragmentationModNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count > 0).First();

            var myModFragments = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragments);
            HashSet<int> neutralMasses = new HashSet<int>(myModFragments.Select(m => (int)m.NeutralMass.ToMz(1)).ToList());
            HashSet<int> expectedMasses = new HashSet<int> { 98,227, 324, 407, 520, 635, 505, 618, 733, //b-ions with and without neutral loss
            148, 263, 376, 459, 556, 685, 557, 654, 783, //y-ions with and without neutral loss
            782}; //molecular ion with neutral loss

            CollectionAssert.AreEquivalent(neutralMasses, expectedMasses);
        }

        [Test]
        public void Test_FragmentationTwoModNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with 2 neutral loss mods
            ModificationMotif.TryGetMotif("Q", out ModificationMotif motifone);
            Modification modone = new Modification(_originalId: "ammonia", _modificationType: "testModType", _target: motifone, _monoisotopicMass: 0, _neutralLosses: new Dictionary<DissociationType,
                List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");

            ModificationMotif.TryGetMotif("T", out ModificationMotif motiftwo);
            Modification modtwo = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motiftwo, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: new Dictionary<DissociationType,
                List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");

            List<Modification> modlistone = new List<Modification> { modone };
            List<Modification> modlisttwo = new List<Modification> { modtwo };

            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEQTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 3, modlistone }, { 4, modlisttwo } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count == 2).First();

            var myModFragments = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragments);
            HashSet<int> neutralMasses = new HashSet<int>(myModFragments.Select(m => (int)m.NeutralMass.ToMz(1)).ToList());
            HashSet<int> expectedMasses = new HashSet<int> { 98, 227, 355, 536, 649, 764, 438, 551, 666, 338, 519, 632, 747, // b-ions with and without neutral losses
                                                             148, 263, 376, 557, 685, 814, 668, 797, 459, 587, 716, //y ions with and without neutral losses
                                                               813, 894, }; //molecular ion with neutral losses (phospho and ammonia respectively)

            CollectionAssert.AreEquivalent(neutralMasses, expectedMasses);
        }

        [Test]
        public void Test_FragmentationTwoModNeutralLossTwoFragTypes()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

            Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
            {
                { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
                { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
            };

            Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count == 1).First();

            var myModFragmentsHCD = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragmentsHCD);

            var neutralMassesHCD = myModFragmentsHCD.Select(m => (int)m.NeutralMass.ToMz(1));
            var expectedMassesHCD = new HashSet<int> { 98, 227, 324, 407, 520, 635, 505, 618, 733,// b-ions with and without neutral loss
                                                        148, 263, 376, 459, 556, 685, 557, 654, 783,//y-ions with and without neutral loss
                                                        782};// molecular ion with neutral loss

            CollectionAssert.AreEquivalent(expectedMassesHCD, neutralMassesHCD);

            //Now try the other half
            var myModFragmentsETD = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.ETD, FragmentationTerminus.Both, myModFragmentsETD);

            var neutralMassesETD = myModFragmentsETD.Select(m => (int)m.NeutralMass.ToMz(1));
            var expectedMassesETD = new HashSet<int> { 115, 244, 341, 505, 618, 733, 522, 635, 750,  // c-ions and c-17 ions
            148, 263, 376, 540, 637, 766, 557, 654, 783,          // y and y-17 ions
            133, 248, 361, 525, 622, 751, 542, 639, 768,         // z+1 and z+1-17 ions
            863 };//Molecular ions minus ammonia

            CollectionAssert.AreEquivalent(expectedMassesHCD, neutralMassesHCD);
        }
        
        [Test]
        public static void TestSerializationPeptideFromString()
        {
            // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
            // has the same properties as before it was serialized. This peptide is unmodified and generated from reading in a string
            string sequence = "PEPTIDE";
            PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>(), 0, null, null, 1, 7, 0);
            PeptideWithSetModifications deserializedPeptide = null;

            string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromString");
            System.IO.Directory.CreateDirectory(dir);
            string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

            var messageTypes = typeof(PeptideWithSetModifications);
            var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

            using (var file = System.IO.File.Create(path))
            {
                ser.Serialize(file, peptide);
            }

            using (var file = System.IO.File.OpenRead(path))
            {
                deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
            }

            deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, Protein>(), null);

            // not asserting any protein properties - since the peptide was created from a sequence string it didn't have a protein to begin with

            Assert.That(peptide.Equals(deserializedPeptide));
            Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
            Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

            var products = new List<Product>();

            deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

            Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        }

        [Test]
        public static void TestSerializationPeptideFromProtein()
        {
            // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
            // has the same properties as before it was serialized. This peptide is unmodified and generated from digesting a protein
            Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein");

            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            PeptideWithSetModifications deserializedPeptide = null;

            string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProtein");
            System.IO.Directory.CreateDirectory(dir);
            string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

            var messageTypes = typeof(PeptideWithSetModifications);
            var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

            using (var file = System.IO.File.Create(path))
            {
                ser.Serialize(file, peptide);
            }

            using (var file = System.IO.File.OpenRead(path))
            {
                deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
            }

            deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, Protein> { { protein.Accession, protein } }, peptide.DigestionParams);

            Assert.That(peptide.DigestionParams.Equals(deserializedPeptide.DigestionParams));
            Assert.That(peptide.Equals(deserializedPeptide));
            Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
            Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
            Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

            var products = new List<Product>();

            deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

            Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        }

        [Test]
        public static void TestSerializationPeptideFromProteinWithMod()
        {
            // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
            // has the same properties as before it was serialized. This peptide is modified with a phosphorylation

            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

            Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
            {
                { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
                { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
            };

            Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mod } } };

            Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein", oneBasedModifications: mods);

            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(v => v.AllModsOneIsNterminus.Count == 1).First();
            PeptideWithSetModifications deserializedPeptide = null;

            string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProteinWithMod");
            System.IO.Directory.CreateDirectory(dir);
            string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

            var messageTypes = typeof(PeptideWithSetModifications);
            var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

            using (var file = System.IO.File.Create(path))
            {
                ser.Serialize(file, peptide);
            }

            using (var file = System.IO.File.OpenRead(path))
            {
                deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
            }

            Dictionary<string, Modification> stringToMod = new Dictionary<string, Modification> { { mods.Values.First().First().IdWithMotif, mods.Values.First().First() } };

            deserializedPeptide.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, Protein> { { protein.Accession, protein } }, peptide.DigestionParams);

            Assert.That(peptide.Equals(deserializedPeptide));
            Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
            Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
            Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

            var products = new List<Product>();

            deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

            Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        }

        [Test]
        public static void TestFragmentNterminalModifiedPeptide()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification nTermMod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "N-terminal.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 1, new List<Modification> { nTermMod } } };

            Protein protein = new Protein("PEPTIDE", "", oneBasedModifications: mods);
            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(p => p.AllModsOneIsNterminus.Count == 1).First();
            Assert.That(peptide.FullSequence == "[testModType:acetylation on P]PEPTIDE");

            var fragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);

            var roundedFragments = fragments.Select(f => (int)f.NeutralMass).ToList();
            CollectionAssert.AreEquivalent(roundedFragments, new int[] { 139, 268, 365, 466, 579, 694, 147, 262, 375, 476, 573, 702 });
        }

        [Test]
        public static void TestFragmentCTerminalModifiedPeptide()
        {
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification cTermMod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "C-terminal.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 7, new List<Modification> { cTermMod } } };

            Protein protein = new Protein("PEPTIDE", "", oneBasedModifications: mods);
            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(p => p.AllModsOneIsNterminus.Count == 1).First();
            Assert.That(peptide.FullSequence == "PEPTIDE-[testModType:acetylation on E]");

            var fragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);

            var roundedFragments = fragments.Select(f => (int)f.NeutralMass).ToList();
            CollectionAssert.AreEquivalent(roundedFragments, new int[] { 97, 226, 323, 424, 537, 652, 189, 304, 417, 518, 615, 744 });
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
        public void CompareTo_SameModification_ReturnsZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_DifferentIdWithMotif_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod2", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_DifferentModificationType_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type2", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_DifferentTarget_ReturnsNonZero()
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
        public void CompareTo_DifferentLocationRestriction_ReturnsNonZero()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            var chemicalFormula = new ChemicalFormula();
            var mod1 = new Modification("mod1", "accession1", "type1", "feature1", motif, "C-terminal.", chemicalFormula, 100.0);
            var mod2 = new Modification("mod1", "accession1", "type1", "feature1", motif, "N-terminal.", chemicalFormula, 100.0);

            NUnit.Framework.Assert.That(mod1.CompareTo(mod2), Is.LessThan(0));
            NUnit.Framework.Assert.That(mod2.CompareTo(mod1), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_DifferentMonoisotopicMass_ReturnsNonZero()
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