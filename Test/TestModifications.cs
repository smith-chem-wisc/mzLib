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
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;
using MassSpectrometry;

namespace Test
{
    [TestFixture]
    public sealed class TestModifications
    {
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
            var mod1 = new Modification(_id: "mod", _modificationType: "type", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: ChemicalFormula.ParseFormula("H").MonoisotopicMass);
            var mod1string = mod1.ToString();
            Assert.IsTrue(mod1string.Contains("MM"));
            var modAfterWriteRead = PtmListLoader.ReadModsFromString(mod1string + Environment.NewLine + "//").First() as Modification;

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
        public void ModificationSitesTest55()
        {
            Assert.IsTrue(ModificationSites.E.ContainsSites(ModificationSites.Any));
            Assert.IsFalse(ModificationSites.E.ContainsSites(ModificationSites.None));
            Assert.IsTrue(ModificationSites.None.ContainsSites(ModificationSites.None));
        }

        [Test]
        public void ChemicalFormulaModificaiton()
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
            Modification m1 = new Modification("23", null, "unknown", null, null, null, null, null, null, null, null, null, null, null);
            Modification m2 = new Modification("23", null, "unknown", null, null, null, null, null, null, null, null, null, null, null);
            HashSet<Modification> mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.AreEqual(1, mods.Count);
        }

        [Test]
        public void Test_modification2_hash_set()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m1 = new Modification("id1", null, "modificationType", null, motif, "Anywhere.", null, null, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            Modification m2 = new Modification("id1", null, "modificationType", null, motif, "Anywhere.", null, null, new Dictionary<string, IList<string>>(), null, null, null, null, null);
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
            Modification m1 = new Modification(_id: "id1", _modificationType: "modificationType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 1.11111d, _databaseReference: new Dictionary<string, IList<string>>(), _neutralLosses: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 2.222222 } } }, _diagnosticIons: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 1.2233 } } });
            Modification m2 = new Modification(_id: "id1", _modificationType: "modificationType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 1.11111d - 1e-10, _databaseReference: new Dictionary<string, IList<string>>(), _neutralLosses: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 2.222222 + 1e-10 } } }, _diagnosticIons: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 1.2233 } } });
            m1.DatabaseReference.Add("key", new List<string> { "value" });
            m2.DatabaseReference.Add("key", new List<string> { "value" });

            HashSet<Modification> mods = new HashSet<Modification>(new Modification[] { m1, m2 });
            Assert.AreEqual(1, mods.Count);
            Assert.IsTrue(m1.Equals(m2));
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
            var myUnmodFragments = unmodPeptide.GetTheoreticalFragments(DissociationType.HCD, FragmentationTerminus.Both);
            List<double> neutralMasses = new List<double>();
            neutralMasses.AddRange(myUnmodFragments.Select(m => m.Mass).ToList());
            List<double> expectedMasses = new List<double> { 226, 323, 424, 537, 652, 147, 262, 375, 476, 573, 702 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            var firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            var secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            Assert.AreEqual(11, myUnmodFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);
        }

        [Test]
        public void TestFragmentationModNoNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification mod = new Modification(_id: "oxidation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("O1"), _locationRestriction: "Anywhere.");
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
            var myUnmodFragments = unmodPeptide.GetTheoreticalFragments(DissociationType.HCD, FragmentationTerminus.Both);
            var neutralMasses = new List<double>();
            neutralMasses.AddRange(myUnmodFragments.Select(m => m.Mass).ToList());
            var expectedMasses = new List<double> { 226, 323, 424, 537, 652, 147, 262, 375, 476, 573, 702 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            var firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            var secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set without oxidation
            Assert.AreEqual(11, myUnmodFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);

            // with oxidation, no neutral loss
            var modPeptide = ye.Where(p => p.AllModsOneIsNterminus.Count == 1).First();

            var compactPeptide = modPeptide.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);



            var myModFragments = modPeptide.GetTheoreticalFragments(DissociationType.HCD, FragmentationTerminus.Both);
            neutralMasses = new List<double>();
            neutralMasses.AddRange(myModFragments.Select(m => m.Mass).ToList());
            expectedMasses = new List<double> { 226, 323, 440, 553, 668, 147, 262, 375, 492, 589, 718 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set with oxidation
            Assert.AreEqual(11, myUnmodFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);

        }

        [Test]
        public void Test_FragmentationModNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification mod = new Modification(_id: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count > 0).First();

            var myModFragments = peptideWithNeutralMassMod.GetTheoreticalFragments(DissociationType.HCD, FragmentationTerminus.Both);
            var neutralMasses = new List<double>();
            neutralMasses.AddRange(myModFragments.Select(m => m.Mass).ToList());
            var expectedMasses = new List<double> { 226, 323, 504, 617, 732, 406, 519, 634, 147, 262, 375, 556, 653, 782, 458, 555, 684 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            var firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            var secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set with oxidation
            Assert.AreEqual(17, myModFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);

        }

        [Test]
        public void Test_FragmentationTwoModNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification mod = new Modification(_id: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PETTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 3, modlist }, { 4, modlist} });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count == 2).First();

            var myModFragments = peptideWithNeutralMassMod.GetTheoreticalFragments(DissociationType.HCD, FragmentationTerminus.Both);
            var neutralMasses = new List<double>();
            neutralMasses.AddRange(myModFragments.Select(m => m.Mass).ToList());
            var expectedMasses = new List<double> { 226, 407, 588, 701, 816, 309, 490, 603, 718, 392, 505, 620, 147, 262, 375, 556, 737, 866, 458, 639, 768, 541, 670 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }
            neutralMasses = neutralMasses.Distinct().ToList();
            var firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            var secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set with oxidation
            Assert.AreEqual(23, neutralMasses.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);

        }


        [Test]
        public void Test_FragmentationTwoModNeutralLossTwoFragTypes()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

            Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
            {
                { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
                { DissociationType.ETD, new List<double>() { 0 } }
            };

            Modification mod = new Modification(_id: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count == 1).First();

            var myModFragmentsHCD = peptideWithNeutralMassMod.GetTheoreticalFragments(DissociationType.HCD, FragmentationTerminus.Both);



            var neutralMassesHCD = new List<double>();
            neutralMassesHCD.AddRange(myModFragmentsHCD.Select(m => m.Mass).ToList());
            var expectedMassesHCD = new List<double> {226, 323, 504, 617, 732, 406, 519, 634, 147, 262, 375, 556, 653, 782, 458, 555, 684 };
            for (int i = 0; i < neutralMassesHCD.Count; i++)
            {
                neutralMassesHCD[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMassesHCD[i], 0).Value;
            }
            neutralMassesHCD = neutralMassesHCD.Distinct().ToList();
            var firstNotSecond = neutralMassesHCD.Except(expectedMassesHCD).ToList();
            var secondNotFirst = expectedMassesHCD.Except(neutralMassesHCD).ToList();

            //this is the set with oxidation
            Assert.AreEqual(17, neutralMassesHCD.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);



            var myModFragmentsEtd = peptideWithNeutralMassMod.GetTheoreticalFragments(DissociationType.ETD, FragmentationTerminus.Both);


            var neutralMassesEtd = new List<double>();
            neutralMassesEtd.AddRange(myModFragmentsEtd.Select(m => m.Mass).ToList());
            var expectedMassesEtd = new List<double> { 114, 243, 340, 521, 634, 749, 147, 262, 375, 556, 653, 782, 131, 246, 359, 540, 637, 766 };
            for (int i = 0; i < neutralMassesEtd.Count; i++)
            {
                neutralMassesEtd[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMassesEtd[i], 0).Value;
            }
            neutralMassesEtd = neutralMassesEtd.Distinct().ToList();
            firstNotSecond = neutralMassesEtd.Except(expectedMassesEtd).ToList();
            secondNotFirst = expectedMassesEtd.Except(neutralMassesEtd).ToList();

            //this is the set with oxidation
            Assert.AreEqual(18, neutralMassesEtd.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);


        }

    }
}