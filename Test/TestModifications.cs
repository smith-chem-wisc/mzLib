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
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public sealed class TestModifications
    {
        [Test]
        public static void Test_modificationsHashCode()
        {
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var mod1 = new ModificationGeneral("mod", null, "type", null, motif, "Anywhere.", null, 1, null, null, null, null, null, null);
            var mod2 = new ModificationGeneral("mod2", null, "type", null, motif, "Anywhere.", null, 10, null, null, null, null, null, null);

            Assert.AreNotEqual(mod1.GetHashCode(), mod2.GetHashCode());
            Assert.AreNotEqual(mod1, mod2);
            HashSet<ModificationGeneral> myHashSet = new HashSet<ModificationGeneral>
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
            var mod1 = new ModificationGeneral(_Id: "mod", _ModificationType: "type", _Target: motif, _Position: "Anywhere.", _ChemicalFormula: ChemicalFormula.ParseFormula("H"), _MonoisotopicMass: ChemicalFormula.ParseFormula("H").MonoisotopicMass);
            var mod1string = mod1.ToString();
            Assert.IsTrue(mod1string.Contains("MM"));
            var modAfterWriteRead = PtmListLoaderGeneral.ReadModsFromString(mod1string + Environment.NewLine + "//").First() as ModificationGeneral;

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
            ModificationGeneral m1 = new ModificationGeneral("23", null, "unknown", null, null, null, null, null, null, null, null, null, null, null);
            ModificationGeneral m2 = new ModificationGeneral("23", null, "unknown", null, null, null, null, null, null, null, null, null, null, null);
            HashSet<ModificationGeneral> mods = new HashSet<ModificationGeneral>(new ModificationGeneral[] { m1, m2 });
            Assert.AreEqual(1, mods.Count);
        }

        [Test]
        public void Test_modification2_hash_set()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            ModificationGeneral m1 = new ModificationGeneral("id1", null, "modificationType", null, motif, "Anywhere.", null, null, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            ModificationGeneral m2 = new ModificationGeneral("id1", null, "modificationType", null, motif, "Anywhere.", null, null, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            m1.DatabaseReference.Add("key", new List<string> { "value" });
            m2.DatabaseReference.Add("key", new List<string> { "value" });
            HashSet<ModificationGeneral> mods = new HashSet<ModificationGeneral>(new ModificationGeneral[] { m1, m2 });
            Assert.True(m1.Equals(m2));
            Assert.AreEqual(1, mods.Count);
        }

        [Test]
        public void Test_modification3_hash_set() // numerical tolerance is 1e-9 so these two mods need to evaluate as identical
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            ModificationGeneral m1 = new ModificationGeneral(_Id: "id1", _ModificationType: "modificationType", _Target: motif, _Position: "Anywhere.", _MonoisotopicMass: 1.11111d, _DatabaseReference: new Dictionary<string, IList<string>>(), _NeutralLosses: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 2.222222 } } }, _DiagnosticIons: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 1.2233 } } });
            ModificationGeneral m2 = new ModificationGeneral(_Id: "id1", _ModificationType: "modificationType", _Target: motif, _Position: "Anywhere.", _MonoisotopicMass: 1.11111d - 1e-10, _DatabaseReference: new Dictionary<string, IList<string>>(), _NeutralLosses: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 2.222222 + 1e-10 } } }, _DiagnosticIons: new Dictionary<MassSpectrometry.DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 1.2233 } } });
            m1.DatabaseReference.Add("key", new List<string> { "value" });
            m2.DatabaseReference.Add("key", new List<string> { "value" });

            HashSet<ModificationGeneral> mods = new HashSet<ModificationGeneral>(new ModificationGeneral[] { m1, m2 });
            Assert.AreEqual(1, mods.Count);
            Assert.IsTrue(m1.Equals(m2));
        }
    }
}