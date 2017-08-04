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
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestFragments
    {

        #region Private Fields

        private Peptide _mockPeptideEveryAminoAcid;

        #endregion Private Fields

        #region Public Methods

        [SetUp]
        public void SetUp()
        {
            _mockPeptideEveryAminoAcid = new Peptide("ACDEFGHIKLMNPQRSTVWY");
        }

        [Test]
        public void FragmentNumberToHigh()
        {
            Assert.Throws<IndexOutOfRangeException>(() => _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, 25).ToList());
        }

        [Test]
        public void FragmentName()
        {
            Fragment fragment = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.a, 1).ToArray()[0];

            Assert.AreEqual("a1", fragment.ToString());
        }

        [Test]
        public void FragmentAllBIons()
        {
            List<Fragment> fragments = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b).ToList();

            Assert.AreEqual(19, fragments.Count);
        }

        [Test]
        public void FragmentAnotherTest()
        {
            List<Fragment> fragments = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, 1, 2).ToList();

            Assert.AreEqual(2, fragments.Count);
        }

        [Test]
        public void FragmentNTermModTest()
        {
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("O"), "lala", ModificationSites.NTerminus));
            Fragment fragment = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, 1).First();
            Assert.IsTrue(fragment.Modifications.SequenceEqual(new List<OldSchoolModification> { new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("O"), "lala", ModificationSites.NTerminus) }));
        }

        [Test]
        public void FragmentModifications()
        {
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(1, "mod1", ModificationSites.C));
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(2, "mod2", ModificationSites.D));
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(3, "mod3", ModificationSites.A));
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(4, "mod4", ModificationSites.Y));
            Fragment fragment = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, 1).First();
            Fragment fragmentEnd = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.y, 1).Last();

            Assert.IsTrue(fragment.Modifications.SequenceEqual(new List<OldSchoolModification> { new OldSchoolModification(3, "mod3", ModificationSites.A) }));

            Assert.IsTrue(fragmentEnd.Modifications.SequenceEqual(new List<OldSchoolModification> { new OldSchoolModification(4, "mod4", ModificationSites.Y) }));
        }

        [Test]
        public void ChemicalFormulaFragment()
        {
            var a = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, true);
            // Can break in 19 places
            Assert.AreEqual(19, a.Count());
            Assert.IsTrue(a.Select(b => b.Sequence).Contains("ACDEFG"));

            var y = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.y, true);
            // Can break in 19 places
            Assert.IsTrue(y.Select(b => b.Sequence).Contains("TVWY"));

            var c = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, true);

            Assert.AreEqual(a.First(), c.First());
        }

        [Test]
        public void TestGetSiteDeterminingFragments()
        {
            var pep1 = new Peptide("ACDEFG");
            var pep2 = new Peptide("ACTVWY");
            var ok = pep1.GetSiteDeterminingFragments(pep2, FragmentTypes.b);
            Assert.AreEqual(6, ok.Count());
            Assert.Contains("ACT", ok.Select(b => b.Sequence).ToArray());
        }

        [Test]
        public void TestFormulaTerminusMods()
        {
            var pep1 = new Peptide("ACDEFG");
            pep1.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H"), ModificationSites.NTerminus));

            Assert.IsTrue(pep1.Fragment(FragmentTypes.b, true).First() is IHasChemicalFormula);

            var pep2 = new Peptide("ACDEFG");
            pep2.AddModification(new OldSchoolModification(2, "haha", ModificationSites.NTerminus));
            Assert.IsFalse(pep2.Fragment(FragmentTypes.b, true).First() is IHasChemicalFormula);

            var pep3 = new Peptide("ACDEFG");
            pep3.AddModification(new OldSchoolModification(3, "haha", ModificationSites.D));

            var list = pep3.Fragment(FragmentTypes.b, true).ToList();

            Assert.IsTrue(list[0] is IHasChemicalFormula);
            Assert.IsFalse(list[2] is IHasChemicalFormula);
        }

        [Test]
        public void CleavageIndicesTest()
        {
            IEnumerable<IProtease> proteases = new List<TestProtease> { new TestProtease() };
            var ok1 = AminoAcidPolymer.GetCleavageIndexes("ACDEFG", proteases, true).ToList();
            var ok2 = AminoAcidPolymer.GetCleavageIndexes("ACDEFG", proteases, false).ToList();
            var ok3 = AminoAcidPolymer.GetCleavageIndexes("ACDE", proteases, true).ToList();
            var ok4 = AminoAcidPolymer.GetCleavageIndexes("ACDE", proteases, false).ToList();
            Assert.AreEqual(3, ok1.Count());
            Assert.AreEqual(2, ok2.Count());
            Assert.AreEqual(4, ok3.Count());
            Assert.AreEqual(2, ok4.Count());
        }

        [Test]
        public void TestGetIonCapFailFail()
        {
            FragmentTypes f = FragmentTypes.All;
            Assert.That(() => f.GetIonCap(), Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("Fragment Type must be a single value to determine the ion cap"));
        }

        [Test]
        public void TestGetTerminusFail()
        {
            FragmentTypes f = FragmentTypes.a | FragmentTypes.adot;
            Assert.That(() => f.GetTerminus(), Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("Fragment Type must be a single value to determine the terminus"));
        }

        #endregion Public Methods

    }
}