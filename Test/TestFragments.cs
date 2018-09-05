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
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestFragments
    {
        private Peptide _mockPeptideEveryAminoAcid;

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

        [Test]
        public static void Test_GetTheoreticalFragments_UnmodifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));

            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();

            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;226.095356938-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.058243153-0", "y2;248.100836242-0" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_nTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.019094739-0", "b2;306.061687827-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.058243153-0", "y2;248.100836242-0" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_cTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;226.095356938-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 199, 328 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;199.024574042-0", "y2;328.067167131-0" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_internallyModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 306 };
            HashSet<int> foundNTerminalMasses = new HashSet<int>(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());

            Assert.That(expectedNTerminalMasses.SetEquals(foundNTerminalMasses));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;306.061687827-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 328 };
            HashSet<int> foundCTerminalMasses = new HashSet<int>(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());

            Assert.That(expectedCTerminalMasses.SetEquals(foundCTerminalMasses));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.058243153-0", "y2;328.067167131-0" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_nTerminalModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306, 79, 208 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.019094739-0", "b2;306.061687827-0", "b1;79.04219916561-97.97689557339", "b2;208.08479225361-97.97689557339" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.058243153-0", "y2;248.100836242-0" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_cTerminalModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            //var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> foundNTerminalMasses = new HashSet<int>(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };

            Assert.That(expectedNTerminalMasses.SetEquals(foundNTerminalMasses));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;226.095356938-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> foundCTerminalMasses = new HashSet<int>(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 199, 328, 101, 230 };

            Assert.That(expectedCTerminalMasses.SetEquals(foundCTerminalMasses));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;199.024574042-0", "y2;328.067167131-0", "y1;101.04767846861-97.97689557339", "y2;230.09027155761-97.97689557339" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_internallyModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var n = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 306, 208 };
            Assert.That(expectedNTerminalMasses.SetEquals(n.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;306.061687827-0", "b2;208.08479225361-97.97689557339" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var c = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 328, 230 };
            Assert.That(expectedCTerminalMasses.SetEquals(c.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.058243153-0", "y2;328.067167131-0", "y2;230.09027155761-97.97689557339" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_nTerminalModifiedPeptide_NeutralLoss_DissociationTypes_AnyActivationType_and_HCD()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.AnyActivationType, FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306, 79, 208 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.019094739-0", "b2;306.061687827-0", "b1;79.04219916561-97.97689557339", "b2;208.08479225361-97.97689557339" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.058243153-0", "y2;248.100836242-0" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_nTerminalModifiedPeptide_NeutralLoss_DissociationTypes_CID_and_HCD()//there should be no added neutral losses in this case becuase the allowed dissociation type doesn't match the dissociation type used in the experiment
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.019094739-0", "b2;306.061687827-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.058243153-0", "y2;248.100836242-0" };
            Assert.That(expectedCTerminalMassesLabels.SetEquals(cTerminalMassesLabels));
        }

        [Test]
        public static void Test_GetTheoreticalFragments_ProductTypeLabel()
        {
            Protein p = new Protein("PET", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.N);
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;226.095356938-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.AnyActivationType, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;226.095356938-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.CID, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b2;226.095356938-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.ECD, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "c1;114.079312951-0", "c2;243.121906039-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.ETD, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "c1;114.079312951-0", "c2;243.121906039-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.EThcD, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;226.095356938-0", "c1;114.079312951-0", "c2;243.121906039-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.ISCID, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = new List<ProductType> { };
            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.Custom, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.IRMPD, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276385-0", "b2;226.095356938-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.PQD, FragmentationTerminus.N);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));
        }

        [Test]
        public static void Test_Fragment_DiagnosticIons()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var diagnosticIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.D).ToList();
            Assert.AreEqual("D0;97.976895573-0", diagnosticIons.First().ToString());
        }

        [Test]
        public static void Test_Fragment_MolecularIon_NeutralLoss()
        {
            Protein p = new Protein("PTE", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var molecularIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.M).ToList();
            Assert.AreEqual("M0;327.14303540761-97.97689557339", molecularIons.First().ToString());
        }

        [Test]
        public static void Test_Fragment_DiagnosticIons_unmatchedDissociationType()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var diagnosticIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.D).ToList();
            Assert.AreEqual(0, diagnosticIons.Count());
        }

        [Test]
        public static void Test_Fragment_MolecularIon_NeutralLoss_unmatchedDissociationType()
        {
            Protein p = new Protein("PTE", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var molecularIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.M).ToList();
            Assert.AreEqual(0, molecularIons.Count());
        }

        [Test]
        public static void Test_NeutralMassShiftFromProductType()
        {
            foreach (ProductType p in Enum.GetValues(typeof(ProductType)))
            {
                double mass = Chemistry.ClassExtensions.RoundedDouble(DissociationTypeCollection.ProductTypeSpecificFragmentNeutralMass(0, p)).Value;
                switch (p)
                {
                    case ProductType.a:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("C-1O-1").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.aDegree:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("C-1O-2H-2").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.aStar:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("C-1O-1N-1H-3").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.b:
                        Assert.AreEqual(0, mass);
                        break;

                    case ProductType.bDegree:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-2O-1").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.bStar:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("N-1H-3").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.c:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("N1H3").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.D:
                        Assert.AreEqual(0, mass);
                        break;

                    case ProductType.M:
                        Assert.AreEqual(0, mass);
                        break;

                    case ProductType.x:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("C1O2").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.y:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H2O1").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.yDegree:
                        Assert.AreEqual(0, mass);
                        break;

                    case ProductType.yStar:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O1H-1N-1").MonoisotopicMass).Value, mass);
                        break;

                    case ProductType.zPlusOne:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O1H1N-1").MonoisotopicMass).Value, mass);
                        break;
                }
            }
        }

        [Test]
        public static void Test_NeutralMassShiftFromProductType_Exceptions()
        {
            ProductType undefinedProduct = (ProductType)99;

            Assert.Throws<MzLibException>(() => DissociationTypeCollection.ProductTypeSpecificFragmentNeutralMass(0, undefinedProduct), "Unknown product type!");
        }

        [Test]
        public static void Test_CustomDissociationType()
        {
            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom].Add(ProductType.b);
            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom].Add(ProductType.y);
            Assert.IsTrue(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom].Contains(ProductType.b));

            var productCollection = TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.N].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom]);
            Assert.IsTrue(productCollection.Contains(ProductType.b));

            productCollection = TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.C].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom]);
            Assert.IsTrue(productCollection.Contains(ProductType.y));
        }

        [Test]
        public static void Test_TerminusSpecificProductTypes()
        {
            Assert.AreEqual(new List<ProductType> { ProductType.b, ProductType.y }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.Both].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD]));
            Assert.AreEqual(new List<ProductType> { ProductType.b }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.N].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD]));
            Assert.AreEqual(new List<ProductType> { ProductType.y }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.C].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD]));
            Assert.AreEqual(new List<ProductType> { }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.None].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD]));

            Assert.AreEqual(new List<ProductType> { ProductType.c, ProductType.y, ProductType.zPlusOne }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.Both].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD]));
            Assert.AreEqual(new List<ProductType> { ProductType.c }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.N].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD]));
            Assert.AreEqual(new List<ProductType> { ProductType.y, ProductType.zPlusOne }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.C].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD]));
            Assert.AreEqual(new List<ProductType> { }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.None].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD]));

            Assert.AreEqual(new List<ProductType> { ProductType.b, ProductType.y }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.Both].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.CID]));
            Assert.AreEqual(new List<ProductType> { ProductType.b }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.N].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.CID]));
            Assert.AreEqual(new List<ProductType> { ProductType.y }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.C].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.CID]));
            Assert.AreEqual(new List<ProductType> { }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.None].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.CID]));
        }

        [Test]
        public static void Test_TerminusSpecificProductTypesFromPeptideWithSetMods()
        {
            Protein protein = new Protein("PEPTIDE", "accession");
            PeptideWithSetModifications p = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 7, "", 0, new Dictionary<int, Modification>(), 0);

            Assert.AreEqual(new List<ProductType> { ProductType.b, ProductType.y }, p.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { ProductType.b }, p.Fragment(DissociationType.HCD, FragmentationTerminus.N).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { ProductType.y }, p.Fragment(DissociationType.HCD, FragmentationTerminus.C).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { }, p.Fragment(DissociationType.HCD, FragmentationTerminus.None).Select(b => b.ProductType).Distinct().ToList());

            Assert.AreEqual(new List<ProductType> { ProductType.c, ProductType.y, ProductType.zPlusOne }, p.Fragment(DissociationType.ETD, FragmentationTerminus.Both).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { ProductType.c }, p.Fragment(DissociationType.ETD, FragmentationTerminus.N).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { ProductType.y, ProductType.zPlusOne }, p.Fragment(DissociationType.ETD, FragmentationTerminus.C).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { }, p.Fragment(DissociationType.ETD, FragmentationTerminus.None).Select(b => b.ProductType).Distinct().ToList());

            Assert.AreEqual(new List<ProductType> { ProductType.b, ProductType.y }, p.Fragment(DissociationType.CID, FragmentationTerminus.Both).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { ProductType.b }, p.Fragment(DissociationType.CID, FragmentationTerminus.N).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { ProductType.y }, p.Fragment(DissociationType.CID, FragmentationTerminus.C).Select(b => b.ProductType).Distinct().ToList());
            Assert.AreEqual(new List<ProductType> { }, p.Fragment(DissociationType.CID, FragmentationTerminus.None).Select(b => b.ProductType).Distinct().ToList());
        }

        [Test]
        public static void Test_MatchedFragmentIonToString()
        {
            Product P = new Product(ProductType.b, new NeutralTerminusFragment(FragmentationTerminus.N, 1, 1, 1), 0);
            MatchedFragmentIon m = new MatchedFragmentIon(P, 1, 1, 1);
            Assert.AreEqual("b1+1\t;1", m.ToString());
        }

        [Test]
        public static void Test_CID_Fragmentation_No_Unmodified_B1_ions()
        {
            //FOR CID B1 ions should always be missing whether or not there is a modification on first amino acid or not.

            Protein protein = new Protein("PEPTIDE", "accession");
            PeptideWithSetModifications p = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 7, "", 0, new Dictionary<int, Modification>(), 0);

            var f = p.Fragment(DissociationType.CID, FragmentationTerminus.Both);
            Assert.AreEqual(11, f.Count());

            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification m = new Modification(_originalId: "myId", _modificationType: "myModType", _target: motif, _monoisotopicMass: 10, _locationRestriction: "Anywhere.");
            List<Modification> modList = new List<Modification>() { m };
            Dictionary<int, List<Modification>> i = new Dictionary<int, List<Modification>>
            {
                { 1, modList }
            };

            protein = new Protein(sequence: "PEPTIDE", accession: "accession", oneBasedModifications: i);
            IEnumerable<PeptideWithSetModifications> pwsmList = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>());

            PeptideWithSetModifications modifiedPwsm = pwsmList.Where(z => z.AllModsOneIsNterminus.Count == 1).First();
            PeptideWithSetModifications unmodifiedPwsm = pwsmList.Where(z => z.AllModsOneIsNterminus.Count == 0).First();

            IEnumerable<Product> modifiedPwsmFragments = modifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both);
            IEnumerable<Product> unmodifiedPwsmFragments = unmodifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both);
            Assert.AreEqual(11, modifiedPwsmFragments.Count());
            Assert.AreEqual(11, unmodifiedPwsmFragments.Count());

            i = new Dictionary<int, List<Modification>>
            {
                { 2, modList }
            };

            protein = new Protein(sequence: "PPPTIDE", accession: "accession", oneBasedModifications: i);
            pwsmList = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>());

            modifiedPwsm = pwsmList.Where(z => z.AllModsOneIsNterminus.Count == 1).First();
            unmodifiedPwsm = pwsmList.Where(z => z.AllModsOneIsNterminus.Count == 0).First();

            modifiedPwsmFragments = modifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both);
            unmodifiedPwsmFragments = unmodifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both);
            Assert.AreEqual(11, modifiedPwsmFragments.Count());
            Assert.AreEqual(11, unmodifiedPwsmFragments.Count());
        }

        [Test]
        [TestCase(DissociationType.HCD, 12)]//the first part is the test case, the latter part is ther result of the assertion
        [TestCase(DissociationType.ETD, 17)]//the first part is the test case, the latter part is ther result of the assertion
        [TestCase(DissociationType.ECD, 17)]//the first part is the test case, the latter part is ther result of the assertion
        [TestCase(DissociationType.EThcD, 23)]//the first part is the test case, the latter part is ther result of the assertion
        public static void Test_ETD_ECD_EThcD_Fragmentation_No_FragmentsAtProline(DissociationType dissociationType, int fragmentCount)
        {
            Protein protein = new Protein(sequence: "PEPTIDE", accession: "accession");
            IEnumerable<PeptideWithSetModifications> pwsmList = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>());
            IEnumerable<PeptideWithSetModifications> digestionProducts = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>());
            PeptideWithSetModifications myPeptide = digestionProducts.First();
            IEnumerable<Product> myFragments = myPeptide.Fragment(dissociationType, FragmentationTerminus.Both);
            Assert.AreEqual(fragmentCount, myFragments.Count());
        }

        [Test]
        public static void CheckProlineFragments()
        {
            PeptideWithSetModifications p = new PeptideWithSetModifications("MPEPTIDE", new Dictionary<string, Modification>());
            var fragments = p.Fragment(DissociationType.ETD, FragmentationTerminus.Both);

            var z = fragments.Where(f => f.ProductType == ProductType.zPlusOne).ToList();

            var ionNums = z.Select(f => f.TerminusFragment.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 6, 8 };

            Assert.That(expected.SequenceEqual(ionNums));
        }

        [Test]
        public static void CheckProlineFragments2()
        {
            PeptideWithSetModifications p = new PeptideWithSetModifications("MTETTIDE", new Dictionary<string, Modification>());
            var fragments = p.Fragment(DissociationType.ETD, FragmentationTerminus.Both);

            var z = fragments.Where(f => f.ProductType == ProductType.zPlusOne).ToList();

            var ionNums = z.Select(f => f.TerminusFragment.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 5, 6, 7, 8 };

            Assert.That(expected.SequenceEqual(ionNums));
        }

        [Test]
        public static void CheckProlineFragments3()
        {
            PeptideWithSetModifications p = new PeptideWithSetModifications("METPIPEEEE", new Dictionary<string, Modification>());
            var fragments = p.Fragment(DissociationType.ETD, FragmentationTerminus.Both);

            var z = fragments.Where(f => f.ProductType == ProductType.zPlusOne).ToList();

            var ionNums = z.Select(f => f.TerminusFragment.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 6, 8, 9, 10 };

            Assert.That(expected.SequenceEqual(ionNums));
        }

        [Test]
        public static void CheckProlineFragments4()
        {
            ModificationMotif.TryGetMotif("P", out var motif);
            Modification m = new Modification("TEST", "", "OK", null, motif, "Anywhere.", null, 20);
            PeptideWithSetModifications p = new PeptideWithSetModifications("METP[OK:TEST]IPEEEE", new Dictionary<string, Modification> { { "TEST", m } });
            var fragments = p.Fragment(DissociationType.ETD, FragmentationTerminus.Both);

            var z = fragments.Where(f => f.ProductType == ProductType.zPlusOne).ToList();

            var ionNums = z.Select(f => f.TerminusFragment.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 6, 8, 9, 10 };

            Assert.That(expected.SequenceEqual(ionNums));
        }
    }
}