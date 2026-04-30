// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (FragmentModificationTests.cs) is part of Proteomics.
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

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.Omics.FragmentationTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class FragmentModificationTests
    {
        private static Peptide _mockPeptideEveryAminoAcid;
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
            _mockPeptideEveryAminoAcid = new Peptide("ACDEFGHIKLMNPQRSTVWY");
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void FragmentNTermModTest()
        {
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("O"), "lala", ModificationSites.NTerminus));
            global::Proteomics.AminoAcidPolymer.Fragment fragment = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, 1).First();
            Assert.IsTrue(fragment.Modifications.SequenceEqual(new List<OldSchoolModification> { new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("O"), "lala", ModificationSites.NTerminus) }));
        }

        [Test]
        public static void FragmentModifications()
        {
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(1, "mod1", ModificationSites.C));
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(2, "mod2", ModificationSites.D));
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(3, "mod3", ModificationSites.A));
            _mockPeptideEveryAminoAcid.AddModification(new OldSchoolModification(4, "mod4", ModificationSites.Y));
            global::Proteomics.AminoAcidPolymer.Fragment fragment = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.b, 1).First();
            global::Proteomics.AminoAcidPolymer.Fragment fragmentEnd = _mockPeptideEveryAminoAcid.Fragment(FragmentTypes.y, 1).Last();

            Assert.IsTrue(fragment.Modifications.SequenceEqual(new List<OldSchoolModification> { new OldSchoolModification(3, "mod3", ModificationSites.A) }));

            Assert.IsTrue(fragmentEnd.Modifications.SequenceEqual(new List<OldSchoolModification> { new OldSchoolModification(4, "mod4", ModificationSites.Y) }));
        }

        [Test]
        public static void TestFormulaTerminusMods()
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
        public static void Test_GetTheoreticalFragments_nTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.01909-0", "b2;306.06169-0" };
            CollectionAssert.AreEquivalent(expectedNTerminalMassesLabels, nTerminalMassesLabels);

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.05824-0", "y2;248.10084-0" };
            CollectionAssert.AreEquivalent(expectedCTerminalMassesLabels, cTerminalMassesLabels);
        }

        [Test]
        public static void Test_GetTheoreticalFragments_cTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;226.09536-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 199, 328 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;199.02457-0", "y2;328.06717-0" };
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

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 306 };
            HashSet<int> foundNTerminalMasses = new HashSet<int>(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());

            Assert.That(expectedNTerminalMasses.SetEquals(foundNTerminalMasses));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;306.06169-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 328 };
            HashSet<int> foundCTerminalMasses = new HashSet<int>(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());

            Assert.That(expectedCTerminalMasses.SetEquals(foundCTerminalMasses));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.05824-0", "y2;328.06717-0" };
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

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306, 79, 208 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.01909-0", "b2;306.06169-0", "b1;79.04220-97.98", "b2;208.08479-97.98" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.05824-0", "y2;248.10084-0" };
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

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);

            //var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> foundNTerminalMasses = new HashSet<int>(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };

            Assert.That(expectedNTerminalMasses.SetEquals(foundNTerminalMasses));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;226.09536-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> foundCTerminalMasses = new HashSet<int>(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 199, 328, 101, 230 };

            CollectionAssert.AreEquivalent(expectedCTerminalMasses, foundCTerminalMasses);
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;199.02457-0", "y2;328.06717-0", "y1;101.04768-97.98", "y2;230.09027-97.98" };
            CollectionAssert.AreEquivalent(expectedCTerminalMassesLabels, cTerminalMassesLabels);
        }

        [Test]
        public static void Test_GetTheoreticalFragments_internallyModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);

            //evaluate N-terminal masses
            var n = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 306, 208 };
            Assert.That(expectedNTerminalMasses.SetEquals(n.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;306.06169-0", "b2;208.08479-97.98" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var c = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 328, 230 };
            Assert.That(expectedCTerminalMasses.SetEquals(c.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.05824-0", "y2;328.06717-0", "y2;230.09027-97.98" };
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

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.AnyActivationType, FragmentationTerminus.Both, theseTheoreticalFragments);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306, 79, 208 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.01909-0", "b2;306.06169-0", "b1;79.04220-97.98", "b2;208.08479-97.98" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.05824-0", "y2;248.10084-0" };
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

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;177.01909-0", "b2;306.06169-0" };
            CollectionAssert.AreEquivalent(expectedNTerminalMassesLabels, nTerminalMassesLabels);

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.05824-0", "y2;248.10084-0" };
            CollectionAssert.AreEquivalent(expectedCTerminalMassesLabels, cTerminalMassesLabels);
        }

        [Test]
        public static void Test_Fragment_DiagnosticIons()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var diagnosticIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.D).ToList();
            Assert.AreEqual("D99;97.97690-0", diagnosticIons.First().ToString());
        }

        [Test]
        public static void Test_Fragment_DiagnosticIons_unmatchedDissociationType()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var diagnosticIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.D).ToList();
            Assert.AreEqual(0, diagnosticIons.Count());
        }

        [Test]
        public static void Test_Fragment_MolecularIon_NeutralLoss()
        {
            Protein p = new Protein("PTE", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(
                _originalId: "phospho",
                _modificationType: "CommonBiological",
                _target: motif,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"),
                _neutralLosses: new Dictionary<DissociationType, List<double>>
                {
                    { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } }
                },
                _diagnosticIons: new Dictionary<DissociationType, List<double>>
                {
                    { MassSpectrometry.DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } }
                });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var molecularIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.M).ToList();
            Assert.AreEqual("M0;327.14304-97.98", molecularIons.First().ToString());
        }

        [Test]
        public static void Test_Fragment_MolecularIon_NeutralLoss_unmatchedDissociationType()
        {
            Protein p = new Protein("PTE", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var molecularIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.M).ToList();
            Assert.AreEqual(0, molecularIons.Count());
        }

        [Test]
        public static void TestThatDiagnosticIonsDontDuplicate()
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            Modification modWithDiagnosticIons = new Modification(
                _originalId: "Test",
                _modificationType: "TestType",
                _target: motif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 1,
                _diagnosticIons: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { 4.0 } } });

            PeptideWithSetModifications p = new PeptideWithSetModifications("P[TestType:Test]E[TestType:Test]P[TestType:Test]TIDE",
                new Dictionary<string, Modification> { { "Test", modWithDiagnosticIons } });

            var fragments = new List<Product>();
            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);

            var diagnosticIons = fragments.Where(f => f.ProductType == ProductType.D).ToList();

            // if there are 3 diagnostic ions (one for each mod on the peptide),
            // diagnostic ions are being incorrectly duplicated!
            Assert.That(diagnosticIons.Count == 1);
        }

        [Test]
        // This unit test tests two things:
        // 1. that neutral loss fragments are correctly generated when a mod w/ annotated neutral losses is on the C or N terminus of a peptide
        // 2. that neutral loss fragments are correctly generated when a mod w/ annotated neutral losses with a "Any activation type" dissociation type is present
        public static void TestNeutralLossModOnTerminus()
        {
            ModificationMotif.TryGetMotif("X", out var motif);

            var modsDictionary = new Dictionary<string, Modification>
            {
                { @"NeutralLoss-N on X", new Modification(
                    _originalId: @"NeutralLoss-N",
                    _modificationType: "Test Category",
                    _target: motif,
                    _locationRestriction: "Peptide N-terminal.",
                    _monoisotopicMass: 225,
                    _neutralLosses: new Dictionary<DissociationType, List<double>>{ { DissociationType.HCD, new List<double> { 126 } } })
                },

                { @"NeutralLoss-C on X", new Modification(
                    _originalId: @"NeutralLoss-C",
                    _modificationType: "Test Category",
                    _target: motif,
                    _locationRestriction: "Peptide C-terminal.",
                    _monoisotopicMass: 225,
                    _neutralLosses: new Dictionary<DissociationType, List<double>>{ { DissociationType.AnyActivationType, new List<double> { 126 } } })
                },
            };

            List<Product> products = new List<Product>();

            // N-term mod
            PeptideWithSetModifications pep = new PeptideWithSetModifications(@"[Test Category:NeutralLoss-N on X]AAAAAAAAAK", modsDictionary);
            pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            Assert.That(products.Count(p => p.NeutralLoss == 126) == 10);

            // C-term mod mod
            pep = new PeptideWithSetModifications(@"AAAAAAAAAK[Test Category:NeutralLoss-C on X]", modsDictionary);
            pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            Assert.That(products.Count(p => p.NeutralLoss == 126) == 10);
        }

        [Test]
        // This unit test tests:
        // that diagnostic ionsare correctly generated when "Any activation type" is used
        public static void TestDiagnosticIonsAnyActivationType()
        {
            ModificationMotif.TryGetMotif("X", out var motif);

            var modsDictionary = new Dictionary<string, Modification>
            {
                { @"DiagnosticIon-N on X", new Modification(
                    _originalId: @"DiagnosticIon-N",
                    _modificationType: "Test Category",
                    _target: motif,
                    _locationRestriction: "Peptide N-terminal.",
                    _monoisotopicMass: 225,
                    _diagnosticIons: new Dictionary<DissociationType, List<double>>{ { DissociationType.AnyActivationType, new List<double> { 126 } } })
                },

                { @"DiagnosticIon-Anywhere on X", new Modification(
                    _originalId: @"DiagnosticIon-Anywhere",
                    _modificationType: "Test Category",
                    _target: motif,
                    _locationRestriction: "Anywhere.",
                    _monoisotopicMass: 225,
                    _diagnosticIons: new Dictionary<DissociationType, List<double>>{ { DissociationType.HCD, new List<double> { 126 } } })
                },
            };

            List<Product> products = new List<Product>();

            // 2 mods w/ the same diagnostic ion, should end up w/ only 1 diagnostic ion (D127)
            PeptideWithSetModifications pep = new PeptideWithSetModifications(
                @"[Test Category:DiagnosticIon-N on X]AAAA[Test Category:DiagnosticIon-Anywhere on X]AAAAAK", modsDictionary);

            pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            Assert.That(products.Count(p => p.ProductType == ProductType.D) == 1);
            Assert.That(products.Count(p => p.Annotation == "D127") == 1);
        }

        [Test]
        // This unit test tests:
        // ETD generates an extra zDot ion. this ion should have an annotated neutral loss if there is a mod w/ neutral loss on the C-term
        public static void TestCTermEtd()
        {
            ModificationMotif.TryGetMotif("X", out var motif);

            var modsDictionary = new Dictionary<string, Modification>
            {
                { @"NeutralLoss-C on X", new Modification(
                    _originalId: @"NeutralLoss-C",
                    _modificationType: "Test Category",
                    _target: motif,
                    _locationRestriction: "Peptide C-terminal.",
                    _monoisotopicMass: 225,
                    _neutralLosses: new Dictionary<DissociationType, List<double>>{ { DissociationType.AnyActivationType, new List<double> { 126 } } })
                },
            };

            List<Product> products = new List<Product>();

            PeptideWithSetModifications pep = new PeptideWithSetModifications(@"AAAAAAAAAK[Test Category:NeutralLoss-C on X]", modsDictionary);

            pep.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            Assert.That(products.Count(p => p.NeutralLoss == 126 && p.ProductType == ProductType.zDot) == 10);
            Assert.That(products.Count(p => p.NeutralLoss == 126 && p.ProductType == ProductType.c) == 0);
            Assert.That(products.Count(p => p.NeutralLoss == 126 && p.ProductType == ProductType.y) == 9);
            Assert.That(products.Count(p => p.NeutralLoss == 126 && p.ProductType == ProductType.M) == 1);
            Assert.That(products.Count(p => p.NeutralLoss == 126) == 20);
        }
    }
}
