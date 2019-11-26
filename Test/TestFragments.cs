﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
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
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    public sealed class TestFragments
    {
        private Peptide _mockPeptideEveryAminoAcid;
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
        [TestCase(DissociationType.HCD, new[] { ProductType.b, ProductType.y })]
        [TestCase(DissociationType.ECD, new[] { ProductType.c, ProductType.y, ProductType.zDot })]
        [TestCase(DissociationType.ETD, new[] { ProductType.c, ProductType.y, ProductType.zDot })]
        [TestCase(DissociationType.EThcD, new[] { ProductType.b, ProductType.y })]
        [TestCase(DissociationType.AnyActivationType, new[] { ProductType.b, ProductType.y })]
        public void TestDissociationProductTypes(DissociationType dissociationType, IEnumerable<ProductType> expectedProductTypes)
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[dissociationType];
            Assert.IsTrue(expectedProductTypes.All(productType => d.Contains(productType)));
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

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);

            //evaluate N-terminal masses
            var nTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));

            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();

            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;226.09536-0" };
            CollectionAssert.AreEquivalent(expectedNTerminalMassesLabels, nTerminalMassesLabels);

            //evaluate C-terminal masses
            var cTerminalMasses = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 248 };
            CollectionAssert.AreEquivalent(expectedCTerminalMasses, cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)));
            var cTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.C).Select(f => f.ToString()).ToList();
            HashSet<string> expectedCTerminalMassesLabels = new HashSet<string> { "y1;119.05824-0", "y2;248.10084-0" };
            CollectionAssert.AreEquivalent(expectedCTerminalMassesLabels, cTerminalMassesLabels);
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
        public static void Test_GetTheoreticalFragments_ProductTypeLabel()
        {
            Protein p = new Protein("PET", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.N, theseTheoreticalFragments);
            var nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            HashSet<string> expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;226.09536-0" };
            CollectionAssert.AreEquivalent(expectedNTerminalMassesLabels, nTerminalMassesLabels);

            aPeptideWithSetModifications.Fragment(DissociationType.AnyActivationType, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;226.09536-0" };
            CollectionAssert.AreEquivalent(expectedNTerminalMassesLabels, nTerminalMassesLabels);

            aPeptideWithSetModifications.Fragment(DissociationType.CID, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b2;226.09536-0" };
            CollectionAssert.AreEquivalent(expectedNTerminalMassesLabels, nTerminalMassesLabels);

            aPeptideWithSetModifications.Fragment(DissociationType.ECD, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "c1;114.07931-0", "c2;243.12191-0" };
            CollectionAssert.AreEquivalent(expectedNTerminalMassesLabels, nTerminalMassesLabels);

            aPeptideWithSetModifications.Fragment(DissociationType.ETD, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "c1;114.07931-0", "c2;243.12191-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            aPeptideWithSetModifications.Fragment(DissociationType.EThcD, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;226.09536-0", "c1;114.07931-0", "c2;243.12191-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            aPeptideWithSetModifications.Fragment(DissociationType.ISCID, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = new List<ProductType> { };
            aPeptideWithSetModifications.Fragment(DissociationType.Custom, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            aPeptideWithSetModifications.Fragment(DissociationType.IRMPD, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
            expectedNTerminalMassesLabels = new HashSet<string> { "b1;97.05276-0", "b2;226.09536-0" };
            Assert.That(expectedNTerminalMassesLabels.SetEquals(nTerminalMassesLabels));

            aPeptideWithSetModifications.Fragment(DissociationType.PQD, FragmentationTerminus.N, theseTheoreticalFragments);
            nTerminalMassesLabels = theseTheoreticalFragments.Where(f => f.Terminus == FragmentationTerminus.N).Select(f => f.ToString()).ToList();
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

            var theseTheoreticalFragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theseTheoreticalFragments);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            //evaluate N-terminal masses
            var diagnosticIons = theseTheoreticalFragments.Where(f => f.ProductType == ProductType.D).ToList();
            Assert.AreEqual("D99;97.97690-0", diagnosticIons.First().ToString());
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
        [TestCase("AAAAAAAAAA", DissociationType.CID, 0, 0, 0, 0, 0, 0, 17)]
        [TestCase("AAAAAAAAAA", DissociationType.LowCID, 0, 0, 0, 0, 0, 0, 17)]
        [TestCase("RAAAAAAAAA", DissociationType.LowCID, 8, 8, 0, 0, 0, 0, 33)]
        [TestCase("KAAAAAAAAA", DissociationType.LowCID, 8, 8, 0, 0, 0, 0, 33)]
        [TestCase("NAAAAAAAAA", DissociationType.LowCID, 8, 8, 0, 0, 0, 0, 33)]
        [TestCase("QAAAAAAAAA", DissociationType.LowCID, 8, 8, 0, 0, 0, 0, 33)]
        [TestCase("AAAAAAAAAR", DissociationType.LowCID, 0, 0, 0, 0, 9, 0, 26)]
        [TestCase("AAAAAAAAAK", DissociationType.LowCID, 0, 0, 0, 0, 9, 0, 26)]
        [TestCase("AAAAAAAAAN", DissociationType.LowCID, 0, 0, 0, 0, 9, 0, 26)]
        [TestCase("AAAAAAAAAQ", DissociationType.LowCID, 0, 0, 0, 0, 9, 0, 26)]
        [TestCase("AARAAAAAAA", DissociationType.LowCID, 7, 7, 0, 0, 2, 0, 33)]
        [TestCase("AAKAAAAAAA", DissociationType.LowCID, 7, 7, 0, 0, 2, 0, 33)]
        [TestCase("AANAAAAAAA", DissociationType.LowCID, 7, 7, 0, 0, 2, 0, 33)]
        [TestCase("AAQAAAAAAA", DissociationType.LowCID, 7, 7, 0, 0, 2, 0, 33)]
        [TestCase("SAAAAAAAAA", DissociationType.LowCID, 0, 0, 8, 8, 0, 0, 33)]
        [TestCase("TAAAAAAAAA", DissociationType.LowCID, 0, 0, 8, 8, 0, 0, 33)]
        [TestCase("EAAAAAAAAA", DissociationType.LowCID, 0, 0, 8, 8, 0, 0, 33)]
        [TestCase("DAAAAAAAAA", DissociationType.LowCID, 0, 0, 8, 8, 0, 0, 33)]
        [TestCase("AAAAAAAAAS", DissociationType.LowCID, 0, 0, 0, 0, 0, 9, 26)]
        [TestCase("AAAAAAAAAT", DissociationType.LowCID, 0, 0, 0, 0, 0, 9, 26)]
        [TestCase("AAAAAAAAAE", DissociationType.LowCID, 0, 0, 0, 0, 0, 9, 26)]
        [TestCase("AAAAAAAAAE", DissociationType.LowCID, 0, 0, 0, 0, 0, 9, 26)]
        [TestCase("AASAAAAAAA", DissociationType.LowCID, 0, 0, 7, 7, 0, 2, 33)]
        [TestCase("AATAAAAAAA", DissociationType.LowCID, 0, 0, 7, 7, 0, 2, 33)]
        [TestCase("AAEAAAAAAA", DissociationType.LowCID, 0, 0, 7, 7, 0, 2, 33)]
        [TestCase("AADAAAAAAA", DissociationType.LowCID, 0, 0, 7, 7, 0, 2, 33)]
        [TestCase("AARAAAASAA", DissociationType.LowCID, 7, 7, 2, 2, 2, 7, 44)]
        public static void Test_Fragment_ProductTypesWithAminoAcidSpecificities(string fullSequence, DissociationType dissociationType, int aStarCount, int bStarCount, int aDegreeCount, int bDegreeCount, int yStarCount, int yDegreeCount, int totalFragmentCount)
        {
            PeptideWithSetModifications myPeptide = new PeptideWithSetModifications(fullSequence, new Dictionary<string, Modification>());

            var theseTheoreticalFragments = new List<Product>();
            myPeptide.Fragment(dissociationType, FragmentationTerminus.Both, theseTheoreticalFragments);//Note that dissociation type here intentionally mismatched to dissociation type in modification constructor

            Assert.AreEqual(aStarCount, theseTheoreticalFragments.Where(f => f.ProductType == ProductType.aStar).Count());
            Assert.AreEqual(bStarCount, theseTheoreticalFragments.Where(f => f.ProductType == ProductType.bStar).Count());
            Assert.AreEqual(aDegreeCount, theseTheoreticalFragments.Where(f => f.ProductType == ProductType.aDegree).Count());
            Assert.AreEqual(bDegreeCount, theseTheoreticalFragments.Where(f => f.ProductType == ProductType.bDegree).Count());
            Assert.AreEqual(yStarCount, theseTheoreticalFragments.Where(f => f.ProductType == ProductType.yStar).Count());
            Assert.AreEqual(yDegreeCount, theseTheoreticalFragments.Where(f => f.ProductType == ProductType.yDegree).Count());
            Assert.AreEqual(totalFragmentCount, theseTheoreticalFragments.Count());
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

                    case ProductType.zDot:
                        Assert.AreEqual(Chemistry.ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O1N-1H-1").MonoisotopicMass + Constants.ProtonMass + Constants.ElectronMass).Value, mass);
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

            Assert.AreEqual(new List<ProductType> { ProductType.c, ProductType.y, ProductType.zDot }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.Both].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD]));
            Assert.AreEqual(new List<ProductType> { ProductType.c }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.N].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD]));
            Assert.AreEqual(new List<ProductType> { ProductType.y, ProductType.zDot }, TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.C].Intersect(DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD]));
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
            PeptideWithSetModifications p = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 7, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            var fragments = new List<Product>();

            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.b, ProductType.y }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.HCD, FragmentationTerminus.N, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.b }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.HCD, FragmentationTerminus.C, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.y }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.HCD, FragmentationTerminus.None, fragments);
            Assert.AreEqual(new List<ProductType> { }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.ETD, FragmentationTerminus.Both, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.c, ProductType.y, ProductType.zDot }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.ETD, FragmentationTerminus.N, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.c }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.ETD, FragmentationTerminus.C, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.y, ProductType.zDot }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.ETD, FragmentationTerminus.None, fragments);
            Assert.AreEqual(new List<ProductType> { }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.CID, FragmentationTerminus.Both, fragments);
            CollectionAssert.AreEquivalent(new List<ProductType> {  ProductType.b, ProductType.y }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.CID, FragmentationTerminus.N, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.b }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.CID, FragmentationTerminus.C, fragments);
            Assert.AreEqual(new List<ProductType> { ProductType.y }, fragments.Select(b => b.ProductType).Distinct().ToList());

            p.Fragment(DissociationType.CID, FragmentationTerminus.None, fragments);
            Assert.AreEqual(new List<ProductType> { }, fragments.Select(b => b.ProductType).Distinct().ToList());
        }

        [Test]
        public static void Test_MatchedFragmentIonToString()
        {
            Product P = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            MatchedFragmentIon m = new MatchedFragmentIon(ref P, 1, 1, 1);
            Assert.AreEqual("b1+1\t;1", m.ToString());
        }

        [Test]
        public static void Test_CID_Fragmentation_No_Unmodified_B1_ions()
        {
            //FOR CID B1 ions should always be missing whether or not there is a modification on first amino acid or not.

            Protein protein = new Protein("PEPTIDE", "accession");
            PeptideWithSetModifications p = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 7, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            var f = new List<Product>();
            p.Fragment(DissociationType.CID, FragmentationTerminus.Both, f);
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

            List<Product> modifiedPwsmFragments = new List<Product>();
            modifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedPwsmFragments);
            List<Product> unmodifiedPwsmFragments = new List<Product>();
            unmodifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedPwsmFragments);
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

            modifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedPwsmFragments);
            unmodifiedPwsm.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedPwsmFragments);
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
            List<Product> myFragments = new List<Product>();
            myPeptide.Fragment(dissociationType, FragmentationTerminus.Both, myFragments);
            Assert.AreEqual(fragmentCount, myFragments.Count());
        }

        [Test]
        public static void CheckProlineFragments()
        {
            PeptideWithSetModifications p = new PeptideWithSetModifications("MPEPTIDE", new Dictionary<string, Modification>());
            var fragments = new List<Product>();
            p.Fragment(DissociationType.ETD, FragmentationTerminus.Both, fragments);

            var z = fragments.Where(f => f.ProductType == ProductType.zDot).ToList();

            var ionNums = z.Select(f => f.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 6, 8 };

            Assert.That(expected.SequenceEqual(ionNums));
        }

        [Test]
        public static void CheckProlineFragments2()
        {
            PeptideWithSetModifications p = new PeptideWithSetModifications("MTETTIDE", new Dictionary<string, Modification>());
            var fragments = new List<Product>();
            p.Fragment(DissociationType.ETD, FragmentationTerminus.Both, fragments);

            var z = fragments.Where(f => f.ProductType == ProductType.zDot).ToList();

            var ionNums = z.Select(f => f.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 5, 6, 7, 8 };

            Assert.That(expected.SequenceEqual(ionNums));
        }

        [Test]
        public static void CheckProlineFragments3()
        {
            PeptideWithSetModifications p = new PeptideWithSetModifications("METPIPEEEE", new Dictionary<string, Modification>());
            var fragments = new List<Product>();
            p.Fragment(DissociationType.ETD, FragmentationTerminus.Both, fragments);

            var z = fragments.Where(f => f.ProductType == ProductType.zDot).ToList();

            var ionNums = z.Select(f => f.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 6, 8, 9, 10 };

            Assert.That(expected.SequenceEqual(ionNums));
        }

        [Test]
        public static void CheckProlineFragments4()
        {
            ModificationMotif.TryGetMotif("P", out var motif);
            Modification m = new Modification("TEST", "", "OK", null, motif, "Anywhere.", null, 20);
            PeptideWithSetModifications p = new PeptideWithSetModifications("METP[OK:TEST]IPEEEE", new Dictionary<string, Modification> { { "TEST", m } });
            var fragments = new List<Product>();
            p.Fragment(DissociationType.ETD, FragmentationTerminus.Both, fragments);

            var z = fragments.Where(f => f.ProductType == ProductType.zDot).ToList();

            var ionNums = z.Select(f => f.FragmentNumber).ToArray();
            var expected = new[] { 1, 2, 3, 4, 6, 8, 9, 10 };

            Assert.That(expected.SequenceEqual(ionNums));
        }

        [Test]
        public static void TestFragmentAnnotations()
        {
            Product p = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 30.3);
            MatchedFragmentIon f = new MatchedFragmentIon(ref p, 400.0, 1000.0, 3);

            Assert.That(p.Annotation == "b2-30.30");
            Assert.That(f.Annotation == "(b2-30.30)+3");

            p = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 0);
            f = new MatchedFragmentIon(ref p, 400.0, 1000.0, 3);

            Assert.That(p.Annotation == "b2");
            Assert.That(f.Annotation == "b2+3");
        }

        [Test]
        public static void TestFragmentErrors()
        {
            Product p = new Product(ProductType.b, FragmentationTerminus.N, 475.205, 2, 3, 30.3);
            MatchedFragmentIon f = new MatchedFragmentIon(ref p, 159.5, 1000.0, 3);

            double experMass = f.Mz.ToMass(f.Charge);
            double theorMass = p.NeutralMass;

            Assert.AreEqual(0.2732, Math.Round(experMass - theorMass, 4));
            Assert.AreEqual(574.85, Math.Round(f.MassErrorPpm, 2));
            Assert.AreEqual(0.2732, Math.Round(f.MassErrorDa, 4));
        }

        [Test]
        public static void TestFragmentEquality()
        {
            Product p1 = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 30.3);
            MatchedFragmentIon f1 = new MatchedFragmentIon(ref p1, 400.0, 1000.0, 3);

            Product p2 = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 30.3);
            MatchedFragmentIon f2 = new MatchedFragmentIon(ref p2, 400.0, 1000.0, 3);

            Assert.AreEqual(p1, p2);
            Assert.AreEqual(f1, f2);
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
    }
}