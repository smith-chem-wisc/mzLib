// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (PeptideConstructionTests.cs) is part of Proteomics.
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
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using Test.ChemistryTests;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class PeptideConstructionTests
    {
        private Peptide _mockPeptideEveryAminoAcid;
        private Peptide _mockTrypticPeptide;
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

        [SetUp]
        public void SetUp()
        {
            _mockPeptideEveryAminoAcid = new Peptide("ACDEFGHIKLMNPQRSTVWY");
            _mockTrypticPeptide = new Peptide("TTGSSSSSSSK");
        }

        [Test]
        public void PeptideTestReal()
        {
            _mockPeptideEveryAminoAcid = new Peptide("LDNLQQEIDFLTALYQAELSQMQTQISETNVILSMDNNR");
        }

        [Test]
        public void PeptideMassGlycine()
        {
            Peptide pep = new Peptide("G");
            ChemicalFormula formula = new ChemicalFormula(ChemicalFormula.ParseFormula("C2H5NO2"));
            ChemicalFormula formula2;
            formula2 = pep.GetChemicalFormula();

            Assert.AreEqual(formula, formula2);
        }

        [Test]
        public void PeptideMassTryptic()
        {
            ChemicalFormula formula = new ChemicalFormula(ChemicalFormula.ParseFormula("C37H66N12O21"));
            ChemicalFormula formula2;
            formula2 = _mockTrypticPeptide.GetChemicalFormula();
            Assert.AreEqual(formula, formula2);
        }

        [Test]
        public void PeptideAminoAcidCount()
        {
            Assert.AreEqual(20, _mockPeptideEveryAminoAcid.Length);
        }

        [Test]
        public void EmptyStringPeptideConstructorLength()
        {
            Peptide peptide = new Peptide();

            Assert.AreEqual(0, peptide.Length);
        }

        [Test]
        public void EmptyStringPeptideConstructorToString()
        {
            Peptide peptide = new Peptide();

            Assert.AreEqual(string.Empty, peptide.ToString());
        }

        [Test]
        public void ParseNTerminalChemicalFormula()
        {
            Peptide peptide = new Peptide("[C2H3NO]-TTGSSSSSSSK");
            ChemicalFormula formulaA = new ChemicalFormula(ChemicalFormula.ParseFormula("C39H69N13O22"));
            ChemicalFormula formulaB;
            formulaB = peptide.GetChemicalFormula();

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParseCTerminalChemicalFormula()
        {
            Peptide peptide = new Peptide("TTGSSSSSSSK-[C2H3NO]");
            ChemicalFormula formulaA = new ChemicalFormula(ChemicalFormula.ParseFormula("C39H69N13O22"));
            ChemicalFormula formulaB;
            formulaB = peptide.GetChemicalFormula();

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParseCTerminalChemicalFormulaWithLastResidueMod()
        {
            Peptide peptide = new Peptide("TTGSSSSSSSK[H2O]-[C2H3NO]");
            ChemicalFormula formulaA = new ChemicalFormula(ChemicalFormula.ParseFormula("C39H71N13O23"));
            ChemicalFormula formulaB;
            formulaB = peptide.GetChemicalFormula();

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParseCTerminalChemicalFormulaWithLastResidueModStringRepresentation()
        {
            Peptide peptide = new Peptide("TTGSSSSSSSK[H2O]-[C2H3NO]");

            Assert.AreEqual("TTGSSSSSSSK[H2O]-[C2H3NO]", peptide.GetSequenceWithModifications());

            peptide.NTerminus = new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("N"));

            Assert.AreEqual("TTGSSSSSSSK[H2O]-[C2H3NO]", peptide.GetSequenceWithModifications());

            ChemicalFormula formulaA = new ChemicalFormula(ChemicalFormula.ParseFormula("C39H70N14O23"));
            var formulaB = peptide.GetChemicalFormula();
            Assert.AreEqual(formulaA, formulaB);

            peptide.AddModification(new ObjectWithMass100(), 0);

            Assert.AreEqual("[mass: 100]-TTGSSSSSSSK[H2O]-[C2H3NO]", peptide.GetSequenceWithModifications());

            Assert.AreEqual(1, peptide.AddModification(new ObjectWithMass100(), Terminus.C));

            Assert.AreEqual(3, peptide.ModificationCount());

            Assert.AreEqual(0, peptide.ReplaceModification(new ObjectWithMass100(), new ObjectWithMass100()));

            Assert.That(() => peptide.ReplaceModification(null, new ObjectWithMass100()),
            Throws.TypeOf<MzLibException>()
            .With.Property("Message")
            .EqualTo("Cannot replace a null modification"));

            peptide.SetModification(new ObjectWithMass100(), new int[] { 1, 11 });
            Assert.AreEqual(4, peptide.ModificationCount());

            OldSchoolModification mod1 = new OldSchoolModification(5, "mass 5 on T", ModificationSites.T);
            peptide.SetModifications(new List<OldSchoolModification> { mod1 });
            Assert.AreEqual(5, peptide.ModificationCount());
        }

        [Test]
        public void ParseNAndCTerminalChemicalFormula()
        {
            Peptide peptide = new Peptide("[C2H3NO]-TTGSSSSSSSK-[C2H3NO]");
            ChemicalFormula formulaA = new ChemicalFormula(ChemicalFormula.ParseFormula("C41H72N14O23"));
            ChemicalFormula formulaB;
            formulaB = peptide.GetChemicalFormula();

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParseDoubleModificationToString()
        {
            Peptide peptide = new Peptide("THGEAK[25.132]K");

            Assert.AreEqual("THGEAK[25.132]K", peptide.ToString());
        }

        [Test]
        public void ParseNamedChemicalModificationInvalidName()
        {
            Assert.That(() => new Peptide("T[TMT 7-plex]HGEAK[Acetyl]K"),
                        Throws.TypeOf<MzLibException>());
        }

        [Test]
        public void EmptyPeptideLengthIsZero()
        {
            Peptide pepA = new Peptide();

            Assert.AreEqual(0, pepA.Length);
        }

        [Test]
        public void EmptyPeptideSequenceIsEmpty()
        {
            Peptide pepA = new Peptide();

            Assert.AreEqual(string.Empty, pepA.BaseSequence);
        }

        [Test]
        public void EmptyPeptideFormulaIsH2O()
        {
            Peptide pepA = new Peptide();
            ChemicalFormula h2O = new ChemicalFormula(ChemicalFormula.ParseFormula("H2O"));
            ChemicalFormula formulaB;
            formulaB = pepA.GetChemicalFormula();

            Assert.AreEqual(h2O, formulaB);
        }

        [Test]
        public void BadSeqeunce()
        {
            Assert.That(() => new Peptide("ABC"), Throws.TypeOf<MzLibException>()
            .With.Property("Message")
            .EqualTo("Amino Acid Letter B does not exist in the Amino Acid Dictionary. B is also not a valid character"));

            Assert.That(() => new Peptide("A["), Throws.TypeOf<MzLibException>()
            .With.Property("Message")
            .EqualTo("Couldn't find the closing ] for a modification in this sequence: A["));
        }

        [Test]
        public static void TestPyrrolysine()
        {
            Peptide mom = new Peptide("MOM");
            double momMass = mom.MonoisotopicMass;
            Assert.That(517.23926172577, Is.EqualTo(momMass).Within(0.001));

            Peptide mm = new Peptide("MM");
            double mmMass = mm.MonoisotopicMass;

            double deltaMass = momMass - mmMass;
            Assert.That(237.1477268648, Is.EqualTo(deltaMass).Within(0.001));
        }

        [TestCase("PEPTIDEK", ExpectedResult = true)]
        [TestCase("PEPTUDEK", ExpectedResult = true)] // U is selenocysteine
        [TestCase("R", ExpectedResult = true)]
        [TestCase("PEPTJDEK", ExpectedResult = false)]
        [TestCase("peptidek", ExpectedResult = false)]
        [TestCase("", ExpectedResult = false)]
        [TestCase("PEPTIDEK ", ExpectedResult = false)]
        [TestCase("P3PT1D3K", ExpectedResult = false)]
        [TestCase("PEP-TIDEK", ExpectedResult = false)]
        [TestCase("PEP_TIDEK", ExpectedResult = false)]
        [TestCase(".PEPTIDEK", ExpectedResult = false)]
        [TestCase("PEPT[]IDEK", ExpectedResult = false)]
        public bool TestValidBaseSequence(string sequence)
        {
            return sequence.AllSequenceResiduesAreValid();
        }

        [Test]
        public void TestValidBaseSequenceWithResidueAddedToResidueDictionary()
        {
            string testSequenceForThisTest = "PEPTIDEa";

            Residue x = new Residue("a", 'a', "a", new ChemicalFormula(), ModificationSites.All);
            Residue.AddNewResiduesToDictionary(new List<Residue> { x });

            Assert.IsTrue(testSequenceForThisTest.AllSequenceResiduesAreValid());
        }

        [Test]
        public void GetLeucineSequence()
        {
            Peptide pepA = new Peptide("DERIEK");
            string leuSeq = pepA.BaseLeucineSequence;

            Assert.AreEqual("DERLEK", leuSeq);
        }

        [Test]
        public void GetLeucineSequenceNoReplacement()
        {
            Peptide pepA = new Peptide("DERLEK");

            string leuSeq = pepA.BaseLeucineSequence;

            Assert.AreEqual("DERLEK", leuSeq);
        }

        [Test]
        public void GetSequenceCoverage()
        {
            Peptide pepA = new Peptide("DERLEK");
            Peptide pepAa = new Peptide("ER");
            Peptide pepAb = new Peptide("RL");
            Peptide pepAc = new Peptide("LEK");
            List<Peptide> myList = new List<Peptide>
            {
                pepAa,
                pepAb,
                pepAc
            };
            Assert.IsTrue(pepA.GetSequenceCoverage(myList).SequenceEqual(new List<int> { 0, 1, 2, 2, 1, 1 }));
        }

        [Test]
        public void GetSequenceCoverageFraction()
        {
            Peptide pepA = new Peptide("DERLEK");
            Peptide pepAa = new Peptide("ER");
            Peptide pepAb = new Peptide("RL");
            List<Peptide> myList = new List<Peptide>
            {
                pepAa,
                pepAb
            };
            Assert.AreEqual(0.5, pepA.GetSequenceCoverageFraction(myList));
        }

        [Test]
        public void TestChemicalFormula()
        {
            Residue.GetResidue('A');

            Peptide A = new Peptide("A");

            Residue.GetResidue('A');

            ChemicalFormula ok = new ChemicalFormula(Residue.GetResidue('A').ThisChemicalFormula);
            ok.Add(new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("OH")));
            ok.Add(new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("H")));

            Residue.GetResidue('A');

            Residue.GetResidue('A');

            Assert.AreEqual(ok, A.GetChemicalFormula());
        }

        [Test]
        public void TestChemicalFormula2()
        {
            Peptide A = new Peptide("A");
            OldSchoolModification a = new OldSchoolModification(1, "Modification without chemical formula", ModificationSites.A);
            A.AddModification(a);
            Assert.Throws<MzLibException>(() => { A.GetChemicalFormula(); }, "Modification Modification without chemical formula does not have a chemical formula!");
        }

        [Test]
        public void TestLeucineSequence()
        {
            Assert.AreEqual("ACDEFGHLKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.GetSequenceWithModifications(true));
            Assert.AreEqual(20, _mockPeptideEveryAminoAcid.ResidueCount());
            Assert.AreEqual(7, _mockTrypticPeptide.ResidueCount('S'));
            Assert.AreEqual(7, _mockTrypticPeptide.ResidueCount(Residue.GetResidue('S')));
            Assert.AreEqual(2, _mockTrypticPeptide.ResidueCount(Residue.GetResidue('S'), 2, 3));
            Assert.AreEqual(3, _mockTrypticPeptide.ResidueCount('S', 2, 4));

            Peptide peptide = new Peptide("III-[C2H3NO]");
            Assert.AreEqual("LLL-[C2H3NO]", peptide.GetSequenceWithModifications(true));
        }

        [Test]
        public void GenerateIsotopologues()
        {
            Peptide pep = new Peptide("DERLEK");
            var a = pep.GenerateAllModificationCombinations().ToArray();
            Assert.AreEqual(0, a.Count());
            var i = new ModificationWithMultiplePossibilitiesCollection("My Iso Mod", ModificationSites.E);
            i.AddModification(new OldSchoolModification(1, "My Mod1a", ModificationSites.E));
            i.AddModification(new OldSchoolModification(2, "My Mod2b", ModificationSites.E));
            pep.SetModification(i);
            var i2 = new ModificationWithMultiplePossibilitiesCollection("My Iso Mod2", ModificationSites.R);
            i2.AddModification(new OldSchoolModification(1, "My Mod2a", ModificationSites.R));
            i2.AddModification(new OldSchoolModification(2, "My Mod2b", ModificationSites.R));
            i2.AddModification(new OldSchoolModification(3, "My Mod2c", ModificationSites.R));
            pep.SetModification(i2);
            a = pep.GenerateAllModificationCombinations().ToArray();
            // Only 6 and not 12, because in the first modification, it is one choice that is substituted across all modification sites
            Assert.AreEqual(6, a.Count());
        }
    }
}
