// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (TestPeptides.cs) is part of Proteomics.
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
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestPeptides
    {
        private Peptide _mockPeptideEveryAminoAcid;
        private Peptide _mockTrypticPeptide;

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
        public void AApolymerNullEquals()
        {
            Peptide pep = new Peptide("G");
            Assert.IsFalse(pep.Equals(null));
        }

        [Test]
        public void PeptideCountElements()
        {
            Peptide pep = new Peptide("G");
            pep.AddModification(new OldSchoolModification(1));
            Assert.AreEqual(5, pep.ElementCountWithIsotopes("H"));

            pep.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H{1}")));
            Assert.AreEqual(5, pep.ElementCountWithIsotopes("H")); // NOTHING HAS BEEN ADDED!

            pep.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H{1}"), ModificationSites.G));
            Assert.AreEqual(6, pep.ElementCountWithIsotopes("H"));

            Isotope isotope = PeriodicTable.GetElement("H").PrincipalIsotope;
            Assert.AreEqual(1, pep.SpecificIsotopeCount(isotope));
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
        public void SetAminoAcidModification()
        {
            var Asparagine = Residue.GetResidue("N");
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), Asparagine);

            Assert.AreEqual("ACDEFGHIKLMN[Fe]PQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTyped()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), ModificationSites.N);

            Assert.AreEqual("ACDEFGHIKLMN[Fe]PQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedMultipleLocations()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), ModificationSites.N | ModificationSites.F | ModificationSites.V);

            Assert.AreEqual("ACDEF[Fe]GHIKLMN[Fe]PQRSTV[Fe]WY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedAny()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), ModificationSites.Any);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedAll()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), ModificationSites.All);

            Assert.AreEqual("[Fe]-A[Fe]C[Fe]D[Fe]E[Fe]F[Fe]G[Fe]H[Fe]I[Fe]K[Fe]L[Fe]M[Fe]N[Fe]P[Fe]Q[Fe]R[Fe]S[Fe]T[Fe]V[Fe]W[Fe]Y[Fe]-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedNone()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), ModificationSites.None);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedTermini()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), ModificationSites.NPep | ModificationSites.PepC);

            Assert.AreEqual("[Fe]-ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidCharacterModification()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), 'D');

            Assert.AreEqual("ACD[Fe]EFGHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetResiduePositionModification()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), 5);

            Assert.AreEqual("ACDEF[Fe]GHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetResiduePositionModificationOutOfRangeUpper()
        {
            Assert.That(() => _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), 25),
                        Throws.TypeOf<MzLibException>());
        }

        [Test]
        public void SetResiduePositionModificationOutOfRangeLower()
        {
            Assert.That(() =>
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), 0),
                        Throws.TypeOf<MzLibException>());
        }

        [Test]
        public void SetCTerminusModStringRepresentation()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), Terminus.C);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetCTerminusModStringRepresentationofChemicalModification()
        {
            IHasChemicalFormula formula = new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe"), "Test");
            _mockPeptideEveryAminoAcid.SetModification(formula, Terminus.C);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY-[Test]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetNAndCTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), Terminus.C);
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H2NO")), Terminus.N);

            Assert.AreEqual("[H2NO]-ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetSameNAndCTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), Terminus.C | Terminus.N);

            Assert.AreEqual("[Fe]-ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void ClearNTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), Terminus.N);

            _mockPeptideEveryAminoAcid.ClearModifications(Terminus.N);

            Assert.IsNull(_mockPeptideEveryAminoAcid.NTerminusModification);
        }

        [Test]
        public void ClearCTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), Terminus.C);

            _mockPeptideEveryAminoAcid.ClearModifications(Terminus.C);

            Assert.IsNull(_mockPeptideEveryAminoAcid.CTerminusModification);
        }

        [Test]
        public void ClearAllMods()
        {
            _mockPeptideEveryAminoAcid.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Fe")), ModificationSites.All);

            _mockPeptideEveryAminoAcid.ClearModifications();

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void ClearModificationsBySites()
        {
            var peptide = new Peptide("AC[Fe]DEFGHIKLMNP[Fe]QRSTV[Fe]WY");

            peptide.ClearModifications(ModificationSites.C | ModificationSites.V);

            Assert.AreEqual("ACDEFGHIKLMNP[Fe]QRSTVWY", peptide.ToString());
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
        public void PeptideEquality()
        {
            Peptide pepA = new Peptide("DEREK");
            Peptide pepB = new Peptide("DEREK");
            Assert.AreEqual(pepA, pepB);

            Peptide pepC = new Peptide("DEREKK");
            Assert.AreNotEqual(pepA, pepC);
        }

        [Test]
        public void PeptideInEqualityAminoAcidSwitch()
        {
            Peptide pepA = new Peptide("DEREK");
            Peptide pepB = new Peptide("DEERK");
            Assert.AreNotEqual(pepA, pepB);
        }

        [Test]
        public void PeptideInEqualityAminoAcidModification()
        {
            Peptide pepA = new Peptide("DEREK");
            Peptide pepB = new Peptide("DEREK");
            pepB.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H2O")), 'R');
            Assert.AreNotEqual(pepA, pepB);
            pepA.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H2O2")), 'R');
            Assert.AreNotEqual(pepA, pepB);
        }

        [Test]
        public void PeptideCloneEquality()
        {
            Peptide pepA = new Peptide("DEREK");
            Peptide pepB = new Peptide(pepA);
            Assert.AreEqual(pepA, pepB);
        }

        [Test]
        public void PeptideCloneNotSameReference()
        {
            Peptide pepA = new Peptide("DEREK");
            Peptide pepB = new Peptide(pepA);

            Assert.AreNotSame(pepA, pepB);
        }

        [Test]
        public void PeptideCloneWithModifications()
        {
            Peptide pepA = new Peptide("DER[Fe]EK");
            Peptide pepB = new Peptide(pepA);
            Assert.AreEqual("DER[Fe]EK", pepB.ToString());
        }

        [Test]
        public void PeptideCloneWithoutModifications()
        {
            Peptide pepA = new Peptide("DER[Fe]EK");
            Peptide pepB = new Peptide(pepA, false);
            Assert.AreEqual("DEREK", pepB.ToString());
        }

        [Test]
        public void PeptideCloneWithModification()
        {
            Peptide pepA = new Peptide("DEREK");
            pepA.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H2O")), 'R');
            Peptide pepB = new Peptide(pepA);

            Assert.AreEqual(pepB, pepA);
        }

        [Test]
        public void PeptideParitalCloneInternal()
        {
            Peptide pepA = new Peptide("DEREK");
            Peptide pepB = new Peptide(pepA, 1, 3);
            Peptide pepC = new Peptide("ERE");
            Assert.AreEqual(pepB, pepC);
        }

        [Test]
        public void PeptideParitalClonelWithInternalModification()
        {
            Peptide pepA = new Peptide("DER[Fe]EK");
            Peptide pepB = new Peptide(pepA, 2, 3);
            Peptide pepC = new Peptide("R[Fe]EK");
            Assert.AreEqual(pepB, pepC);
        }

        [Test]
        public void PeptideHashing()
        {
            Peptide pep1 = new Peptide("DEREK");
            Peptide pep2 = new Peptide("DEREKN");
            Peptide pep3 = new Peptide("DEREKM");
            Peptide pep4 = new Peptide("DEREKM");
            HashSet<Peptide> uu = new HashSet<Peptide> { pep1, pep2, pep3, pep4 };
            uu.Add(new Peptide("DEREKN"));
            Assert.AreEqual(3, uu.Count);
        }

        [Test]
        public void ClearMods()
        {
            Peptide pepA = new Peptide("DE[Al]R[Fe]EK");
            pepA.ClearModifications(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("Al")));
            Assert.AreEqual("DER[Fe]EK", pepA.ToString());
            pepA.ClearModifications(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("C")));
            Assert.AreEqual("DER[Fe]EK", pepA.ToString());
            pepA.ClearModifications();
            Assert.AreEqual("DEREK", pepA.ToString());
            pepA.ClearModifications();
            Assert.AreEqual("DEREK", pepA.ToString());
            pepA.ClearModifications(ModificationSites.Any);
            Assert.AreEqual("DEREK", pepA.ToString());
            pepA.ClearModifications(Terminus.C);
            Assert.AreEqual("DEREK", pepA.ToString());
        }

        [Test]
        public void PeptideParitalClonelWithInternalModificationTwoMods()
        {
            Peptide pepA = new Peptide("DE[Al]R[Fe]EK");
            Peptide pepB = new Peptide(pepA, 2, 3);
            Peptide pepC = new Peptide("R[Fe]EK");
            Assert.AreEqual(pepB, pepC);
        }

        [Test]
        public void PeptideParitalCloneInternalWithCTerminusModification()
        {
            Peptide pepA = new Peptide("DEREK");
            pepA.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H2O")), Terminus.C);
            Peptide pepB = new Peptide(pepA, 2, 3);

            Peptide pepC = new Peptide("REK");
            pepC.SetModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H2O")), Terminus.C);

            Assert.AreEqual(pepC, pepB);
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
        public void TerminusModification()
        {
            Peptide pepA = new Peptide("DERLEK");
            pepA.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("SO")), Terminus.N);
            Assert.AreEqual("[OS]-DERLEK", pepA.ToString());
        }

        [Test]
        public void DigestionTest()
        {
            IProtease protease = new TestProtease();
            Assert.AreEqual(6, AminoAcidPolymer.Digest(_mockPeptideEveryAminoAcid, protease).Count());
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
        public void TestMultipleModificationsAtSingleResidue()
        {
            Peptide a = new Peptide("ACDEFGHIKLMNPQRSTVWY");
            a.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("O"), ModificationSites.D));
            a.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H"), ModificationSites.D));
            var products = a.Fragment(FragmentTypes.b | FragmentTypes.y, true);
            foreach (Fragment fragment in products)
            {
            }
            foreach (IHasChemicalFormula fragment in products)
            {
            }
        }

        [Test]
        public void TestAApolymerContains()
        {
            Assert.IsFalse(_mockTrypticPeptide.Contains('A'));
            Assert.IsTrue(_mockTrypticPeptide.Contains(Residue.GetResidue('T')));
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
        public void TestClearModifications()
        {
            Peptide a = new Peptide("ACDEFGHIKLMNPQRSTVWY");
            a.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("O"), ModificationSites.D));
            a.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H"), ModificationSites.E));
            Assert.AreEqual(2, a.ModificationCount());
            a.ClearModifications();
            Assert.AreEqual(0, a.ModificationCount());
            a.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("O"), ModificationSites.NTerminus));
            a.AddModification(new OldSchoolModification(1), ModificationSites.TerminusC);
            Assert.AreEqual(2, a.ModificationCount());
            a.Fragment(FragmentTypes.y);

            Peptide peptide = new Peptide("[C2H3NO]-LLL-[C2H3NO]");
            ModificationSites ff = ModificationSites.NPep | ModificationSites.PepC;
            peptide.ClearModifications(ff);
            Assert.AreEqual("LLL", peptide.GetSequenceWithModifications());
        }

        [Test]
        public void TestGetSubPeptide()
        {
            Peptide pep = new Peptide("DERLEK");
            Assert.AreEqual(new Peptide("LE"), pep.GetSubPeptide(3, 2));
        }

        [Test]
        public void TestRealPeptideWithModifications()
        {
            Peptide a = new Peptide("LDNLQQEIDFLTALYQAELSQM[O]QTQISETNVILSM[O]DNNR");
            Assert.AreEqual(2, a.ModificationCount());
        }

        [Test]
        public void TestGetDigestionPointsWithMethionine()
        {
            var ok = AminoAcidPolymer.GetDigestionPointsAndLengths("MDERLEKDERLE", new List<TestProtease> { new TestProtease() }, 0, 0, 10000, true, false).ToList();
            Assert.AreEqual(1, ok[0].Index); // Methionine cleaved, digestion is at 1
            Assert.AreEqual(4, ok[0].Length); // The test protease cleaves at index 4, so after L.
            Assert.AreEqual(0, ok[1].Index); // Regular digestion 1
            Assert.AreEqual(5, ok[1].Length); // Regular digestion 1
            Assert.AreEqual(5, ok[2].Index); // Regular digestion 2
            Assert.AreEqual(1, ok[2].Length); // Regular digestion 2
            Assert.AreEqual(6, ok[3].Index); // Regular digestion 3
            Assert.AreEqual(6, ok[3].Length); // Regular digestion 3
        }

        [Test]
        public void TestGetDigestionPointsWithMethionineAndSemiDigestion()
        {
            var ok = AminoAcidPolymer.GetDigestionPointsAndLengths("MDERLEK", new List<TestProtease> { new TestProtease() }, 0, 0, 10000, true, true).ToList();

            IEqualityComparer<DigestionPointAndLength> jj = new OkComparer();
            var yee = new HashSet<DigestionPointAndLength>(ok, jj);

            Assert.AreEqual(1 + 3 + 1 + (8 - 1) + 1 + 1, yee.Count);
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

        private class OkComparer : IEqualityComparer<DigestionPointAndLength>
        {
            public bool Equals(DigestionPointAndLength x, DigestionPointAndLength y)
            {
                return x.Index.Equals(y.Index) && x.Length.Equals(y.Length);
            }

            public int GetHashCode(DigestionPointAndLength obj)
            {
                return obj.Length + obj.Index * 256;
            }
        }
    }

    internal class TestProtease : IProtease
    {
        public IEnumerable<int> GetDigestionSites(AminoAcidPolymer aminoAcidSequence)
        {
            return GetDigestionSites(aminoAcidSequence.BaseSequence);
        }

        public IEnumerable<int> GetDigestionSites(string aminoAcidSequence)
        {
            yield return 4;
            yield return 5;
        }

        public int MissedCleavages(AminoAcidPolymer aminoAcidSequence)
        {
            throw new NotImplementedException();
        }

        public int MissedCleavages(string sequence)
        {
            throw new NotImplementedException();
        }
    }
}