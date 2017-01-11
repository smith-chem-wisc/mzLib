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
using NUnit.Framework;
using Proteomics;
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
            ChemicalFormula formula = new ChemicalFormula("C2H5NO2");
            ChemicalFormula formula2;
            formula2 = pep.GetChemicalFormula();

            Assert.AreEqual(formula, formula2);
        }

        [Test]
        public void PeptideCountElements()
        {
            Peptide pep = new Peptide("G");
            pep.AddModification(new Modification(1));
            Assert.AreEqual(5, pep.ElementCountWithIsotopes("H"));
            Isotope isotope = PeriodicTable.GetElement("H").PrincipalIsotope;
            Assert.AreEqual(0, pep.SpecificIsotopeCount(isotope));

        }

        [Test]
        public void PeptideMassTryptic()
        {
            ChemicalFormula formula = new ChemicalFormula("C37H66N12O21");
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
            ChemicalFormula formulaA = new ChemicalFormula("C39H69N13O22");
            ChemicalFormula formulaB;
            formulaB = peptide.GetChemicalFormula();

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParseCTerminalChemicalFormula()
        {
            Peptide peptide = new Peptide("TTGSSSSSSSK-[C2H3NO]");
            ChemicalFormula formulaA = new ChemicalFormula("C39H69N13O22");
            ChemicalFormula formulaB;
            formulaB = peptide.GetChemicalFormula();

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParseCTerminalChemicalFormulaWithLastResidueMod()
        {
            Peptide peptide = new Peptide("TTGSSSSSSSK[H2O]-[C2H3NO]");
            ChemicalFormula formulaA = new ChemicalFormula("C39H71N13O23");
            ChemicalFormula formulaB;
            formulaB = peptide.GetChemicalFormula();

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParseCTerminalChemicalFormulaWithLastResidueModStringRepresentation()
        {
            Peptide peptide = new Peptide("TTGSSSSSSSK[H2O]-[C2H3NO]");

            Assert.AreEqual("TTGSSSSSSSK[H2O]-[C2H3NO]", peptide.GetSequenceWithModifications());
        }

        [Test]
        public void ParseNAndCTerminalChemicalFormula()
        {
            Peptide peptide = new Peptide("[C2H3NO]-TTGSSSSSSSK-[C2H3NO]");
            ChemicalFormula formulaA = new ChemicalFormula("C41H72N14O23");
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
                        Throws.TypeOf<ArgumentException>());
        }

        [Test]
        public void SetAminoAcidModification()
        {
            var Asparagine = Residue.GetResidue("N");
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), Asparagine);

            Assert.AreEqual("ACDEFGHIKLMN[Fe]PQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTyped()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), ModificationSites.N);

            Assert.AreEqual("ACDEFGHIKLMN[Fe]PQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedMultipleLocations()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), ModificationSites.N | ModificationSites.F | ModificationSites.V);

            Assert.AreEqual("ACDEF[Fe]GHIKLMN[Fe]PQRSTV[Fe]WY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedAny()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), ModificationSites.Any);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedAll()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), ModificationSites.All);

            Assert.AreEqual("[Fe]-A[Fe]C[Fe]D[Fe]E[Fe]F[Fe]G[Fe]H[Fe]I[Fe]K[Fe]L[Fe]M[Fe]N[Fe]P[Fe]Q[Fe]R[Fe]S[Fe]T[Fe]V[Fe]W[Fe]Y[Fe]-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedNone()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), ModificationSites.None);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidModificationStronglyTypedTermini()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), ModificationSites.NPep | ModificationSites.PepC);

            Assert.AreEqual("[Fe]-ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetAminoAcidCharacterModification()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), 'D');

            Assert.AreEqual("ACD[Fe]EFGHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetResiduePositionModification()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), 5);

            Assert.AreEqual("ACDEF[Fe]GHIKLMNPQRSTVWY", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetResiduePositionModificationOutOfRangeUpper()
        {
            Assert.That(() => _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), 25),
                        Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        [Test]
        public void SetResiduePositionModificationOutOfRangeLower()
        {
            Assert.That(() =>
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), 0),
                        Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        [Test]
        public void SetCTerminusModStringRepresentation()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), Terminus.C);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetCTerminusModStringRepresentationofChemicalModification()
        {
            IHasChemicalFormula formula = new ChemicalFormulaModification("Fe", "Test");
            _mockPeptideEveryAminoAcid.SetModification(formula, Terminus.C);

            Assert.AreEqual("ACDEFGHIKLMNPQRSTVWY-[Test]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetNAndCTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), Terminus.C);
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("H2NO"), Terminus.N);

            Assert.AreEqual("[H2NO]-ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void SetSameNAndCTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), Terminus.C | Terminus.N);

            Assert.AreEqual("[Fe]-ACDEFGHIKLMNPQRSTVWY-[Fe]", _mockPeptideEveryAminoAcid.ToString());
        }

        [Test]
        public void ClearNTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), Terminus.N);

            _mockPeptideEveryAminoAcid.ClearModifications(Terminus.N);

            Assert.IsNull(_mockPeptideEveryAminoAcid.NTerminusModification);
        }

        [Test]
        public void ClearCTerminusMod()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), Terminus.C);

            _mockPeptideEveryAminoAcid.ClearModifications(Terminus.C);

            Assert.IsNull(_mockPeptideEveryAminoAcid.CTerminusModification);
        }

        [Test]
        public void ClearAllMods()
        {
            _mockPeptideEveryAminoAcid.SetModification(new ChemicalFormulaModification("Fe"), ModificationSites.All);

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

            Assert.AreEqual(String.Empty, pepA.BaseSequence);
        }

        [Test]
        public void EmptyPeptideFormulaIsH2O()
        {
            Peptide pepA = new Peptide();
            ChemicalFormula h2O = new ChemicalFormula("H2O");
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
            pepB.SetModification(new ChemicalFormulaModification("H2O"), 'R');

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
            pepA.SetModification(new ChemicalFormulaModification("H2O"), 'R');
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
            pepA.SetModification(new ChemicalFormulaModification("H2O"), Terminus.C);
            Peptide pepB = new Peptide(pepA, 2, 3);

            Peptide pepC = new Peptide("REK");
            pepC.SetModification(new ChemicalFormulaModification("H2O"), Terminus.C);

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
            List<Peptide> myList = new List<Peptide>();
            myList.Add(pepAa);
            myList.Add(pepAb);
            myList.Add(pepAc);
            Assert.IsTrue(pepA.GetSequenceCoverage(myList).SequenceEqual(new List<int> { 0, 1, 2, 2, 1, 1 }));
        }

        [Test]
        public void GenerateIsotopologues()
        {
            Peptide pep = new Peptide("DERLEK");
            var a = pep.GenerateAllModificationCombinations().ToArray();
            Assert.AreEqual(0, a.Count());
            var i = new ModificationWithMultiplePossibilitiesCollection("My Iso Mod", ModificationSites.E);
            i.AddModification(new Modification(1, "My Mod1a", ModificationSites.E));
            i.AddModification(new Modification(2, "My Mod2b", ModificationSites.E));
            pep.SetModification(i);
            var i2 = new ModificationWithMultiplePossibilitiesCollection("My Iso Mod2", ModificationSites.R);
            i2.AddModification(new Modification(1, "My Mod2a", ModificationSites.R));
            i2.AddModification(new Modification(2, "My Mod2b", ModificationSites.R));
            i2.AddModification(new Modification(3, "My Mod2c", ModificationSites.R));
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
            List<Peptide> myList = new List<Peptide>();
            myList.Add(pepAa);
            myList.Add(pepAb);
            Assert.AreEqual(0.5, pepA.GetSequenceCoverageFraction(myList));
        }

        [Test]
        public void TerminusModification()
        {
            Peptide pepA = new Peptide("DERLEK");
            pepA.AddModification(new ChemicalFormulaModification("SO"), Terminus.N);
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
            ok.Add(new ChemicalFormulaTerminus("OH"));
            ok.Add(new ChemicalFormulaTerminus("H"));

            Residue.GetResidue('A');

            Residue.GetResidue('A');

            Assert.AreEqual(ok, A.GetChemicalFormula());
        }

        [Test]
        public void TestChemicalFormula2()
        {
            Peptide A = new Peptide("A");
            Modification a = new Modification(1, "Modification without chemical formula", ModificationSites.A);
            A.AddModification(a);
            Assert.Throws<InvalidCastException>(() => { A.GetChemicalFormula(); }, "Modification Modification without chemical formula does not have a chemical formula!");
        }

        [Test]
        public void TestMultipleModificationsAtSingleResidue()
        {
            Peptide a = new Peptide("ACDEFGHIKLMNPQRSTVWY");
            a.AddModification(new ChemicalFormulaModification("O", ModificationSites.D));
            a.AddModification(new ChemicalFormulaModification("H", ModificationSites.D));
            var products = a.Fragment(FragmentTypes.b | FragmentTypes.y, true);
            foreach (Fragment fragment in products)
            {
            }
            foreach (IHasChemicalFormula fragment in products)
            {
            }
        }

        [Test]
        public void TestAApolymerOutOfRangeInitialization()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => { new Peptide(_mockPeptideEveryAminoAcid, -1, 0, false); }, "The first residue index is outside the valid range [0-21]");
            Assert.Throws<ArgumentOutOfRangeException>(() => { new Peptide(_mockPeptideEveryAminoAcid, 0, 100, false); }, "The length + firstResidue value is too large");
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
            a.AddModification(new ChemicalFormulaModification("O", ModificationSites.D));
            a.AddModification(new ChemicalFormulaModification("H", ModificationSites.E));
            Assert.AreEqual(2, a.ModificationCount());
            a.ClearModifications();
            Assert.AreEqual(0, a.ModificationCount());
            a.AddModification(new ChemicalFormulaModification("O", ModificationSites.NTerminus));
            a.AddModification(new Modification(1), ModificationSites.TerminusC);
            Assert.AreEqual(2, a.ModificationCount());
            a.Fragment(FragmentTypes.y);
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