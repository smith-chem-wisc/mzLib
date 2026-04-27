// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (PeptideModificationTests.cs) is part of Proteomics.
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
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Fragmentation;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class PeptideModificationTests
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
        public void TerminusModification()
        {
            Peptide pepA = new Peptide("DERLEK");
            pepA.AddModification(new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("SO")), Terminus.N);
            Assert.AreEqual("[OS]-DERLEK", pepA.ToString());
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
        public void TestRealPeptideWithModifications()
        {
            Peptide a = new Peptide("LDNLQQEIDFLTALYQAELSQM[O]QTQISETNVILSM[O]DNNR");
            Assert.AreEqual(2, a.ModificationCount());
        }
    }
}
