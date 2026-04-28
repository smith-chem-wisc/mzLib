// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (PeptideDigestionTests.cs) is part of Proteomics.
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
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class PeptideDigestionTests
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
        public void DigestionTest()
        {
            IProtease protease = new TestProtease();
            Assert.AreEqual(6, AminoAcidPolymer.Digest(_mockPeptideEveryAminoAcid, protease).Count());
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
        public void TestNonSpecificOverride()
        {
            string trypsin = "trypsin";
            DigestionParams digestionParams = new DigestionParams(trypsin);
            Assert.AreEqual(digestionParams.Protease.Name, trypsin);

            digestionParams = new DigestionParams(trypsin, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            Assert.AreEqual(digestionParams.Protease.Name, "singleN");
            Assert.AreEqual(digestionParams.SpecificProtease.Name, trypsin);

            digestionParams = new DigestionParams(trypsin, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            Assert.AreEqual(digestionParams.Protease.Name, "singleC");
            Assert.AreEqual(digestionParams.SpecificProtease.Name, trypsin);
        }

        [Test]
        public void TestAApolymerContains()
        {
            Assert.IsFalse(_mockTrypticPeptide.Contains('A'));
            Assert.IsTrue(_mockTrypticPeptide.Contains(Residue.GetResidue('T')));
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
}
