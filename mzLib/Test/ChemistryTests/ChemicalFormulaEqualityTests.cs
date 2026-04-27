// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaEqualityTests.cs) is part of Chemistry Library.
//
// Chemistry Library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chemistry Library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Chemistry Library. If not, see <http://www.gnu.org/licenses/>

using System;
using System.Diagnostics.CodeAnalysis;
using Chemistry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ChemistryTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ChemicalFormulaEqualityTests
    {
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

        [Test]
        public static void FormulaValueInequality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("NC1OH3");

            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public static void FormulaValueInequalityHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public static void FormulaValueEqualityItself()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public static void FormulaValueEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("NC2OH3");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void FormulaEquality_AreEqual_SameReference()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public static void FormulaEquality_AreEqual_DifferentReference()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void FormulaEquality_Rearranged()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3NOC2");
            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void FormulaEquality_null()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.AreNotEqual(formulaA, null);
        }

        [Test]
        public static void FormulaEquality_OtherType()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.AreNotEqual(formulaA, "C2H3NO");
        }

        [Test]
        public static void FormulaAlmostEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C{12}2H3NO");
            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public static void FormulaEquality_AreEqual_DifferentReference_Interface()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            IHasChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.IsTrue(formulaA.Equals(formulaB));
        }

        [Test]
        public static void FormulaEquality_Rearranged_Interface()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            IHasChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3NOC2");
            Assert.IsTrue(formulaA.Equals(formulaB));
        }

        [Test]
        public static void FormulaAlmostEquality_Interface()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            IHasChemicalFormula formulaB = ChemicalFormula.ParseFormula("C{12}2H3NO");
            Assert.IsFalse(formulaA.Equals(formulaB));
        }

        [Test]
        public static void HashCodeEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3C2NO");

            Assert.AreEqual(formulaA.GetHashCode(), formulaB.GetHashCode());
        }

        [Test]
        public static void HashCodeImmutableEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.AreEqual(formulaA.GetHashCode(), formulaA.GetHashCode());
        }

        [Test]
        public static void HashCodeCheck()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("Al");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("Al{27}");
            Assert.AreNotEqual(formulaA.GetHashCode(), formulaB.GetHashCode());
        }

        [Test]
        public static void EqualsFalse()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("OCHHCHN");
            Assert.IsFalse(formulaA.Equals("C"));
        }

        [Test]
        public static void Equals()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("OCHHCHN");
            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public static void NotEqual()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C15O15H15S15N15");
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("N15S15H15O15C15");
            Assert.AreEqual(formulaA, formulaB);
            Assert.IsTrue(Math.Abs(formulaA.MonoisotopicMass - formulaB.MonoisotopicMass) < 1e-9);
        }

        [Test]
        public static void TestEquality()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CO");
            Assert.AreEqual(formulaB, formulaB);
        }

        [Test]
        public static void TestToChemicalFormula()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CO");
            Assert.AreEqual(ChemicalFormula.ParseFormula("CO"), formulaB);
        }

        [Test]
        public static void IsSuperSetOf()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CH3NO{17}C");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CHNO{16}");

            Assert.IsFalse(formulaA.IsSupersetOf(formulaB));
        }

        [Test]
        public static void ChemicalForulaIsSubSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H3NO");

            Assert.IsTrue(formulaA.IsSubsetOf(formulaB));
        }

        [Test]
        public static void ChemicalFormulaIsNotSubSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H2NO");

            Assert.IsFalse(formulaA.IsSubsetOf(formulaB));
        }

        [Test]
        public static void ChemicalFormulaIsSuperSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H3NO");

            Assert.IsTrue(formulaB.IsSupersetOf(formulaA));
        }

        [Test]
        public static void ChemicalFormulaIsNotSuperSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO2");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H3NO");

            Assert.IsFalse(formulaB.IsSupersetOf(formulaA));
        }

        private class PhysicalObjectWithChemicalFormula : IHasChemicalFormula
        {
            public PhysicalObjectWithChemicalFormula(string v)
            {
                ThisChemicalFormula = ChemicalFormula.ParseFormula(v);
            }

            public double MonoisotopicMass
            {
                get
                {
                    return ThisChemicalFormula.MonoisotopicMass;
                }
            }

            public ChemicalFormula ThisChemicalFormula
            {
                get; private set;
            }
        }
    }
}
