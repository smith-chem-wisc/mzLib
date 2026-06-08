// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaArithmeticTests.cs) is part of Chemistry Library.
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
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Chemistry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ChemistryTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ChemicalFormulaArithmeticTests
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
        public static void Multiply()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            formulaA.Multiply(2);
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C4H6N2O2");

            Assert.AreEqual(formulaA, formulaB);

            Assert.IsFalse(formulaA.Equals(null));
        }

        [Test]
        public static void AddFormulaToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H2O");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H5NO2");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public static void AddFormulasWithNegativeIsotopeValues()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H-1N-1O");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H-1N-1O");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("H-2N-2O2");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public static void AddFormulaToItself()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C4H6N2O2");

            formulaA.Add(new ChemicalFormula(formulaA));

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void AddIChemicalFormulaToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            IHasChemicalFormula formulaB = new PhysicalObjectWithChemicalFormula("H2O");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H5NO2");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public static void AddNegativeFormulaToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C-1H-2");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("CHNO");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public static void RemoveEmptyFormulaFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.Remove(new ChemicalFormula());

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void RemoveFormulaFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H5NOO{16}");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H2O{16}");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.Remove(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public static void RemoveNonExistantIsotopeFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H5NO2");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("Fe");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H5Fe-1NO2");

            formulaA.Remove(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public static void TestRemoveObjectFromChemicalFormula()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CO");
            var ok = new PhysicalObjectWithChemicalFormula("C");
            formulaB.Remove(ok);

            Assert.AreEqual("O", formulaB.Formula);
        }

        [Test]
        public static void CombineTest()
        {
            List<IHasChemicalFormula> theList = new List<IHasChemicalFormula>
            {
                new PhysicalObjectWithChemicalFormula("C2H3NO"),
                new PhysicalObjectWithChemicalFormula("CO")
            };
            var c = ChemicalFormula.Combine(theList);

            Assert.AreEqual("C3H3NO2", c.Formula);
        }

        [Test]
        public static void TestReplaceIsotopes()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}2H3NO");

            formulaA.Replace(PeriodicTable.GetElement("C")[13], PeriodicTable.GetElement("C")[12]);
            Assert.AreEqual("CC{12}2H3NO", formulaA.Formula);
        }

        [Test]
        public static void TestAddChemicalFormula()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C");
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{12}");
            formulaB.Add(formulaA);
            Assert.AreEqual("CC{12}", formulaB.Formula);
        }

        [Test]
        public static void TestAddChemicalFormulaOperator()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C");
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{12}");

            var addedFormula = formulaA + formulaB;
            formulaB.Add(formulaA);

            Assert.AreEqual("CC{12}", formulaB.Formula);
            Assert.AreEqual("CC{12}", addedFormula.Formula);

            var leftNull = null + formulaB;
            Assert.AreEqual(formulaB, leftNull);

            var rightNull = formulaB + null;
            Assert.AreEqual(formulaB, rightNull);

            ChemicalFormula nullFormula = null;
            var bothNull = nullFormula + nullFormula;
            Assert.AreEqual(null, bothNull);
        }

        [Test]
        public static void TestMultiplyChemicalFormulaOperator()
        {
            ChemicalFormula formula = ChemicalFormula.ParseFormula("C3");

            var multipliedFormula = formula * 3;
            Assert.AreEqual("C9", multipliedFormula.Formula);

            multipliedFormula = 3 * formula;
            Assert.AreEqual("C9", multipliedFormula.Formula);

            ChemicalFormula nullFormula = null;
            var leftNull = nullFormula * 3;
            Assert.AreEqual(null, leftNull);

            var rightNull = 3 * nullFormula;
            Assert.AreEqual(null, rightNull);
        }

        [Test]
        [TestCase("C", "N", "CN-1")]
        [TestCase(null, "N", "N-1")]
        [TestCase("C", null, "C")]
        public static void TestSubtractChemicalFormulaOperator(string formA, string formB, string expected)
        {
            ChemicalFormula formulaA = formA == null ? null : ChemicalFormula.ParseFormula(formA);
            ChemicalFormula formulaB = formB == null ? null : ChemicalFormula.ParseFormula(formB);

            var subtractedFormula = formulaA - formulaB;
            Assert.AreEqual(expected, subtractedFormula.Formula);
        }

        [Test]
        public static void ChemicalFormulaMyTest()
        {
            ChemicalFormula formula = new ChemicalFormula();
            formula.Add(ChemicalFormula.ParseFormula("C3H5NO"));
            Assert.AreEqual(PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass * 3 + PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 5 + PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass, formula.MonoisotopicMass);
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
