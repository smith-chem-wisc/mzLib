// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaConstructionTests.cs) is part of Chemistry Library.
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
using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ChemistryTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ChemicalFormulaConstructionTests
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
        public static void AddIsotopeWithExistingMassNumber()
        {
            Element al = PeriodicTable.GetElement("Al");
            Assert.Throws<MzLibException>(() =>
            {
                al.AddIsotope(27, 28, 1);
            }, "Isotope with mass number " + 28 + " already exists");
        }

        [Test]
        public static void AddElementToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3N2O");

            Element n = PeriodicTable.GetElement(7);

            formulaA.Add(n, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void AddIsotopeToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3H{1}NO");

            Isotope h1 = PeriodicTable.GetElement("H")[1];

            formulaA.Add(h1, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        /// <summary>
        /// This tests that the array for chemical formula properly expands
        /// </summary>
        [Test]
        public static void AddLargeElementToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NOFe");

            Element fe = PeriodicTable.GetElement("Fe");

            formulaA.Add(fe, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void AddNegativeIsotopeToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2HH{1}2NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2HNO");

            Isotope h1 = PeriodicTable.GetElement("H")[1];

            formulaA.Add(h1, -2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void AddELementByAtomicNumber()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H2NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2HNO");

            formulaB.Add(1, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void AddNonExistentSymbolToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.Throws<MzLibException>(() => { formulaA.AddPrincipalIsotopesOf("Faa", 1); });
        }

        [Test]
        public static void AddNonExistentAtomicNumberToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.Throws<MzLibException>(() => { formulaA.AddPrincipalIsotopesOf(1000000, 1); });
        }

        [Test]
        public static void AddZeroElementToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Element n = PeriodicTable.GetElement("N");

            formulaA.Add(n, 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void AddZeroIsotopeToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Isotope h1 = PeriodicTable.GetElement("H")[1];

            formulaA.Add(h1, 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void AddZeroSymbolToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.AddPrincipalIsotopesOf("H", 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ClearFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            formulaA.Clear();
            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public static void ConstructorBlankStringEqualsEmptyFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula();

            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public static void ConstructorDefaultEqualsEmptyFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula();

            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public static void CopyConstructorValueEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula(formulaA);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void CopyConstructorReferenceInequality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula(formulaA);

            Assert.AreNotSame(formulaA, formulaB);
        }

        [Test]
        public static void InexistingElement()
        {
            Assert.Throws<KeyNotFoundException>(() => { ChemicalFormula.ParseFormula("Q"); });
        }

        [Test]
        public static void BadFormula()
        {
            Assert.Throws<MzLibException>(() => { ChemicalFormula.ParseFormula("!@#$"); }, "Input string for chemical formula was in an incorrect format");
        }

        [Test]
        public static void InvalidChemicalElement()
        {
            Assert.Throws<KeyNotFoundException>(() => { PeriodicTable.GetElement("Faa"); });
        }

        [Test]
        public static void InvalidElementIsotope()
        {
            Assert.IsNull(PeriodicTable.GetElement("C")[100]);
        }

        [Test]
        public static void ParsingFormulaNoNumbers()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CCHHHNO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ParsingFormulaWithInternalSpaces()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2 H3 N O");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ParsingFormulaWithSpacesAtEnd()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO  ");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ParsingFormulaWithSpacesAtBeginning()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("    C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ParsingFormulaWithSpaces()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("  C2 H3 N O  ");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ParsingFormulaNoNumbersRandomOrder()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("OCHHCHN");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ParsingFormulaRepeatedElements()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CH3NOC");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ParsingFormulaRepeatedElementsCancelEachOther()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NOC-2");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void RemoveElementCompletelyFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2NO");

            formulaA.RemoveIsotopesOf(PeriodicTable.GetElement("H"));

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public static void RemoveElementCompletelyFromFromulaBySymbol()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2NO");

            formulaA.RemoveIsotopesOf("H");

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public static void RemoveElementCompletelyFromFromulaWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2C{13}H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3NO");

            formulaA.RemoveIsotopesOf("C");

            Assert.AreEqual(formulaA, formulaB);
            Assert.AreEqual(formulaA.MonoisotopicMass, formulaB.MonoisotopicMass);
        }

        [Test]
        public static void RemoveIsotopeCompletelyFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2NO");

            formulaA.RemoveIsotopesOf(PeriodicTable.GetElement("H"));

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void RemoveElementFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2HNO");

            formulaA.Remove(PeriodicTable.GetElement("H"), 2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void RemoveIsotopeFromFromulaEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3O");

            formulaA.Remove("N", 1);

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public static void RemoveNegativeElementFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H5NO");

            formulaA.Remove(PeriodicTable.GetElement("H"), -2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void RemoveZeroIsotopeFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.Remove(PeriodicTable.GetElement("H")[1], 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public static void ContainsSpecificIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H5NOO{16}");

            Assert.IsTrue(formulaA.ContainsSpecificIsotope(PeriodicTable.GetElement("O")[16]));
        }

        [Test]
        public static void ContainsIsotopesOf()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("O{16}");
            Assert.IsTrue(formulaA.ContainsIsotopesOf("O"));
            Assert.IsTrue(formulaA.ContainsSpecificIsotope("O", 16));
            Assert.AreEqual(1, formulaA.CountSpecificIsotopes("O", 16));
        }

        [Test]
        public static void ContainsIsotopesOfYe()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.IsTrue(formulaA.ContainsIsotopesOf(PeriodicTable.GetElement("C")));
        }

        [Test]
        public static void TestAddLargeNegative()
        {
            ChemicalFormula f = ChemicalFormula.ParseFormula("CO");
            f.Add("O", -10);
            Assert.That(f.Formula == "CO-9");
            Assert.That(f.NumberOfUniqueElementsByAtomicNumber == 2);
            Assert.That(f.MonoisotopicMass, Is.EqualTo(-131.95423157613).Within(0.001));
        }
    }
}
