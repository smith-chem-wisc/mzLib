// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaTestFixture.cs) is part of Chemistry Library.
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

using Chemistry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class ChemicalFormulaTestFixture
    {
        #region Public Methods

        [Test]
        public void AddIsotopeWithExistingMassNumber()
        {
            Element al = PeriodicTable.GetElement("Al");
            Assert.Throws<MzLibException>(() =>
            {
                al.AddIsotope(27, 28, 1);
            }, "Isotope with mass number " + 28 + " already exists");
        }

        [Test]
        public void AddElementToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3N2O");

            Element n = PeriodicTable.GetElement(7);

            formulaA.Add(n, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void Multiply()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            formulaA.Multiply(2);
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C4H6N2O2");

            Assert.AreEqual(formulaA, formulaB);

            Assert.IsFalse(formulaA.Equals(null));
        }

        [Test]
        public void CheckToStringOfElements()
        {
            Element n = PeriodicTable.GetElement("N");
            Assert.AreEqual("" + n, "N");
        }

        [Test]
        public void AddFormulaToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H2O");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H5NO2");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void AddFormulaToItself()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C4H6N2O2");

            formulaA.Add(new ChemicalFormula(formulaA));

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddIChemicalFormulaToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            IHasChemicalFormula formulaB = new PhysicalObjectWithChemicalFormula("H2O");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H5NO2");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void AddIsotopeToFormula()
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
        public void AddLargeElementToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NOFe");

            Element fe = PeriodicTable.GetElement("Fe");

            formulaA.Add(fe, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddNegativeFormulaToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C-1H-2");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("CHNO");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void AddNegativeIsotopeToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2HH{1}2NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2HNO");

            Isotope h1 = PeriodicTable.GetElement("H")[1];

            formulaA.Add(h1, -2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddELementByAtomicNumber()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H2NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2HNO");

            formulaB.Add(1, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddNonExistentSymbolToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.Throws<KeyNotFoundException>(() => { formulaA.AddPrincipalIsotopesOf("Faa", 1); });
        }

        [Test]
        public void InexistingElement()
        {
            Assert.Throws<KeyNotFoundException>(() => { ChemicalFormula.ParseFormula("Q"); });
        }

        [Test]
        public void AddZeroElementToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Element n = PeriodicTable.GetElement("N");

            formulaA.Add(n, 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddZeroIsotopeToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Isotope h1 = PeriodicTable.GetElement("H")[1];

            formulaA.Add(h1, 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddZeroSymbolToFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.AddPrincipalIsotopesOf("H", 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ClearFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            formulaA.Clear();
            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public void ConstructorBlankStringEqualsEmptyFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula();

            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public void ConstructorDefaultEqualsEmptyFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula();

            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public void CopyConstructorValueEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula(formulaA);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void CopyConstructorReferenceInequality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula(formulaA);

            Assert.AreNotSame(formulaA, formulaB);
        }

        [Test]
        public void EmptyMonoisotopicMassIsZero()
        {
            Assert.AreEqual(0.0, new ChemicalFormula().MonoisotopicMass);
        }

        [Test]
        public void EmptyAverageMassIsZero()
        {
            Assert.AreEqual(0.0, new ChemicalFormula().AverageMass);
        }

        [Test]
        public void EmptyStringIsBlank()
        {
            Assert.IsEmpty(new ChemicalFormula().Formula);
        }

        [Test]
        public void EmptyAtomCountIsZero()
        {
            Assert.AreEqual(0, new ChemicalFormula().AtomCount);
        }

        [Test]
        public void EmptyElementCountIsZero()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public void EmptyIsotopeCountIsZero()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueIsotopes);
        }

        [Test]
        public void FormulaValueInequality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("NC1OH3");

            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public void FormulaValueInequalityHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public void FormulaValueEqualityItself()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public void FormulaValueEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("NC2OH3");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void FormulaEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public void FormulaAlmostEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C{12}2H3NO");
            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public void HashCodeEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3C2NO");

            Assert.AreEqual(formulaA.GetHashCode(), formulaB.GetHashCode());
        }

        [Test]
        public void HashCodeImmutableEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            Assert.AreEqual(formulaA.GetHashCode(), formulaA.GetHashCode());
        }

        [Test]
        public void HashCodeCheck()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("Al");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("Al{27}");
            Assert.AreNotEqual(formulaA.GetHashCode(), formulaB.GetHashCode());
        }

        [Test]
        public void HillNotation()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC2O");

            Assert.AreEqual("C2H3NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationNoCarbon()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("BrH");

            Assert.AreEqual("HBr", formulaA.Formula);
        }

        [Test]
        public void HillNotationNoCarbonNoHydrogen()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("Ca5O14Br6");

            Assert.AreEqual("Br6Ca5O14", formulaA.Formula);
        }

        [Test]
        public void HillNotationNoHydrogen()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("NC2O");

            Assert.AreEqual("C2NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC2C{13}2O");

            Assert.AreEqual("C2C{13}2H3NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationWithNegativeCount()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC-2O");

            Assert.AreEqual("C-2H3NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationWithHeavyIsotopeNegativeCount()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC2C{13}-2O");

            Assert.AreEqual("C2C{13}-2H3NO", formulaA.Formula);
        }

        [Test]
        public void BadFormula()
        {
            Assert.Throws<MzLibException>(() => { ChemicalFormula.ParseFormula("!@#$"); }, "Input string for chemical formula was in an incorrect format");
        }

        [Test]
        public void InvalidChemicalElement()
        {
            Assert.Throws<KeyNotFoundException>(() => { PeriodicTable.GetElement("Faa"); });
        }

        [Test]
        public void InvalidElementIsotope()
        {
            Assert.IsNull(PeriodicTable.GetElement("C")[100]);
        }

        [Test]
        public void NumberOfAtoms()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(7, formulaA.AtomCount);
        }

        [Test]
        public void NumberOfAtomsOfEmptyFormula()
        {
            Assert.AreEqual(0, new ChemicalFormula().AtomCount);
        }

        [Test]
        public void NumberOfAtomsOfNegativeFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C-2H-3N-1O-1");

            Assert.AreEqual(-7, formulaA.AtomCount);
        }

        [Test]
        public void ParsingFormulaNoNumbers()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CCHHHNO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithInternalSpaces()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2 H3 N O");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithSpacesAtEnd()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO  ");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithSpacesAtBeginning()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("    C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithSpaces()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("  C2 H3 N O  ");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaNoNumbersRandomOrder()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("OCHHCHN");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void EqualsFalse()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("OCHHCHN");
            Assert.IsFalse(formulaA.Equals("C"));
        }

        [Test]
        public void Equals()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("OCHHCHN");
            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public void ParsingFormulaRepeatedElements()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CH3NOC");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void IsSuperSetOf()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CH3NO{17}C");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CHNO{16}");

            Assert.IsFalse(formulaA.IsSupersetOf(formulaB));
        }

        [Test]
        public void ParsingFormulaRepeatedElementsCancelEachOther()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NOC-2");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveElementCompletelyFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2NO");

            formulaA.RemoveIsotopesOf(PeriodicTable.GetElement("H"));

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public void RemoveElementCompletelyFromFromulaBySymbol()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2NO");

            formulaA.RemoveIsotopesOf("H");

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public void RemoveElementCompletelyFromFromulaWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2C{13}H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H3NO");

            formulaA.RemoveIsotopesOf("C");

            Assert.AreEqual(formulaA, formulaB);
            Assert.AreEqual(formulaA.MonoisotopicMass, formulaB.MonoisotopicMass);
        }

        [Test]
        public void RemoveEmptyFormulaFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.Remove(new ChemicalFormula());

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveFormulaFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H5NOO{16}");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("H2O{16}");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.Remove(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void ContainsSpecificIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H5NOO{16}");

            Assert.IsTrue(formulaA.ContainsSpecificIsotope(PeriodicTable.GetElement("O")[16]));
        }

        [Test]
        public void ContainsIsotopesOf()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("O{16}");
            Assert.IsTrue(formulaA.ContainsIsotopesOf("O"));
            Assert.IsTrue(formulaA.ContainsSpecificIsotope("O", 16));
            Assert.AreEqual(1, formulaA.CountSpecificIsotopes("O", 16));
        }

        [Test]
        public void HydrogenCarbonRatio()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C8H4");
            Assert.AreEqual(0.5, formulaA.HydrogenCarbonRatio);
        }

        [Test]
        public void RemoveIsotopeCompletelyFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2NO");

            formulaA.RemoveIsotopesOf(PeriodicTable.GetElement("H"));

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveElementFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2HNO");

            formulaA.Remove(PeriodicTable.GetElement("H"), 2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveIsotopeFromFromulaEquality()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3O");

            formulaA.Remove("N", 1);

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public void RemoveNegativeElementFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H5NO");

            formulaA.Remove(PeriodicTable.GetElement("H"), -2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveNonExistantIsotopeFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H5NO2");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("Fe");
            ChemicalFormula formulaC = ChemicalFormula.ParseFormula("C2H5Fe-1NO2");

            formulaA.Remove(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void RemoveZeroIsotopeFromFromula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C2H3NO");

            formulaA.Remove(PeriodicTable.GetElement("H")[1], 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void TotalProtons()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{13}2H3NO");

            Assert.AreEqual(30, formulaA.ProtonCount);
        }

        [Test]
        public void TotalProtons2()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{12}2H3NO");

            Assert.AreEqual(30, formulaA.ProtonCount);
        }

        [Test]
        public void AverageMass()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C");

            Assert.AreEqual(PeriodicTable.GetElement("C").AverageMass, formulaA.AverageMass);
        }

        [Test]
        public void UniqueElements()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(4, formulaA.NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public void UniqueElementsOfEmptyFormula()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public void UniqueElementsWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.AreEqual(4, formulaA.NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public void UniqueIsotopes()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(0, formulaA.NumberOfUniqueIsotopes);
        }

        [Test]
        public void UniqueIsotopesOfEmptyFormula()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueIsotopes);
        }

        [Test]
        public void UniqueIsotopesWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.AreEqual(1, formulaA.NumberOfUniqueIsotopes);
        }

        [Test]
        public void ContainsIsotopesOfYe()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.IsTrue(formulaA.ContainsIsotopesOf(PeriodicTable.GetElement("C")));
        }

        [Test]
        public void TestReplaceIsotopes()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}2H3NO");

            formulaA.Replace(PeriodicTable.GetElement("C")[13], PeriodicTable.GetElement("C")[12]);
            Assert.AreEqual("CC{12}2H3NO", formulaA.Formula);
        }

        [Test]
        public void ChemicalForulaIsSubSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H3NO");

            Assert.IsTrue(formulaA.IsSubsetOf(formulaB));
        }

        [Test]
        public void ChemicalForulaIsNotSubSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H2NO");

            Assert.IsFalse(formulaA.IsSubsetOf(formulaB));
        }

        [Test]
        public void ChemicalForulaIsSuperSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H3NO");

            Assert.IsTrue(formulaB.IsSupersetOf(formulaA));
        }

        [Test]
        public void ChemicalForulaIsNotSuperSet()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO2");
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C3H3NO");

            Assert.IsFalse(formulaB.IsSupersetOf(formulaA));
        }

        [Test]
        public void ChemicalForulaMyTest()
        {
            ChemicalFormula formula = new ChemicalFormula();
            formula.Add(ChemicalFormula.ParseFormula("C3H5NO"));
            Assert.AreEqual(PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass * 3 + PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 5 + PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass, formula.MonoisotopicMass);
        }

        [Test]
        public void TestIsotopicDistribution()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            var a = IsotopicDistribution.GetDistribution(formulaA);

            Assert.True(Math.Abs(formulaA.MonoisotopicMass - a.Masses.ToArray()[Array.IndexOf(a.Intensities.ToArray(), a.Intensities.Max())]) < 1e-9);
        }

        [Test]
        public void TestIsotopicDistribution2()
        {
            IsotopicDistribution.GetDistribution(ChemicalFormula.ParseFormula("AlO{16}"));
        }

        [Test]
        public void TestIsotopicDistribution3()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CO");

            // Distinguish between O and C isotope masses
            var a = IsotopicDistribution.GetDistribution(formulaA, 0.0001);
            Assert.AreEqual(6, a.Masses.Count());

            // Do not distinguish between O and C isotope masses
            IsotopicDistribution.GetDistribution(formulaA, 0.001);

            // Do not distinguish between O and C isotope masses
            var b = IsotopicDistribution.GetDistribution(formulaA);
            Assert.AreEqual(4, b.Masses.Count());

            IsotopicDistribution.GetDistribution(formulaA, 0.1);

            PhysicalObjectWithChemicalFormula formulaB = new PhysicalObjectWithChemicalFormula("CO");
            IsotopicDistribution.GetDistribution(formulaB.ThisChemicalFormula, 1);
        }

        [Test]
        public void CatchIsotopicDistributionStuff()
        {
            ChemicalFormula formula = (ChemicalFormula.ParseFormula("C500O50H250N50"));
            IsotopicDistribution.GetDistribution(formula, 0.001, 1e-1, 1e-15);
        }

        [Test]
        public void CatchProbStuff()
        {
            ChemicalFormula formula = (ChemicalFormula.ParseFormula("C50O50"));
            IsotopicDistribution.GetDistribution(formula, 0.001, 1e-50, 1e-15);
        }

        [Test]
        public void I0j1()
        {
            ChemicalFormula formula = (ChemicalFormula.ParseFormula("C50O50"));
            IsotopicDistribution.GetDistribution(formula, 0.01, 0.1);

            IsotopicDistribution.GetDistribution(formula, 0.01, 0.5);

            IsotopicDistribution.GetDistribution(formula, 0.01, 0.75);
        }

        [Test]
        public void ThresholdProbability()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CO");

            // Only the principal isotopes have joint probability of 0.5! So one result when calcuating isotopic distribution
            var a = IsotopicDistribution.GetDistribution(formulaA, 0.0001, 0.5);
            Assert.AreEqual(1, a.Masses.Count());
            Assert.IsTrue(Math.Abs((PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass - a.Masses.First())) < 1e-9);
        }

        [Test]
        public void TestAnotherFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H{1}CC{13}2H3NO{16}");
            Assert.AreEqual("CC{13}2H3H{1}NO{16}", formulaA.Formula);
        }

        [Test]
        public void NeutronCount()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{12}O{16}");
            Assert.AreEqual(14, formulaA.NeutronCount());
        }

        [Test]
        public void NeutronCountFail()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CO");
            Assert.Throws<MzLibException>(() => { formulaA.NeutronCount(); }, "Cannot know for sure what the number of neutrons is!");
        }

        [Test]
        public void CombineTest()
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
        public void ValidatePeriodicTable()
        {
            Assert.IsTrue(PeriodicTable.ValidateAverageMasses(1e-2));
            Assert.IsFalse(PeriodicTable.ValidateAverageMasses(1e-3));
            Assert.IsTrue(PeriodicTable.ValidateAbundances(1e-15));
            Assert.IsFalse(PeriodicTable.ValidateAbundances(0));
        }

        [Test]
        public void TestAddChemicalFormula()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C");
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{12}");
            formulaB.Add(formulaA);
            Assert.AreEqual("CC{12}", formulaB.Formula);
        }

        [Test]
        public void NotEqual()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("C15O15H15S15N15");
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("N15S15H15O15C15");
            Assert.AreEqual(formulaA, formulaB);
            Assert.IsTrue(Math.Abs(formulaA.MonoisotopicMass - formulaB.MonoisotopicMass) < 1e-9);
        }

        [Test]
        public void TestRemoveObjectFromChemicalFormula()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CO");
            var ok = new PhysicalObjectWithChemicalFormula("C");
            formulaB.Remove(ok);

            Assert.AreEqual("O", formulaB.Formula);
        }

        [Test]
        public void TestEquality()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CO");
            Assert.AreEqual(formulaB, formulaB);
        }

        [Test]
        public void TestToChemicalFormula()
        {
            ChemicalFormula formulaB = ChemicalFormula.ParseFormula("CO");
            Assert.AreEqual(ChemicalFormula.ParseFormula("CO"), formulaB);
        }

        [Test]
        public void IsoTest()
        {
            ChemicalFormula formula = ChemicalFormula.ParseFormula("C5H8NO");

            IsotopicDistribution d = IsotopicDistribution.GetDistribution(formula, Math.Pow(2, -14));

            Assert.AreEqual(324, d.Intensities.Count());

            d = IsotopicDistribution.GetDistribution(formula, Math.Pow(2, -1));

            Assert.AreEqual(17, d.Intensities.Count());
        }

        #endregion Public Methods

        #region Private Classes

        private class PhysicalObjectWithChemicalFormula : IHasChemicalFormula
        {
            #region Public Constructors

            public PhysicalObjectWithChemicalFormula(string v)
            {
                ThisChemicalFormula = ChemicalFormula.ParseFormula(v);
            }

            #endregion Public Constructors

            #region Public Properties

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

            #endregion Public Properties
        }

        #endregion Private Classes
    }
}