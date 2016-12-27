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
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class ChemicalFormulaTestFixture
    {
        private class PhysicalObjectWithChemicalFormula : IHasChemicalFormula
        {
            public PhysicalObjectWithChemicalFormula(string v)
            {
                ThisChemicalFormula = new ChemicalFormula(v);
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

        [Test]
        public void AddIsotopeWithExistingMassNumber()
        {
            Element al = PeriodicTable.GetElement("Al");
            Assert.Throws<ArgumentException>(() =>
            {
                al.AddIsotope(27, 28, 1);
            }, "Isotope with mass number " + 28 + " already exists");
        }

        [Test]
        public void AddElementToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3N2O");

            Element n = PeriodicTable.GetElement(7);

            formulaA.Add(n, 1);

            Assert.AreEqual(formulaA, formulaB);
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
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("H2O");
            ChemicalFormula formulaC = new ChemicalFormula("C2H5NO2");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void AddFormulaToItself()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C4H6N2O2");

            formulaA.Add(new ChemicalFormula(formulaA));

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddIChemicalFormulaToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            IHasChemicalFormula formulaB = new PhysicalObjectWithChemicalFormula("H2O");
            ChemicalFormula formulaC = new ChemicalFormula("C2H5NO2");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void AddIsotopeToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3H{1}NO");

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
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NOFe");

            Element fe = PeriodicTable.GetElement("Fe");

            formulaA.Add(fe, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddNegativeFormulaToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C-1H-2");
            ChemicalFormula formulaC = new ChemicalFormula("CHNO");

            formulaA.Add(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void AddNegativeIsotopeToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2HH{1}2NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2HNO");

            Isotope h1 = PeriodicTable.GetElement("H")[1];

            formulaA.Add(h1, -2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddELementByAtomicNumber()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H2NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2HNO");

            formulaB.Add(1, 1);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddNonExistentSymbolToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");

            Assert.Throws<KeyNotFoundException>(() => { formulaA.AddPrincipalIsotopesOf("Faa", 1); });
        }

        [Test]
        public void InexistingElement()
        {
            Assert.Throws<KeyNotFoundException>(() => { new ChemicalFormula("Q"); });
        }

        [Test]
        public void AddZeroElementToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Element n = PeriodicTable.GetElement("N");

            formulaA.Add(n, 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddZeroIsotopeToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Isotope h1 = PeriodicTable.GetElement("H")[1];

            formulaA.Add(h1, 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void AddZeroSymbolToFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            formulaA.AddPrincipalIsotopesOf("H", 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ClearFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            formulaA.Clear();
            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public void ConstructorBlankStringEqualsEmptyFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("");

            Assert.AreEqual(formulaA, new ChemicalFormula());
        }

        [Test]
        public void ConstructorEmptyStringEqualsEmptyFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula(string.Empty);

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
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula(formulaA);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void CopyConstructorReferenceInequality()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
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
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("NC1OH3");

            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public void FormulaValueInequalityHeavyIsotope()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("CC{13}H3NO");

            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public void FormulaValueEqualityItself()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public void FormulaValueEquality()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("NC2OH3");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void FormulaEquality()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public void FormulaAlmostEquality()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C{12}2H3NO");
            Assert.AreNotEqual(formulaA, formulaB);
        }

        [Test]
        public void HashCodeEquality()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("H3C2NO");

            Assert.AreEqual(formulaA.GetHashCode(), formulaB.GetHashCode());
        }

        [Test]
        public void HashCodeImmutableEquality()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            Assert.AreEqual(formulaA.GetHashCode(), formulaA.GetHashCode());
        }

        [Test]
        public void HashCodeCheck()
        {
            ChemicalFormula formulaA = new ChemicalFormula("Al");
            ChemicalFormula formulaB = new ChemicalFormula("Al{27}");
            Assert.AreNotEqual(formulaA.GetHashCode(), formulaB.GetHashCode());
        }

        [Test]
        public void HillNotation()
        {
            ChemicalFormula formulaA = new ChemicalFormula("H3NC2O");

            Assert.AreEqual("C2H3NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationNoCarbon()
        {
            ChemicalFormula formulaA = new ChemicalFormula("BrH");

            Assert.AreEqual("HBr", formulaA.Formula);
        }

        [Test]
        public void HillNotationNoCarbonNoHydrogen()
        {
            ChemicalFormula formulaA = new ChemicalFormula("Ca5O14Br6");

            Assert.AreEqual("Br6Ca5O14", formulaA.Formula);
        }

        [Test]
        public void HillNotationNoHydrogen()
        {
            ChemicalFormula formulaA = new ChemicalFormula("NC2O");

            Assert.AreEqual("C2NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationWithHeavyIsotope()
        {
            ChemicalFormula formulaA = new ChemicalFormula("H3NC2C{13}2O");

            Assert.AreEqual("C2C{13}2H3NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationWithNegativeCount()
        {
            ChemicalFormula formulaA = new ChemicalFormula("H3NC-2O");

            Assert.AreEqual("C-2H3NO", formulaA.Formula);
        }

        [Test]
        public void HillNotationWithHeavyIsotopeNegativeCount()
        {
            ChemicalFormula formulaA = new ChemicalFormula("H3NC2C{13}-2O");

            Assert.AreEqual("C2C{13}-2H3NO", formulaA.Formula);
        }

        [Test]
        public void ImplicitConstructor()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = "C2H3NO";

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void BadFormula()
        {
            Assert.Throws<FormatException>(() => { new ChemicalFormula("!@#$"); }, "Input string for chemical formula was in an incorrect format");
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
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");

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
            ChemicalFormula formulaA = new ChemicalFormula("C-2H-3N-1O-1");

            Assert.AreEqual(-7, formulaA.AtomCount);
        }

        [Test]
        public void ParsingFormulaNoNumbers()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CCHHHNO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithInternalSpaces()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2 H3 N O");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithSpacesAtEnd()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO  ");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithSpacesAtBeginning()
        {
            ChemicalFormula formulaA = new ChemicalFormula("    C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaWithSpaces()
        {
            ChemicalFormula formulaA = new ChemicalFormula("  C2 H3 N O  ");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void ParsingFormulaNoNumbersRandomOrder()
        {
            ChemicalFormula formulaA = new ChemicalFormula("OCHHCHN");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void EqualsFalse()
        {
            ChemicalFormula formulaA = new ChemicalFormula("OCHHCHN");
            Assert.IsFalse(formulaA.Equals("C"));
        }

        [Test]
        public void Equals()
        {
            ChemicalFormula formulaA = new ChemicalFormula("OCHHCHN");
            Assert.AreEqual(formulaA, formulaA);
        }

        [Test]
        public void ParsingFormulaRepeatedElements()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CH3NOC");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void IsSuperSetOf()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CH3NO{17}C");
            ChemicalFormula formulaB = new ChemicalFormula("CHNO{16}");

            Assert.IsFalse(formulaA.IsSupersetOf(formulaB));
        }

        [Test]
        public void ParsingFormulaRepeatedElementsCancelEachOther()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NOC-2");
            ChemicalFormula formulaB = new ChemicalFormula("H3NO");

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveElementCompletelyFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2NO");

            formulaA.RemoveIsotopesOf(PeriodicTable.GetElement("H"));

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public void RemoveElementCompletelyFromFromulaBySymbol()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2NO");

            formulaA.RemoveIsotopesOf("H");

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public void RemoveElementCompletelyFromFromulaWithHeavyIsotope()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2C{13}H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("H3NO");

            formulaA.RemoveIsotopesOf("C");

            Assert.AreEqual(formulaA, formulaB);
            Assert.AreEqual(formulaA.MonoisotopicMass, formulaB.MonoisotopicMass);
        }

        [Test]
        public void RemoveEmptyFormulaFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            formulaA.Remove(new ChemicalFormula());

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveFormulaFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H5NOO{16}");
            ChemicalFormula formulaB = new ChemicalFormula("H2O{16}");
            ChemicalFormula formulaC = new ChemicalFormula("C2H3NO");

            formulaA.Remove(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void ContainsSpecificIsotope()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H5NOO{16}");

            Assert.IsTrue(formulaA.ContainsSpecificIsotope(PeriodicTable.GetElement("O")[16]));
        }

        [Test]
        public void ContainsIsotopesOf()
        {
            ChemicalFormula formulaA = new ChemicalFormula("O{16}");
            Assert.IsTrue(formulaA.ContainsIsotopesOf("O"));
            Assert.IsTrue(formulaA.ContainsSpecificIsotope("O", 16));
            Assert.AreEqual(1, formulaA.CountSpecificIsotopes("O", 16));
        }

        [Test]
        public void HydrogenCarbonRatio()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C8H4");
            Assert.AreEqual(0.5, formulaA.HydrogenCarbonRatio);
        }

        [Test]
        public void RemoveIsotopeCompletelyFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2NO");

            formulaA.RemoveIsotopesOf(PeriodicTable.GetElement("H"));

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveElementFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2HNO");

            formulaA.Remove(PeriodicTable.GetElement("H"), 2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveIsotopeFromFromulaEquality()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3O");

            formulaA.Remove("N", 1);

            Assert.AreEqual(formulaB, formulaA);
        }

        [Test]
        public void RemoveNegativeElementFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H5NO");

            formulaA.Remove(PeriodicTable.GetElement("H"), -2);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void RemoveNonExistantIsotopeFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H5NO2");
            ChemicalFormula formulaB = new ChemicalFormula("Fe");
            ChemicalFormula formulaC = new ChemicalFormula("C2H5Fe-1NO2");

            formulaA.Remove(formulaB);

            Assert.AreEqual(formulaA, formulaC);
        }

        [Test]
        public void RemoveZeroIsotopeFromFromula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C2H3NO");

            formulaA.Remove(PeriodicTable.GetElement("H")[1], 0);

            Assert.AreEqual(formulaA, formulaB);
        }

        [Test]
        public void TotalProtons()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C{13}2H3NO");

            Assert.AreEqual(30, formulaA.ProtonCount);
        }

        [Test]
        public void TotalProtons2()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C{12}2H3NO");

            Assert.AreEqual(30, formulaA.ProtonCount);
        }

        [Test]
        public void AverageMass()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C");

            Assert.AreEqual(PeriodicTable.GetElement("C").AverageMass, formulaA.AverageMass);
        }

        [Test]
        public void UniqueElements()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");

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
            ChemicalFormula formulaA = new ChemicalFormula("CC{13}H3NO");

            Assert.AreEqual(4, formulaA.NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public void UniqueIsotopes()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");

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
            ChemicalFormula formulaA = new ChemicalFormula("CC{13}H3NO");

            Assert.AreEqual(1, formulaA.NumberOfUniqueIsotopes);
        }

        [Test]
        public void ContainsIsotopesOfYe()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CC{13}H3NO");

            Assert.IsTrue(formulaA.ContainsIsotopesOf(PeriodicTable.GetElement("C")));
        }

        [Test]
        public void TestReplaceIsotopes()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CC{13}2H3NO");

            formulaA.Replace(PeriodicTable.GetElement("C")[13], PeriodicTable.GetElement("C")[12]);
            Assert.AreEqual("CC{12}2H3NO", formulaA.Formula);
        }

        [Test]
        public void ChemicalForulaIsSubSet()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C3H3NO");

            Assert.IsTrue(formulaA.IsSubsetOf(formulaB));
        }

        [Test]
        public void ChemicalForulaIsNotSubSet()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C3H2NO");

            Assert.IsFalse(formulaA.IsSubsetOf(formulaB));
        }

        [Test]
        public void ChemicalForulaIsSuperSet()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");
            ChemicalFormula formulaB = new ChemicalFormula("C3H3NO");

            Assert.IsTrue(formulaB.IsSupersetOf(formulaA));
        }

        [Test]
        public void ChemicalForulaIsNotSuperSet()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO2");
            ChemicalFormula formulaB = new ChemicalFormula("C3H3NO");

            Assert.IsFalse(formulaB.IsSupersetOf(formulaA));
        }

        [Test]
        public void ChemicalForulaMyTest()
        {
            ChemicalFormula formula = new ChemicalFormula();
            formula.Add(new ChemicalFormula("C3H5NO"));
            Assert.AreEqual(PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass * 3 + PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 5 + PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass, formula.MonoisotopicMass);
        }

        [Test]
        public void TestIsotopicDistribution()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C2H3NO");

            var a = new IsotopicDistribution(formulaA);

            Assert.True(formulaA.MonoisotopicMass.MassEquals(a.Masses[a.Intensities.IndexOf(a.Intensities.Max())]));
        }

        [Test]
        public void TestIsotopicDistribution2()
        {
            new IsotopicDistribution("AlO{16}");
        }

        [Test]
        public void TestIsotopicDistribution3()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CO");

            // Distinguish between O and C isotope masses
            var a = new IsotopicDistribution(formulaA, 0.0001);
            Assert.AreEqual(6, a.Masses.Count());

            // Do not distinguish between O and C isotope masses
            new IsotopicDistribution(formulaA, 0.001);

            // Do not distinguish between O and C isotope masses
            var b = new IsotopicDistribution(formulaA);
            Assert.AreEqual(4, b.Masses.Count());

            new IsotopicDistribution(formulaA, 0.1);

            PhysicalObjectWithChemicalFormula formulaB = new PhysicalObjectWithChemicalFormula("CO");
            new IsotopicDistribution(formulaB.ThisChemicalFormula, 1);
        }

        [Test]
        public void CatchIsotopicDistributionStuff()
        {
            ChemicalFormula formula = (new ChemicalFormula("C500O50H250N50"));
            new IsotopicDistribution(formula, 0.001, 1e-1, 1e-15);
            Console.WriteLine("");
        }

        [Test]
        public void catchProbStuff()
        {
            ChemicalFormula formula = (new ChemicalFormula("C50O50"));
            new IsotopicDistribution(formula, 0.001, 1e-50, 1e-15);
        }

        [Test]
        public void i0j1()
        {
            ChemicalFormula formula = (new ChemicalFormula("C50O50"));
            new IsotopicDistribution(formula, 0.01, 0.1);
            //Console.WriteLine(String.Join(", ", masses));
            //Console.WriteLine(String.Join(", ", intensities));

            new IsotopicDistribution(formula, 0.01, 0.5);
            //Console.WriteLine(String.Join(", ", masses));
            //Console.WriteLine(String.Join(", ", intensities));

            new IsotopicDistribution(formula, 0.01, 0.75);
            //Console.WriteLine(String.Join(", ", masses));
            //Console.WriteLine(String.Join(", ", intensities));
        }

        [Test]
        public void ThresholdProbability()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CO");

            // Only the principal isotopes have joint probability of 0.5! So one result when calcuating isotopic distribution
            var a = new IsotopicDistribution(formulaA, 0.0001, 0.5);
            Assert.AreEqual(1, a.Masses.Count());
            Assert.IsTrue((PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass).MassEquals(a.Masses[0]));
        }

        [Test]
        public void TestAnotherFormula()
        {
            ChemicalFormula formulaA = new ChemicalFormula("H{1}CC{13}2H3NO{16}");
            Assert.AreEqual("CC{13}2H3H{1}NO{16}", formulaA.Formula);
        }

        [Test]
        public void NeutronCount()
        {
            ChemicalFormula formulaA = new ChemicalFormula("C{12}O{16}");
            Assert.AreEqual(14, formulaA.NeutronCount);
        }

        [Test]
        public void NeutronCountFail()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CO");
            Assert.Throws<NotSupportedException>(() => { var a = formulaA.NeutronCount; }, "Cannot know for sure what the number of neutrons is!");
        }

        [Test]
        public void CombineTest()
        {
            List<IHasChemicalFormula> theList = new List<IHasChemicalFormula>();
            theList.Add(new PhysicalObjectWithChemicalFormula("C2H3NO"));
            theList.Add(new PhysicalObjectWithChemicalFormula("CO"));
            var c = ChemicalFormula.Combine(theList);

            Assert.AreEqual("C3H3NO2", c.Formula);
        }

        [Test]
        public void ValidatePeriodicTable()
        {
            Assert.AreEqual("Validation passed", PeriodicTable.ValidateAverageMasses(1e-2).Message);
            Assert.AreEqual(ValidationResult.PassedAverageMassValidation, PeriodicTable.ValidateAverageMasses(1e-2).ThisValidationResult);
            Assert.AreEqual("Average mass of Li is 6.94003660291572 instead of 6.9675", PeriodicTable.ValidateAverageMasses(1e-3).Message);
            Assert.AreEqual(ValidationResult.FailedAverageMassValidation, PeriodicTable.ValidateAverageMasses(1e-3).ThisValidationResult);

            Assert.AreEqual("Validation passed", PeriodicTable.ValidateAbundances(1e-15).Message);
            Assert.AreEqual(ValidationResult.PassedAbundanceValidation, PeriodicTable.ValidateAbundances(1e-15).ThisValidationResult);
        }

        [Test]
        public void TestNullArguments()
        {
            ChemicalFormula formulaA = new ChemicalFormula("CO");
            IHasChemicalFormula ok = null;
            Assert.AreEqual("item", Assert.Throws<ArgumentNullException>(() => { formulaA.Add(ok); }).ParamName);
            ChemicalFormula ok2 = null;
            Assert.AreEqual("formula", Assert.Throws<ArgumentNullException>(() => { formulaA.Add(ok2); }).ParamName);
            Assert.AreEqual("other", Assert.Throws<ArgumentNullException>(() => { new ChemicalFormula(ok2); }).ParamName);
            Element ok3 = null;
            Assert.AreEqual("element", Assert.Throws<ArgumentNullException>(() => { formulaA.AddPrincipalIsotopesOf(ok3, 0); }).ParamName);
            Assert.AreEqual("item", Assert.Throws<ArgumentNullException>(() => { formulaA.Remove(ok); }).ParamName);
            Assert.AreEqual("formula", Assert.Throws<ArgumentNullException>(() => { formulaA.Remove(ok2); }).ParamName);
            Assert.AreEqual("formula", Assert.Throws<ArgumentNullException>(() => { formulaA.IsSubsetOf(ok2); }).ParamName);
            Assert.AreEqual("formula", Assert.Throws<ArgumentNullException>(() => { formulaA.IsSupersetOf(ok2); }).ParamName);
            Assert.AreEqual("element", Assert.Throws<ArgumentNullException>(() => { formulaA.CountWithIsotopes(ok3); }).ParamName);
            Assert.AreEqual("element", Assert.Throws<ArgumentNullException>(() => { formulaA.CountSpecificIsotopes(ok3, 0); }).ParamName);
            Assert.IsFalse(formulaA.Equals(ok2));
            IEnumerable<IHasChemicalFormula> ok4 = null;
            Assert.AreEqual("formulas", Assert.Throws<ArgumentNullException>(() => { ChemicalFormula.Combine(ok4); }).ParamName);
            Assert.AreEqual("element", Assert.Throws<ArgumentNullException>(() => { PeriodicTable.Add(ok3); }).ParamName);

            Assert.AreEqual("formula", Assert.Throws<ArgumentNullException>(() => { new IsotopicDistribution(ok2); }).ParamName);

            IHasMass ok5 = null;

            Assert.AreEqual("objectWithMass", Assert.Throws<ArgumentNullException>(() => { ok5.ToMZ(0); }).ParamName);

            new PhysicalObjectWithChemicalFormula("C");
        }

        [Test]
        public void TestAddChemicalFormula()
        {
            ChemicalFormula formulaB = new ChemicalFormula("C");
            ChemicalFormula formulaA = new ChemicalFormula("C{12}");
            formulaB.Add(formulaA);
            Assert.AreEqual("CC{12}", formulaB.Formula);
        }

        [Test]
        public void NotEqual()
        {
            ChemicalFormula formulaB = new ChemicalFormula("COHSN");
            ChemicalFormula formulaA = new ChemicalFormula("NSHOC");
            Assert.AreEqual(formulaA, formulaB);
            Assert.IsTrue(formulaA.MonoisotopicMass.MassEquals(formulaB.MonoisotopicMass));
            Assert.IsFalse(formulaA.MonoisotopicMass.MassEquals(formulaB.MonoisotopicMass, 0));
        }

        [Test]
        public void TestRemoveObjectFromChemicalFormula()
        {
            ChemicalFormula formulaB = new ChemicalFormula("CO");
            var ok = new PhysicalObjectWithChemicalFormula("C");
            formulaB.Remove(ok);

            Assert.AreEqual("O", formulaB.Formula);
        }

        [Test]
        public void TestEquality()
        {
            ChemicalFormula formulaB = new ChemicalFormula("CO");
            Assert.IsTrue(formulaB.Equals(formulaB));
        }

        [Test]
        public void TestToChemicalFormula()
        {
            ChemicalFormula formulaB = ChemicalFormula.ToChemicalFormula("CO");
            Assert.AreEqual(new ChemicalFormula("CO"), formulaB);
        }

        [Test]
        public void IsoTest()
        {
            ChemicalFormula formula = new ChemicalFormula("C5H8NO");

            IsotopicDistribution d = new IsotopicDistribution(formula, Math.Pow(2, -14));

            Assert.AreEqual(324, d.Intensities.Count);

            d = new IsotopicDistribution(formula, Math.Pow(2, -1));

            Assert.AreEqual(17, d.Intensities.Count);
        }
    }
}