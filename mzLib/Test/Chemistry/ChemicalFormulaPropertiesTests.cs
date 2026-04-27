// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaPropertiesTests.cs) is part of Chemistry Library.
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
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ChemicalFormulaPropertiesTests
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
        public static void EmptyMonoisotopicMassIsZero()
        {
            Assert.AreEqual(0.0, new ChemicalFormula().MonoisotopicMass);
        }

        [Test]
        public static void EmptyAverageMassIsZero()
        {
            Assert.AreEqual(0.0, new ChemicalFormula().AverageMass);
        }

        [Test]
        public static void EmptyStringIsBlank()
        {
            Assert.IsEmpty(new ChemicalFormula().Formula);
        }

        [Test]
        public static void EmptyAtomCountIsZero()
        {
            Assert.AreEqual(0, new ChemicalFormula().AtomCount);
        }

        [Test]
        public static void EmptyElementCountIsZero()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public static void EmptyIsotopeCountIsZero()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueIsotopes);
        }

        [Test]
        public static void NumberOfAtoms()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(7, formulaA.AtomCount);
        }

        [Test]
        public static void NumberOfAtomsOfEmptyFormula()
        {
            Assert.AreEqual(0, new ChemicalFormula().AtomCount);
        }

        [Test]
        public static void NumberOfAtomsOfNegativeFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C-2H-3N-1O-1");

            Assert.AreEqual(-7, formulaA.AtomCount);
        }

        [Test]
        public static void TotalProtons()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{13}2H3NO");

            Assert.AreEqual(30, formulaA.ProtonCount);
        }

        [Test]
        public static void TotalProtons2()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{12}2H3NO");

            Assert.AreEqual(30, formulaA.ProtonCount);
        }

        [Test]
        public static void AverageMass()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C");

            Assert.AreEqual(PeriodicTable.GetElement("C").AverageMass, formulaA.AverageMass);
        }

        [Test]
        public static void UniqueElements()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(4, formulaA.NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public static void UniqueElementsOfEmptyFormula()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public static void UniqueElementsWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.AreEqual(4, formulaA.NumberOfUniqueElementsByAtomicNumber);
        }

        [Test]
        public static void UniqueIsotopes()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            Assert.AreEqual(0, formulaA.NumberOfUniqueIsotopes);
        }

        [Test]
        public static void UniqueIsotopesOfEmptyFormula()
        {
            Assert.AreEqual(0, new ChemicalFormula().NumberOfUniqueIsotopes);
        }

        [Test]
        public static void UniqueIsotopesWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CC{13}H3NO");

            Assert.AreEqual(1, formulaA.NumberOfUniqueIsotopes);
        }

        [Test]
        public static void HillNotation()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC2O");

            Assert.AreEqual("C2H3NO", formulaA.Formula);
        }

        [Test]
        public static void HillNotationNoCarbon()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("BrH");

            Assert.AreEqual("HBr", formulaA.Formula);
        }

        [Test]
        public static void HillNotationNoCarbonNoHydrogen()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("Ca5O14Br6");

            Assert.AreEqual("Br6Ca5O14", formulaA.Formula);
        }

        [Test]
        public static void HillNotationNoHydrogen()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("NC2O");

            Assert.AreEqual("C2NO", formulaA.Formula);
        }

        [Test]
        public static void HillNotationWithHeavyIsotope()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC2C{13}2O");

            Assert.AreEqual("C2C{13}2H3NO", formulaA.Formula);
        }

        [Test]
        public static void HillNotationWithNegativeCount()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC-2O");

            Assert.AreEqual("C-2H3NO", formulaA.Formula);
        }

        [Test]
        public static void HillNotationWithHeavyIsotopeNegativeCount()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H3NC2C{13}-2O");

            Assert.AreEqual("C2C{13}-2H3NO", formulaA.Formula);
        }

        [Test]
        public static void TestIsotopicDistribution()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C2H3NO");

            var a = IsotopicDistribution.GetDistribution(formulaA);

            Assert.True(Math.Abs(formulaA.MonoisotopicMass - a.Masses.ToArray()[Array.IndexOf(a.Intensities.ToArray(), a.Intensities.Max())]) < 1e-9);
        }

        [Test]
        public static void TestIsotopicDistribution2()
        {
            IsotopicDistribution.GetDistribution(ChemicalFormula.ParseFormula("AlO{16}"));
        }

        [Test]
        public static void TestIsotopicDistribution3()
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
        public static void CatchIsotopicDistributionStuff()
        {
            ChemicalFormula formula = (ChemicalFormula.ParseFormula("C500O50H250N50"));
            IsotopicDistribution.GetDistribution(formula, 0.001, 1e-1, 1e-15);
        }

        [Test]
        public static void CatchProbStuff()
        {
            ChemicalFormula formula = (ChemicalFormula.ParseFormula("C50O50"));
            IsotopicDistribution.GetDistribution(formula, 0.001, 1e-50, 1e-15);
        }

        [Test]
        public static void I0j1()
        {
            ChemicalFormula formula = (ChemicalFormula.ParseFormula("C50O50"));
            IsotopicDistribution.GetDistribution(formula, 0.01, 0.1);

            IsotopicDistribution.GetDistribution(formula, 0.01, 0.5);

            IsotopicDistribution.GetDistribution(formula, 0.01, 0.75);
        }

        [Test]
        public static void ThresholdProbability()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CO");

            // Only the principal isotopes have joint probability of 0.5! So one result when calcuating isotopic distribution
            var a = IsotopicDistribution.GetDistribution(formulaA, 0.0001, 0.5);
            Assert.AreEqual(1, a.Masses.Count());
            Assert.IsTrue(Math.Abs((PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass - a.Masses.First())) < 1e-9);
        }

        [Test]
        public static void TestAnotherFormula()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("H{1}CC{13}2H3NO{16}");
            Assert.AreEqual("CC{13}2H3H{1}NO{16}", formulaA.Formula);
        }

        [Test]
        public static void NeutronCount()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C{12}O{16}");
            Assert.AreEqual(14, formulaA.NeutronCount());
        }

        [Test]
        public static void NeutronCountFail()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("CO");
            Assert.Throws<MzLibException>(() => { formulaA.NeutronCount(); }, "Cannot know for sure what the number of neutrons is!");
        }

        [Test]
        public static void IsoTest()
        {
            ChemicalFormula formula = ChemicalFormula.ParseFormula("C5H8NO");

            IsotopicDistribution d = IsotopicDistribution.GetDistribution(formula, Math.Pow(2, -14));

            Assert.AreEqual(324, d.Intensities.Count());

            d = IsotopicDistribution.GetDistribution(formula, Math.Pow(2, -1));

            Assert.AreEqual(17, d.Intensities.Count());
        }

        [Test]
        public static void ValidatePeriodicTable()
        {
            Assert.IsTrue(PeriodicTable.ValidateAverageMasses(1e-2));
            Assert.IsFalse(PeriodicTable.ValidateAverageMasses(1e-3));
            Assert.IsTrue(PeriodicTable.ValidateAbundances(1e-15));
            Assert.IsFalse(PeriodicTable.ValidateAbundances(0));
        }

        [Test]
        public static void HydrogenCarbonRatio()
        {
            ChemicalFormula formulaA = ChemicalFormula.ParseFormula("C8H4");
            Assert.AreEqual(0.5, formulaA.HydrogenCarbonRatio);
        }

        [Test]
        public static void CheckToStringOfElements()
        {
            Element n = PeriodicTable.GetElement("N");
            Assert.AreEqual("" + n, "N");
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
