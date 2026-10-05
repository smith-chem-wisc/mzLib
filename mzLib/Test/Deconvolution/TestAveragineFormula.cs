using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.Deconvolution
{
    /// <summary>
    /// Tests for <see cref="AveragineFormula"/>, which fills the part of a species' mass that has no known formula with
    /// averagine atoms and returns an isotopic envelope anchored at the species' exact monoisotopic mass.
    /// </summary>
    [TestFixture]
    public static class TestAveragineFormula
    {
        private const double TmtMass = 229.162932;

        private static readonly ChemicalFormula KnownPart = ChemicalFormula.ParseFormula("C50H80N14O16S");

        private static double AveragineMass(IReadOnlyDictionary<char, double> composition) =>
            composition.Sum(kvp => PeriodicTable.GetElement(kvp.Key.ToString()).AverageMass * kvp.Value);

        private static (double[] Masses, double[] Intensities) ShiftedDistribution(ChemicalFormula formula, double monoisotopicMass)
        {
            var distribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);
            double[] masses = distribution.Masses.ToArray();
            for (int i = 0; i < masses.Length; i++)
            {
                masses[i] += (monoisotopicMass - formula.MonoisotopicMass);
            }
            return (masses, distribution.Intensities.ToArray());
        }

        [Test]
        public static void AddAveragine_RoundsEachElementToWholeAtoms()
        {
            var composition = new Averagine().GetAverageChemicalFormula();
            var formula = new ChemicalFormula();

            // ten averagine residues: C49.384 H77.583 N13.577 O14.773 S0.417
            AveragineFormula.AddAveragine(formula, 10 * AveragineMass(composition), composition);

            Assert.AreEqual("C49H78N14O15", formula.Formula);
        }

        [Test]
        public static void AddAveragine_NegativeMassSubtractsAtoms()
        {
            var composition = new Averagine().GetAverageChemicalFormula();
            var formula = ChemicalFormula.ParseFormula("C100H200N30O40S2");

            AveragineFormula.AddAveragine(formula, -10 * AveragineMass(composition), composition);

            Assert.AreEqual("C51H122N16O25S2", formula.Formula);
        }

        [Test]
        public static void GetAnchoredDistribution_AnchorsMonoisotopicPeakAndLeavesFormulaUnchanged()
        {
            string before = KnownPart.Formula;
            double monoisotopicMass = KnownPart.MonoisotopicMass + TmtMass;

            var (masses, intensities) = AveragineFormula.GetAnchoredDistribution(KnownPart, monoisotopicMass,
                new Averagine().GetAverageChemicalFormula(), 0.125, 1e-8);

            Assert.AreEqual(monoisotopicMass, masses[0], 1e-9);
            Assert.AreEqual(masses.Length, intensities.Length);
            Assert.AreEqual(before, KnownPart.Formula);
        }

        [Test]
        public static void GetAnchoredDistribution_FillsAGapLargerThanTheThreshold()
        {
            var composition = new Averagine().GetAverageChemicalFormula();
            double monoisotopicMass = KnownPart.MonoisotopicMass + TmtMass;
            var filled = new ChemicalFormula(KnownPart);
            AveragineFormula.AddAveragine(filled, TmtMass, composition);
            var expected = ShiftedDistribution(filled, monoisotopicMass);

            var actual = AveragineFormula.GetAnchoredDistribution(KnownPart, monoisotopicMass, composition, 0.125, 1e-8);

            Assert.AreEqual(expected.Masses, actual.Masses);
            Assert.AreEqual(expected.Intensities, actual.Intensities);
            Assert.AreNotEqual(ShiftedDistribution(KnownPart, monoisotopicMass).Intensities, actual.Intensities);
        }

        [Test]
        [TestCase(15.9949, 20)]
        [TestCase(TmtMass, double.PositiveInfinity)]
        public static void GetAnchoredDistribution_LeavesAGapWithinTheThresholdUnfilled(double gap, double threshold)
        {
            double monoisotopicMass = KnownPart.MonoisotopicMass + gap;
            var expected = ShiftedDistribution(KnownPart, monoisotopicMass);

            var actual = AveragineFormula.GetAnchoredDistribution(KnownPart, monoisotopicMass,
                new Averagine().GetAverageChemicalFormula(), 0.125, 1e-8, threshold);

            Assert.AreEqual(expected.Masses, actual.Masses);
            Assert.AreEqual(expected.Intensities, actual.Intensities);
        }

        [Test]
        public static void GetAnchoredDistribution_FillsWithTheCompositionItIsGiven()
        {
            double monoisotopicMass = KnownPart.MonoisotopicMass + 1000;

            var amino = AveragineFormula.GetAnchoredDistribution(KnownPart, monoisotopicMass,
                new Averagine().GetAverageChemicalFormula(), 0.125, 1e-8);
            var ribo = AveragineFormula.GetAnchoredDistribution(KnownPart, monoisotopicMass,
                new OxyriboAveragine().GetAverageChemicalFormula(), 0.125, 1e-8);

            Assert.AreEqual(monoisotopicMass, ribo.Masses[0], 1e-9);
            Assert.AreNotEqual(amino.Intensities, ribo.Intensities);
        }

        [Test]
        public static void NullArgumentsThrow()
        {
            var composition = new Averagine().GetAverageChemicalFormula();

            Assert.Throws<ArgumentNullException>(() => AveragineFormula.AddAveragine(null, 100, composition));
            Assert.Throws<ArgumentNullException>(() => AveragineFormula.AddAveragine(new ChemicalFormula(), 100, null));
            Assert.Throws<ArgumentNullException>(() => AveragineFormula.GetAnchoredDistribution(null, 100, composition, 0.125, 1e-8));
            Assert.Throws<ArgumentNullException>(() => AveragineFormula.GetAnchoredDistribution(KnownPart, 100, null, 0.125, 1e-8));
        }
    }
}
