using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public static class TestWeighting
    {
        [Test]
        public static void TestWeightedAverage()
        {
            double[] values = new double[] { 10, 0 };
            double[] weights = new double[] { 8, 2 };
            double average = SpectralMerging.MergePeakValuesToAverage(values, weights);
            Assert.That(average, Is.EqualTo(8));

            values = new double[] { 10, 2, 0 };
            weights = new double[] { 9, 1, 0 };
            average = SpectralMerging.MergePeakValuesToAverage(values, weights);
            Assert.That(Math.Round(average, 4), Is.EqualTo(9.200));
        }

        [Test]
        public static void TestWeightByNormalDistribution()
        {
            double[] test = new double[] { 10, 8, 6, 5, 4, 3, 2, 1 };
            double[] weights = new double[test.Length];
            BinWeighting.WeightByNormalDistribution(test, ref weights);
            double weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(4.5749));
        }

        [Test]
        public static void TestWeightByCauchyDistribution()
        {
            double[] test = new double[] { 10, 8, 6, 5, 4, 3, 2, 1 };
            double[] weights = new double[test.Length];
            BinWeighting.WeightByCauchyDistribution(test, ref weights);
            double weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(4.6449));
        }

        [Test]
        public static void TestWeightByPoissonDistribution()
        {
            double[] test = new double[] { 10, 8, 6, 5, 4, 3, 2, 1 };
            double[] weights = new double[test.Length];
            BinWeighting.WeightByPoissonDistribution(test, ref weights);
            double weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(5.0244));
        }

        [Test]
        public static void TestWeightByGammaDistribution()
        {
            double[] test = new double[] { 10, 8, 6, 5, 4, 3, 2, 1 };
            double[] weights = new double[test.Length];
            BinWeighting.WeightByGammaDistribution(test, ref weights);
            double weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(4.6598));
        }

        [Test]
        public static void TestWeightingSwitch()
        {
            double[] test = new double[] { 10, 8, 6, 5, 4, 3, 2, 1 };
            double[] weights = new double[test.Length];
            weights = BinWeighting.CalculateWeights(test, WeightingType.NormalDistribution);
            double weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(4.5749));

            weights = BinWeighting.CalculateWeights(test, WeightingType.CauchyDistribution);
            weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(4.6449));

            weights = BinWeighting.CalculateWeights(test, WeightingType.PoissonDistribution);
            weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(5.0244));

            weights = BinWeighting.CalculateWeights(test, WeightingType.GammaDistribution);
            weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(4.6598));

            weights = BinWeighting.CalculateWeights(test, WeightingType.NoWeight);
            weightedAverage = SpectralMerging.MergePeakValuesToAverage(test, weights);
            Assert.That(Math.Round(weightedAverage, 4), Is.EqualTo(Math.Round(test.Average(), 4)));
        }
	}
}
