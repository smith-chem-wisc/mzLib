using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibSpectralAveraging;
using MzLibUtil;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public static class TestWeighting
    {
        public static BinnedSpectra BinnedSpectra { get; set; }

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            SpectralAveragingOptions options = new();
            var xArray = new[] { 1.0, 2, 3, 4, 5 };
            var yArray1 = new [] { 10.0, 10, 10, 10, 10 };
            var yArray2 = new [] { 20.0, 20, 20, 20, 20 };
            var yArray3 = new [] { 30.0, 30, 30, 30, 30 };
            var spec1 = new MzSpectrum(xArray, yArray1, true);
            var spec2 = new MzSpectrum(xArray, yArray2, true);
            var spec3 = new MzSpectrum(xArray, yArray3, true);
            List<MzSpectrum> spectra = new() { spec1, spec2, spec3 };

            BinnedSpectra = new(spectra, options.BinSize);
        }

        [Test]
        public static void TestWeightingSwitchError()
        {
            var exception = Assert.Throws<MzLibException>(() =>
            {
                SpectralWeighting.CalculateSpectraWeights(BinnedSpectra, new SpectralAveragingOptions()
                    { SpectralWeightingType = (SpectraWeightingType)(-1) });
            });
            Assert.That(exception.Message == "Spectra Weighting Type Not Implemented");
        }

        [Test]
        public static void TestWeighEvenly()
        {
            SpectralWeighting.CalculateSpectraWeights(BinnedSpectra, new SpectralAveragingOptions() 
                {SpectralWeightingType = SpectraWeightingType.WeightEvenly});
            var expected = new[] { 1.0, 1, 1 };
            Assert.That(expected.SequenceEqual(BinnedSpectra.Weights.Select(p => p.Value)));
        }

        [Test]
        public static void TestWeightByTicValue()
        {
            SpectralWeighting.CalculateSpectraWeights(BinnedSpectra, new SpectralAveragingOptions()
                { SpectralWeightingType = SpectraWeightingType.TicValue });
            var expected = new[] { (1.0 / 3), (2.0 / 3), 1 };
            Assert.That(expected.SequenceEqual(BinnedSpectra.Weights.Select(p => p.Value)));
        }

        [Test]
        public static void TestWeightByMrsNoiseEstimation()
        {
            var xArrays = new[]
            {
                new double[] { 0, 1, 2, 3, 3.49, 4 },
                new double[] { 0, 1, 2, 3, 4 },
                new double[] { 0.1, 1.1, 2.1, 3.1, 4.1}
            };
            var yArrays = new[]
            {
                new double[] { 10, 11, 12, 12, 13, 14 },
                new double[] { 11, 12, 13, 14, 15 },
                new double[] { 20, 25, 30, 35, 40 }
            };
            var binSize = 1.0;

            BinnedSpectra bs = new(xArrays, yArrays, binSize);
            bs.PerformNormalization();
            SpectralWeighting.CalculateSpectraWeights(bs, new SpectralAveragingOptions() { SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation });
            double[] expectedWeights = new[]
            {
                //1539.23913, 
                //1636.63560,
                //755.37045
                2073.60000,
                3785.99353127,
                889.4099838
            };
            Assert.That(bs.Weights.Values.ToArray(),
                Is.EqualTo(expectedWeights).Within(0.01));
        }
    }
}
