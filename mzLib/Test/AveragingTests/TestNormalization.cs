using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public static class TestNormalization
    {
        [Test]
        public static void TestBaseNormalizationFunctions()
        {
            double[] sampleData = new double[] { 100, 80, 70, 60, 50, 40 };
            double[] expected = new double[] { 0.25, 0.2, 0.175, 0.15, 0.125, 0.1 }; 
            var copy = new double[sampleData.Length];

            SpectrumNormalization.NormalizeSpectrumToTic(ref sampleData, 400);
            Assert.That(Math.Abs(sampleData.Sum() - 1) < 0.001);
            Assert.That(sampleData.SequenceEqual(expected));


            sampleData.CopyTo(copy, 0);
            SpectrumNormalization.NormalizeSpectrumToTic(ref sampleData, 0);
            Assert.That(sampleData.SequenceEqual(copy));

            sampleData = new double[] { 100, 80, 70, 60, 50, 40 };
            SpectrumNormalization.NormalizeSpectrumToTic(sampleData, 400);
            Assert.That(Math.Abs(sampleData.Sum() - 1) < 0.001);
            Assert.That(sampleData.SequenceEqual(expected));

            sampleData.CopyTo(copy, 0);
            SpectrumNormalization.NormalizeSpectrumToTic(sampleData, 0);
            Assert.That(sampleData.SequenceEqual(copy));
        }

        [Test]
        public static void TestAvgTicNormalizedSpectrum()
        {
            double[] sampleData = new double[] { 100, 80, 70, 60, 50, 40 };
            var copy = new double[sampleData.Length];

            double[] expected = new double[] { 25, 20, 17.5, 15, 12.5, 10 };
            SpectrumNormalization.NormalizeSpectrumToTic(sampleData, 400, 100);
            Assert.That(Math.Abs(sampleData.Sum() - 100) < 0.001);
            Assert.That(sampleData.SequenceEqual(expected));


            sampleData.CopyTo(copy, 0);
            SpectrumNormalization.NormalizeSpectrumToTic(sampleData, 0, 100);
            Assert.That(sampleData.SequenceEqual(copy));
        }

        
    }
}
