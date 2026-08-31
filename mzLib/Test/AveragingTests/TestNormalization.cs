using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Printing;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;
using NUnit.Framework;
using SpectralAveraging;

namespace Test.AveragingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestNormalization
    {

        #region Set up test data
        public struct TestCase
        {
            internal NormalizationType normalizationType { get; init; }
            internal double[][] yArray { get; init; }
            internal double[][] expected { get; init; }

            internal TestCase(NormalizationType type, double[][] yArr)
            {
                normalizationType = type;
                yArray = yArr.SubArray(0, yArr.Length);

                expected = new double[yArr.Length][];
                switch (normalizationType)
                {
                    case NormalizationType.NoNormalization:
                        expected = yArr;
                        break;

                    case NormalizationType.RelativeToTics:
                        expected = yArr.Select(p => p.Select(m => m / p.Sum() * yArr.Average(n => n.Sum())).ToArray()).ToArray();
                        break;

                    case NormalizationType.AbsoluteToTic:
                        expected = yArr.Select(p => p.Select(m => m / p.Sum()).ToArray()).ToArray();
                        break;

                    case NormalizationType.RelativeIntensity:
                        expected = yArr.Select(p => p.Select(m => m / p.Max()).ToArray()).ToArray();
                        break;

                    default:
                        Debugger.Break();
                        break;
                }
            }
        }

        public static double[][] singleYArray =
        {
            new double[] { 100, 80, 70, 60, 50, 40 }
        };

        public static double[][] twoYArrays =
        {
            new double[] { 100, 80, 70, 60, 50, 40 },
            new double[] { 50, 40, 35, 30, 25, 20 },
        };

        public static double[][] threeYArrays =
        {
            new double[] { 100, 80, 70, 60, 50, 40 },
            new double[] { 100, 80, 70, 60, 50, 40 },
            new double[] { 50, 40, 35, 30, 25, 20 },
        };

        public static readonly object[] testCases =
        {
            // single y array
            new object[] { new TestCase(NormalizationType.NoNormalization, singleYArray) },
            new object[] { new TestCase(NormalizationType.RelativeToTics, singleYArray)},
            new object[] { new TestCase(NormalizationType.AbsoluteToTic, singleYArray) },
            new object[] { new TestCase(NormalizationType.RelativeIntensity, singleYArray) },

            //// two y arrays
            new object[] { new TestCase(NormalizationType.NoNormalization, twoYArrays) },
            new object[] { new TestCase(NormalizationType.RelativeToTics, twoYArrays) },
            new object[] { new TestCase(NormalizationType.AbsoluteToTic, twoYArrays) },
            new object[] { new TestCase(NormalizationType.RelativeIntensity, twoYArrays) },

            //// three y arrays
            new object[] { new TestCase(NormalizationType.NoNormalization, threeYArrays) },
            new object[] { new TestCase(NormalizationType.RelativeToTics, threeYArrays) },
            new object[] { new TestCase(NormalizationType.AbsoluteToTic, threeYArrays) },
            new object[] { new TestCase(NormalizationType.RelativeIntensity, threeYArrays) },
        };

        #endregion

        [Test]
        public static void TestNormalizationSwitchError()
        {
            var exception = Assert.Throws<MzLibException>(() =>
            {
                SpectraNormalization.NormalizeSpectra( new double[1][], (NormalizationType)(-1));
            });
            Assert.That(exception.Message == "Normalization Type not yet implemented");
        }

        [Test]
        [TestCaseSource(nameof(testCases))]
        public static void TestNormalizationFunctions(TestCase testCase)
        {
            var arraysToNormalize = testCase.yArray;

            SpectraNormalization.NormalizeSpectra(arraysToNormalize, testCase.normalizationType);
            for (var index = 0; index < arraysToNormalize.Length; index++)
            {
                var expected = testCase.expected[index].Select(p => Math.Round(p, 7));
                var result = arraysToNormalize[index].Select(p => Math.Round(p, 7));

                Assert.That(expected.SequenceEqual(result));
            }
        }

    }
}
