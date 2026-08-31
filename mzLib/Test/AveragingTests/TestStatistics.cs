using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using SpectralAveraging;

namespace Test.AveragingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestStatistics
    {
        [Test]
        public static void TestMedianAbsoluteDeviationFromMedian()
        {
            double[] arr1 = new double[] { 0, 0, 0, 1, 2, 3, 4, 5, 6, 0, 0, 0 };
            double expectedArr1 = 0.5;
            double results = BasicStatistics.MedianAbsoluteDeviationFromMedian(arr1);
            
            Assert.That(expectedArr1, Is.EqualTo(results).Within(0.01));
        }

        [Test]
        public static void TestBiweightMidvariance()
        {
            double[] arr1 = new double[] { 1, 2, 3, 4, 5, 6 };
            double expectedArr1 = 1.5;
            double results = BasicStatistics.MedianAbsoluteDeviationFromMedian(arr1);
            Assert.That(expectedArr1, Is.EqualTo(results).Within(0.01));
        }
    }
}
