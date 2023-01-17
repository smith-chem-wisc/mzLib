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
    internal class TestStatistics
    {
        
        [Test]
        public static void TestNonZeroMedian()
        {
            double[] arr = new double[] { 1, 2, 3, 4, 5 };
            double[] arr2 = new double[] { 1, 2, 3, 4, 5, 6 };
            double[] arr3 = new double[] { 1, 2, 3, 4, 5, 0 };
            double[] arr4 = new double[] { 0, 0, 0, 1, 2, 3, 4, 5, 6, 0, 0, 0 };
            double[] arr5 = new double[] { 0, 0, 0, 0, 0, 0 };

            double med1 = BasicStatistics.CalculateNonZeroMedian(arr);
            double med2 = BasicStatistics.CalculateNonZeroMedian(arr2);
            double med3 = BasicStatistics.CalculateNonZeroMedian(arr3);
            double med4 = BasicStatistics.CalculateNonZeroMedian(arr4);
            double med5 = BasicStatistics.CalculateNonZeroMedian(arr5);

            Assert.That(Math.Abs(med1 - 3) < 0.001);
            Assert.That(Math.Abs(med2 - 3.5) < 0.001);
            Assert.That(Math.Abs(med3 - 3) < 0.001);
            Assert.That(Math.Abs(med4 - 3.5) < 0.001);
            Assert.That(Math.Abs(med5 - 0) < 0.001);
        }

        [Test]
        public static void TestStandardDeviation()
        {
            double[] arr = new double[] { 1, 2, 3, 4, 5 };
            double[] arr2 = new double[] { 1, 2, 3, 4, 5, 6 };
            double[] arr3 = new double[] { 4, 1, 3, 5, 2 };
            double[] arr4 = new double[] { 5, 2, 6, 3, 4, 1 };

            double std1 = BasicStatistics.CalculateStandardDeviation(arr);
            double std1Avg = BasicStatistics.CalculateStandardDeviation(arr, arr.Average());
            Assert.That(Math.Abs(std1 - arr.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std1Avg - arr.StandardDeviation()) < 0.001);

            std1 = BasicStatistics.CalculateStandardDeviation((IEnumerable<double>)arr);
            std1Avg = BasicStatistics.CalculateStandardDeviation((IEnumerable<double>)arr, arr.Average());
            Assert.That(Math.Abs(std1 - arr.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std1Avg - arr.StandardDeviation()) < 0.001);
            
            double std2 = BasicStatistics.CalculateStandardDeviation(arr2);
            double std2Avg = BasicStatistics.CalculateStandardDeviation(arr2, arr2.Average());
            Assert.That(Math.Abs(std2 - arr2.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std2Avg - arr2.StandardDeviation()) < 0.001);

            double std3 = BasicStatistics.CalculateStandardDeviation(arr3);
            double std3Avg = BasicStatistics.CalculateStandardDeviation(arr3, arr3.Average());
            Assert.That(Math.Abs(std3 - arr3.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std3Avg - arr3.StandardDeviation()) < 0.001);
            
            double std4 = BasicStatistics.CalculateStandardDeviation(arr4);
            double std4Avg = BasicStatistics.CalculateStandardDeviation(arr4, arr4.Average());
            Assert.That(Math.Abs(std4 - arr4.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std4Avg - arr4.StandardDeviation()) < 0.001);
        }

        [Test]
        public static void TestNonZeroStandardDeviation()
        {
            double[] arr = new double[] { 1, 2, 3, 4, 5 };
            double[] arr2 = new double[] { 1, 2, 3, 4, 5, 6 };
            double[] arr3 = new double[] { 4, 1, 3, 5, 2 };
            double[] arr4 = new double[] { 5, 2, 6, 3, 4, 1 };
            double[] arr5 = new double[] { 1, 2, 3, 4, 5, 0 };
            double[] arr6 = new double[] { 0, 0, 0, 1, 2, 3, 4, 5, 6, 0, 0, 0 };
            double[] arr7 = new double[] { 0, 0, 0, 0, 0, 0 };

            double std1 = BasicStatistics.CalculateNonZeroStandardDeviation(arr);
            double std1Avg = BasicStatistics.CalculateNonZeroStandardDeviation(arr, arr.Average());
            Assert.That(Math.Abs(std1 - arr.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std1Avg - arr.StandardDeviation()) < 0.001);

            double std2 = BasicStatistics.CalculateNonZeroStandardDeviation(arr2);
            double std2Avg = BasicStatistics.CalculateNonZeroStandardDeviation(arr2, arr2.Average());
            Assert.That(Math.Abs(std2 - arr2.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std2Avg - arr2.StandardDeviation()) < 0.001);

            double std3 = BasicStatistics.CalculateNonZeroStandardDeviation(arr3);
            double std3Avg = BasicStatistics.CalculateNonZeroStandardDeviation(arr3, arr3.Average());
            Assert.That(Math.Abs(std3 - arr3.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std3Avg - arr3.StandardDeviation()) < 0.001);

            double std4 = BasicStatistics.CalculateNonZeroStandardDeviation(arr4);
            double std4Avg = BasicStatistics.CalculateNonZeroStandardDeviation(arr4, arr4.Average());
            Assert.That(Math.Abs(std4 - arr4.StandardDeviation()) < 0.001);
            Assert.That(Math.Abs(std4Avg - arr4.StandardDeviation()) < 0.001);

            double std5 = BasicStatistics.CalculateNonZeroStandardDeviation(arr5);
            Assert.That(Math.Abs(std5 - arr.StandardDeviation()) < 0.001);

            double std6 = BasicStatistics.CalculateNonZeroStandardDeviation(arr6);
            Assert.That(Math.Abs(std6 - arr2.StandardDeviation()) < 0.001);

            double std7 = BasicStatistics.CalculateNonZeroStandardDeviation(arr7);
            Assert.That(Math.Abs(std7 - 0) < 0.001);
        }
    }
}
