using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using NUnit.Framework.Internal;
using SpectralAveraging;

namespace Test.AveragingTests
{



    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestOutlierRejection
    {
        private SpectralAveragingParameters _parameters;
        
        private static readonly double[] MinMaxTest = { 10, 9, 8, 7, 6, 5 };
        private static readonly double[] MinMaxExpected = { 9, 8, 7, 6 };

        private static readonly double[] PercentileTest = new double[]
            { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
        private static readonly double[] PercentileExpected = new double[]
            { 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5 };

        private static readonly double[] SigmaTest = { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
        private static readonly double[] SigmaExpected = { 60, 50, 40, 30, 20, 10, 0 };

        private static readonly double[] WinsorizedTest = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
        private static readonly double[] WinsorizedExpected = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };

        private static readonly double[] TestAveragedSigma =
            { 120, 65, 64, 63, 62, 61, 60, 59, 59, 58, 57, 56, 30, 15 };
        private static readonly double[] ExpectedAveragedSigma = TestAveragedSigma[1..^1];

        private static readonly double[] TestThreshold = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
        private static readonly double[] ExpectedThreshold = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };

        private static readonly object[] Arguments =
        {
            new object[] { new TestCase(OutlierRejectionType.PercentileClipping, PercentileTest, PercentileExpected) },
            new object[]
                { new TestCase(OutlierRejectionType.WinsorizedSigmaClipping, WinsorizedTest, WinsorizedExpected) },
            new object[]
                { new TestCase(OutlierRejectionType.AveragedSigmaClipping, TestAveragedSigma, ExpectedAveragedSigma) },
            new object[] { new TestCase(OutlierRejectionType.MinMaxClipping, MinMaxTest, MinMaxExpected) },
            new object[] { new TestCase(OutlierRejectionType.SigmaClipping, SigmaTest, SigmaExpected) },
            new object[] { new TestCase(OutlierRejectionType.NoRejection, MinMaxTest, MinMaxTest) },
            new object[]
                { new TestCase(OutlierRejectionType.BelowThresholdRejection, TestThreshold, ExpectedThreshold) }
        };

        public class TestCase
        {
            public OutlierRejectionType OutlierRejectionType { get; set; }
            public double[] TestArray { get; set; }
            public double[] ExpectedArray { get; set; }

            public TestCase(OutlierRejectionType outlierRejection, double[] test, double[] expected)
            {
                OutlierRejectionType = outlierRejection;
                TestArray = test;
                ExpectedArray = expected;
            }
        }

        [OneTimeSetUp]
        public void OneTimeSetup()
        {
            _parameters = new();
            _parameters.SetDefaultValues();
            _parameters.MaxSigmaValue = 1.5;
            _parameters.MinSigmaValue = 1.5;
            _parameters.BinSize = 1.0d;
            _parameters.NormalizationType = NormalizationType.RelativeToTics;
            _parameters.Percentile = 0.9;
            _parameters.SpectralAveragingType = SpectralAveragingType.MzBinning;
        }

        [Test]
        [TestCaseSource(nameof(Arguments))]
        public void TestRejectionMethods(TestCase testCase)
        {
            // test pixel stack entry point
            var testArray = testCase.TestArray;
            var expectedArray = testCase.ExpectedArray;
            _parameters.OutlierRejectionType = testCase.OutlierRejectionType;

            double[] testXarray = Enumerable.Repeat(1.0, testArray.Length).ToArray();

            // test double array entry point
            var results = OutlierRejection.RejectOutliers(testArray, _parameters);
            Assert.That(results, Is.EqualTo(expectedArray));

        }

        [Test]
        public static void TestMinMaxClipping()
        {
            // test array entry point
            SpectralAveragingParameters parameters = new()
                { OutlierRejectionType = OutlierRejectionType.MinMaxClipping };
            double[] test = { 10, 9, 8, 7, 6, 5 };
            double[] expected = { 9, 8, 7, 6 };
            double[] minMaxClipped = OutlierRejection.RejectOutliers(test, parameters);
            Assert.That(minMaxClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));
        }

        [Test]
        public static void TestPercentileClipping()
        {
            // test array entry point
            SpectralAveragingParameters parameters = new()
                { OutlierRejectionType = OutlierRejectionType.PercentileClipping, Percentile = 0.9 };
            double[] test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
            double[] expected = new double[] { 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5 };
            double[] percentileClipped = OutlierRejection.RejectOutliers(test, parameters);
            Assert.That(percentileClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            // test array entry point
            test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
            expected = new double[] { 60, 50, 40, 30, 20, 10 };
            percentileClipped = OutlierRejection.RejectOutliers(test, parameters);
            Assert.That(percentileClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));
        }

        [Test]
        public static void TestSigmaClipping()
        {
            SpectralAveragingParameters parameters = new()
                { OutlierRejectionType = OutlierRejectionType.SigmaClipping, MaxSigmaValue = 1.5, MinSigmaValue = 1.5 };
            var test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
            double[] sigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            var expected = new double[] { 60, 50, 40, 30, 20, 10, 0 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1;
            parameters.MaxSigmaValue = 1;
            test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
            sigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 60, 50, 40 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1.3;
            parameters.MaxSigmaValue = 1.3;
            sigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 70, 60, 50, 40, 30 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1;
            parameters.MaxSigmaValue = 1.3;
            sigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 80, 70, 60 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1.5;
            parameters.MaxSigmaValue = 1.5;
            sigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));
        }

        [Test]
        public static void TestWinsorizedSigmaClipping()
        {
            // test array entry point
            SpectralAveragingParameters parameters = new()
            {
                OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping,
                MinSigmaValue = 1.5,
                MaxSigmaValue = 1.5
            };
            var test = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
            double[] windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            var expected = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
            Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            test = new double[] { 15, 30, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 100 };
            expected = new double[] { 56d, 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d, 65d };
            windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1.3;
            parameters.MaxSigmaValue = 1.3;
            windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d };
            Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MaxSigmaValue = 1.5;
            windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d, 65d };
            Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));
        }

        [Test]
        public static void TestAveragedSigmaClipping()
        {
            SpectralAveragingParameters parameters = new()
            {
                OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping,
                MaxSigmaValue = 3,
                MinSigmaValue = 3
            };
            var test = new double[] { 120, 65, 64, 63, 62, 61, 60, 59, 59, 58, 57, 56, 30, 15 };
            double[] averagedSigmaClipping = OutlierRejection.RejectOutliers(test, parameters);
            var expected = test[..];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1;
            parameters.MaxSigmaValue = 1;
            averagedSigmaClipping = OutlierRejection.RejectOutliers(test, parameters);
            expected = test[1..^2];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1.3;
            parameters.MaxSigmaValue = 1.3;
            averagedSigmaClipping = OutlierRejection.RejectOutliers(test, parameters);
            expected = test[1..^2];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));

            parameters.MinSigmaValue = 1.5;
            parameters.MaxSigmaValue = 1.5;
            averagedSigmaClipping = OutlierRejection.RejectOutliers(test, parameters);
            expected = test[1..^1];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));
            Assert.That(RejectAsBinnedPeak(test, parameters).SequenceEqual(expected));
        }

        [Test]
        public static void TestSetValuesAndRejectOutliersSwitch()
        {
            SpectralAveragingParameters parameters = new SpectralAveragingParameters();

            parameters.SetDefaultValues();
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.NoRejection);
            Assert.That(parameters.SpectralWeightingType == SpectraWeightingType.WeightEvenly);
            Assert.That(0.1, Is.EqualTo(parameters.Percentile));
            Assert.That(1.5, Is.EqualTo(parameters.MinSigmaValue));
            Assert.That(1.5, Is.EqualTo(parameters.MaxSigmaValue));

            parameters.SetValues(OutlierRejectionType.MinMaxClipping, SpectraWeightingType.WeightEvenly, SpectralAveragingType.MzBinning,
                normalizationType: NormalizationType.RelativeToTics, percentile: .8, minSigma: 2, maxSigma: 4);
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.MinMaxClipping);
            Assert.That(0.8, Is.EqualTo(parameters.Percentile));
            Assert.That(2, Is.EqualTo(parameters.MinSigmaValue));
            Assert.That(4, Is.EqualTo(parameters.MaxSigmaValue));

            // no rejection
            parameters.SetDefaultValues();
            double[] test = new double[] { 10, 8, 6, 5, 4, 2, 0 };
            double[] output = OutlierRejection.RejectOutliers(test, parameters);
            double[] expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
            Assert.That(output, Is.EqualTo(expected));

            // min max
            parameters.SetValues(outlierRejectionType: OutlierRejectionType.MinMaxClipping);
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.MinMaxClipping);
            output = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 8, 6, 5, 4, 2 };
            Assert.That(output, Is.EqualTo(expected));

            // percentile
            parameters.SetValues(outlierRejectionType: OutlierRejectionType.PercentileClipping);
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.PercentileClipping);
            output = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 6, 5, 4 };
            Assert.That(output, Is.EqualTo(expected));

            // sigma
            parameters.SetValues(outlierRejectionType: OutlierRejectionType.SigmaClipping);
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.SigmaClipping);
            output = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
            Assert.That(output, Is.EqualTo(expected));

            // winsorized sigma
            parameters.SetValues(outlierRejectionType: OutlierRejectionType.WinsorizedSigmaClipping);
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.WinsorizedSigmaClipping);
            output = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
            Assert.That(output, Is.EqualTo(expected));
            output = OutlierRejection.RejectOutliers(new double[] { }, parameters);
            Assert.That(output.SequenceEqual(new double[] { }));

            // averaged sigma
            parameters.SetValues(outlierRejectionType: OutlierRejectionType.AveragedSigmaClipping, minSigma: 1, maxSigma:1);
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.AveragedSigmaClipping);
            output = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 8, 6, 5, 4, 2 };
            Assert.That(output, Is.EqualTo(expected));

            // below threshold
            parameters.SetValues(outlierRejectionType: OutlierRejectionType.BelowThresholdRejection);
            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.BelowThresholdRejection);
            output = OutlierRejection.RejectOutliers(test, parameters);
            expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
            Assert.That(output, Is.EqualTo(expected));
        }

        [Test]
        public static void TestBelowThresholdRejection()
        {
            SpectralAveragingParameters parameters = new()
            {
                OutlierRejectionType = OutlierRejectionType.BelowThresholdRejection
            };
            var arr = new double[] { 0, 10, 0, 0, 0, 0, 0 };
            var arr2 = new double[] { 0, 10, 0, 0, 0 };
            var arr3 = new double[] { 0, 10, 0, 0 };
            var copy3 = new double[arr3.Length];
            Array.Copy(arr3, copy3, arr3.Length);

            // test array entry point
            var result = OutlierRejection.RejectOutliers(arr, parameters);
            Assert.That(result.All(p => p == 0));
            Assert.That(RejectAsBinnedPeak(arr, parameters).All(p => p == 0));

            result = OutlierRejection.RejectOutliers(arr2, parameters);
            Assert.That(result.All(p => p == 0));
            Assert.That(RejectAsBinnedPeak(arr2, parameters).All(p => p == 0));

            result = OutlierRejection.RejectOutliers(arr3, parameters);
            Assert.That(result.SequenceEqual(copy3));
            Assert.That(RejectAsBinnedPeak(arr3, parameters).SequenceEqual(copy3));
        }

        private static double[] RejectAsBinnedPeak(double[] testArray, SpectralAveragingParameters parameters)
        {
            List<BinnedPeak> peaks = new();
            for (int i = 0; i < testArray.Length; i++)
            {
                peaks.Add(new BinnedPeak(i, 1, testArray[i], 1));
            }

            return OutlierRejection.RejectOutliers(peaks, parameters).Select(p => p.Intensity).ToArray();
        }
    }
}