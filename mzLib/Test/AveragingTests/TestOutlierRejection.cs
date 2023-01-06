using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public static class TestOutlierRejection
    {
        [Test]
		public static void TestMinMaxClipping()
        {
            double[] test = { 10, 9, 8, 7, 6, 5 };
			double[] expected = { 9, 8, 7, 6 };
			double[] minMaxClipped = OutlierRejection.MinMaxClipping(test);
			Assert.That(minMaxClipped, Is.EqualTo(expected));
        }

        [Test]
        public static void TestPercentileClipping()
        {
            double[] test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
            double[] expected = new double[] { 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5 };
            double[] percentileClipped = OutlierRejection.PercentileClipping(test, 0.9);
            Assert.That(percentileClipped, Is.EqualTo(expected));

            test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
            expected = new double[] { 60, 50, 40, 30, 20, 10 };
            percentileClipped = OutlierRejection.PercentileClipping(test, 0.9);
            Assert.That(percentileClipped, Is.EqualTo(expected));
        }

        [Test]
        public static void TestSigmaClipping()
        {
            var test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
            double[] sigmaClipped = OutlierRejection.SigmaClipping(test, 1.5, 1.5);
            var expected = new double[] { 60, 50, 40, 30, 20, 10, 0 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));

            test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
            sigmaClipped = OutlierRejection.SigmaClipping(test, 1, 1);
            expected = new double[] { 60, 50, 40 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));

            sigmaClipped = OutlierRejection.SigmaClipping(test, 1.3, 1.3);
            expected = new double[] { 70, 60, 50, 40, 30 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));

            sigmaClipped = OutlierRejection.SigmaClipping(test, 1, 1.3);
            expected = new double[] { 80, 70, 60 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));

            sigmaClipped = OutlierRejection.SigmaClipping(test, 1.5, 1.5);
            expected = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
            Assert.That(sigmaClipped, Is.EqualTo(expected));
		}

        [Test]
        public static void TestWinsorizedSigmaClipping()
        {
            var test = new double[] { 0, 10, 20, 30, 40 ,50, 60, 70, 80, 90, 100  };
            double[] windsorizedSigmaClipped = OutlierRejection.WinsorizedSigmaClipping(test, 1.5, 1.5);
            var expected = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
            Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));

            test = new double[] { 100, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 30, 15 };
            test = new double[] { 15, 30, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 100 };
            expected = new double[] { 56d, 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d, 65d };
            windsorizedSigmaClipped = OutlierRejection.WinsorizedSigmaClipping(test, 1.5, 1.5);
            Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));

            windsorizedSigmaClipped = OutlierRejection.WinsorizedSigmaClipping(test, 1.3, 1.3);
			expected = new double[] { 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d };
			Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));

            windsorizedSigmaClipped = OutlierRejection.WinsorizedSigmaClipping(test, 1.3, 1.5);
			expected = new double[] { 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d, 65d };
			Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));
        }

        [Test]
        public static void TestAveragedSigmaClipping()
        {
            var test = new double[] { 120, 65, 64, 63, 62, 61, 60, 59, 59, 58, 57, 56, 30, 15 };
            double[] averagedSigmaClipping = OutlierRejection.AveragedSigmaClipping(test, 3, 3);
            var expected = test[1..test.Length];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));

            averagedSigmaClipping = OutlierRejection.AveragedSigmaClipping(test, 1, 1);
            expected = test[1..^2];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));

            averagedSigmaClipping = OutlierRejection.AveragedSigmaClipping(test, 1.3, 1.3);
            expected = test[1..^2];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));

            averagedSigmaClipping = OutlierRejection.AveragedSigmaClipping(test, 1.5, 1.5);
            expected = test[1..^2];
            Assert.That(averagedSigmaClipping, Is.EqualTo(expected));
        }

        [Test]
        public static void TestBelowThresholdRejection()
        {
            var arr = new double[] { 0, 10, 0, 0, 0, 0, 0};
            var arr2 = new double[] { 0, 10, 0, 0, 0};
            var arr3 = new double[] { 0, 10, 0, 0};
            var copy3 = new double[arr3.Length];
            Array.Copy(arr3, copy3, arr3.Length);

            var result = OutlierRejection.BelowThresholdRejection(arr);
            Assert.That(result.All(p => p == 0));
            
            var result2 = OutlierRejection.BelowThresholdRejection(arr2);
            Assert.That(result2.All(p => p == 0));

            var result3 = OutlierRejection.BelowThresholdRejection(arr3);
            Assert.That(result3.SequenceEqual(copy3));
        }


		[Test]
		public static void TestSetValuesAndRejectOutliersSwitch()
		{
            SpectralAveragingOptions options = new SpectralAveragingOptions();

			options.SetDefaultValues();
			Assert.That(options.RejectionType == RejectionType.NoRejection);
			Assert.That(options.WeightingType == WeightingType.NoWeight);
			Assert.That(0.1, Is.EqualTo(options.Percentile));
			Assert.That(1.5, Is.EqualTo(options.MinSigmaValue));
			Assert.That(1.5, Is.EqualTo(options.MaxSigmaValue));

			options.SetValues(RejectionType.MinMaxClipping, WeightingType.NoWeight, SpectrumMergingType.MzBinning, true, .8, 2, 4);
			Assert.That(options.RejectionType == RejectionType.MinMaxClipping);
			Assert.That(0.8, Is.EqualTo(options.Percentile));
			Assert.That(2, Is.EqualTo(options.MinSigmaValue));
			Assert.That(4, Is.EqualTo(options.MaxSigmaValue));

            // no rejection
			options.SetDefaultValues();
			double[] test = new double[] { 10, 8, 6, 5, 4, 2, 0 };
			double[] output = OutlierRejection.RejectOutliers(test, options);
			double[] expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
			Assert.That(output, Is.EqualTo(expected));

            // min max
			options.SetValues(rejectionType: RejectionType.MinMaxClipping);
			Assert.That(options.RejectionType == RejectionType.MinMaxClipping);
			output = OutlierRejection.RejectOutliers(test, options);
			expected = new double[] { 8, 6, 5, 4, 2 };
			Assert.That(output, Is.EqualTo(expected));

            // percentile
			options.SetValues(rejectionType: RejectionType.PercentileClipping);
			Assert.That(options.RejectionType == RejectionType.PercentileClipping);
            output = OutlierRejection.RejectOutliers(test, options);
            expected = new double[] { 6, 5, 4 };
            Assert.That(output, Is.EqualTo(expected));

            // sigma
            options.SetValues(rejectionType: RejectionType.SigmaClipping);
			Assert.That(options.RejectionType == RejectionType.SigmaClipping);
            output = OutlierRejection.RejectOutliers(test, options);
            expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
            Assert.That(output, Is.EqualTo(expected));

            // winsorized sigma
            options.SetValues(rejectionType: RejectionType.WinsorizedSigmaClipping);
			Assert.That(options.RejectionType == RejectionType.WinsorizedSigmaClipping);
            output = OutlierRejection.RejectOutliers(test, options);
            expected = new double[] { 10, 8, 6, 5, 4, 2 };
            Assert.That(output, Is.EqualTo(expected));
            output = OutlierRejection.RejectOutliers(new double[]{}, options);
            Assert.That(output.SequenceEqual(new double[] {}));

            // averaged sigma
            options.SetValues(rejectionType: RejectionType.AveragedSigmaClipping);
			Assert.That(options.RejectionType == RejectionType.AveragedSigmaClipping);
            output = OutlierRejection.RejectOutliers(test, options);
            expected = new double[] { 6, 5 };
            Assert.That(output, Is.EqualTo(expected));

            // below threshold
            options.SetValues(rejectionType: RejectionType.BelowThresholdRejection);
            Assert.That(options.RejectionType == RejectionType.BelowThresholdRejection);
            output = OutlierRejection.RejectOutliers(test, options);
            expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
            Assert.That(output, Is.EqualTo(expected));
        }
    }

}
