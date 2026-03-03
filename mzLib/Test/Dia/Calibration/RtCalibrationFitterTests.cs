// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using MassSpectrometry.Dia;
using NUnit.Framework;

namespace Test.Dia.Calibration
{
    [TestFixture]
    public class RtCalibrationFitterTests
    {
        /// <summary>
        /// Helper to disambiguate overload resolution: double[] → ReadOnlySpan&lt;double&gt;.
        /// </summary>
        private static RtCalibrationModel Fit(double[] xs, double[] ys,
            RtCalibrationFitter.FitOptions options = null)
            => RtCalibrationFitter.Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, options);
        // ── Perfect linear data ─────────────────────────────────────────────

        [Test]
        public void Fit_PerfectLinearData_RecoversSlopeAndIntercept()
        {
            // y = 0.5x + 2.0, no noise
            int n = 50;
            var xs = new double[n];
            var ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xs[i] = i * 2.0; // iRT: 0, 2, 4, ..., 98
                ys[i] = 0.5 * xs[i] + 2.0;
            }

            var model = Fit(xs, ys);

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Slope, Is.EqualTo(0.5).Within(1e-8));
            Assert.That(model.Intercept, Is.EqualTo(2.0).Within(1e-8));
            Assert.That(model.SigmaMinutes, Is.LessThan(1e-8));
            Assert.That(model.RSquared, Is.EqualTo(1.0).Within(1e-8));
            Assert.That(model.AnchorCount, Is.EqualTo(n));
            Assert.That(model.IsReliable, Is.True);
        }

        [Test]
        public void Fit_PerfectLinearData_OlsOnly()
        {
            int n = 30;
            var xs = new double[n];
            var ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xs[i] = i * 3.0;
                ys[i] = 0.8 * xs[i] + 1.5;
            }

            var options = new RtCalibrationFitter.FitOptions { UseRansac = false };
            var model = Fit(xs, ys, options);

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Slope, Is.EqualTo(0.8).Within(1e-8));
            Assert.That(model.Intercept, Is.EqualTo(1.5).Within(1e-8));
        }

        // ── Noisy data ──────────────────────────────────────────────────────

        [Test]
        public void Fit_NoisyLinearData_RecoversApproximateFit()
        {
            // y = 0.6x + 3.0 + noise
            var rng = new Random(42);
            int n = 100;
            var xs = new double[n];
            var ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xs[i] = i * 1.0;
                ys[i] = 0.6 * xs[i] + 3.0 + (rng.NextDouble() - 0.5) * 0.4;
            }

            var model = Fit(xs, ys);

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Slope, Is.EqualTo(0.6).Within(0.05));
            Assert.That(model.Intercept, Is.EqualTo(3.0).Within(0.5));
            Assert.That(model.RSquared, Is.GreaterThan(0.99));
            Assert.That(model.IsReliable, Is.True);
        }

        // ── RANSAC outlier robustness ───────────────────────────────────────

        [Test]
        public void Fit_WithOutliers_RansacRejectsOutliers()
        {
            // True model: y = 0.5x + 2.0
            // Add 20% gross outliers
            var rng = new Random(123);
            int n = 100;
            int outlierCount = 20;
            var xs = new double[n];
            var ys = new double[n];

            for (int i = 0; i < n; i++)
            {
                xs[i] = i * 1.0;
                ys[i] = 0.5 * xs[i] + 2.0 + (rng.NextDouble() - 0.5) * 0.2;
            }

            // Inject gross outliers
            for (int i = 0; i < outlierCount; i++)
            {
                int idx = rng.Next(n);
                ys[idx] += (rng.NextDouble() > 0.5 ? 1 : -1) * 20.0;
            }

            var model = Fit(xs, ys);

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Slope, Is.EqualTo(0.5).Within(0.05));
            Assert.That(model.Intercept, Is.EqualTo(2.0).Within(1.0));
            Assert.That(model.AnchorCount, Is.LessThan(n)); // Some points rejected
            Assert.That(model.RSquared, Is.GreaterThan(0.95));
        }

        [Test]
        public void Fit_WithOutliers_OlsOnlyIsWorse()
        {
            // Same data as above but without RANSAC — OLS result should be worse
            var rng = new Random(123);
            int n = 100;
            var xs = new double[n];
            var ys = new double[n];

            for (int i = 0; i < n; i++)
            {
                xs[i] = i * 1.0;
                ys[i] = 0.5 * xs[i] + 2.0 + (rng.NextDouble() - 0.5) * 0.2;
            }

            for (int i = 0; i < 20; i++)
            {
                int idx = rng.Next(n);
                ys[idx] += (rng.NextDouble() > 0.5 ? 1 : -1) * 20.0;
            }

            var ransacModel = Fit(xs, ys,
                new RtCalibrationFitter.FitOptions { UseRansac = true });
            var olsModel = Fit(xs, ys,
                new RtCalibrationFitter.FitOptions { UseRansac = false, OutlierRejectionPasses = 0 });

            Assert.That(ransacModel, Is.Not.Null);
            Assert.That(olsModel, Is.Not.Null);

            // RANSAC slope should be closer to true 0.5
            double ransacError = Math.Abs(ransacModel.Slope - 0.5);
            double olsError = Math.Abs(olsModel.Slope - 0.5);
            Assert.That(ransacError, Is.LessThan(olsError));
        }

        // ── Edge cases ──────────────────────────────────────────────────────

        [Test]
        public void Fit_TooFewPoints_ReturnsNull()
        {
            var xs = new double[] { 1.0, 2.0 };
            var ys = new double[] { 3.0, 4.0 };

            // Default MinAnchors = 5
            var model = Fit(xs, ys);
            Assert.That(model, Is.Null);
        }

        [Test]
        public void Fit_EmptyArrays_ReturnsNull()
        {
            var model = Fit(
                Array.Empty<double>(), Array.Empty<double>());
            Assert.That(model, Is.Null);
        }

        [Test]
        public void Fit_MismatchedLengths_Throws()
        {
            Assert.Throws<ArgumentException>(() =>
                Fit(new double[] { 1, 2 }, new double[] { 1 }));
        }

        [Test]
        public void Fit_AllSameX_ReturnsNull()
        {
            // Degenerate: all x-values identical
            int n = 20;
            var xs = new double[n];
            var ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xs[i] = 42.0;
                ys[i] = i * 1.0;
            }

            var model = Fit(xs, ys);
            Assert.That(model, Is.Null);
        }

        [Test]
        public void Fit_MinimumPointsExactly_Succeeds()
        {
            // Exactly 5 points (default MinAnchors)
            var xs = new double[] { 0, 25, 50, 75, 100 };
            var ys = new double[] { 2, 14.5, 27, 39.5, 52 }; // y = 0.5x + 2

            var model = Fit(xs, ys);
            Assert.That(model, Is.Not.Null);
            Assert.That(model.Slope, Is.EqualTo(0.5).Within(1e-6));
        }

        // ── Convenience overload with Lists ─────────────────────────────────

        [Test]
        public void Fit_ListOverload_SameResultAsSpan()
        {
            int n = 30;
            var xsArray = new double[n];
            var ysArray = new double[n];
            var xsList = new List<double>(n);
            var ysList = new List<double>(n);

            for (int i = 0; i < n; i++)
            {
                double x = i * 2.0;
                double y = 0.7 * x + 1.0;
                xsArray[i] = x;
                ysArray[i] = y;
                xsList.Add(x);
                ysList.Add(y);
            }

            var options = new RtCalibrationFitter.FitOptions { UseRansac = false };
            var modelSpan = RtCalibrationFitter.Fit(xsArray.AsSpan(), ysArray.AsSpan(), options);
            var modelList = RtCalibrationFitter.Fit(xsList, ysList, options);

            Assert.That(modelSpan.Slope, Is.EqualTo(modelList.Slope).Within(1e-10));
            Assert.That(modelSpan.Intercept, Is.EqualTo(modelList.Intercept).Within(1e-10));
        }

        // ── FitOptions configuration ────────────────────────────────────────

        [Test]
        public void FitOptions_CustomMinAnchors_Respected()
        {
            var xs = new double[] { 0, 10, 20 };
            var ys = new double[] { 0, 5, 10 };

            // Default MinAnchors=5 → null
            Assert.That(Fit(xs, ys), Is.Null);

            // Custom MinAnchors=3 → succeeds
            var options = new RtCalibrationFitter.FitOptions { MinAnchors = 3, UseRansac = false };
            var model = Fit(xs, ys, options);
            Assert.That(model, Is.Not.Null);
        }

        [Test]
        public void Fit_Deterministic_SameResultsWithSameSeed()
        {
            var rng = new Random(999);
            int n = 80;
            var xs = new double[n];
            var ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xs[i] = i * 1.0;
                ys[i] = 0.5 * xs[i] + 2.0 + (rng.NextDouble() - 0.5) * 0.5;
            }
            // Add a few outliers
            ys[10] += 15.0;
            ys[50] -= 12.0;

            var options = new RtCalibrationFitter.FitOptions { RandomSeed = 42 };
            var model1 = Fit(xs, ys, options);
            var model2 = Fit(xs, ys, options);

            Assert.That(model1.Slope, Is.EqualTo(model2.Slope).Within(1e-10));
            Assert.That(model1.Intercept, Is.EqualTo(model2.Intercept).Within(1e-10));
            Assert.That(model1.AnchorCount, Is.EqualTo(model2.AnchorCount));
        }

        // ── Realistic proteomics-scale data ─────────────────────────────────

        [Test]
        public void Fit_RealisticProteomicsData_GoodCalibration()
        {
            // Simulate: iRT range [0, 120], run 5-65 min, ~200 anchors, small noise
            // True model: RT = 0.5 * iRT + 5.0
            var rng = new Random(42);
            int n = 200;
            var xs = new double[n];
            var ys = new double[n];

            for (int i = 0; i < n; i++)
            {
                xs[i] = rng.NextDouble() * 120.0; // random iRT in [0, 120]
                ys[i] = 0.5 * xs[i] + 5.0 + rng.NextGaussian() * 0.3;
            }

            var model = Fit(xs, ys);

            Assert.That(model, Is.Not.Null);
            Assert.That(model.IsReliable, Is.True);
            Assert.That(model.Slope, Is.EqualTo(0.5).Within(0.02));
            Assert.That(model.Intercept, Is.EqualTo(5.0).Within(0.5));
            Assert.That(model.SigmaMinutes, Is.LessThan(0.5));
            Assert.That(model.RSquared, Is.GreaterThan(0.99));
        }

        [Test]
        public void Fit_NegativeSlope_Handled()
        {
            // Some LC setups might have reverse correlation
            int n = 50;
            var xs = new double[n];
            var ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xs[i] = i * 2.0;
                ys[i] = -0.3 * xs[i] + 60.0;
            }

            var model = Fit(xs, ys);
            Assert.That(model, Is.Not.Null);
            Assert.That(model.Slope, Is.EqualTo(-0.3).Within(1e-6));
        }
    }

    /// <summary>
    /// Extension for generating Gaussian random numbers in tests.
    /// </summary>
    internal static class RandomExtensions
    {
        public static double NextGaussian(this Random rng)
        {
            // Box-Muller transform
            double u1 = 1.0 - rng.NextDouble();
            double u2 = rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
    }
}