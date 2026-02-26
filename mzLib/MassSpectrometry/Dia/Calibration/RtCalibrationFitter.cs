// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Fits a linear calibration model (iRT → observed RT minutes) from anchor matches.
    /// 
    /// Supports:
    ///   - Ordinary Least Squares (OLS) for well-behaved data
    ///   - RANSAC + OLS refit for robustness against outliers (default)
    ///   - Iterative outlier rejection (post-fit 3σ pruning)
    /// 
    /// All methods are static and allocation-minimal (work on arrays, not LINQ).
    /// Thread-safe: no mutable state.
    /// </summary>
    public static class RtCalibrationFitter
    {
        /// <summary>
        /// Configuration for the calibration fitting process.
        /// </summary>
        public sealed class FitOptions
        {
            /// <summary>Whether to use RANSAC for initial outlier robustness. Default: true.</summary>
            public bool UseRansac { get; set; } = true;

            /// <summary>Number of RANSAC iterations. Default: 200.</summary>
            public int RansacIterations { get; set; } = 200;

            /// <summary>RANSAC inlier threshold in minutes. Default: 2.0.</summary>
            public double RansacInlierThresholdMinutes { get; set; } = 2.0;

            /// <summary>Minimum fraction of points that must be inliers for a valid RANSAC model. Default: 0.5.</summary>
            public double RansacMinInlierFraction { get; set; } = 0.5;

            /// <summary>Number of post-RANSAC outlier rejection passes (3σ pruning). Default: 2.</summary>
            public int OutlierRejectionPasses { get; set; } = 2;

            /// <summary>Outlier rejection threshold in standard deviations. Default: 3.0.</summary>
            public double OutlierSigmaThreshold { get; set; } = 3.0;

            /// <summary>Minimum anchor count required to attempt fitting. Default: 5.</summary>
            public int MinAnchors { get; set; } = 5;

            /// <summary>Random seed for RANSAC reproducibility. Default: 42.</summary>
            public int RandomSeed { get; set; } = 42;

            /// <summary>Default options suitable for most DIA runs.</summary>
            public static FitOptions Default => new();
        }

        /// <summary>
        /// Fits a calibration model from anchor pairs (library_iRT, observed_RT_minutes).
        /// 
        /// Pipeline:
        ///   1. If UseRansac: find consensus inlier set via RANSAC
        ///   2. OLS fit on inliers (or all points if RANSAC disabled)
        ///   3. Iterative 3σ outlier rejection
        ///   4. Final OLS fit → RtCalibrationModel
        /// 
        /// Returns null if insufficient anchors or degenerate data.
        /// </summary>
        /// <param name="libraryIrts">Library iRT values (x-axis).</param>
        /// <param name="observedRtMinutes">Observed RT in minutes (y-axis). Parallel to libraryIrts.</param>
        /// <param name="options">Fitting configuration. Null uses defaults.</param>
        public static RtCalibrationModel Fit(
            ReadOnlySpan<double> libraryIrts,
            ReadOnlySpan<double> observedRtMinutes,
            FitOptions options = null)
        {
            options ??= FitOptions.Default;

            if (libraryIrts.Length != observedRtMinutes.Length)
                throw new ArgumentException("Input arrays must have the same length.");

            if (libraryIrts.Length < options.MinAnchors)
                return null;

            // Work with mutable arrays for outlier rejection
            int n = libraryIrts.Length;
            var xs = new double[n];
            var ys = new double[n];
            libraryIrts.CopyTo(xs);
            observedRtMinutes.CopyTo(ys);

            // Track which points are active (not rejected)
            var active = new bool[n];
            int activeCount = n;
            for (int i = 0; i < n; i++) active[i] = true;

            // Step 1: RANSAC to find inlier consensus
            if (options.UseRansac && n >= 2)
            {
                var inliers = RansacFindInliers(xs, ys, options);
                if (inliers != null)
                {
                    active = inliers;
                    activeCount = CountTrue(active);
                }
                // If RANSAC fails (no good model found), proceed with all points
            }

            if (activeCount < options.MinAnchors)
                return null;

            // Step 2: Iterative outlier rejection
            for (int pass = 0; pass < options.OutlierRejectionPasses; pass++)
            {
                var fit = FitOls(xs, ys, active, activeCount);
                if (fit == null) return null;

                // Compute residuals and reject outliers > threshold * σ
                double threshold = options.OutlierSigmaThreshold * fit.Value.Sigma;
                if (threshold < 1e-12) break; // σ ≈ 0, no outliers possible

                bool changed = false;
                for (int i = 0; i < n; i++)
                {
                    if (!active[i]) continue;
                    double predicted = fit.Value.Slope * xs[i] + fit.Value.Intercept;
                    double residual = Math.Abs(ys[i] - predicted);
                    if (residual > threshold)
                    {
                        active[i] = false;
                        activeCount--;
                        changed = true;
                    }
                }

                if (!changed || activeCount < options.MinAnchors)
                    break;
            }

            if (activeCount < options.MinAnchors)
                return null;

            // Step 3: Final fit on surviving points
            var finalFit = FitOls(xs, ys, active, activeCount);
            if (finalFit == null) return null;

            // Compute R²
            double rSquared = ComputeRSquared(xs, ys, active, activeCount, finalFit.Value);

            return new RtCalibrationModel(
                slope: finalFit.Value.Slope,
                intercept: finalFit.Value.Intercept,
                sigmaMinutes: finalFit.Value.Sigma,
                rSquared: rSquared,
                anchorCount: activeCount);
        }

        /// <summary>
        /// Convenience overload accepting List&lt;double&gt;.
        /// Uses IList&lt;T&gt; to avoid ambiguity with the ReadOnlySpan overload
        /// (double[] is implicitly convertible to both Span and IReadOnlyList).
        /// </summary>
        public static RtCalibrationModel Fit(
            IList<double> libraryIrts,
            IList<double> observedRtMinutes,
            FitOptions options = null)
        {
            if (libraryIrts.Count != observedRtMinutes.Count)
                throw new ArgumentException("Input lists must have the same length.");

            var xs = new double[libraryIrts.Count];
            var ys = new double[observedRtMinutes.Count];
            for (int i = 0; i < xs.Length; i++)
            {
                xs[i] = libraryIrts[i];
                ys[i] = observedRtMinutes[i];
            }
            return Fit((ReadOnlySpan<double>)xs, (ReadOnlySpan<double>)ys, options);
        }

        // ── Internal OLS fit ────────────────────────────────────────────────

        private readonly struct OlsResult
        {
            public readonly double Slope;
            public readonly double Intercept;
            public readonly double Sigma;

            public OlsResult(double slope, double intercept, double sigma)
            {
                Slope = slope;
                Intercept = intercept;
                Sigma = sigma;
            }
        }

        /// <summary>
        /// Ordinary Least Squares fit on active subset.
        /// Returns null if degenerate (zero variance in x, or too few points).
        /// </summary>
        private static OlsResult? FitOls(double[] xs, double[] ys, bool[] active, int activeCount)
        {
            if (activeCount < 2)
                return null;

            // Compute means
            double sumX = 0, sumY = 0;
            for (int i = 0; i < xs.Length; i++)
            {
                if (!active[i]) continue;
                sumX += xs[i];
                sumY += ys[i];
            }
            double meanX = sumX / activeCount;
            double meanY = sumY / activeCount;

            // Compute slope and intercept
            double num = 0, den = 0;
            for (int i = 0; i < xs.Length; i++)
            {
                if (!active[i]) continue;
                double dx = xs[i] - meanX;
                num += dx * (ys[i] - meanY);
                den += dx * dx;
            }

            if (Math.Abs(den) < 1e-15)
                return null; // All x values identical

            double slope = num / den;
            double intercept = meanY - slope * meanX;

            // Compute residual standard deviation
            double sumSqResidual = 0;
            for (int i = 0; i < xs.Length; i++)
            {
                if (!active[i]) continue;
                double predicted = slope * xs[i] + intercept;
                double residual = ys[i] - predicted;
                sumSqResidual += residual * residual;
            }

            // Use N-2 for unbiased estimate (2 parameters: slope + intercept)
            double sigma = activeCount > 2
                ? Math.Sqrt(sumSqResidual / (activeCount - 2))
                : Math.Sqrt(sumSqResidual / activeCount);

            return new OlsResult(slope, intercept, sigma);
        }

        // ── RANSAC ──────────────────────────────────────────────────────────

        /// <summary>
        /// RANSAC: repeatedly sample 2 points, fit line, count inliers.
        /// Returns the best inlier mask, or null if no valid model found.
        /// </summary>
        private static bool[] RansacFindInliers(double[] xs, double[] ys, FitOptions options)
        {
            int n = xs.Length;
            if (n < 2) return null;

            var rng = new Random(options.RandomSeed);
            int bestInlierCount = 0;
            bool[] bestInliers = null;
            int minRequired = Math.Max(2, (int)(n * options.RansacMinInlierFraction));
            double threshold = options.RansacInlierThresholdMinutes;

            var candidateInliers = new bool[n];

            for (int iter = 0; iter < options.RansacIterations; iter++)
            {
                // Sample 2 distinct random points
                int i1 = rng.Next(n);
                int i2;
                do { i2 = rng.Next(n); } while (i2 == i1);

                // Fit line through the two points
                double dx = xs[i2] - xs[i1];
                if (Math.Abs(dx) < 1e-15) continue; // Degenerate: same x

                double slope = (ys[i2] - ys[i1]) / dx;
                double intercept = ys[i1] - slope * xs[i1];

                // Count inliers
                int inlierCount = 0;
                for (int i = 0; i < n; i++)
                {
                    double predicted = slope * xs[i] + intercept;
                    double residual = Math.Abs(ys[i] - predicted);
                    candidateInliers[i] = residual <= threshold;
                    if (candidateInliers[i]) inlierCount++;
                }

                if (inlierCount > bestInlierCount)
                {
                    bestInlierCount = inlierCount;
                    bestInliers = new bool[n];
                    Array.Copy(candidateInliers, bestInliers, n);
                }
            }

            if (bestInlierCount < minRequired)
                return null;

            return bestInliers;
        }

        // ── R² computation ──────────────────────────────────────────────────

        private static double ComputeRSquared(double[] xs, double[] ys, bool[] active, int activeCount,
            OlsResult fit)
        {
            if (activeCount < 2) return 0.0;

            double sumY = 0;
            for (int i = 0; i < ys.Length; i++)
                if (active[i]) sumY += ys[i];
            double meanY = sumY / activeCount;

            double ssTot = 0, ssRes = 0;
            for (int i = 0; i < ys.Length; i++)
            {
                if (!active[i]) continue;
                double dy = ys[i] - meanY;
                ssTot += dy * dy;

                double predicted = fit.Slope * xs[i] + fit.Intercept;
                double residual = ys[i] - predicted;
                ssRes += residual * residual;
            }

            if (ssTot < 1e-15) return 1.0; // All y values identical → perfect fit
            return 1.0 - (ssRes / ssTot);
        }

        // ── Utility ─────────────────────────────────────────────────────────

        private static int CountTrue(bool[] arr)
        {
            int count = 0;
            for (int i = 0; i < arr.Length; i++)
                if (arr[i]) count++;
            return count;
        }
    }
}