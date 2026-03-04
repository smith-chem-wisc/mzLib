// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompt 3: LOWESS non-linear calibration model
// Placement: MassSpectrometry/Dia/Calibration/LowessRtModel.cs

using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Dia.Calibration
{
    /// <summary>
    /// LOWESS (Locally Weighted Scatterplot Smoothing) RT calibration model.
    ///
    /// For each anchor point, a local weighted linear regression is fit using a tricube kernel,
    /// producing a smooth non-parametric mapping from library RT to observed RT. After fitting,
    /// an optional isotonic regression pass enforces monotonicity (pool-adjacent-violators).
    ///
    /// Prediction at arbitrary library RTs uses linear interpolation between fitted anchor values.
    /// Local σ is computed from a rolling window of ±10 anchors around the query point.
    ///
    /// Because <see cref="RtCalibrationModel"/> is sealed, this class implements
    /// <see cref="IRtCalibrationModel"/> instead.
    /// </summary>
    public sealed class LowessRtModel : IRtCalibrationModel
    {
        // ── Stored fitted data (sorted by library RT) ──────────────────────

        /// <summary>Anchor library RTs, sorted ascending.</summary>
        public double[] FittedLibraryRts { get; }

        /// <summary>Fitted (smoothed, monotonic-enforced) observed RTs at each anchor.</summary>
        public double[] FittedObservedRts { get; }

        /// <summary>Residuals (original observed - fitted) at each anchor point.</summary>
        private readonly double[] _residuals;

        /// <summary>Bandwidth parameter α used during fitting.</summary>
        public double Bandwidth { get; }

        /// <summary>Whether monotonicity was enforced via isotonic regression.</summary>
        public bool MonotonicEnforced { get; }

        // ── IRtCalibrationModel ────────────────────────────────────────────

        public double SigmaMinutes { get; }
        public double RSquared { get; }
        public int AnchorCount { get; }
        public bool IsReliable => RSquared >= RtCalibrationModel.MinReliableRSquared
                               && AnchorCount >= RtCalibrationModel.MinReliableAnchors;
        public RtCalibrationModelType ModelType => RtCalibrationModelType.Lowess;
        public double Slope { get; }
        public double Intercept { get; }

        // ── Construction ───────────────────────────────────────────────────

        private LowessRtModel(
            double[] fittedLibraryRts,
            double[] fittedObservedRts,
            double[] residuals,
            double bandwidth,
            bool monotonicEnforced,
            double globalSigma,
            double globalRSquared,
            double globalSlope,
            double globalIntercept,
            int anchorCount)
        {
            FittedLibraryRts = fittedLibraryRts;
            FittedObservedRts = fittedObservedRts;
            _residuals = residuals;
            Bandwidth = bandwidth;
            MonotonicEnforced = monotonicEnforced;
            SigmaMinutes = globalSigma;
            RSquared = globalRSquared;
            Slope = globalSlope;
            Intercept = globalIntercept;
            AnchorCount = anchorCount;
        }

        // ── Prediction ─────────────────────────────────────────────────────

        public double ToMinutes(double libraryRtOrIrt)
        {
            int n = FittedLibraryRts.Length;

            // Before first anchor — extrapolate
            if (libraryRtOrIrt <= FittedLibraryRts[0])
            {
                if (n < 2) return FittedObservedRts[0];
                double slope = (FittedObservedRts[1] - FittedObservedRts[0])
                             / (FittedLibraryRts[1] - FittedLibraryRts[0]);
                return FittedObservedRts[0] + slope * (libraryRtOrIrt - FittedLibraryRts[0]);
            }

            // After last anchor — extrapolate
            if (libraryRtOrIrt >= FittedLibraryRts[n - 1])
            {
                if (n < 2) return FittedObservedRts[n - 1];
                double slope = (FittedObservedRts[n - 1] - FittedObservedRts[n - 2])
                             / (FittedLibraryRts[n - 1] - FittedLibraryRts[n - 2]);
                return FittedObservedRts[n - 1] + slope * (libraryRtOrIrt - FittedLibraryRts[n - 1]);
            }

            // Binary search for bracketing anchors
            int lo = 0, hi = n - 2;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (libraryRtOrIrt < FittedLibraryRts[mid + 1]) hi = mid;
                else lo = mid + 1;
            }

            double leftLib = FittedLibraryRts[lo];
            double rightLib = FittedLibraryRts[lo + 1];
            double leftObs = FittedObservedRts[lo];
            double rightObs = FittedObservedRts[lo + 1];

            double fraction = (libraryRtOrIrt - leftLib) / (rightLib - leftLib);
            return leftObs + fraction * (rightObs - leftObs);
        }

        /// <summary>
        /// Returns the local σ computed from a rolling window of ±10 anchors
        /// around the nearest anchor to the given library RT.
        /// </summary>
        public double GetLocalSigma(double libraryRtOrIrt)
        {
            int n = FittedLibraryRts.Length;
            if (n == 0) return SigmaMinutes;

            int centerIdx = BinarySearchClosest(FittedLibraryRts, libraryRtOrIrt);

            int windowRadius = 10;
            int lo = Math.Max(0, centerIdx - windowRadius);
            int hi = Math.Min(n - 1, centerIdx + windowRadius);
            int count = hi - lo + 1;
            if (count < 3) return SigmaMinutes;

            double sumSq = 0;
            for (int i = lo; i <= hi; i++)
                sumSq += _residuals[i] * _residuals[i];

            return Math.Sqrt(sumSq / count);
        }

        public RtCalibrationModel ToRtCalibrationModel()
        {
            double slope = Math.Abs(Slope) < 1e-12 ? 1.0 : Slope;
            return new RtCalibrationModel(
                slope: slope,
                intercept: Intercept,
                sigmaMinutes: SigmaMinutes,
                rSquared: RSquared,
                anchorCount: AnchorCount);
        }

        // ── Fitting ────────────────────────────────────────────────────────

        /// <summary>
        /// Fits a LOWESS model to anchor pairs. Returns null if fewer than 50 anchors.
        ///
        /// Algorithm for each anchor point x_i:
        /// 1. Compute distances from x_i to all other library RTs
        /// 2. Select the nearest (bandwidth × N) points
        /// 3. Apply tricube weights: w_j = (1 - (|x_j - x_i| / max_dist)³)³
        /// 4. Fit weighted OLS: y = a + b·x with weights w_j
        /// 5. Store fitted value ŷ_i
        ///
        /// If enforceMonotonic is true, runs pool-adjacent-violators on the fitted
        /// values to ensure ŷ is non-decreasing.
        /// </summary>
        public static LowessRtModel Fit(
            double[] libraryRts,
            double[] observedRts,
            double bandwidth = 0.3,
            bool enforceMonotonic = true)
        {
            int n = libraryRts.Length;
            if (n < 50) return null;

            // Sort by library RT
            var indices = Enumerable.Range(0, n).ToArray();
            var libSorted = new double[n];
            var obsSorted = new double[n];
            Array.Copy(libraryRts, libSorted, n);
            Array.Sort(libSorted, indices);
            for (int i = 0; i < n; i++)
                obsSorted[i] = observedRts[indices[i]];

            // Neighborhood size
            int neighborhoodSize = Math.Max(3, (int)Math.Ceiling(bandwidth * n));

            // LOWESS fitting
            var fitted = new double[n];
            var distances = new double[n];

            for (int i = 0; i < n; i++)
            {
                double xi = libSorted[i];

                for (int j = 0; j < n; j++)
                    distances[j] = Math.Abs(libSorted[j] - xi);

                // Find the k-th distance to define neighborhood
                var sortedDist = new double[n];
                Array.Copy(distances, sortedDist, n);
                Array.Sort(sortedDist);
                double maxDist = sortedDist[Math.Min(neighborhoodSize - 1, n - 1)];
                if (maxDist < 1e-10) maxDist = 1e-10;

                // Weighted OLS with tricube kernel
                double sumW = 0, sumWX = 0, sumWY = 0, sumWXX = 0, sumWXY = 0;

                for (int j = 0; j < n; j++)
                {
                    if (distances[j] > maxDist) continue;

                    double u = distances[j] / maxDist;
                    double u3 = u * u * u;
                    double w = 1.0 - u3;
                    w = w * w * w; // tricube

                    sumW += w;
                    sumWX += w * libSorted[j];
                    sumWY += w * obsSorted[j];
                    sumWXX += w * libSorted[j] * libSorted[j];
                    sumWXY += w * libSorted[j] * obsSorted[j];
                }

                double denom = sumW * sumWXX - sumWX * sumWX;
                if (Math.Abs(denom) < 1e-12)
                {
                    fitted[i] = sumW > 0 ? sumWY / sumW : obsSorted[i];
                }
                else
                {
                    double b = (sumW * sumWXY - sumWX * sumWY) / denom;
                    double a = (sumWY - b * sumWX) / sumW;
                    fitted[i] = a + b * xi;
                }
            }

            // Isotonic regression
            if (enforceMonotonic)
                fitted = PoolAdjacentViolators(fitted);

            // Compute residuals and global statistics
            var residuals = new double[n];
            double ssRes = 0, ssTot = 0, meanObs = 0;
            for (int i = 0; i < n; i++) meanObs += obsSorted[i];
            meanObs /= n;

            for (int i = 0; i < n; i++)
            {
                residuals[i] = obsSorted[i] - fitted[i];
                ssRes += residuals[i] * residuals[i];
                double dev = obsSorted[i] - meanObs;
                ssTot += dev * dev;
            }

            double globalSigma = Math.Sqrt(ssRes / n);
            double globalRSq = ssTot > 0 ? 1.0 - ssRes / ssTot : 0;

            // Global linear fit for backward-compat Slope/Intercept
            FitGlobalOLS(libSorted, obsSorted, n, out double gSlope, out double gIntercept);

            return new LowessRtModel(
                libSorted, fitted, residuals,
                bandwidth, enforceMonotonic,
                globalSigma, globalRSq, gSlope, gIntercept, n);
        }

        // ── Pool-Adjacent-Violators ────────────────────────────────────────

        private static double[] PoolAdjacentViolators(double[] values)
        {
            int n = values.Length;
            var result = new double[n];
            Array.Copy(values, result, n);

            var blockStarts = new List<int>();
            var blockEnds = new List<int>();
            var blockValues = new List<double>();

            blockStarts.Add(0);
            blockEnds.Add(0);
            blockValues.Add(result[0]);

            for (int i = 1; i < n; i++)
            {
                blockStarts.Add(i);
                blockEnds.Add(i);
                blockValues.Add(result[i]);

                while (blockValues.Count >= 2 &&
                       blockValues[blockValues.Count - 1] < blockValues[blockValues.Count - 2])
                {
                    int last = blockValues.Count - 1;
                    int prev = last - 1;

                    int countPrev = blockEnds[prev] - blockStarts[prev] + 1;
                    int countLast = blockEnds[last] - blockStarts[last] + 1;
                    double pooled = (blockValues[prev] * countPrev + blockValues[last] * countLast)
                                  / (countPrev + countLast);

                    blockEnds[prev] = blockEnds[last];
                    blockValues[prev] = pooled;

                    blockStarts.RemoveAt(last);
                    blockEnds.RemoveAt(last);
                    blockValues.RemoveAt(last);
                }
            }

            for (int b = 0; b < blockStarts.Count; b++)
                for (int i = blockStarts[b]; i <= blockEnds[b]; i++)
                    result[i] = blockValues[b];

            return result;
        }

        // ── Helpers ────────────────────────────────────────────────────────

        private static int BinarySearchClosest(double[] sorted, double target)
        {
            int lo = 0, hi = sorted.Length - 1;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (sorted[mid] < target) lo = mid + 1;
                else hi = mid;
            }
            if (lo > 0 && Math.Abs(sorted[lo - 1] - target) < Math.Abs(sorted[lo] - target))
                return lo - 1;
            return lo;
        }

        private static void FitGlobalOLS(double[] xs, double[] ys, int n, out double slope, out double intercept)
        {
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
            double d = n * sxx - sx * sx;
            if (Math.Abs(d) < 1e-12) { slope = 1.0; intercept = (sy - sx) / n; return; }
            slope = (n * sxy - sx * sy) / d;
            intercept = (sy - slope * sx) / n;
        }
    }
}