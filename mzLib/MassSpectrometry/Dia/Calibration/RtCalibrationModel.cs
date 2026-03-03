// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Linear calibration model mapping between library iRT values and observed RT in minutes.
    /// 
    /// Model: observed_RT_minutes = Slope * library_iRT + Intercept
    /// Inverse: calibrated_iRT = (observed_RT_minutes - Intercept) / Slope
    /// 
    /// This model is fit per-run from high-confidence anchor matches (Phase 1 broad search),
    /// then used to:
    ///   1. Convert scan RT to iRT space for library candidate selection
    ///   2. Compute RT windows as k * SigmaIrt (data-driven, not fixed)
    ///   3. Compute Gaussian RT scoring: rtScore = -(residual^2) / (2 * σ_iRT^2)
    /// 
    /// Immutable after construction. Thread-safe for concurrent reads.
    /// </summary>
    public sealed class RtCalibrationModel
    {
        /// <summary>Slope of the linear fit: minutes per iRT unit.</summary>
        public double Slope { get; }

        /// <summary>Intercept of the linear fit (minutes).</summary>
        public double Intercept { get; }

        /// <summary>Standard deviation of residuals in minutes.</summary>
        public double SigmaMinutes { get; }

        /// <summary>Standard deviation of residuals in iRT space: σ_minutes / |slope|.</summary>
        public double SigmaIrt { get; }

        /// <summary>R² goodness-of-fit statistic (0–1). Values below 0.9 suggest poor calibration.</summary>
        public double RSquared { get; }

        /// <summary>Number of anchor points used in the final fit.</summary>
        public int AnchorCount { get; }

        /// <summary>Whether this model is considered reliable for refined search.</summary>
        public bool IsReliable => RSquared >= MinReliableRSquared && AnchorCount >= MinReliableAnchors;

        /// <summary>Minimum R² threshold for reliable calibration.</summary>
        public const double MinReliableRSquared = 0.90;

        /// <summary>Minimum anchor count for reliable calibration.</summary>
        public const int MinReliableAnchors = 10;

        public RtCalibrationModel(
            double slope,
            double intercept,
            double sigmaMinutes,
            double rSquared,
            int anchorCount)
        {
            if (Math.Abs(slope) < 1e-12)
                throw new ArgumentException("Slope must be non-zero.", nameof(slope));

            Slope = slope;
            Intercept = intercept;
            SigmaMinutes = sigmaMinutes;
            SigmaIrt = sigmaMinutes / Math.Abs(slope);
            RSquared = rSquared;
            AnchorCount = anchorCount;
        }

        /// <summary>
        /// Converts observed RT (minutes) to calibrated iRT.
        /// </summary>
        public double ToIrt(double observedRtMinutes)
            => (observedRtMinutes - Intercept) / Slope;

        /// <summary>
        /// Converts library iRT to predicted observed RT (minutes).
        /// </summary>
        public double ToMinutes(double irt)
            => Slope * irt + Intercept;

        /// <summary>
        /// Computes the RT window half-width in iRT space: k * σ_iRT.
        /// </summary>
        /// <param name="k">Number of standard deviations (typically 3).</param>
        public double GetIrtWindowHalfWidth(double k = 3.0)
            => k * SigmaIrt;

        /// <summary>
        /// Computes the RT window half-width in minutes: k * σ_minutes.
        /// </summary>
        /// <param name="k">Number of standard deviations (typically 3).</param>
        public double GetMinutesWindowHalfWidth(double k = 3.0)
            => k * SigmaMinutes;

        /// <summary>
        /// Gaussian log-likelihood RT score for a given residual in iRT space.
        /// rtScore = -(residual_iRT^2) / (2 * σ_iRT^2)
        /// More negative = worse fit. Zero = perfect match.
        /// </summary>
        public double ComputeRtScore(double residualIrt)
        {
            if (SigmaIrt < 1e-12) return 0.0;
            return -(residualIrt * residualIrt) / (2.0 * SigmaIrt * SigmaIrt);
        }

        /// <summary>
        /// Computes the RT score given a library iRT and an observed scan iRT.
        /// Convenience method combining residual computation + scoring.
        /// </summary>
        public double ComputeRtScore(double libraryIrt, double calibratedScanIrt)
        {
            double residual = calibratedScanIrt - libraryIrt;
            return ComputeRtScore(residual);
        }

        /// <summary>
        /// Creates an uncalibrated (provisional) model that uses a simple linear scaling
        /// from iRT range to observed RT range. Used for Phase 1 broad search before
        /// any anchor matches are available.
        /// 
        /// provisional_iRT = observed_RT * (iRT_range / run_time_range)
        /// </summary>
        /// <param name="runRtMin">Minimum observed RT in the run (minutes).</param>
        /// <param name="runRtMax">Maximum observed RT in the run (minutes).</param>
        /// <param name="libraryIrtMin">Minimum iRT in the library.</param>
        /// <param name="libraryIrtMax">Maximum iRT in the library.</param>
        /// <param name="initialWindowIrt">Initial broad window half-width in iRT units (e.g., 20).</param>
        public static RtCalibrationModel CreateProvisional(
            double runRtMin, double runRtMax,
            double libraryIrtMin, double libraryIrtMax,
            double initialWindowIrt = 20.0)
        {
            double runRange = runRtMax - runRtMin;
            double irtRange = libraryIrtMax - libraryIrtMin;

            if (runRange < 1e-6 || irtRange < 1e-6)
            {
                // Degenerate case: identity mapping with large window
                return new RtCalibrationModel(
                    slope: 1.0,
                    intercept: 0.0,
                    sigmaMinutes: initialWindowIrt / 3.0,
                    rSquared: 0.0,
                    anchorCount: 0);
            }

            // observed_RT = slope * iRT + intercept
            // Map iRT range [irtMin, irtMax] → [rtMin, rtMax]
            double slope = runRange / irtRange;
            double intercept = runRtMin - slope * libraryIrtMin;

            // σ chosen so that k=3 window = initialWindowIrt in iRT space
            double sigmaMinutes = (initialWindowIrt / 3.0) * Math.Abs(slope);

            return new RtCalibrationModel(
                slope: slope,
                intercept: intercept,
                sigmaMinutes: sigmaMinutes,
                rSquared: 0.0,
                anchorCount: 0);
        }

        public override string ToString()
        {
            return $"RtCalibration: RT = {Slope:F6} * iRT + {Intercept:F4}, " +
                   $"σ={SigmaMinutes:F4} min ({SigmaIrt:F4} iRT), " +
                   $"R²={RSquared:F4}, anchors={AnchorCount}" +
                   (IsReliable ? "" : " [UNRELIABLE]");
        }
    }
}
