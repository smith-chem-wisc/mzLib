using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace StatisticalModels
{
    /// <summary>Whether <see cref="EffectScale.Estimate"/> could estimate the scale.</summary>
    public enum EffectScaleStatus
    {
        /// <summary>The scale and its standard error were estimated.</summary>
        Estimated,

        /// <summary>Fewer than <see cref="EffectScale.MinimumFeatures"/> features, or none with a non-zero reference
        /// effect: a scale and an overdispersion cannot both be estimated. Every estimate is NaN.</summary>
        TooFewFeatures,

        /// <summary>Newton's method on the estimating equations did not settle within <see cref="EffectScale.MaxIterations"/>
        /// steps from the alternating fit's values; those values are reported.</summary>
        NotConverged,
    }

    /// <summary>The scale between two sets of effects of the same features, and the estimators that frame it.</summary>
    public sealed class EffectScaleResult
    {
        internal EffectScaleResult() { }

        /// <summary>Features used.</summary>
        public int Features { get; internal set; }
        /// <summary>Whether the scale could be estimated.</summary>
        public EffectScaleStatus Status { get; internal set; }
        /// <summary>The scale c in scaled ≈ c × reference, by the robust adjusted profile score without its robust loss
        /// (Zhao et al. 2020): it allows for the reference effects' own error and for overdispersion τ².</summary>
        public double Scale { get; internal set; }
        /// <summary>Its standard error (sandwich).</summary>
        public double StandardError { get; internal set; }
        /// <summary>Lower bound of the confidence interval (normal reference).</summary>
        public double ConfidenceLow { get; internal set; }
        /// <summary>Upper bound of the confidence interval (normal reference).</summary>
        public double ConfidenceHigh { get; internal set; }
        /// <summary>Overdispersion τ²: variance of the scaled effects beyond their standard errors, in their units.</summary>
        public double Tau2 { get; internal set; }
        /// <summary>Standard error of <see cref="Tau2"/> (sandwich).</summary>
        public double Tau2StandardError { get; internal set; }
        /// <summary>The inverse-variance-weighted slope through the origin, Σ x y / s_y² ÷ Σ x² / s_y². It ignores the
        /// reference effects' error, so it is shrunk toward zero (regression dilution); shown so the shrinkage is visible.</summary>
        public double InverseVarianceScale { get; internal set; }
        /// <summary>The method-of-moments slope through the origin (Akritas and Bershady 1996), Σ x y ÷ (Σ x² − Σ s_x²).
        /// NaN when the reference effects' error variance is not smaller than their spread.</summary>
        public double BcesScale { get; internal set; }
        /// <summary>Bowden et al.'s (2016) I²_GX, in [0, 1]: the share of the reference effects' spread that is not
        /// their measurement error. The inverse-variance slope is shrunk by about this factor; below 0.9 the reference
        /// effects are too noisy to trust it.</summary>
        public double ISquaredGx { get; internal set; }
        /// <summary>Rounds of the alternating fit that gave Newton's method its start.</summary>
        public int Iterations { get; internal set; }
    }

    /// <summary>
    /// The multiplicative scale c between two sets of effect estimates of the same features, scaled ≈ c × reference,
    /// when both carry measurement error: for example a TMT dataset's effects against label-free effects of the same
    /// proteins, where co-isolation shrinks the TMT effects toward zero.
    /// </summary>
    /// <remarks>
    /// The model is that of Zhao et al.'s (2020, Ann. Stat. 48:1742) overdispersed two-sample estimator:
    /// x_j ~ N(θ_j, s_xj²) and y_j ~ N(c θ_j, s_yj² + τ²). With D_j = c² s_xj² + s_yj² + τ², c and τ² solve
    /// Σ (y − c x)(y s_x² c + x (s_y² + τ²)) / D² = 0 and Σ s_x² ((y − c x)² − D) / D² = 0. They are found as mr.raps 0.2's
    /// mr.raps.overdispersed finds them, by alternately maximising the two profile log-likelihoods whose derivatives
    /// these are, over the same bounded intervals and from the same start (its mr.raps.simple); the equations are then
    /// solved by Newton's method from there, because a maximum found by comparing values is placed only to about √ε. The
    /// standard errors are mr.raps's sandwich. Unlike mr.raps, the line goes through the origin by definition (no intercept) and the
    /// reference effects are not selected: callers must not choose features by the significance of these effects,
    /// which biases the scale. The robust-loss variant is not implemented.
    /// </remarks>
    public static class EffectScale
    {
        /// <summary>The fewest features for which a scale and an overdispersion are estimated.</summary>
        public const int MinimumFeatures = 3;

        /// <summary>The most rounds of the alternating fit.</summary>
        public const int MaxIterations = 100;

        /// <summary>The fit has settled when the relative changes in c and τ² over a round sum to at most this: the square
        /// root of machine epsilon, mr.raps's own tolerance. A maximum located by comparing a smooth function's values
        /// cannot be placed more finely than that (the function is flat at its peak), so a tighter tolerance never settles.</summary>
        public const double ConvergenceTolerance = 1.4901161193847656e-8;

        /// <summary>The relative tolerance of each one-dimensional maximisation.</summary>
        private const double LineTolerance = 1e-12;

        /// <summary>
        /// Estimates the scale c in scaled ≈ c × reference, with its standard error and interval, the overdispersion,
        /// and the inverse-variance and method-of-moments slopes and I²_GX beside it.
        /// </summary>
        /// <param name="reference">The reference effects (x), one per feature.</param>
        /// <param name="referenceStandardErrors">Their standard errors; finite and positive.</param>
        /// <param name="scaled">The scaled effects (y), the same features in the same order.</param>
        /// <param name="scaledStandardErrors">Their standard errors; finite and positive.</param>
        /// <param name="confidenceLevel">Coverage of the reported interval, e.g. 0.95.</param>
        /// <exception cref="ArgumentNullException">A list is null.</exception>
        /// <exception cref="ArgumentException">The lists differ in length, an effect is not finite, or a standard error
        /// is not finite and positive.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="confidenceLevel"/> is not in (0, 1).</exception>
        public static EffectScaleResult Estimate(IReadOnlyList<double> reference, IReadOnlyList<double> referenceStandardErrors,
            IReadOnlyList<double> scaled, IReadOnlyList<double> scaledStandardErrors, double confidenceLevel = 0.95)
        {
            ArgumentNullException.ThrowIfNull(reference);
            ArgumentNullException.ThrowIfNull(referenceStandardErrors);
            ArgumentNullException.ThrowIfNull(scaled);
            ArgumentNullException.ThrowIfNull(scaledStandardErrors);
            int k = reference.Count;
            if (referenceStandardErrors.Count != k || scaled.Count != k || scaledStandardErrors.Count != k)
                throw new ArgumentException("reference, scaled and their standard errors differ in length.", nameof(scaled));
            if (!(confidenceLevel > 0 && confidenceLevel < 1))
                throw new ArgumentOutOfRangeException(nameof(confidenceLevel), confidenceLevel, "The confidence level must be in (0, 1).");
            for (int j = 0; j < k; j++)
            {
                if (!double.IsFinite(reference[j]) || !double.IsFinite(scaled[j]))
                    throw new ArgumentException($"Feature {j} has an effect that is not finite.", nameof(reference));
                if (!(referenceStandardErrors[j] > 0) || !double.IsFinite(referenceStandardErrors[j])
                    || !(scaledStandardErrors[j] > 0) || !double.IsFinite(scaledStandardErrors[j]))
                    throw new ArgumentException($"Feature {j} has a standard error that is not finite and positive.", nameof(referenceStandardErrors));
            }

            double[] x = reference.ToArray(), sx2 = referenceStandardErrors.Select(s => s * s).ToArray();
            double[] y = scaled.ToArray(), sy2 = scaledStandardErrors.Select(s => s * s).ToArray();
            var result = new EffectScaleResult { Features = k };
            double[] ratios = Enumerable.Range(0, k).Where(j => x[j] != 0).Select(j => Math.Abs(y[j] / x[j])).ToArray();
            if (k < MinimumFeatures || ratios.Length == 0)
            {
                SetNotEstimable(result);
                return result;
            }

            result.InverseVarianceScale = Sum(k, j => x[j] * y[j] / sy2[j]) / Sum(k, j => x[j] * x[j] / sy2[j]);
            double bcesDenominator = Sum(k, j => x[j] * x[j]) - Sum(k, j => sx2[j]);
            result.BcesScale = bcesDenominator > 0 ? Sum(k, j => x[j] * y[j]) / bcesDenominator : double.NaN;
            double weightedMean = Sum(k, j => x[j] / sx2[j]) / Sum(k, j => 1 / sx2[j]);
            double q = Sum(k, j => (x[j] - weightedMean) * (x[j] - weightedMean) / sx2[j]);
            result.ISquaredGx = q > 0 ? Math.Max(0, (q - (k - 1)) / q) : 0;

            // Start as mr.raps.overdispersed does, from mr.raps.simple: the profile likelihood with no overdispersion.
            double ratioQuantile = Quantile7(ratios, 0.95);
            double scale = MaximizeScaleWithoutOverdispersion(x, y, sx2, sy2, 2 * ratioQuantile);
            double tauBound = 10 * Quantile7(sy2, 0.95);
            double scaleBound = 10 * ratioQuantile;
            double tau2 = 0;
            int iteration = 0;
            while (iteration < MaxIterations)
            {
                iteration++;
                double scaleBefore = scale, tauBefore = tau2;
                double s = scale;
                tau2 = MixedModel.Brent(t => -TauProfile(x, y, sx2, sy2, s, t), 0, tauBound, LineTolerance);
                (scale, scaleBound) = MaximizeScaleExtending(x, y, sx2, sy2, tau2, scaleBound);
                if (Math.Abs(scaleBefore - scale) / Math.Abs(scale + 1e-10) + Math.Abs(tauBefore - tau2) / Math.Abs(tau2 + 1e-10)
                    <= ConvergenceTolerance)
                {
                    break;
                }
            }

            // The alternation places the maxima only to about √ε and can creep for many rounds when c and τ² are coupled,
            // so the estimating equations themselves are then solved by Newton's method from where it stopped.
            bool converged = SolveEstimatingEquations(x, y, sx2, sy2, ref scale, ref tau2);

            // The sandwich of mr.raps.overdispersed: the scores' variance between the inverse Jacobian and its transpose.
            double v11 = 0, v22 = 0, i11 = 0, i12 = 0, i22 = 0;
            for (int j = 0; j < k; j++)
            {
                double ty = tau2 + sy2[j];
                double d2 = Math.Pow(ty + sx2[j] * scale * scale, 2);
                double core = (x[j] * x[j] - sx2[j]) * ty + (y[j] * y[j] - ty) * sx2[j];
                v11 += (core + sx2[j] * ty) / d2;
                v22 += 2 * sx2[j] * sx2[j] / d2;
                i11 -= core / d2;
                i12 -= sx2[j] * scale / d2;
                i22 -= sx2[j] / d2;
            }
            // Jacobian [[i11, i12], [0, i22]] (R fills column-major: I[1,2] is the derivative of the scale score in τ²).
            double det = i11 * i22;
            double a11 = i22 / det, a12 = -i12 / det, a22 = i11 / det;
            double var11 = a11 * a11 * v11 + a12 * a12 * v22;
            double var22 = a22 * a22 * v22;

            result.Status = converged ? EffectScaleStatus.Estimated : EffectScaleStatus.NotConverged;
            result.Iterations = iteration;
            result.Scale = scale;
            result.Tau2 = tau2;
            result.StandardError = Math.Sqrt(var11);
            result.Tau2StandardError = Math.Sqrt(var22);
            double z = Normal.InvCDF(0, 1, 0.5 + confidenceLevel / 2);
            result.ConfidenceLow = scale - z * result.StandardError;
            result.ConfidenceHigh = scale + z * result.StandardError;
            return result;
        }

        /// <summary>
        /// An effect on the reference scale: <paramref name="estimate"/> / c, with the standard error from the delta
        /// method, √(se² / c² + estimate² Var(c) / c⁴). The second term is shared by every effect rescaled by the same c,
        /// so rescaled effects of one dataset are not independent of each other. NaN for both when the scale was not
        /// estimated or is zero.
        /// <para>The formula has no covariance term: it assumes <paramref name="estimate"/> and c are uncorrelated. That holds
        /// for an effect that was not one of the features c was fitted on. For a calibration feature itself, its estimate
        /// helped fit c, so the two do covary, and the returned standard error leaves that covariance out.</para>
        /// </summary>
        public static (double Estimate, double StandardError) Rescale(double estimate, double standardError, EffectScaleResult scale)
        {
            ArgumentNullException.ThrowIfNull(scale);
            double c = scale.Scale;
            if (scale.Status != EffectScaleStatus.Estimated || !double.IsFinite(c) || c == 0)
                return (double.NaN, double.NaN);
            double varC = scale.StandardError * scale.StandardError;
            return (estimate / c, Math.Sqrt(standardError * standardError / (c * c) + estimate * estimate * varC / Math.Pow(c, 4)));
        }

        private static void SetNotEstimable(EffectScaleResult result)
        {
            result.Status = EffectScaleStatus.TooFewFeatures;
            result.Scale = result.StandardError = result.ConfidenceLow = result.ConfidenceHigh = double.NaN;
            result.Tau2 = result.Tau2StandardError = double.NaN;
            result.InverseVarianceScale = result.BcesScale = result.ISquaredGx = double.NaN;
        }

        private static double Sum(int k, Func<int, double> term)
        {
            double total = 0;
            for (int j = 0; j < k; j++) total += term(j);
            return total;
        }

        /// <summary>mr.raps's profile log-likelihood in τ² at a fixed scale, each feature weighted by s_x².</summary>
        private static double TauProfile(double[] x, double[] y, double[] sx2, double[] sy2, double scale, double tau2)
        {
            double total = 0;
            for (int j = 0; j < x.Length; j++)
            {
                double d = tau2 + sy2[j] + sx2[j] * scale * scale;
                double r = y[j] - x[j] * scale;
                total -= 0.5 * sx2[j] * Math.Log(d) + 0.5 * sx2[j] * r * r / d;
            }
            return total;
        }

        /// <summary>mr.raps's profile log-likelihood in the scale at a fixed τ² (τ² = 0 and no log term: mr.raps.simple).</summary>
        private static double ScaleProfile(double[] x, double[] y, double[] sx2, double[] sy2, double scale, double tau2)
        {
            double total = 0;
            for (int j = 0; j < x.Length; j++)
            {
                double r = y[j] - x[j] * scale;
                total -= 0.5 * r * r / (tau2 + sy2[j] + sx2[j] * scale * scale);
            }
            return total;
        }

        /// <summary>mr.raps.simple: the scale maximising the profile with no overdispersion on [−bound, bound], doubling
        /// the bound while the maximum is within 5% of it.</summary>
        private static double MaximizeScaleWithoutOverdispersion(double[] x, double[] y, double[] sx2, double[] sy2, double bound) =>
            MaximizeScaleExtending(x, y, sx2, sy2, 0, bound).Scale;

        /// <summary>The scale maximising the profile at τ² on [−bound, bound], doubling the bound while the maximum is
        /// within 5% of it; the doubled bound is kept for the next round, as mr.raps.overdispersed keeps it.</summary>
        private static (double Scale, double Bound) MaximizeScaleExtending(double[] x, double[] y, double[] sx2, double[] sy2,
            double tau2, double bound)
        {
            double scale = MixedModel.Brent(c => -ScaleProfile(x, y, sx2, sy2, c, tau2), -bound, bound, LineTolerance);
            for (int extension = 0; extension < 60 && Math.Abs(scale) > 0.95 * bound; extension++)
            {
                bound *= 2;
                scale = MixedModel.Brent(c => -ScaleProfile(x, y, sx2, sy2, c, tau2), -bound, bound, LineTolerance);
            }
            return (scale, bound);
        }

        /// <summary>The scale's estimating equation, the derivative of <see cref="ScaleProfile"/> in the scale.</summary>
        private static double ScaleScore(double[] x, double[] y, double[] sx2, double[] sy2, double scale, double tau2)
        {
            double total = 0;
            for (int j = 0; j < x.Length; j++)
            {
                double ty = tau2 + sy2[j];
                double d = ty + sx2[j] * scale * scale;
                total += (y[j] - scale * x[j]) * (y[j] * sx2[j] * scale + x[j] * ty) / (d * d);
            }
            return total;
        }

        /// <summary>The overdispersion's estimating equation, twice the derivative of <see cref="TauProfile"/> in τ².</summary>
        private static double TauScore(double[] x, double[] y, double[] sx2, double[] sy2, double scale, double tau2)
        {
            double total = 0;
            for (int j = 0; j < x.Length; j++)
            {
                double d = tau2 + sy2[j] + sx2[j] * scale * scale;
                double r = y[j] - scale * x[j];
                total += sx2[j] * (r * r - d) / (d * d);
            }
            return total;
        }

        /// <summary>
        /// Solves the two estimating equations by Newton's method (central-difference Jacobian) from the alternation's
        /// values. Where the overdispersion's equation is not positive at τ² = 0 for the solved scale, τ² stays on its
        /// boundary, 0, and only the scale's equation is solved, as the bounded maximisation would place it there.
        /// True when the steps fell below 1e-12 relative within <see cref="MaxIterations"/>; otherwise the values are
        /// left as the alternation gave them.
        /// </summary>
        private static bool SolveEstimatingEquations(double[] x, double[] y, double[] sx2, double[] sy2, ref double scale, ref double tau2)
        {
            double c = scale, t = tau2;
            for (int iteration = 0; iteration < MaxIterations; iteration++)
            {
                double hc = 1e-6 * Math.Max(Math.Abs(c), 1e-3);
                bool boundary = t <= 0 || TauScore(x, y, sx2, sy2, c, 0) <= 0;
                double stepC, stepT;
                if (boundary)
                {
                    t = 0;
                    double f = ScaleScore(x, y, sx2, sy2, c, 0);
                    double df = (ScaleScore(x, y, sx2, sy2, c + hc, 0) - ScaleScore(x, y, sx2, sy2, c - hc, 0)) / (2 * hc);
                    if (!(df != 0) || !double.IsFinite(df)) return false;
                    stepC = -f / df;
                    stepT = 0;
                }
                else
                {
                    double ht = 1e-6 * Math.Max(t, 1e-8);
                    double f1 = ScaleScore(x, y, sx2, sy2, c, t), f2 = TauScore(x, y, sx2, sy2, c, t);
                    double j11 = (ScaleScore(x, y, sx2, sy2, c + hc, t) - ScaleScore(x, y, sx2, sy2, c - hc, t)) / (2 * hc);
                    double j12 = (ScaleScore(x, y, sx2, sy2, c, t + ht) - ScaleScore(x, y, sx2, sy2, c, Math.Max(0, t - ht))) / (t + ht - Math.Max(0, t - ht));
                    double j21 = (TauScore(x, y, sx2, sy2, c + hc, t) - TauScore(x, y, sx2, sy2, c - hc, t)) / (2 * hc);
                    double j22 = (TauScore(x, y, sx2, sy2, c, t + ht) - TauScore(x, y, sx2, sy2, c, Math.Max(0, t - ht))) / (t + ht - Math.Max(0, t - ht));
                    double det = j11 * j22 - j12 * j21;
                    if (!(det != 0) || !double.IsFinite(det)) return false;
                    stepC = -(j22 * f1 - j12 * f2) / det;
                    stepT = -(-j21 * f1 + j11 * f2) / det;
                }
                c += stepC;
                t = Math.Max(0, t + stepT);
                if (!double.IsFinite(c) || !double.IsFinite(t)) return false;
                if (Math.Abs(stepC) <= 1e-12 * Math.Max(Math.Abs(c), 1e-10) && Math.Abs(stepT) <= 1e-12 * Math.Max(t, 1e-10))
                {
                    scale = c;
                    tau2 = t;
                    return true;
                }
            }
            return false;
        }

        /// <summary>R's default (type 7) sample quantile.</summary>
        private static double Quantile7(IReadOnlyList<double> values, double p)
        {
            double[] sorted = values.OrderBy(v => v).ToArray();
            double h = (sorted.Length - 1) * p;
            int lo = (int)Math.Floor(h);
            int hi = Math.Min(lo + 1, sorted.Length - 1);
            return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);
        }
    }
}
