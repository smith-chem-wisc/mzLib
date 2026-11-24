using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace StatisticalModels
{
    /// <summary>One pooled estimate across studies.</summary>
    public sealed class MetaAnalysisResult
    {
        internal MetaAnalysisResult() { }

        /// <summary>Number of studies pooled.</summary>
        public int Studies { get; internal set; }
        /// <summary>Random-effects pooled estimate.</summary>
        public double Estimate { get; internal set; }
        /// <summary>Standard error of the pooled estimate.</summary>
        public double StandardError { get; internal set; }
        /// <summary>Lower bound of the confidence interval.</summary>
        public double ConfidenceLow { get; internal set; }
        /// <summary>Upper bound of the confidence interval.</summary>
        public double ConfidenceHigh { get; internal set; }
        /// <summary>Two-sided p-value of the pooled estimate against zero, from the reference <see cref="Inference"/> names.</summary>
        public double PValue { get; internal set; }
        /// <summary>How the standard error, interval and p were formed. The estimate and τ² are DerSimonian–Laird's either way.</summary>
        public RandomEffectsInference Inference { get; internal set; }
        /// <summary>Degrees of freedom of the reference: +∞ for <see cref="RandomEffectsInference.Normal"/>, k − 1 for
        /// <see cref="RandomEffectsInference.HartungKnapp"/> and <see cref="RandomEffectsInference.HartungKnappAdHoc"/>.</summary>
        public double DegreesOfFreedom { get; internal set; }
        /// <summary>Between-study variance τ² (DerSimonian–Laird). Zero for a single study.</summary>
        public double Tau2 { get; internal set; }
        /// <summary>Cochran's Q. NaN for a single study.</summary>
        public double Q { get; internal set; }
        /// <summary>I², the share of variation attributed to heterogeneity, in [0, 1]. NaN for a single study,
        /// which is not a defect: one study has no heterogeneity to measure.</summary>
        public double ISquared { get; internal set; }
        /// <summary>Studies whose estimate has the same sign as the pooled one.</summary>
        public int DirectionAgree { get; internal set; }
        /// <summary>Studies whose estimate has the opposite sign.</summary>
        public int DirectionDisagree { get; internal set; }
        /// <summary>Largest |change| in the pooled estimate from dropping any one study. NaN for one study.</summary>
        public double LeaveOneOutMaxDelta { get; internal set; }
    }

    /// <summary>How a random-effects pool forms its standard error, interval and p.</summary>
    public enum RandomEffectsInference
    {
        /// <summary>The DerSimonian–Laird standard error, √(1/Σwᵢ*), with a normal reference.</summary>
        Normal,
        /// <summary>
        /// Hartung–Knapp (2001) and Sidik–Jonkman (2002): the variance Σwᵢ*(yᵢ − μ̂)² / ((k − 1) Σwᵢ*), with a t reference
        /// on k − 1 degrees of freedom (IntHout et al. 2014). Needs at least two studies. When every study gives the same
        /// estimate the standard error is 0 and the interval is that estimate.
        /// </summary>
        HartungKnapp,
        /// <summary>
        /// <see cref="HartungKnapp"/> with Knapp and Hartung's (2003) ad hoc variance correction: the standard error is the
        /// larger of the Hartung–Knapp and the DerSimonian–Laird ones, with the same t reference on k − 1 degrees of
        /// freedom (Röver, Knapp and Friede 2015; metafor's <c>test = "adhoc"</c>). Plain Hartung–Knapp's standard error
        /// shrinks towards 0 as the studies agree, whatever their own standard errors; this one never falls below the
        /// DerSimonian–Laird random-effects value, so it is never 0. Needs at least two studies.
        /// </summary>
        HartungKnappAdHoc,
    }

    /// <summary>
    /// Random-effects meta-analysis of per-study effect sizes by the DerSimonian–Laird (1986) moment
    /// estimator of the between-study variance.
    /// </summary>
    public static class RandomEffectsMeta
    {
        /// <summary>
        /// Pools one effect size per study into a random-effects estimate, with heterogeneity (τ², Q, I²),
        /// direction agreement and the largest leave-one-out change.
        /// </summary>
        /// <param name="estimates">One effect size per study.</param>
        /// <param name="standardErrors">Its standard error; must be finite and positive.</param>
        /// <param name="confidenceLevel">Coverage of the reported interval, e.g. 0.95.</param>
        /// <param name="inference">How the standard error, interval and p are formed; the normal reference by default.</param>
        /// <returns>
        /// The pooled result. Standard errors so small or so large (below about 1e-154, above about 1e154) that the weights
        /// 1/se² cannot be formed in double precision give a NaN estimate or a standard error of 0; a caller checks them.
        /// </returns>
        /// <exception cref="ArgumentException"><see cref="RandomEffectsInference.HartungKnapp"/> or
        /// <see cref="RandomEffectsInference.HartungKnappAdHoc"/> with one study, which has no degrees of freedom; and the
        /// input checks.</exception>
        public static MetaAnalysisResult Pool(IReadOnlyList<double> estimates, IReadOnlyList<double> standardErrors,
            double confidenceLevel = 0.95, RandomEffectsInference inference = RandomEffectsInference.Normal)
        {
            ArgumentNullException.ThrowIfNull(estimates);
            ArgumentNullException.ThrowIfNull(standardErrors);
            int k = estimates.Count;
            if (k == 0) throw new ArgumentException("No studies to pool.", nameof(estimates));
            if (standardErrors.Count != k) throw new ArgumentException("estimates and standardErrors differ in length.", nameof(standardErrors));
            for (int i = 0; i < k; i++)
            {
                if (!double.IsFinite(estimates[i])) throw new ArgumentException($"Estimate {i} is not finite.", nameof(estimates));
                if (!(standardErrors[i] > 0) || !double.IsFinite(standardErrors[i]))
                    throw new ArgumentException($"Standard error {i} must be finite and positive.", nameof(standardErrors));
            }
            if (!(confidenceLevel > 0 && confidenceLevel < 1)) throw new ArgumentOutOfRangeException(nameof(confidenceLevel));
            if (!Enum.IsDefined(inference)) throw new ArgumentOutOfRangeException(nameof(inference));
            bool hartungKnapp = inference is RandomEffectsInference.HartungKnapp or RandomEffectsInference.HartungKnappAdHoc;
            if (hartungKnapp && k < 2)
                throw new ArgumentException("Hartung–Knapp needs at least two studies: one study has no degrees of freedom.", nameof(estimates));

            var (estimate, se, tau2, q) = DerSimonianLaird(estimates, standardErrors);
            double df = double.PositiveInfinity, quantile, p;
            if (hartungKnapp)
            {
                double seDerSimonianLaird = se;
                df = k - 1;
                double[] ws = standardErrors.Select(s => 1 / (s * s + tau2)).ToArray();
                double sws = ws.Sum();
                // Residuals about the pooled estimate, taken from the first study so that equal estimates give exactly 0.
                double shift = ws.Zip(estimates, (wi, bi) => wi * (bi - estimates[0])).Sum() / sws;
                double spread = ws.Zip(estimates, (wi, bi) => wi * Math.Pow(bi - estimates[0] - shift, 2)).Sum();
                se = Math.Sqrt(spread / (df * sws));
                if (inference == RandomEffectsInference.HartungKnappAdHoc) se = Math.Max(se, seDerSimonianLaird);
                quantile = StudentT.InvCDF(0, 1, df, 0.5 + confidenceLevel / 2);
                p = 2 * StudentT.CDF(0, 1, df, -Math.Abs(estimate / se));
            }
            else
            {
                quantile = Normal.InvCDF(0, 1, 0.5 + confidenceLevel / 2);
                p = 2 * Normal.CDF(0, 1, -Math.Abs(estimate / se));
            }
            // Standard errors far outside any measured range (below about 1e-154, above about 1e154) make the weights 1/se²
            // overflow or underflow, and the estimate NaN; it then agrees and disagrees with no study.
            int sign = double.IsNaN(estimate) ? 0 : Math.Sign(estimate);
            var result = new MetaAnalysisResult
            {
                Studies = k,
                Estimate = estimate,
                StandardError = se,
                ConfidenceLow = estimate - quantile * se,
                ConfidenceHigh = estimate + quantile * se,
                PValue = p,
                Inference = inference,
                DegreesOfFreedom = df,
                Tau2 = tau2,
                Q = k > 1 ? q : double.NaN,
                ISquared = k > 1 ? (q > 0 ? Math.Max(0, (q - (k - 1)) / q) : 0) : double.NaN,
                DirectionAgree = estimates.Count(b => Math.Sign(b) == sign && b != 0),
                DirectionDisagree = estimates.Count(b => Math.Sign(b) == -sign && b != 0),
                LeaveOneOutMaxDelta = double.NaN,
            };
            if (k > 1)
            {
                double max = 0;
                for (int drop = 0; drop < k; drop++)
                {
                    var b = estimates.Where((_, i) => i != drop).ToList();
                    var s = standardErrors.Where((_, i) => i != drop).ToList();
                    max = Math.Max(max, Math.Abs(DerSimonianLaird(b, s).estimate - estimate));
                }
                result.LeaveOneOutMaxDelta = max;
            }
            return result;
        }

        private static (double estimate, double se, double tau2, double q) DerSimonianLaird(
            IReadOnlyList<double> b, IReadOnlyList<double> se)
        {
            int k = b.Count;
            double[] w = se.Select(s => 1 / (s * s)).ToArray();
            double sw = w.Sum();
            double fixedEffect = w.Zip(b, (wi, bi) => wi * bi).Sum() / sw;
            double q = w.Zip(b, (wi, bi) => wi * (bi - fixedEffect) * (bi - fixedEffect)).Sum();
            double tau2 = 0;
            if (k > 1)
            {
                double c = sw - w.Sum(wi => wi * wi) / sw;
                tau2 = c > 0 ? Math.Max(0, (q - (k - 1)) / c) : 0;
            }
            double[] ws = se.Select(s => 1 / (s * s + tau2)).ToArray();
            double sws = ws.Sum();
            double estimate = ws.Zip(b, (wi, bi) => wi * bi).Sum() / sws;
            return (estimate, Math.Sqrt(1 / sws), tau2, q);
        }
    }
}
