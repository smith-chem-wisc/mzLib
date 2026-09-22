using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace Quantification.DifferentialAbundance
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
        /// <summary>Two-sided p-value of the pooled estimate against zero (normal reference).</summary>
        public double PValue { get; internal set; }
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
        public static MetaAnalysisResult Pool(IReadOnlyList<double> estimates, IReadOnlyList<double> standardErrors,
            double confidenceLevel = 0.95)
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

            var (estimate, se, tau2, q) = DerSimonianLaird(estimates, standardErrors);
            double z = Normal.InvCDF(0, 1, 0.5 + confidenceLevel / 2);
            var result = new MetaAnalysisResult
            {
                Studies = k,
                Estimate = estimate,
                StandardError = se,
                ConfidenceLow = estimate - z * se,
                ConfidenceHigh = estimate + z * se,
                PValue = 2 * Normal.CDF(0, 1, -Math.Abs(estimate / se)),
                Tau2 = tau2,
                Q = k > 1 ? q : double.NaN,
                ISquared = k > 1 ? (q > 0 ? Math.Max(0, (q - (k - 1)) / q) : 0) : double.NaN,
                DirectionAgree = estimates.Count(b => Math.Sign(b) == Math.Sign(estimate) && b != 0),
                DirectionDisagree = estimates.Count(b => Math.Sign(b) == -Math.Sign(estimate) && b != 0),
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
