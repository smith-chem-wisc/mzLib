using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;

namespace Quantification.DifferentialAbundance
{
    /// <summary>
    /// The prior distribution of per-feature residual variances: s²_g ~ s0² · χ²(d0)/d0, fitted across
    /// all features. With a trend, s0² depends on a covariate (the feature's average response).
    /// </summary>
    public sealed class VariancePrior
    {
        /// <summary>
        /// A variance of exactly zero has no logarithm, so variances are floored at this fraction of the
        /// median variance before the prior is fitted.
        /// </summary>
        public const double ZeroVarianceFloor = 1e-5;

        internal VariancePrior(double df, double[] scale, bool trended, int splineBasisCount)
        {
            Df = df;
            Scale = scale;
            Trended = trended;
            SplineBasisCount = splineBasisCount;
        }

        /// <summary>Prior degrees of freedom d0. Positive infinity when the observed variances are no more
        /// dispersed than sampling alone explains, in which case every feature takes the prior variance.</summary>
        public double Df { get; }

        /// <summary>Prior variance s0², one per feature (all equal when <see cref="Trended"/> is false).
        /// NaN for a feature that did not enter the fit.</summary>
        public IReadOnlyList<double> Scale { get; }

        /// <summary>Whether s0² was allowed to vary with the covariate.</summary>
        public bool Trended { get; }

        /// <summary>Basis functions (intercept included) of the natural cubic spline the trend used; 1 without a trend.</summary>
        public int SplineBasisCount { get; }
    }

    /// <summary>Moderated t-statistics for one coefficient across all features. Read-only to callers.</summary>
    public sealed class ModeratedTest
    {
        internal ModeratedTest(string coefficient, VariancePrior prior, IReadOnlyList<FeatureFitStatus> status)
        {
            int n = status.Count;
            double[] Nan() => Enumerable.Repeat(double.NaN, n).ToArray();
            Coefficient = coefficient;
            Prior = prior;
            Status = status.ToArray();
            EstimateValues = Nan(); StandardErrorValues = Nan(); TValues = Nan(); DfTotalValues = Nan();
            PValues = Nan(); AdjustedValues = Nan(); PosteriorVarianceValues = Nan();
        }

        internal double[] EstimateValues { get; }
        internal double[] StandardErrorValues { get; }
        internal double[] TValues { get; }
        internal double[] DfTotalValues { get; }
        internal double[] PValues { get; }
        internal double[] AdjustedValues { get; }
        internal double[] PosteriorVarianceValues { get; }

        /// <summary>The coefficient tested.</summary>
        public string Coefficient { get; }
        /// <summary>The fitted prior.</summary>
        public VariancePrior Prior { get; }
        /// <summary>Per-feature status carried over from the fit.</summary>
        public IReadOnlyList<FeatureFitStatus> Status { get; }
        /// <summary>Least-squares estimate of the coefficient (moderation does not change it).</summary>
        public IReadOnlyList<double> Estimate => EstimateValues;
        /// <summary>Moderated standard error: sqrt(posterior variance) × unscaled standard deviation.</summary>
        public IReadOnlyList<double> StandardError => StandardErrorValues;
        /// <summary>Moderated t-statistic.</summary>
        public IReadOnlyList<double> T => TValues;
        /// <summary>Degrees of freedom of the moderated t: residual df plus prior df, capped at the pooled residual df.</summary>
        public IReadOnlyList<double> DfTotal => DfTotalValues;
        /// <summary>Two-sided p-value.</summary>
        public IReadOnlyList<double> PValue => PValues;
        /// <summary>Benjamini–Hochberg adjusted p-values over the features that were tested. Not a target-decoy q-value.</summary>
        public IReadOnlyList<double> BenjaminiHochbergAdjusted => AdjustedValues;
        /// <summary>Posterior (moderated) residual variance.</summary>
        public IReadOnlyList<double> PosteriorVariance => PosteriorVarianceValues;
    }

    /// <summary>
    /// Empirical-Bayes moderation of per-feature variances (Smyth, 2004, Stat. Appl. Genet. Mol. Biol.
    /// 3:3), optionally with a prior variance that trends with average intensity.
    /// </summary>
    /// <remarks>
    /// Implemented from the paper's method-of-moments estimator on log variances, not from any existing
    /// implementation's source. Validated against reference outputs rather than by construction.
    /// </remarks>
    public static class EmpiricalBayes
    {
        /// <summary>Default number of spline basis functions (intercept included) for an intensity trend.</summary>
        public const int DefaultSplineBasisCount = 4;

        /// <summary>
        /// Fits the variance prior by matching moments of e_g = log s²_g − ψ(d_g/2) + log(d_g/2), whose
        /// variance is ψ′(d_g/2) + ψ′(d0/2).
        /// </summary>
        /// <param name="variances">Residual variances s²_g; non-finite or negative entries are ignored.</param>
        /// <param name="df">Residual degrees of freedom d_g; entries &lt;= 0 are ignored.</param>
        /// <param name="covariate">If given, the prior mean of e_g is a natural cubic spline in this covariate.</param>
        /// <param name="splineBasisCount">Basis functions for the trend, intercept included.</param>
        public static VariancePrior FitPrior(IReadOnlyList<double> variances, IReadOnlyList<double> df,
            IReadOnlyList<double>? covariate = null, int splineBasisCount = DefaultSplineBasisCount)
        {
            ArgumentNullException.ThrowIfNull(variances);
            ArgumentNullException.ThrowIfNull(df);
            int n = variances.Count;
            if (df.Count != n) throw new ArgumentException("variances and df differ in length.", nameof(df));
            if (covariate != null && covariate.Count != n) throw new ArgumentException("covariate and variances differ in length.", nameof(covariate));

            var use = Enumerable.Range(0, n).Where(i =>
                double.IsFinite(variances[i]) && variances[i] >= 0 && df[i] > 0
                && (covariate == null || double.IsFinite(covariate[i]))).ToArray();
            var scale = Enumerable.Repeat(double.NaN, n).ToArray();
            if (use.Length < 2)
                throw new ArgumentException($"A variance prior needs at least 2 usable features; {use.Length} were usable.", nameof(variances));

            // A variance of exactly zero has no logarithm. Offset them away from zero relative to the median.
            double median = use.Select(i => variances[i]).Median();
            double floor = median > 0 ? VariancePrior.ZeroVarianceFloor * median : VariancePrior.ZeroVarianceFloor;
            double[] e = use.Select(i =>
                Math.Log(Math.Max(variances[i], floor)) - SpecialFunctions.DiGamma(df[i] / 2) + Math.Log(df[i] / 2)).ToArray();

            double[] mean;
            double residualVariance;
            int basis = 1;
            if (covariate == null)
            {
                double m = e.Average();
                mean = Enumerable.Repeat(m, e.Length).ToArray();
                residualVariance = e.Sum(v => (v - m) * (v - m)) / (e.Length - 1);
            }
            else
            {
                double[] x = use.Select(i => covariate[i]).ToArray();
                basis = Math.Min(splineBasisCount, Math.Min(x.Distinct().Count(), e.Length - 1));
                basis = Math.Max(basis, 1);
                var design = NaturalSplineBasis(x, basis);
                var qr = design.QR(MathNet.Numerics.LinearAlgebra.Factorization.QRMethod.Thin);
                var ev = Vector<double>.Build.DenseOfArray(e);
                var fitted = design * qr.Solve(ev);
                mean = fitted.ToArray();
                var r = ev - fitted;
                residualVariance = r.DotProduct(r) / (e.Length - basis);
            }

            double excess = residualVariance - use.Average(i => Polygamma.Trigamma(df[i] / 2));
            double d0 = excess > 0 ? 2 * Polygamma.TrigammaInverse(excess) : double.PositiveInfinity;
            for (int k = 0; k < use.Length; k++)
                scale[use[k]] = double.IsPositiveInfinity(d0)
                    ? Math.Exp(mean[k])
                    : Math.Exp(mean[k] + SpecialFunctions.DiGamma(d0 / 2) - Math.Log(d0 / 2));

            return new VariancePrior(d0, scale, covariate != null, basis);
        }

        /// <summary>
        /// Moderated t-test of one coefficient: each fitted feature's residual variance is shrunk toward the
        /// prior, s̃² = (d0·s0² + d·s²) / (d0 + d), and t = β / (s̃ · unscaled SD) on d + d0 degrees of freedom.
        /// </summary>
        /// <param name="fit">Output of <see cref="LinearModel.Fit"/>.</param>
        /// <param name="coefficient">Name of the coefficient to test, e.g. "age_decades".</param>
        /// <param name="trend">Let the prior variance trend with <see cref="LinearModelFit.AverageResponse"/>.</param>
        /// <param name="splineBasisCount">Basis functions for the trend, intercept included.</param>
        public static ModeratedTest Moderate(LinearModelFit fit, string coefficient, bool trend,
            int splineBasisCount = DefaultSplineBasisCount)
        {
            ArgumentNullException.ThrowIfNull(fit);
            int j = fit.IndexOf(coefficient);
            int n = fit.FeatureCount;
            var fitted = Enumerable.Range(0, n).Where(f => fit.Status[f] == FeatureFitStatus.Fitted).ToArray();
            if (fitted.Length < 2)
                throw new ArgumentException(
                    $"Moderation needs at least 2 fitted features to estimate a variance prior; this fit has {fitted.Length} " +
                    $"of {n} (the rest are too sparse or rank-deficient).", nameof(fit));

            var s2 = new double[n];
            var df = new double[n];
            for (int f = 0; f < n; f++)
            {
                bool ok = fit.Status[f] == FeatureFitStatus.Fitted;
                s2[f] = ok ? fit.Sigma[f] * fit.Sigma[f] : double.NaN;
                df[f] = ok ? fit.DfResidual[f] : 0;
            }
            var prior = FitPrior(s2, df, trend ? fit.AverageResponse : null, splineBasisCount);
            double pooledDf = fitted.Sum(f => (double)fit.DfResidual[f]);

            var test = new ModeratedTest(coefficient, prior, fit.Status);
            foreach (int f in fitted)
            {
                double post = double.IsPositiveInfinity(prior.Df)
                    ? prior.Scale[f]
                    : (prior.Df * prior.Scale[f] + df[f] * s2[f]) / (prior.Df + df[f]);
                double dfTotal = Math.Min(df[f] + prior.Df, pooledDf);
                double se = Math.Sqrt(post) * fit.StdevUnscaled(f, j);
                double t = fit.Coefficient(f, j) / se;
                test.EstimateValues[f] = fit.Coefficient(f, j);
                test.PosteriorVarianceValues[f] = post;
                test.StandardErrorValues[f] = se;
                test.TValues[f] = t;
                test.DfTotalValues[f] = dfTotal;
                test.PValues[f] = TwoSidedP(t, dfTotal);
            }
            var adjusted = MultipleTesting.BenjaminiHochberg(test.PValues);
            Array.Copy(adjusted, test.AdjustedValues, n);
            return test;
        }

        internal static double TwoSidedP(double t, double df)
        {
            if (!double.IsFinite(t)) return double.NaN;
            double a = -Math.Abs(t);
            double tail = double.IsPositiveInfinity(df) ? Normal.CDF(0, 1, a) : StudentT.CDF(0, 1, df, a);
            return Math.Min(1.0, 2 * tail);
        }

        /// <summary>
        /// A basis of natural cubic splines (linear beyond the boundary knots) with <paramref name="count"/>
        /// functions including the intercept. Knots: the range of x, plus count − 2 interior knots at evenly
        /// spaced sample quantiles. Only the spanned space matters to a least-squares fit, so any basis of it
        /// gives the same fitted values; this is the truncated-power form, on x rescaled to [0, 1].
        /// </summary>
        internal static Matrix<double> NaturalSplineBasis(double[] x, int count)
        {
            int n = x.Length;
            double lo = x.Min(), hi = x.Max();
            double span = hi > lo ? hi - lo : 1;
            double[] u = x.Select(v => (v - lo) / span).ToArray();
            if (count <= 1) return Matrix<double>.Build.Dense(n, 1, 1.0);
            if (count == 2) return Matrix<double>.Build.Dense(n, 2, (i, j) => j == 0 ? 1 : u[i]);

            var sorted = u.OrderBy(v => v).ToArray();
            int k = count;                       // total knots = basis functions
            var knots = new double[k];
            knots[0] = 0; knots[k - 1] = 1;
            for (int i = 1; i < k - 1; i++)
                knots[i] = Statistics.QuantileCustom(sorted, (double)i / (k - 1), QuantileDefinition.R7);

            double D(double v, int idx)
            {
                double a = Math.Max(0, v - knots[idx]), b = Math.Max(0, v - knots[k - 1]);
                return (a * a * a - b * b * b) / (knots[k - 1] - knots[idx]);
            }
            return Matrix<double>.Build.Dense(n, k, (i, j) =>
                j == 0 ? 1 : j == 1 ? u[i] : D(u[i], j - 2) - D(u[i], k - 2));
        }
    }
}
