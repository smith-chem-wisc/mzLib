using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;

namespace StatisticalModels
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
            Scale = Array.AsReadOnly(scale);
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
        internal ModeratedTest(string coefficient, VariancePrior prior, IReadOnlyList<FeatureFitStatus> status,
            bool residualDfDiffer)
        {
            int n = status.Count;
            double[] Nan() => Enumerable.Repeat(double.NaN, n).ToArray();
            Coefficient = coefficient;
            Prior = prior;
            Status = status.ToArray();
            ResidualDfDiffer = residualDfDiffer;
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
        /// <summary>
        /// True when the fitted features do not all have the same residual degrees of freedom, as happens
        /// whenever missing values are omitted. For such input this result follows limma's legacy estimator
        /// (<c>eBayes(legacy = TRUE)</c>), and limma 3.61 and later, by default, would give a different
        /// prior and so different moderated statistics. See <see cref="EmpiricalBayes"/>.
        /// </summary>
        public bool ResidualDfDiffer { get; }
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
    /// <para>
    /// Implemented from the paper's method-of-moments estimator on log variances, not from any existing
    /// implementation's source. It reproduces limma's <c>eBayes(legacy = TRUE)</c> to 1e-8 relative
    /// (see the reference comparison tests), with the intensity trend and without.
    /// </para>
    /// <para>
    /// <b>This is limma's legacy estimator, not its current default.</b> Since limma 3.61, when residual
    /// degrees of freedom differ between features, <c>eBayes</c> uses a different prior estimator
    /// (<c>fitFDistUnequalDF1</c>: features with fewer degrees of freedom are down-weighted, the prior df
    /// is estimated by profile likelihood, and the trend is a lowess curve). Omitting missing values per
    /// feature makes residual df differ, so for most label-free proteomics data these results differ from
    /// default limma. <see cref="ModeratedTest.ResidualDfDiffer"/> reports when that applies. With equal
    /// residual df, default limma also uses the legacy estimator unless asked otherwise.
    /// </para>
    /// <para>
    /// Not implemented: <c>robust = TRUE</c>, contrasts, the B-statistic, observation weights.
    /// </para>
    /// </remarks>
    public static class EmpiricalBayes
    {
        /// <summary>
        /// The most spline basis functions (intercept included) the default intensity trend uses. As in
        /// limma's <c>fitFDist</c>, the default is 1 + (n ≥ 3) + (n ≥ 6) + (n ≥ 30) for n usable features,
        /// capped at the number of distinct covariate values, so it reaches this only from 30 features.
        /// </summary>
        public const int DefaultSplineBasisCount = 4;

        /// <summary>
        /// Fits the variance prior by matching moments of e_g = log s²_g − ψ(d_g/2) + log(d_g/2), whose
        /// variance is ψ′(d_g/2) + ψ′(d0/2).
        /// </summary>
        /// <param name="variances">Residual variances s²_g; non-finite or negative entries are ignored.</param>
        /// <param name="df">Residual degrees of freedom d_g; entries &lt;= 0 are ignored.</param>
        /// <param name="covariate">If given, the prior mean of e_g is a natural cubic spline in this covariate.
        /// It must be finite for every usable feature (limma refuses a missing covariate too).</param>
        /// <param name="splineBasisCount">Basis functions for the trend, intercept included. Null (the default)
        /// takes limma's rule for the number of usable features (see <see cref="DefaultSplineBasisCount"/>).
        /// Either way it is capped at the distinct covariate values, and knots that tied covariate values
        /// make coincide are merged, which can leave fewer; <see cref="VariancePrior.SplineBasisCount"/>
        /// reports how many were used.</param>
        public static VariancePrior FitPrior(IReadOnlyList<double> variances, IReadOnlyList<double> df,
            IReadOnlyList<double>? covariate = null, int? splineBasisCount = null)
        {
            ArgumentNullException.ThrowIfNull(variances);
            ArgumentNullException.ThrowIfNull(df);
            int n = variances.Count;
            if (df.Count != n) throw new ArgumentException("variances and df differ in length.", nameof(df));
            if (covariate != null && covariate.Count != n) throw new ArgumentException("covariate and variances differ in length.", nameof(covariate));

            var use = Enumerable.Range(0, n).Where(i =>
                double.IsFinite(variances[i]) && variances[i] >= 0 && df[i] > 0).ToArray();
            if (covariate != null && use.FirstOrDefault(i => !double.IsFinite(covariate[i]), -1) is int bad and >= 0)
                throw new ArgumentException($"The covariate of feature {bad} is {covariate[bad]}; a trend needs a finite covariate for every usable feature.", nameof(covariate));
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
                var design = NaturalSplineBasis(x, TrendBasisCount(e.Length, x.Distinct().Count(), splineBasisCount));
                basis = design.ColumnCount;
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
        /// <param name="splineBasisCount">Basis functions for the trend, intercept included; null takes limma's
        /// rule (see <see cref="FitPrior"/>).</param>
        public static ModeratedTest Moderate(LinearModelFit fit, string coefficient, bool trend,
            int? splineBasisCount = null)
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

            bool dfDiffer = fitted.Any(f => fit.DfResidual[f] != fit.DfResidual[fitted[0]]);
            var test = new ModeratedTest(coefficient, prior, fit.Status, dfDiffer);
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
        /// Basis functions for a trend over <paramref name="usable"/> features whose covariate takes
        /// <paramref name="distinct"/> values. Without a request, limma's <c>fitFDist</c> rule; a request is
        /// capped so the fit keeps at least one residual degree of freedom.
        /// </summary>
        internal static int TrendBasisCount(int usable, int distinct, int? requested)
        {
            int count = requested is int r
                ? Math.Min(r, usable - 1)
                : 1 + (usable >= 3 ? 1 : 0) + (usable >= 6 ? 1 : 0) + (usable >= 30 ? 1 : 0);
            return Math.Max(1, Math.Min(count, distinct));
        }

        /// <summary>
        /// A basis of natural cubic splines (linear beyond the boundary knots) with <paramref name="count"/>
        /// functions including the intercept. Knots: the range of x, plus count − 2 interior knots at evenly
        /// spaced sample quantiles. Only the spanned space matters to a least-squares fit, so any basis of it
        /// gives the same fitted values; this is the truncated-power form, on x rescaled to [0, 1].
        /// Heavily tied x can put two knots on one value. Coinciding knots are merged, so the basis can
        /// come back with fewer columns than asked (down to the straight line, or the intercept alone when
        /// x is constant). R's <c>ns</c> would keep a repeated knot instead, so limma's trend differs there.
        /// </summary>
        internal static Matrix<double> NaturalSplineBasis(double[] x, int count)
        {
            int n = x.Length;
            double lo = x.Min(), hi = x.Max();
            double span = hi > lo ? hi - lo : 1;
            double[] u = x.Select(v => (v - lo) / span).ToArray();
            if (count <= 1) return Matrix<double>.Build.Dense(n, 1, 1.0);
            if (count == 2) return Matrix<double>.Build.Dense(n, 2, (i, j) => j == 0 ? 1 : u[i]);

            if (!(hi > lo)) return Matrix<double>.Build.Dense(n, 1, 1.0);
            var sorted = u.OrderBy(v => v).ToArray();
            const double tie = 1e-12;            // on the [0, 1] scale
            var knotList = new List<double> { 0 };
            for (int i = 1; i < count - 1; i++)
            {
                double q = MathNet.Numerics.Statistics.Statistics.QuantileCustom(sorted, (double)i / (count - 1), QuantileDefinition.R7);
                if (q > knotList[^1] + tie && q < 1 - tie) knotList.Add(q);
            }
            knotList.Add(1);
            var knots = knotList.ToArray();
            int k = knots.Length;                // total knots = basis functions
            if (k == 2) return Matrix<double>.Build.Dense(n, 2, (i, j) => j == 0 ? 1 : u[i]);

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
