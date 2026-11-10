using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;

namespace StatisticalModels
{
    /// <summary>How the variance prior's hyperparameters (d0, s0²) are estimated.</summary>
    public enum VariancePriorEstimator
    {
        /// <summary>
        /// Method of moments on log variances (Smyth 2004), as limma's <c>eBayes(legacy = TRUE)</c>. Exact when
        /// every feature has the same residual df; with unequal df it treats them through an averaged trigamma
        /// term and can return d0 = ∞ where the data do not support it.
        /// </summary>
        MomentsLegacy,
        /// <summary>
        /// Maximum marginal likelihood: s²_g / s0²_g ~ F(d_g, d0) with each feature's OWN residual df, maximized
        /// over d0 (∞ allowed) and s0² (or its trend). Unequal df, as omitting missing values produces, are
        /// handled exactly by the likelihood. An independent estimator, not limma's <c>fitFDistUnequalDF1</c>,
        /// so its numbers differ from default limma's. d0 is searched from about 0.05 upward, so a reported
        /// d0 at e^-3 ≈ 0.05 is that floor, not an interior optimum.
        /// </summary>
        MarginalLikelihood,
    }

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

        internal VariancePrior(double df, double[] scale, bool trended, int splineBasisCount,
            VariancePriorEstimator estimator = VariancePriorEstimator.MomentsLegacy)
        {
            Estimator = estimator;
            Df = df;
            Scale = Array.AsReadOnly(scale);
            Trended = trended;
            SplineBasisCount = splineBasisCount;
        }

        /// <summary>How d0 and s0² were estimated.</summary>
        public VariancePriorEstimator Estimator { get; }

        /// <summary>Prior degrees of freedom d0. Positive infinity when the observed variances are no more
        /// dispersed than sampling alone explains, in which case every feature takes the prior variance. Zero with fewer
        /// than 3 features, as in limma: no prior is fitted and each feature keeps its own variance.</summary>
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
        /// whenever missing values are omitted. With <see cref="VariancePriorEstimator.MomentsLegacy"/> this
        /// result then follows limma's legacy estimator (<c>eBayes(legacy = TRUE)</c>), and limma 3.61 and
        /// later, by default, would give a different prior. <see cref="VariancePriorEstimator.MarginalLikelihood"/>
        /// handles unequal df exactly. The prior's estimator is <see cref="VariancePrior.Estimator"/>.
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

    /// <summary>A named contrast: one weight per coefficient of the fit, in <see cref="LinearModelFit.CoefficientNames"/> order.</summary>
    /// <param name="Name">Unique within one call.</param>
    /// <param name="Weights">Finite, not all zero; e.g. (0, −1, 1, 0) for the third coefficient minus the second.</param>
    public sealed record ContrastWeights(string Name, IReadOnlyList<double> Weights);

    /// <summary>
    /// One linear-model fit and the contrasts tested on it, as a part of a pooled moderation
    /// (<see cref="EmpiricalBayes.ModerateContrasts(IReadOnlyList{LinearModelModerationPart}, double)"/>).
    /// </summary>
    /// <param name="Fit">Output of <see cref="LinearModel.Fit"/>.</param>
    /// <param name="Contrasts">Weights on this fit's coefficients; names unique within the part; may be empty.</param>
    public sealed record LinearModelModerationPart(LinearModelFit Fit, IReadOnlyList<ContrastWeights> Contrasts);

    /// <summary>
    /// One fitted feature's part in a contrast c, as a model gives it to the shared moderation
    /// (<see cref="EmpiricalBayes.ModerateContrastsCore"/>).
    /// </summary>
    /// <param name="Estimate">c′β.</param>
    /// <param name="StdevUnscaled">sqrt(c′Uc), U the coefficients' covariance divided by the residual variance.</param>
    /// <param name="Df">
    /// The degrees of freedom of this contrast's residual variance, before the prior's are added: the residual df of a
    /// linear model; a mixed model's Satterthwaite df for c. Not finite: the feature has no result for this contrast.
    /// </param>
    internal readonly record struct FeatureContrast(double Estimate, double StdevUnscaled, double Df);

    /// <summary>
    /// A moderated test of one contrast across all features, with its confidence interval, from a linear model
    /// (<see cref="EmpiricalBayes.ModerateContrasts"/>) or a nested mixed model (<see cref="NestedMixedModel.ModerateContrasts"/>).
    /// Read-only to callers. Every per-feature value is NaN for a feature with no result: one that was not fitted, or, for
    /// a mixed model, one whose variance df or contrast df is not finite.
    /// </summary>
    public sealed class ModeratedContrast
    {
        internal ModeratedContrast(string name, IReadOnlyList<double> weights, VariancePrior prior,
            IReadOnlyList<FeatureFitStatus> status, bool residualDfDiffer, double confidenceLevel)
        {
            int n = status.Count;
            double[] Nan() => Enumerable.Repeat(double.NaN, n).ToArray();
            Name = name;
            Weights = weights.ToArray();
            Prior = prior;
            Status = status.ToArray();
            ResidualDfDiffer = residualDfDiffer;
            ConfidenceLevel = confidenceLevel;
            EstimateValues = Nan(); StandardErrorValues = Nan(); TValues = Nan(); DfTotalValues = Nan();
            PValues = Nan(); AdjustedValues = Nan(); PosteriorVarianceValues = Nan(); LowValues = Nan(); HighValues = Nan();
        }

        internal double[] EstimateValues { get; }
        internal double[] StandardErrorValues { get; }
        internal double[] TValues { get; }
        internal double[] DfTotalValues { get; }
        internal double[] PValues { get; }
        internal double[] AdjustedValues { get; }
        internal double[] PosteriorVarianceValues { get; }
        internal double[] LowValues { get; }
        internal double[] HighValues { get; }

        /// <summary>The contrast's name.</summary>
        public string Name { get; }
        /// <summary>The contrast's weights, one per coefficient.</summary>
        public IReadOnlyList<double> Weights { get; }
        /// <summary>The fitted prior, shared by every contrast of one call.</summary>
        public VariancePrior Prior { get; }
        /// <summary>Per-feature status carried over from the fit.</summary>
        public IReadOnlyList<FeatureFitStatus> Status { get; }
        /// <summary>
        /// True when the moderated features' variance df d are not all equal. For a linear model d is the residual df, and
        /// this is <see cref="ModeratedTest.ResidualDfDiffer"/>. For a mixed model d is each feature's joint Satterthwaite
        /// df of the condition terms (<see cref="NestedMixedModelFit.JointSatterthwaiteDf"/>), and the prior is always
        /// <see cref="VariancePriorEstimator.MomentsLegacy"/>.
        /// </summary>
        public bool ResidualDfDiffer { get; }
        /// <summary>
        /// The contrast's estimate, c′β: least squares for a linear model, the REML fit's for a mixed model. Moderation does
        /// not change it.
        /// </summary>
        public IReadOnlyList<double> Estimate => EstimateValues;
        /// <summary>
        /// Moderated standard error: sqrt(posterior variance) × sqrt(c′Uc), with U the coefficients' covariance divided by
        /// the residual variance: (XᵀX)⁻¹ for a linear model, (XᵀV⁻¹X)⁻¹ for a mixed model
        /// (<see cref="NestedMixedModelFit.UnscaledCovariance"/>).
        /// </summary>
        public IReadOnlyList<double> StandardError => StandardErrorValues;
        /// <summary>Moderated t-statistic.</summary>
        public IReadOnlyList<double> T => TValues;
        /// <summary>
        /// Degrees of freedom of the moderated t: the contrast's own df plus the prior df, capped at the sum of the moderated
        /// features' variance df d. For a linear model those are the residual df and the pooled residual df; for a mixed
        /// model, the contrast's unmoderated Satterthwaite df (<see cref="NestedMixedModelFit.SatterthwaiteDf"/>) and the sum
        /// of the features' joint Satterthwaite df.
        /// </summary>
        public IReadOnlyList<double> DfTotal => DfTotalValues;
        /// <summary>Two-sided p-value.</summary>
        public IReadOnlyList<double> PValue => PValues;
        /// <summary>Benjamini–Hochberg adjusted p-values over this contrast's tested features. Not a target-decoy q-value.</summary>
        public IReadOnlyList<double> BenjaminiHochbergAdjusted => AdjustedValues;
        /// <summary>Posterior (moderated) residual variance.</summary>
        public IReadOnlyList<double> PosteriorVariance => PosteriorVarianceValues;
        /// <summary>The two-sided level of <see cref="ConfidenceLow"/> and <see cref="ConfidenceHigh"/>, e.g. 0.95.</summary>
        public double ConfidenceLevel { get; }
        /// <summary>Lower bound: estimate − t-quantile((1 + level)/2, <see cref="DfTotal"/>) × <see cref="StandardError"/>.</summary>
        public IReadOnlyList<double> ConfidenceLow => LowValues;
        /// <summary>Upper bound: estimate + the same half-width.</summary>
        public IReadOnlyList<double> ConfidenceHigh => HighValues;
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
    /// <see cref="VariancePriorEstimator.MarginalLikelihood"/> is this library's answer to that gap: the
    /// prior is fitted by maximum marginal likelihood with each feature's own residual df, which handles
    /// unequal df exactly without porting limma's code. It is a different estimator from limma's default,
    /// so its numbers are compared with limma's, not required to equal them.
    /// </para>
    /// <para>
    /// Contrasts, with confidence intervals: <see cref="ModerateContrasts"/>. Not implemented: <c>robust = TRUE</c>,
    /// the B-statistic, observation weights.
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
        /// Fits the variance prior (d0, s0²) with the chosen <paramref name="estimator"/>: by default by matching
        /// moments of e_g = log s²_g − ψ(d_g/2) + log(d_g/2), whose variance is ψ′(d_g/2) + ψ′(d0/2), or by
        /// maximum marginal likelihood.
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
        /// <param name="estimator">
        /// <see cref="VariancePriorEstimator.MomentsLegacy"/> (default, limma legacy parity) or
        /// <see cref="VariancePriorEstimator.MarginalLikelihood"/> (unequal residual df handled exactly).
        /// </param>
        public static VariancePrior FitPrior(IReadOnlyList<double> variances, IReadOnlyList<double> df,
            IReadOnlyList<double>? covariate = null, int? splineBasisCount = null,
            VariancePriorEstimator estimator = VariancePriorEstimator.MomentsLegacy)
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
            // As limma's squeezeVar: with fewer than 3 features there is no prior to fit, so each feature
            // keeps its own variance (d0 = 0) and moderation leaves it unchanged.
            if (n < 3)
            {
                foreach (int i in use) scale[i] = variances[i];
                return new VariancePrior(0, scale, covariate != null, 1, estimator);
            }

            // A variance of exactly zero has no logarithm. Offset them away from zero relative to the median.
            double median = use.Select(i => variances[i]).Median();
            double floor = median > 0 ? VariancePrior.ZeroVarianceFloor * median : VariancePrior.ZeroVarianceFloor;
            if (estimator == VariancePriorEstimator.MarginalLikelihood)
                return FitPriorByLikelihood(variances, df, covariate, splineBasisCount, use, floor, scale);

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

        private static VariancePrior FitPriorByLikelihood(IReadOnlyList<double> variances, IReadOnlyList<double> df,
            IReadOnlyList<double>? covariate, int? splineBasisCount, int[] use, double floor, double[] scale)
        {
            double[] s2 = use.Select(i => Math.Max(variances[i], floor)).ToArray();
            double[] d = use.Select(i => df[i]).ToArray();
            int basis = 1;
            Matrix<double> b;
            if (covariate == null) b = Matrix<double>.Build.Dense(use.Length, 1, 1.0);
            else
            {
                double[] x = use.Select(i => covariate[i]).ToArray();
                b = NaturalSplineBasis(x, TrendBasisCount(use.Length, x.Distinct().Count(), splineBasisCount));
                basis = b.ColumnCount;
            }
            var (d0, logScale) = VariancePriorLikelihood.Fit(s2, d, b);
            for (int k = 0; k < use.Length; k++) scale[use[k]] = Math.Exp(logScale[k]);
            return new VariancePrior(d0, scale, covariate != null, basis, VariancePriorEstimator.MarginalLikelihood);
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
        /// <param name="estimator">How the prior is estimated; see <see cref="VariancePriorEstimator"/>.</param>
        public static ModeratedTest Moderate(LinearModelFit fit, string coefficient, bool trend,
            int? splineBasisCount = null,
            VariancePriorEstimator estimator = VariancePriorEstimator.MomentsLegacy)
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
            var prior = FitPrior(s2, df, trend ? fit.AverageResponse : null, splineBasisCount, estimator);
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

        /// <summary>
        /// Moderated tests of contrasts, with confidence intervals, as limma's <c>contrasts.fit</c> + <c>eBayes</c> +
        /// <c>topTable(confint = …)</c>. The variance prior is fitted once, exactly as <see cref="Moderate"/> fits it,
        /// and shared by every contrast. For each fitted feature and contrast c: estimate c′β; standard error
        /// s̃ · sqrt(c′(XᵀX)⁻¹c) with the feature's own (XᵀX)⁻¹ (<see cref="LinearModelFit.ContrastStdevUnscaled"/>);
        /// t on d + d0 degrees of freedom, capped at the pooled residual df; and the interval
        /// estimate ± t-quantile((1 + level)/2, df) × SE. As in limma, the cap keeps df finite even when d0 is infinite.
        /// </summary>
        /// <remarks>
        /// A contrast of one coefficient (weight 1 there, 0 elsewhere) gives exactly <see cref="Moderate"/>'s numbers.
        /// With missing values the unscaled SD is exact per feature, where limma's <c>contrasts.fit</c> approximates it;
        /// with complete data the two agree.
        /// </remarks>
        /// <param name="fit">Output of <see cref="LinearModel.Fit"/>.</param>
        /// <param name="contrasts">At least one; names unique; one finite weight per coefficient, not all zero.</param>
        /// <param name="trend">Let the prior variance trend with <see cref="LinearModelFit.AverageResponse"/>.</param>
        /// <param name="splineBasisCount">Basis functions for the trend, intercept included; null takes limma's rule.</param>
        /// <param name="estimator">How the prior is estimated; see <see cref="VariancePriorEstimator"/>.</param>
        /// <param name="confidenceLevel">Two-sided interval level, strictly between 0 and 1. Default 0.95.</param>
        public static IReadOnlyList<ModeratedContrast> ModerateContrasts(LinearModelFit fit,
            IReadOnlyList<ContrastWeights> contrasts, bool trend, int? splineBasisCount = null,
            VariancePriorEstimator estimator = VariancePriorEstimator.MomentsLegacy, double confidenceLevel = 0.95)
        {
            ArgumentNullException.ThrowIfNull(fit);
            ArgumentNullException.ThrowIfNull(contrasts);
            // Length and finiteness of the weights, with the coefficient named.
            CheckContrasts(contrasts, confidenceLevel, c => fit.ContrastEstimate(0, c.Weights));

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
            var prior = FitPrior(s2, df, trend ? fit.AverageResponse : null, splineBasisCount, estimator);
            return ModerateContrastsCore(contrasts, prior, s2, df, fitted, fit.Status, confidenceLevel,
                (f, w) => new FeatureContrast(fit.ContrastEstimate(f, w), fit.ContrastStdevUnscaled(f, w), df[f]));
        }

        /// <summary>
        /// Moderated tests of contrasts over several linear-model fits that share one variance prior: limma's
        /// <c>eBayes(legacy = TRUE, trend = FALSE)</c> over every feature of every fit, as if they were one fit. Use it when
        /// features are fitted on different designs (e.g. each on the columns its own observed samples can estimate), so
        /// they cannot share one fit, but their variances still share one prior.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b>Pooled.</b> The prior is fitted once, from every fitted feature of every part in part order, as
        /// <see cref="FitPrior"/> with <see cref="VariancePriorEstimator.MomentsLegacy"/> and no trend. The df cap is the
        /// sum of the residual df of those features, and Benjamini–Hochberg runs per contrast name over the features of
        /// every part that tests it. Each part's contrasts are its own weights on its own coefficients; parts may test
        /// different contrasts, and a part with none still informs the prior.
        /// </para>
        /// <para>
        /// <b>Few features.</b> The prior counts FITTED features only: with fewer than 3, no prior is fitted (d0 = 0), so
        /// each keeps its own variance and df, which is limma's result when it is given the fitted features alone; with
        /// none, every value is NaN. Nothing is thrown for too few features (QuantProject STAT1 ruling 5, 2026-11-10).
        /// </para>
        /// <para>
        /// <b>Same numbers as one fit.</b> With one part and at least 3 fitted features this gives exactly what
        /// <see cref="ModerateContrasts(LinearModelFit, IReadOnlyList{ContrastWeights}, bool, int?, VariancePriorEstimator, double)"/>
        /// gives with no trend; features split over several parts in their original order give the same numbers as one fit
        /// of them all.
        /// </para>
        /// </remarks>
        /// <param name="parts">At least one: each a fit and the contrasts tested on it.</param>
        /// <param name="confidenceLevel">Two-sided interval level, strictly between 0 and 1. Default 0.95.</param>
        /// <returns>Per part, in order, one result per contrast of that part, in its order.</returns>
        public static IReadOnlyList<IReadOnlyList<ModeratedContrast>> ModerateContrasts(
            IReadOnlyList<LinearModelModerationPart> parts, double confidenceLevel = 0.95)
        {
            ArgumentNullException.ThrowIfNull(parts);
            CheckLevel(confidenceLevel);
            if (parts.Count == 0) throw new ArgumentException("No part was given.", nameof(parts));
            var pooled = new List<PooledPart>(parts.Count);
            for (int i = 0; i < parts.Count; i++)
            {
                var part = parts[i] ?? throw new ArgumentException($"Part {i} is null.", nameof(parts));
                var fit = part.Fit ?? throw new ArgumentException($"Part {i} has no fit.", nameof(parts));
                var contrasts = part.Contrasts ?? throw new ArgumentException($"Part {i} has no contrast list.", nameof(parts));
                int p = fit.CoefficientNames.Count;
                if (contrasts.Count > 0)
                    CheckContrasts(contrasts, confidenceLevel, c =>
                    {
                        if (c.Weights is null || c.Weights.Count != p || c.Weights.Any(w => !double.IsFinite(w)))
                            throw new ArgumentException($"Contrast '{c.Name}' of part {i} needs {p} finite weights.", nameof(parts));
                    });
                int n = fit.FeatureCount;
                var s2 = new double[n];
                var df = new double[n];
                var usable = new List<int>();
                for (int f = 0; f < n; f++)
                {
                    s2[f] = double.NaN;
                    if (fit.Status[f] != FeatureFitStatus.Fitted) continue;
                    s2[f] = fit.Sigma[f] * fit.Sigma[f];
                    df[f] = fit.DfResidual[f];
                    usable.Add(f);
                }
                pooled.Add(new PooledPart(contrasts, s2, df, usable, fit.Status,
                    (f, w) => new FeatureContrast(fit.ContrastEstimate(f, w), fit.ContrastStdevUnscaled(f, w), df[f])));
            }
            return ModeratePooled(pooled, confidenceLevel);
        }

        /// <summary>
        /// One fit's share of a pooled moderation: its contrasts; every feature's variance and its df (read for the usable
        /// features only); the usable features, ascending; every feature's status; and the model's contrast.
        /// </summary>
        internal sealed record PooledPart(IReadOnlyList<ContrastWeights> Contrasts, double[] Variances, double[] Df,
            IReadOnlyList<int> Usable, IReadOnlyList<FeatureFitStatus> Status, Func<int, double[], FeatureContrast> ContrastOf);

        /// <summary>
        /// The moderation the pooled <c>ModerateContrasts</c> overloads share. The prior is fitted from the usable features
        /// of every part, concatenated in part order, so <see cref="FitPrior"/> counts usable features only: fewer than 3
        /// fit no prior (d0 = 0); 1 keeps its own variance; none leaves every value NaN. Then each part is moderated by
        /// <see cref="ModerateContrastsCore"/> with that prior, the df cap and the df-differ flag of all parts together, and
        /// Benjamini–Hochberg is redone per contrast name over every part.
        /// </summary>
        internal static IReadOnlyList<IReadOnlyList<ModeratedContrast>> ModeratePooled(IReadOnlyList<PooledPart> parts,
            double confidenceLevel)
        {
            var s2 = new List<double>();
            var d = new List<double>();
            foreach (var part in parts)
                foreach (int f in part.Usable)
                {
                    s2.Add(part.Variances[f]);
                    d.Add(part.Df[f]);
                }
            double d0;
            IReadOnlyList<double> scale;
            if (s2.Count >= 2)
            {
                var prior = FitPrior(s2, d, null, null, VariancePriorEstimator.MomentsLegacy);
                d0 = prior.Df;
                scale = prior.Scale;
            }
            else
            {
                d0 = 0; // 0 or 1 usable feature: no prior, as FitPrior gives for 2
                scale = s2;
            }
            double dfCap = d.Sum();
            bool dfDiffer = d.Any(x => x != d[0]);

            var results = new List<IReadOnlyList<ModeratedContrast>>(parts.Count);
            int next = 0;
            foreach (var part in parts)
            {
                var partScale = Enumerable.Repeat(double.NaN, part.Status.Count).ToArray();
                foreach (int f in part.Usable) partScale[f] = scale[next++];
                var partPrior = new VariancePrior(d0, partScale, false, 1, VariancePriorEstimator.MomentsLegacy);
                results.Add(ModerateContrastsCore(part.Contrasts, partPrior, part.Variances, part.Df, part.Usable, part.Status,
                    confidenceLevel, part.ContrastOf, dfCap, dfDiffer));
            }

            foreach (string name in parts.SelectMany(p => p.Contrasts).Select(c => c.Name).Distinct(StringComparer.Ordinal))
            {
                var members = results.SelectMany(r => r).Where(c => string.Equals(c.Name, name, StringComparison.Ordinal)).ToList();
                var adjusted = MultipleTesting.BenjaminiHochberg(members.SelectMany(c => c.PValues).ToArray());
                int offset = 0;
                foreach (var c in members)
                {
                    Array.Copy(adjusted, offset, c.AdjustedValues, 0, c.AdjustedValues.Length);
                    offset += c.AdjustedValues.Length;
                }
            }
            return results;
        }

        /// <summary>Refuses a confidence level not strictly between 0 and 1.</summary>
        internal static void CheckLevel(double confidenceLevel)
        {
            if (!(confidenceLevel > 0 && confidenceLevel < 1))
                throw new ArgumentOutOfRangeException(nameof(confidenceLevel), confidenceLevel, "The confidence level must be strictly between 0 and 1.");
        }

        /// <summary>
        /// The checks every <c>ModerateContrasts</c> makes of its request: a level strictly between 0 and 1; at least one
        /// contrast; names not blank and unique; weights that <paramref name="checkWeights"/> accepts (it throws otherwise),
        /// not all zero.
        /// </summary>
        internal static void CheckContrasts(IReadOnlyList<ContrastWeights> contrasts, double confidenceLevel,
            Action<ContrastWeights> checkWeights)
        {
            if (!(confidenceLevel > 0 && confidenceLevel < 1))
                throw new ArgumentOutOfRangeException(nameof(confidenceLevel), confidenceLevel, "The confidence level must be strictly between 0 and 1.");
            if (contrasts.Count == 0) throw new ArgumentException("No contrast was given.", nameof(contrasts));
            var seen = new HashSet<string>(StringComparer.Ordinal);
            foreach (var c in contrasts)
            {
                ArgumentNullException.ThrowIfNull(c);
                ArgumentException.ThrowIfNullOrWhiteSpace(c.Name, nameof(contrasts));
                if (!seen.Add(c.Name)) throw new ArgumentException($"Two contrasts are named '{c.Name}'.", nameof(contrasts));
                checkWeights(c);
                if (c.Weights.All(w => w == 0))
                    throw new ArgumentException($"Every weight of contrast '{c.Name}' is zero.", nameof(contrasts));
            }
        }

        /// <summary>
        /// The moderation every <c>ModerateContrasts</c> shares, once the caller has checked the request and fitted the prior
        /// from (s², d). Per contrast c, for each feature in <paramref name="features"/>: the model's
        /// <see cref="FeatureContrast"/> (none when its df is not finite); the posterior variance
        /// s̃² = (d0·s0² + d·s²) / (d0 + d), or s0² when d0 is infinite; SE = s̃ · sqrt(c′Uc); t on the contrast's df + d0,
        /// capped at the sum of d over <paramref name="features"/>; the two-sided p; the interval
        /// estimate ± t-quantile((1 + level)/2, df) × SE; and Benjamini–Hochberg over the features with a result.
        /// </summary>
        /// <param name="contrasts">Already checked (<see cref="CheckContrasts"/>).</param>
        /// <param name="prior">Fitted from <paramref name="variances"/> and <paramref name="df"/>.</param>
        /// <param name="variances">s² per feature; read only for <paramref name="features"/>.</param>
        /// <param name="df">d, the degrees of freedom of s², per feature; read only for <paramref name="features"/>.</param>
        /// <param name="features">The features moderated, in ascending order; at least one.</param>
        /// <param name="status">Every feature's fit status; its length is the feature count.</param>
        /// <param name="confidenceLevel">Already checked.</param>
        /// <param name="contrastOf">The model's estimate, unscaled SD and df of a contrast (its weights) for a feature.</param>
        /// <param name="dfCapOverride">The df cap when the prior was fitted over more features than these (a pooled
        /// moderation, <see cref="ModeratePooled"/>); null caps at the sum of d over <paramref name="features"/>.</param>
        /// <param name="dfDifferOverride">Whether the d of every moderated feature differ, for a pooled moderation; null
        /// judges <paramref name="features"/> alone.</param>
        internal static IReadOnlyList<ModeratedContrast> ModerateContrastsCore(IReadOnlyList<ContrastWeights> contrasts,
            VariancePrior prior, IReadOnlyList<double> variances, IReadOnlyList<double> df, IReadOnlyList<int> features,
            IReadOnlyList<FeatureFitStatus> status, double confidenceLevel, Func<int, double[], FeatureContrast> contrastOf,
            double? dfCapOverride = null, bool? dfDifferOverride = null)
        {
            int n = status.Count;
            double dfCap = dfCapOverride ?? features.Sum(f => df[f]);
            bool dfDiffer = dfDifferOverride ?? features.Any(f => df[f] != df[features[0]]);
            double upper = (1 + confidenceLevel) / 2;

            var results = new List<ModeratedContrast>(contrasts.Count);
            foreach (var c in contrasts)
            {
                var w = c.Weights.ToArray();
                var result = new ModeratedContrast(c.Name, w, prior, status, dfDiffer, confidenceLevel);
                foreach (int f in features)
                {
                    var feature = contrastOf(f, w);
                    if (!double.IsFinite(feature.Df)) continue;
                    double post = double.IsPositiveInfinity(prior.Df)
                        ? prior.Scale[f]
                        : (prior.Df * prior.Scale[f] + df[f] * variances[f]) / (prior.Df + df[f]);
                    double dfTotal = Math.Min(feature.Df + prior.Df, dfCap);
                    double se = Math.Sqrt(post) * feature.StdevUnscaled;
                    double t = feature.Estimate / se;
                    double half = Quantile(upper, dfTotal) * se;
                    result.EstimateValues[f] = feature.Estimate;
                    result.PosteriorVarianceValues[f] = post;
                    result.StandardErrorValues[f] = se;
                    result.TValues[f] = t;
                    result.DfTotalValues[f] = dfTotal;
                    result.PValues[f] = TwoSidedP(t, dfTotal);
                    result.LowValues[f] = feature.Estimate - half;
                    result.HighValues[f] = feature.Estimate + half;
                }
                var adjusted = MultipleTesting.BenjaminiHochberg(result.PValues);
                Array.Copy(adjusted, result.AdjustedValues, n);
                results.Add(result);
            }
            return results;
        }

        /// <summary>The <paramref name="probability"/> quantile of Student's t on <paramref name="df"/> degrees of freedom; normal when df is infinite.</summary>
        internal static double Quantile(double probability, double df) =>
            double.IsPositiveInfinity(df) ? Normal.InvCDF(0, 1, probability) : StudentT.InvCDF(0, 1, df, probability);
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
