using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using NUnit.Framework;
using StatisticalModels;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// Closed-form and simulation checks of the per-feature statistics engine. Reference values from
/// another implementation belong in a separate, fixture-driven test; these establish that each piece
/// is correct against mathematics that can be written down.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialAbundanceTests
{
    // ---- Polygamma -------------------------------------------------------------------------------

    [Test]
    public void Trigamma_MatchesClosedForms()
    {
        Assert.That(Polygamma.Trigamma(1), Is.EqualTo(Math.PI * Math.PI / 6).Within(1e-13));
        Assert.That(Polygamma.Trigamma(0.5), Is.EqualTo(Math.PI * Math.PI / 2).Within(1e-12));
        Assert.That(Polygamma.Trigamma(2), Is.EqualTo(Math.PI * Math.PI / 6 - 1).Within(1e-13));
        Assert.That(Polygamma.Trigamma(1000), Is.EqualTo(1.0 / 1000 + 1.0 / (2 * 1e6) + 1.0 / (6 * 1e9)).Within(1e-15));
        // Between 10 and 20 the series alone was off by up to 2e-13 relative; the recurrence now runs to 20.
        // ψ′(10) = π²/6 − Σ_{k=1..9} 1/k².
        double tri10 = Math.PI * Math.PI / 6 - Enumerable.Range(1, 9).Sum(k => 1.0 / ((double)k * k));
        Assert.That(Polygamma.Trigamma(10), Is.EqualTo(tri10).Within(1e-14 * tri10));
    }

    [Test]
    public void Tetragamma_AtOne_IsMinusTwiceZetaThree() =>
        Assert.That(Polygamma.Tetragamma(1), Is.EqualTo(-2 * 1.2020569031595942).Within(1e-12));

    [TestCase(1e-5)]
    [TestCase(0.01)]
    [TestCase(0.3)]
    [TestCase(1)]
    [TestCase(7)]
    [TestCase(1e4)]
    [TestCase(1e8)]
    public void TrigammaInverse_RoundTrips(double x)
    {
        double y = Polygamma.TrigammaInverse(x);
        Assert.That(Polygamma.Trigamma(y), Is.EqualTo(x).Within(1e-8 * x));
    }

    [Test]
    public void Polygamma_RejectsNonPositiveArguments()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Polygamma.Trigamma(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Polygamma.TrigammaInverse(-1));
    }

    // ---- LinearModel -----------------------------------------------------------------------------

    private static double[,] InterceptAndSlope(double[] x) =>
        ToMatrix(x.Length, 2, (i, j) => j == 0 ? 1 : x[i]);

    private static double[,] ToMatrix(int rows, int cols, Func<int, int, double> f)
    {
        var m = new double[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i, j] = f(i, j);
        return m;
    }

    [Test]
    public void LinearModel_RecoversSimpleRegressionClosedForm()
    {
        double[] x = { 1, 2, 3, 4, 5, 6 };
        double[] y = { 1.1, 1.9, 3.2, 3.8, 5.1, 6.3 };
        var fit = LinearModel.Fit(ToMatrix(1, 6, (_, s) => y[s]), InterceptAndSlope(x), new[] { "intercept", "x" });

        double xm = x.Average(), ym = y.Average();
        double sxx = x.Sum(v => (v - xm) * (v - xm));
        double slope = x.Zip(y, (a, b) => (a - xm) * (b - ym)).Sum() / sxx;
        double intercept = ym - slope * xm;
        double rss = x.Zip(y, (a, b) => Math.Pow(b - intercept - slope * a, 2)).Sum();

        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Coefficient(0, 1), Is.EqualTo(slope).Within(1e-12));
        Assert.That(fit.Coefficient(0, 0), Is.EqualTo(intercept).Within(1e-12));
        Assert.That(fit.DfResidual[0], Is.EqualTo(4));
        Assert.That(fit.Sigma[0], Is.EqualTo(Math.Sqrt(rss / 4)).Within(1e-12));
        Assert.That(fit.StdevUnscaled(0, 1), Is.EqualTo(1 / Math.Sqrt(sxx)).Within(1e-12));
        Assert.That(fit.AverageResponse[0], Is.EqualTo(ym).Within(1e-12));
    }

    /// <summary>A missing value is omitted, so the fit equals the fit on the observed samples alone.</summary>
    [Test]
    public void LinearModel_OmitsMissingValuesPerFeature()
    {
        double[] x = { 1, 2, 3, 4, 5, 6 };
        var withGap = ToMatrix(1, 6, (_, s) => s == 2 ? double.NaN : 0.3 * x[s] + (s % 2 == 0 ? 0.1 : -0.1));
        var full = LinearModel.Fit(withGap, InterceptAndSlope(x));

        double[] xs = x.Where((_, s) => s != 2).ToArray();
        var subset = LinearModel.Fit(ToMatrix(1, 5, (_, s) => withGap[0, s < 2 ? s : s + 1]), InterceptAndSlope(xs));

        Assert.That(full.Observed[0], Is.EqualTo(5));
        Assert.That(full.Coefficient(0, 1), Is.EqualTo(subset.Coefficient(0, 1)).Within(1e-14));
        Assert.That(full.Sigma[0], Is.EqualTo(subset.Sigma[0]).Within(1e-14));
    }

    [Test]
    public void LinearModel_ReportsWhyAFeatureCannotBeFitted()
    {
        // Columns: intercept, age, male. Samples 0-3 are female, 4-7 male.
        var design = ToMatrix(8, 3, (i, j) => j == 0 ? 1 : j == 1 ? 30 + 5 * i : (i >= 4 ? 1 : 0));
        var y = ToMatrix(2, 8, (f, s) =>
            f == 0 ? (s < 3 ? 1.0 * s : double.NaN)                         // 3 observations, 3 coefficients: no residual df
                   : (s >= 4 ? 0.5 * s + (s % 2) * 0.1 : double.NaN));    // 4 observed, all male: sex not estimable
        var fit = LinearModel.Fit(y, design, new[] { "intercept", "age", "male" });

        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.TooFewObservations));
        Assert.That(fit.Status[1], Is.EqualTo(FeatureFitStatus.RankDeficient));
        Assert.That(double.IsNaN(fit.Coefficient(0, 1)) && double.IsNaN(fit.Coefficient(1, 1)));
    }

    [Test]
    public void LinearModel_RejectsAMisuedDesign()
    {
        var y = ToMatrix(1, 4, (_, s) => s);
        Assert.Throws<ArgumentException>(() => LinearModel.Fit(y, ToMatrix(3, 1, (_, _) => 1)));                    // wrong rows
        Assert.Throws<ArgumentException>(() => LinearModel.Fit(y, ToMatrix(4, 2, (_, _) => 1)));                    // redundant column
        Assert.Throws<ArgumentException>(() => LinearModel.Fit(y, ToMatrix(4, 1, (i, _) => i == 0 ? double.NaN : 1)));
    }

    [Test]
    public void LinearModel_RankTestIsIndependentOfColumnUnits()
    {
        // An intercept beside a covariate on a raw-intensity scale (about 1e8, e.g. total ion current). R[0,0] is
        // sqrt(12) = 3.5 while the covariate's diagonal element is about 1.2e8, so a tolerance relative to the
        // largest diagonal element (1e-7 * 1.2e8 = 12) called the intercept redundant and refused the design.
        double[] tic = Enumerable.Range(0, 12).Select(i => 1e8 * (1 + 0.1 * i)).ToArray();
        var y = ToMatrix(1, 12, (_, s) => 2.0 + 3e-8 * tic[s] + (s % 2 == 0 ? 0.05 : -0.05));
        var raw = LinearModel.Fit(y, InterceptAndSlope(tic));
        var rescaled = LinearModel.Fit(y, InterceptAndSlope(tic.Select(t => t / 1e8).ToArray()));

        Assert.That(raw.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(raw.Coefficient(0, 0), Is.EqualTo(rescaled.Coefficient(0, 0)).Within(1e-9));
        Assert.That(raw.Coefficient(0, 1) * 1e8, Is.EqualTo(rescaled.Coefficient(0, 1)).Within(1e-9));
        // A genuinely redundant pair on the same scale is still refused.
        Assert.Throws<ArgumentException>(() => LinearModel.Fit(y, ToMatrix(12, 2, (i, j) => j == 0 ? tic[i] : 2 * tic[i])));
    }

    [Test]
    public void LinearModel_OutputDoesNotDependOnThreadCount()
    {
        var (y, design) = SimulatedAgeStudy(features: 500, samples: 12, trueSlope: 0.2, seed: 7);
        var one = LinearModel.Fit(y, design, maxThreads: 1);
        var many = LinearModel.Fit(y, design, maxThreads: 8);
        Assert.That(many.CoefficientMatrix, Is.EqualTo(one.CoefficientMatrix));
        Assert.That(many.Sigma, Is.EqualTo(one.Sigma));
    }

    /// <summary>-1 means "all cores but one", as in FlashLFQ, not "one thread".</summary>
    [Test]
    public void MaxThreads_MinusOneMeansAllCoresButOne()
    {
        Assert.That(LinearModel.ResolveThreads(-1), Is.EqualTo(Math.Max(1, Environment.ProcessorCount - 1)));
        Assert.That(LinearModel.ResolveThreads(0), Is.EqualTo(1));
        Assert.That(LinearModel.ResolveThreads(1), Is.EqualTo(1));
    }

    /// <summary>A fit with too few usable features says so about the fit, not about an internal argument.</summary>
    [Test]
    public void Moderate_RefusesAFitWithFewerThanTwoFittedFeatures()
    {
        var design = InterceptAndSlope(new double[] { 1, 2, 3, 4 });
        var y = ToMatrix(2, 4, (f, s) => f == 0 ? s + 0.1 * (s % 2) : (s == 0 ? 1 : double.NaN));
        var fit = LinearModel.Fit(y, design, new[] { "intercept", "x" });
        var ex = Assert.Throws<ArgumentException>(() => EmpiricalBayes.Moderate(fit, "x", trend: false));
        Assert.That(ex!.ParamName, Is.EqualTo("fit"));
    }

    // ---- MultipleTesting ------------------------------------------------------------------------

    [Test]
    public void BenjaminiHochberg_WorkedExample()
    {
        // Sorted: 0.005 -> 0.02, 0.01 -> 0.02, 0.03 -> 0.04, 0.04 -> 0.04, after the running minimum.
        var q = MultipleTesting.BenjaminiHochberg(new[] { 0.01, 0.04, 0.03, 0.005 });
        Assert.That(q, Is.EqualTo(new[] { 0.02, 0.04, 0.04, 0.02 }).Within(1e-15));
    }

    [Test]
    public void BenjaminiHochberg_LeavesUntestedFeaturesOutOfTheFamily()
    {
        var q = MultipleTesting.BenjaminiHochberg(new[] { 0.01, double.NaN, 0.02 });
        Assert.That(double.IsNaN(q[1]));
        Assert.That(q[0], Is.EqualTo(0.02).Within(1e-15));   // m = 2, not 3
        Assert.That(q[2], Is.EqualTo(0.02).Within(1e-15));
    }

    /// <summary>
    /// The two worked examples above are already monotone, so they would also pass a naive p·m/rank. Here the
    /// step-up must act. Expected values are R's p.adjust(p, "BH").
    /// </summary>
    [Test]
    public void BenjaminiHochberg_StepUpLowersAnOutOfOrderValue()
    {
        // Raw p·m/rank: 0.9, 0.08, 0.004, 0.0547 (sorted ranks 4, 2, 1, 3). The running minimum from the top
        // lowers rank 2's 0.08 to rank 3's 0.0547.
        var q = MultipleTesting.BenjaminiHochberg(new[] { 0.9, 0.04, 0.001, 0.041 });
        Assert.That(q[0], Is.EqualTo(0.9).Within(1e-15));
        Assert.That(q[1], Is.EqualTo(0.041 * 4 / 3).Within(1e-15));
        Assert.That(q[2], Is.EqualTo(0.004).Within(1e-15));
        Assert.That(q[3], Is.EqualTo(0.041 * 4 / 3).Within(1e-15));
    }

    [Test]
    public void BenjaminiHochberg_RawValueAboveOneIsBroughtDown()
    {
        // Rank 1's raw p·m/rank is 0.6 · 2 = 1.2; the step-up returns rank 2's 0.61 for both, never above 1.
        var q = MultipleTesting.BenjaminiHochberg(new[] { 0.61, 0.6 });
        Assert.That(q[0], Is.EqualTo(0.61).Within(1e-15));
        Assert.That(q[1], Is.EqualTo(0.61).Within(1e-15));
        Assert.That(q.All(v => v <= 1));
    }

    // ---- RandomEffectsMeta ----------------------------------------------------------------------

    [Test]
    public void RandomEffectsMeta_HeterogeneousWorkedExample()
    {
        // Equal weights: fixed effect 2, Q = 8, C = 2, tau2 = (8 - 2) / 2 = 3, I2 = 0.75.
        var r = RandomEffectsMeta.Pool(new[] { 0.0, 2.0, 4.0 }, new[] { 1.0, 1.0, 1.0 });
        Assert.That(r.Estimate, Is.EqualTo(2).Within(1e-14));
        Assert.That(r.Tau2, Is.EqualTo(3).Within(1e-14));
        Assert.That(r.Q, Is.EqualTo(8).Within(1e-14));
        Assert.That(r.ISquared, Is.EqualTo(0.75).Within(1e-14));
        Assert.That(r.StandardError, Is.EqualTo(Math.Sqrt(4.0 / 3)).Within(1e-14));
        Assert.That(r.DirectionAgree, Is.EqualTo(2));
        Assert.That(r.LeaveOneOutMaxDelta, Is.EqualTo(1).Within(1e-14));
    }

    [Test]
    public void RandomEffectsMeta_HomogeneousStudiesHaveNoBetweenStudyVariance()
    {
        var r = RandomEffectsMeta.Pool(new[] { 1.0, 2.0 }, new[] { 1.0, 1.0 });
        Assert.That(r.Tau2, Is.EqualTo(0));
        Assert.That(r.Estimate, Is.EqualTo(1.5).Within(1e-14));
        Assert.That(r.ISquared, Is.EqualTo(0));
    }

    [Test]
    public void RandomEffectsMeta_OneStudyIsAnHonestRowNotAMetaAnalysis()
    {
        var r = RandomEffectsMeta.Pool(new[] { 0.4 }, new[] { 0.1 });
        Assert.That(r.Estimate, Is.EqualTo(0.4).Within(1e-15));
        Assert.That(r.StandardError, Is.EqualTo(0.1).Within(1e-15));
        Assert.That(double.IsNaN(r.ISquared) && double.IsNaN(r.LeaveOneOutMaxDelta));
    }

    // ---- EmpiricalBayes -------------------------------------------------------------------------

    /// <summary>
    /// Simulates the model moderation assumes: σ²_g drawn from a scaled inverse χ² with d0 and s0², then
    /// s²_g ~ σ²_g·χ²(d)/d. With many features the moment estimator must recover d0 and s0².
    /// </summary>
    [Test]
    public void FitPrior_RecoversTheHyperparametersItAssumes()
    {
        const double d0 = 6, s0Sq = 0.05, d = 5;
        var rng = new MersenneTwister(11);
        var chiD0 = new ChiSquared(d0, rng);
        var chiD = new ChiSquared(d, rng);
        int n = 40000;
        var s2 = new double[n];
        for (int g = 0; g < n; g++)
        {
            double sigma2 = d0 * s0Sq / chiD0.Sample();
            s2[g] = sigma2 * chiD.Sample() / d;
        }
        var prior = EmpiricalBayes.FitPrior(s2, Enumerable.Repeat(d, n).ToArray());

        Assert.That(prior.Df, Is.EqualTo(d0).Within(0.1 * d0));
        Assert.That(prior.Scale[0], Is.EqualTo(s0Sq).Within(0.05 * s0Sq));
        Assert.That(prior.Trended, Is.False);
    }

    /// <summary>When every feature shares one true variance, sample variances are no more spread than χ²
    /// sampling explains, the prior df is infinite, and every feature takes the prior variance.</summary>
    [Test]
    public void FitPrior_CommonVarianceGivesInfinitePriorDf()
    {
        var rng = new MersenneTwister(3);
        var chi = new ChiSquared(4, rng);
        var s2 = Enumerable.Range(0, 20000).Select(_ => 0.2 * chi.Sample() / 4).ToArray();
        var prior = EmpiricalBayes.FitPrior(s2, Enumerable.Repeat(4.0, s2.Length).ToArray());
        Assert.That(prior.Df, Is.GreaterThan(200));   // infinite or very large: no detectable dispersion
        Assert.That(prior.Scale[0], Is.EqualTo(0.2).Within(0.01));
    }

    /// <summary>A prior variance that falls with intensity is captured by the trend and missed without it.</summary>
    [Test]
    public void FitPrior_TrendFollowsIntensity()
    {
        var rng = new MersenneTwister(5);
        var chi = new ChiSquared(6, rng);
        int n = 20000;
        var amean = Enumerable.Range(0, n).Select(g => 18 + 12.0 * g / n).ToArray();
        var s2 = amean.Select(a => Math.Exp(-0.3 * (a - 24)) * 0.05 * chi.Sample() / 6).ToArray();
        var df = Enumerable.Repeat(6.0, n).ToArray();

        var trended = EmpiricalBayes.FitPrior(s2, df, amean);
        double low = trended.Scale[n / 10], high = trended.Scale[9 * n / 10];
        double expectedRatio = Math.Exp(-0.3 * (amean[n / 10] - amean[9 * n / 10]));
        Assert.That(trended.Trended);
        Assert.That(low / high, Is.EqualTo(expectedRatio).Within(0.1 * expectedRatio));
    }

    /// <summary>
    /// limma's fitFDist sizes the default trend spline by the number of usable features,
    /// 1 + (n ≥ 3) + (n ≥ 6) + (n ≥ 30), capped at the distinct covariate values; a size the caller
    /// asks for is capped to leave one residual degree of freedom.
    /// </summary>
    [TestCase(2, 2, null, 1)]
    [TestCase(3, 3, null, 2)]
    [TestCase(5, 5, null, 2)]
    [TestCase(6, 6, null, 3)]
    [TestCase(29, 29, null, 3)]
    [TestCase(30, 30, null, 4)]
    [TestCase(5000, 5000, null, 4)]
    [TestCase(30, 2, null, 2)]
    [TestCase(30, 1, null, 1)]
    [TestCase(20, 20, 4, 4)]
    [TestCase(3, 3, 4, 2)]
    [TestCase(100, 3, 6, 3)]
    public void TrendBasisCount_FollowsLimmaWhenNotGiven(int usable, int distinct, int? requested, int expected)
    {
        Assert.That(EmpiricalBayes.TrendBasisCount(usable, distinct, requested), Is.EqualTo(expected));
    }

    /// <summary>
    /// Tied covariate values can put two knots on one value, or an interior knot on the upper boundary.
    /// Either used to give 0/0 or a rank-deficient design and NaN everywhere; the knots are now merged.
    /// </summary>
    [TestCase(0.9, 30.0, VariancePriorEstimator.MomentsLegacy)]    // tied at the top: interior knots land on the upper boundary
    [TestCase(0.6, 24.0, VariancePriorEstimator.MomentsLegacy)]    // tied in the middle: both interior knots coincide
    [TestCase(0.45, 18.0, VariancePriorEstimator.MomentsLegacy)]   // tied at the bottom: an interior knot lands on the lower boundary
    [TestCase(0.9, 30.0, VariancePriorEstimator.MarginalLikelihood)]
    [TestCase(0.6, 24.0, VariancePriorEstimator.MarginalLikelihood)]
    [TestCase(0.45, 18.0, VariancePriorEstimator.MarginalLikelihood)]
    public void FitPrior_TiedCovariatesMergeKnotsInsteadOfGivingNaN(double tiedFraction, double tiedValue,
        VariancePriorEstimator estimator)
    {
        var rng = new MersenneTwister(23);
        var chi = new ChiSquared(5, rng);
        int n = 200, tied = (int)(tiedFraction * n);
        var amean = Enumerable.Range(0, n).Select(g => g < tied ? tiedValue : 18 + 12.0 * rng.NextDouble()).ToArray();
        var s2 = amean.Select(_ => 0.05 * chi.Sample() / 5).ToArray();
        var prior = EmpiricalBayes.FitPrior(s2, Enumerable.Repeat(5.0, n).ToArray(), amean, estimator: estimator);

        Assert.That(prior.Scale.All(double.IsFinite), "every feature gets a finite prior variance");
        Assert.That(double.IsNaN(prior.Df), Is.False);
        Assert.That(prior.SplineBasisCount, Is.InRange(2, EmpiricalBayes.DefaultSplineBasisCount));
    }

    [Test]
    public void NaturalSplineBasis_MergesCoincidingKnotsAndStaysFullRank()
    {
        var x = Enumerable.Repeat(5.0, 8).Concat(new[] { 6.0, 7.0 }).ToArray();   // quantile knots all at 5
        var basis = EmpiricalBayes.NaturalSplineBasis(x, 4);
        Assert.That(basis.ColumnCount, Is.LessThan(4));
        Assert.That(basis.Rank(), Is.EqualTo(basis.ColumnCount));
        Assert.That(basis.Enumerate().All(double.IsFinite));

        Assert.That(EmpiricalBayes.NaturalSplineBasis(Enumerable.Repeat(3.0, 6).ToArray(), 1).ColumnCount, Is.EqualTo(1));
    }

    /// <summary>limma refuses a missing covariate; silently dropping the feature hid it.</summary>
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void FitPrior_RefusesANonFiniteCovariateOnAUsableFeature(double bad)
    {
        var s2 = new[] { 0.1, 0.2, 0.15, 0.3 };
        var df = new[] { 4.0, 4, 4, 4 };
        var ex = Assert.Throws<ArgumentException>(() => EmpiricalBayes.FitPrior(s2, df, new[] { 20, bad, 22, 23 }));
        Assert.That(ex!.Message, Does.Contain("feature 1"));

        // A feature that is not usable anyway (no df) may carry any covariate.
        Assert.DoesNotThrow(() => EmpiricalBayes.FitPrior(s2, new[] { 4.0, 0, 4, 4 }, new[] { 20, bad, 22, 23 }));
    }

    [Test]
    public void VariancePrior_ScaleCannotBeCastBackAndChanged()
    {
        var prior = EmpiricalBayes.FitPrior(new[] { 0.1, 0.2, 0.15, 0.3 }, new[] { 4.0, 4, 4, 4 });
        Assert.That(prior.Scale, Is.Not.InstanceOf<double[]>());
        Assert.Throws<NotSupportedException>(() => ((System.Collections.Generic.IList<double>)prior.Scale)[0] = 1);
    }

    /// <summary>Under the null, moderated p-values are uniform: about 5% fall below 0.05.</summary>
    [Test]
    public void Moderate_NullPValuesAreUniform()
    {
        var (y, design) = SimulatedAgeStudy(features: 6000, samples: 10, trueSlope: 0, seed: 13);
        var fit = LinearModel.Fit(y, design, new[] { "intercept", "age_decades", "male" }, maxThreads: 4);
        var test = EmpiricalBayes.Moderate(fit, "age_decades", trend: true);

        var p = test.PValue.Where(double.IsFinite).ToArray();
        double below = p.Count(v => v < 0.05) / (double)p.Length;
        Assert.That(p.Length, Is.GreaterThan(5500));
        Assert.That(below, Is.EqualTo(0.05).Within(0.01));
        Assert.That(test.BenjaminiHochbergAdjusted.Where(double.IsFinite).Count(q => q < 0.05), Is.LessThan(30));
    }

    /// <summary>A real age effect on a tenth of features is found, and moderation does not move the estimate.</summary>
    [Test]
    public void Moderate_FindsTrueAgeEffectsAndLeavesEstimatesAlone()
    {
        var (y, design) = SimulatedAgeStudy(features: 3000, samples: 12, trueSlope: 0.6, seed: 17, affectedEvery: 10);
        var fit = LinearModel.Fit(y, design, new[] { "intercept", "age_decades", "male" });
        var test = EmpiricalBayes.Moderate(fit, "age_decades", trend: false);

        int truePositives = Enumerable.Range(0, 3000).Count(f => f % 10 == 0 && test.BenjaminiHochbergAdjusted[f] < 0.05);
        int falsePositives = Enumerable.Range(0, 3000).Count(f => f % 10 != 0 && test.BenjaminiHochbergAdjusted[f] < 0.05);
        Assert.That(truePositives, Is.GreaterThan(150));
        Assert.That(falsePositives / (double)Math.Max(1, truePositives + falsePositives), Is.LessThan(0.1));
        Assert.That(test.Estimate[0], Is.EqualTo(fit.Coefficient(0, 1)));
        Assert.That(test.DfTotal[0], Is.GreaterThan(fit.DfResidual[0]));
    }

    /// <summary>
    /// Omitting a missing value lowers that feature's residual df, which is when this result departs from
    /// limma's current default, so the result must say so. Complete data leave every df equal. A feature
    /// that could not be fitted has no residual df and must not raise the flag.
    /// </summary>
    [Test]
    public void Moderate_ReportsWhenResidualDfDiffer()
    {
        var design = InterceptAndSlope(new double[] { 1, 2, 3, 4, 5, 6 });
        double[,] Responses(bool oneMissing) => ToMatrix(5, 6, (f, s) =>
            f == 4 ? (s == 0 ? 1 : double.NaN)                  // too sparse to fit
            : oneMissing && f == 0 && s == 5 ? double.NaN
            : 10 + f + 0.3 * s + 0.1 * ((f + 1) * (s + 2) % 5));
        var names = new[] { "intercept", "x" };

        var complete = EmpiricalBayes.Moderate(LinearModel.Fit(Responses(false), design, names), "x", trend: false);
        Assert.That(complete.ResidualDfDiffer, Is.False);

        var withMissing = EmpiricalBayes.Moderate(LinearModel.Fit(Responses(true), design, names), "x", trend: false);
        Assert.That(withMissing.ResidualDfDiffer, Is.True);
    }

    /// <summary>
    /// log2 abundances for a cross-sectional cohort: ages 30-75 spread across the samples, sex alternating,
    /// per-feature baseline and variance drawn from a scaled inverse χ², about 10% missing at random.
    /// Features with index divisible by <paramref name="affectedEvery"/> change by
    /// <paramref name="trueSlope"/> log2 units per decade.
    /// </summary>
    private static (double[,] y, double[,] design) SimulatedAgeStudy(int features, int samples, double trueSlope,
        int seed, int affectedEvery = 1)
    {
        var rng = new MersenneTwister(seed);
        var normal = new Normal(0, 1, rng);
        var chi = new ChiSquared(5, rng);
        var ageDecades = Enumerable.Range(0, samples).Select(s => (30 + 45.0 * s / (samples - 1) - 50) / 10).ToArray();
        var design = ToMatrix(samples, 3, (i, j) => j == 0 ? 1 : j == 1 ? ageDecades[i] : i % 2);
        var y = new double[features, samples];
        for (int f = 0; f < features; f++)
        {
            double baseline = 20 + 4 * normal.Sample();
            double sd = Math.Sqrt(5 * 0.04 / chi.Sample());
            double slope = f % affectedEvery == 0 ? trueSlope : 0;
            for (int s = 0; s < samples; s++)
                y[f, s] = rng.NextDouble() < 0.1
                    ? double.NaN
                    : baseline + slope * ageDecades[s] + 0.2 * (s % 2) + sd * normal.Sample();
        }
        return (y, design);
    }
}
