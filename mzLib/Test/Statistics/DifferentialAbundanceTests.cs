using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using NUnit.Framework;
using Statistics;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Statistics;

/// <summary>
/// Closed-form and simulation checks of the per-feature statistics engine. Reference values from
/// another implementation belong in a separate, fixture-driven test; these establish that each piece
/// is correct against mathematics that can be written down.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialAbundanceTests
{
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
