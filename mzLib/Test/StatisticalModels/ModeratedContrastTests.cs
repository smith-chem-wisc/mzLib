using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using NUnit.Framework;
using StatisticalModels;

namespace Test.StatisticalModels;

/// <summary>
/// What moderated contrasts promise beyond the limma reference values: a unit contrast reproduces
/// <see cref="EmpiricalBayes.Moderate"/> exactly under both prior estimators; the stored unscaled covariance agrees
/// with <see cref="LinearModelFit.StdevUnscaled"/>; the interval is symmetric about the estimate at the requested
/// level, on the pooled residual df when the prior df is infinite; non-fitted features stay NaN; misuse throws.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ModeratedContrastTests
{
    private static readonly string[] Names = { "intercept", "groupB", "age" };

    /// <summary>40 features, 9 samples, intercept + group + age; a few values missing; one feature too sparse to fit.</summary>
    private static LinearModelFit Fit(int seed = 20261008)
    {
        var rng = new MersenneTwister(seed);
        int g = 40, n = 9;
        var design = new double[n, 3];
        for (int i = 0; i < n; i++) { design[i, 0] = 1; design[i, 1] = i % 2; design[i, 2] = (20 + 7 * i - 50) / 10.0; }
        var y = new double[g, n];
        for (int f = 0; f < g; f++)
        {
            double sd = 0.2 + 0.4 * rng.NextDouble();
            for (int i = 0; i < n; i++)
                y[f, i] = 20 + 0.5 * design[i, 1] - 0.3 * design[i, 2] + sd * Normal.Sample(rng, 0, 1);
        }
        y[3, 1] = double.NaN; y[7, 4] = double.NaN; y[7, 5] = double.NaN; y[11, 8] = double.NaN;
        for (int i = 0; i < n - 2; i++) y[39, i] = double.NaN; // too few observations
        return LinearModel.Fit(y, design, Names);
    }

    [TestCase(VariancePriorEstimator.MomentsLegacy, false)]
    [TestCase(VariancePriorEstimator.MomentsLegacy, true)]
    [TestCase(VariancePriorEstimator.MarginalLikelihood, false)]
    [TestCase(VariancePriorEstimator.MarginalLikelihood, true)]
    public void AUnitContrastReproducesModerate(VariancePriorEstimator estimator, bool trend)
    {
        var fit = Fit();
        var single = EmpiricalBayes.Moderate(fit, "groupB", trend, estimator: estimator);
        var contrast = EmpiricalBayes.ModerateContrasts(fit, new[] { new ContrastWeights("groupB", new[] { 0.0, 1, 0 }) },
            trend, estimator: estimator).Single();
        Assert.That(contrast.Prior.Df, Is.EqualTo(single.Prior.Df));
        Assert.That(contrast.Estimate, Is.EqualTo(single.Estimate));
        Assert.That(contrast.StandardError, Is.EqualTo(single.StandardError));
        Assert.That(contrast.T, Is.EqualTo(single.T));
        Assert.That(contrast.DfTotal, Is.EqualTo(single.DfTotal));
        Assert.That(contrast.PValue, Is.EqualTo(single.PValue));
        Assert.That(contrast.BenjaminiHochbergAdjusted, Is.EqualTo(single.BenjaminiHochbergAdjusted));
        Assert.That(contrast.PosteriorVariance, Is.EqualTo(single.PosteriorVariance));
        Assert.That(contrast.ResidualDfDiffer, Is.EqualTo(single.ResidualDfDiffer));
        Assert.That(contrast.Status, Is.EqualTo(single.Status));
    }

    [Test]
    public void UnscaledCovarianceIsSymmetricAndItsDiagonalIsStdevUnscaledSquared()
    {
        var fit = Fit();
        for (int f = 0; f < 39; f++)
            for (int i = 0; i < 3; i++)
            {
                Assert.That(fit.UnscaledCovariance(f, i, i), Is.EqualTo(fit.StdevUnscaled(f, i) * fit.StdevUnscaled(f, i)).Within(1e-12).Percent);
                for (int j = 0; j < 3; j++)
                    Assert.That(fit.UnscaledCovariance(f, i, j), Is.EqualTo(fit.UnscaledCovariance(f, j, i)));
            }
    }

    [Test]
    public void ContrastStdevUnscaledIsTheQuadraticForm()
    {
        var fit = Fit();
        double[] c = { 0, 1, -2 };
        for (int f = 0; f < 39; f++)
        {
            double q = 0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) q += c[i] * c[j] * fit.UnscaledCovariance(f, i, j);
            Assert.That(fit.ContrastStdevUnscaled(f, c), Is.EqualTo(Math.Sqrt(q)).Within(1e-12));
            Assert.That(fit.ContrastEstimate(f, c), Is.EqualTo(fit.Coefficient(f, 1) - 2 * fit.Coefficient(f, 2)).Within(1e-12));
        }
    }

    [Test]
    public void ANonFittedFeatureIsNaNEverywhere()
    {
        var fit = Fit();
        Assert.That(fit.Status[39], Is.EqualTo(FeatureFitStatus.TooFewObservations));
        Assert.That(fit.UnscaledCovariance(39, 0, 1), Is.NaN);
        Assert.That(fit.ContrastStdevUnscaled(39, new[] { 0.0, 1, 0 }), Is.NaN);
        var r = EmpiricalBayes.ModerateContrasts(fit, new[] { new ContrastWeights("x", new[] { 0.0, 1, 1 }) }, trend: false).Single();
        Assert.That(new[] { r.Estimate[39], r.StandardError[39], r.T[39], r.PValue[39], r.ConfidenceLow[39], r.ConfidenceHigh[39] },
            Is.All.NaN);
        Assert.That(r.Estimate.Take(39), Is.All.Not.NaN);
    }

    [TestCase(0.95)]
    [TestCase(0.80)]
    [TestCase(0.99)]
    public void TheIntervalIsSymmetricAtTheRequestedLevel(double level)
    {
        var fit = Fit();
        var r = EmpiricalBayes.ModerateContrasts(fit, new[] { new ContrastWeights("x", new[] { 0.0, 1, -1 }) }, false,
            confidenceLevel: level).Single();
        Assert.That(r.ConfidenceLevel, Is.EqualTo(level));
        for (int f = 0; f < 39; f++)
        {
            double half = StudentT.InvCDF(0, 1, r.DfTotal[f], (1 + level) / 2) * r.StandardError[f];
            Assert.That((r.ConfidenceLow[f] + r.ConfidenceHigh[f]) / 2, Is.EqualTo(r.Estimate[f]).Within(1e-12));
            Assert.That(r.ConfidenceHigh[f] - r.Estimate[f], Is.EqualTo(half).Within(1e-10));
        }
    }

    [Test]
    public void AnInfinitePriorDfUsesThePooledResidualDf()
    {
        // Every feature has the same residuals up to its mean, so the variances are not over-dispersed: d0 = Inf.
        // As limma, df.total = min(d + d0, pooled df) is then the pooled residual df, and the interval is t on it.
        double[] noise = { 0.3, -0.2, 0.1, -0.4, 0.25, -0.05, 0.15, -0.15, 0.0, 0.2 };
        int g = 12, n = noise.Length;
        var design = new double[n, 2];
        var y = new double[g, n];
        for (int i = 0; i < n; i++) { design[i, 0] = 1; design[i, 1] = i % 2; }
        for (int f = 0; f < g; f++)
            for (int i = 0; i < n; i++) y[f, i] = 10 + f + 0.1 * f * design[i, 1] + noise[i];
        var fit = LinearModel.Fit(y, design, new[] { "intercept", "x" });
        var r = EmpiricalBayes.ModerateContrasts(fit, new[] { new ContrastWeights("x", new[] { 0.0, 1 }) }, false).Single();
        Assert.That(r.Prior.Df, Is.EqualTo(double.PositiveInfinity));
        double pooled = fit.DfResidual.Sum();
        double q = StudentT.InvCDF(0, 1, pooled, 0.975);
        for (int f = 0; f < g; f++)
        {
            Assert.That(r.DfTotal[f], Is.EqualTo(pooled));
            Assert.That(r.ConfidenceHigh[f] - r.Estimate[f], Is.EqualTo(q * r.StandardError[f]).Within(1e-12));
        }
    }

    [Test]
    public void ContrastsShareOnePrior()
    {
        var fit = Fit();
        var rs = EmpiricalBayes.ModerateContrasts(fit, new[]
        {
            new ContrastWeights("a", new[] { 0.0, 1, 0 }), new ContrastWeights("b", new[] { 0.0, 0, 1 }),
        }, trend: true);
        Assert.That(rs[0].Prior, Is.SameAs(rs[1].Prior));
        Assert.That(rs[0].PosteriorVariance, Is.EqualTo(rs[1].PosteriorVariance));
        Assert.That(rs[0].Weights, Is.EqualTo(new[] { 0.0, 1, 0 }));
    }

    [Test]
    public void MisuseThrows()
    {
        var fit = Fit();
        IReadOnlyList<ModeratedContrast> Run(double[] w, double level = 0.95) =>
            EmpiricalBayes.ModerateContrasts(fit, new[] { new ContrastWeights("c", w) }, false, confidenceLevel: level);
        Assert.Throws<ArgumentException>(() => Run(new[] { 0.0, 1 }), "wrong length");
        Assert.Throws<ArgumentException>(() => Run(new[] { 0.0, 0, 0 }), "all zero");
        Assert.Throws<ArgumentException>(() => Run(new[] { 0.0, double.NaN, 1 }), "not finite");
        Assert.Throws<ArgumentOutOfRangeException>(() => Run(new[] { 0.0, 1, 0 }, 1.0), "level 1");
        Assert.Throws<ArgumentOutOfRangeException>(() => Run(new[] { 0.0, 1, 0 }, 0.0), "level 0");
        Assert.Throws<ArgumentException>(() => EmpiricalBayes.ModerateContrasts(fit, Array.Empty<ContrastWeights>(), false));
        Assert.Throws<ArgumentException>(() => EmpiricalBayes.ModerateContrasts(fit, new[]
            { new ContrastWeights("c", new[] { 0.0, 1, 0 }), new ContrastWeights("c", new[] { 0.0, 0, 1 }) }, false), "duplicate name");
        Assert.Throws<ArgumentException>(() => fit.ContrastStdevUnscaled(0, new[] { 1.0 }));
    }
}
