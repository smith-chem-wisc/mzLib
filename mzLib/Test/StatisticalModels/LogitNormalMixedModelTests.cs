using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using StatisticalModels;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// The logit-normal mixed model: its quadrature and count probabilities in closed form, its reduction to
/// <see cref="MixedModel"/> when nothing is counted and to logistic regression when nothing varies between groups,
/// recovery of simulated parameters, separation, the detection-link fit, and thread-count determinism.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class LogitNormalMixedModelTests
{
    private static double Logit(double p) => Math.Log(p / (1 - p));
    private static double Expit(double x) => 1 / (1 + Math.Exp(-x));

    [Test]
    public void GaussHermiteIntegratesPolynomialsExactly()
    {
        // ∫ x^(2k) e^(−x²) dx = Γ(k + ½); odd moments vanish. An n-node rule is exact to degree 2n − 1.
        foreach (int n in new[] { 20, 40 })
        {
            var (x, w) = LogitNormalMixedModel.GaussHermite(n);
            for (int k = 0; k <= 10; k++)
            {
                double even = x.Zip(w, (xi, wi) => wi * Math.Pow(xi, 2 * k)).Sum();
                Assert.That(even, Is.EqualTo(SpecialFunctions.Gamma(k + 0.5)).Within(1e-10 * SpecialFunctions.Gamma(k + 0.5)), $"n {n}, x^{2 * k}");
                Assert.That(x.Zip(w, (xi, wi) => wi * Math.Pow(xi, 2 * k + 1)).Sum(), Is.EqualTo(0).Within(1e-9));
            }
        }
    }

    [Test]
    public void CountProbabilitiesAreBinomialOrBetaBinomial()
    {
        foreach (double g in new[] { 0.02, 0.3, 0.5, 0.9 })
        {
            double logit = Logit(g);
            for (int k = 0; k <= 7; k++)
                Assert.That(LogitNormalMixedModel.LogCountProbability(k, 7, logit, 0), Is.EqualTo(Binomial.PMFLn(g, 7, k)).Within(1e-12));
            foreach (double rho in new[] { 0.05, 0.3, 0.8 })
            {
                double s = (1 - rho) / rho, a = g * s, b = (1 - g) * s;
                double total = 0;
                for (int k = 0; k <= 7; k++)
                {
                    double closed = SpecialFunctions.BinomialLn(7, k) + SpecialFunctions.BetaLn(k + a, 7 - k + b) - SpecialFunctions.BetaLn(a, b);
                    double mine = LogitNormalMixedModel.LogCountProbability(k, 7, logit, rho);
                    Assert.That(mine, Is.EqualTo(closed).Within(1e-10));
                    total += Math.Exp(mine);
                }
                Assert.That(total, Is.EqualTo(1).Within(1e-12));
            }
        }
        Assert.That(LogitNormalMixedModel.LogCountProbability(3, 9, 0.4, 1e-9), Is.EqualTo(LogitNormalMixedModel.LogCountProbability(3, 9, 0.4, 0)).Within(1e-6));

        // The first and second derivatives in the logit, against central differences.
        foreach (double rho in new[] { 0.0, 0.2 })
            foreach (var (k, n) in new[] { (0, 4), (2, 5), (6, 6) })
                foreach (double logit in new[] { -2.0, 0.3, 1.5 })
                {
                    const double h = 1e-4;
                    double P(double l) => LogitNormalMixedModel.LogCountProbability(k, n, l, rho);
                    LogitNormalMixedModel.LogCountProbability(k, n, logit, rho, SpecialFunctions.BinomialLn(n, k), out double d1, out double d2);
                    Assert.That(d1, Is.EqualTo((P(logit + h) - P(logit - h)) / (2 * h)).Within(1e-7), $"d1 {k}/{n} ρ {rho}");
                    Assert.That(d2, Is.EqualTo((P(logit + h) - 2 * P(logit) + P(logit - h)) / (h * h)).Within(1e-5), $"d2 {k}/{n} ρ {rho}");
                }
    }

    /// <summary>8 groups of 3; an intercept, a within-group covariate and a group-level trait. Nothing counted.</summary>
    private static (double[,] y, double[,] design, string[] groups) MeasuredOnly(int seed)
    {
        int G = 8, k = 3, n = G * k;
        var rng = new Random(seed);
        var design = new double[n, 3];
        var y = new double[1, n];
        var groups = new string[n];
        var u = Enumerable.Range(0, G).Select(_ => Normal.Sample(rng, 0, 0.8)).ToArray();
        for (int s = 0; s < n; s++)
        {
            int g = s / k;
            design[s, 0] = 1; design[s, 1] = s % k; design[s, 2] = g % 2;
            groups[s] = $"animal{g}";
            y[0, s] = -1 + 0.3 * design[s, 1] + 0.7 * design[s, 2] + u[g] + Normal.Sample(rng, 0, 0.5);
        }
        return (y, design, groups);
    }

    private static BinomialCount?[,] NoCounts(int features, int samples) => new BinomialCount?[features, samples];

    [Test]
    public void WithNothingCountedItIsMixedModel()
    {
        var (y, design, groups) = MeasuredOnly(11);
        var counts = NoCounts(1, y.GetLength(1));
        var reference = MixedModel.Fit(y, design, groups);
        var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, DetectionLink.Identity);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Estimator, Is.EqualTo(VarianceEstimator.Reml));
        Assert.That(fit.Counted[0], Is.EqualTo(0));
        Assert.That(fit.Measured[0], Is.EqualTo(24));
        for (int j = 0; j < 3; j++)
        {
            Assert.That(fit.Coefficient(0, j), Is.EqualTo(reference.Coefficient(0, j)).Within(1e-6));
            Assert.That(fit.StandardError(0, j), Is.EqualTo(reference.StandardError(0, j)).Within(1e-5 * reference.StandardError(0, j)));
            Assert.That(fit.DegreesOfFreedom(0, j), Is.EqualTo(reference.DegreesOfFreedom(0, j)));
            Assert.That(fit.PValue(0, j), Is.EqualTo(reference.PValue(0, j)).Within(1e-5));
        }
        Assert.That(fit.ResidualVariance[0], Is.EqualTo(reference.ResidualVariance[0]).Within(1e-5 * reference.ResidualVariance[0]));
        Assert.That(fit.GroupVariance[0], Is.EqualTo(reference.GroupVariance[0]).Within(1e-5 * reference.GroupVariance[0]));
        Assert.That(fit.LogLikelihood[0], Is.EqualTo(reference.LogLikelihood[0]).Within(1e-7), "the Laplace REML is exact here");

        var ml = LogitNormalMixedModel.Fit(y, counts, design, groups, DetectionLink.Identity, estimator: VarianceEstimator.MaximumLikelihood);
        var referenceMl = MixedModel.Fit(y, design, groups, estimator: VarianceEstimator.MaximumLikelihood);
        Assert.That(ml.Coefficient(0, 2), Is.EqualTo(referenceMl.Coefficient(0, 2)).Within(1e-6));
        Assert.That(ml.ResidualVariance[0], Is.EqualTo(referenceMl.ResidualVariance[0]).Within(1e-5 * referenceMl.ResidualVariance[0]));
        Assert.That(ml.LogLikelihood[0], Is.EqualTo(referenceMl.LogLikelihood[0]).Within(1e-7));
    }

    [Test]
    public void WithNoGroupVarianceAndTinyScatterCountsAreLogisticRegression()
    {
        // Every group holds the same counts at the same covariate values, so the between-group variance is on its
        // boundary; with σ² given as 1e-6 the model is a binomial GLM on logit g = xᵀβ (identity link).
        int G = 4;
        int[] trials = { 6, 6, 6, 6, 6 };
        int[] successes = { 1, 2, 3, 4, 5 };
        double[] covariate = { -2, -1, 0, 1, 2 };
        int n = G * covariate.Length;
        var design = new double[n, 2];
        var groups = new string[n];
        var counts = new BinomialCount?[1, n];
        var y = new double[1, n];
        var bernoulli = new List<(double x, double outcome)>();
        for (int s = 0; s < n; s++)
        {
            int i = s % covariate.Length;
            design[s, 0] = 1; design[s, 1] = covariate[i];
            groups[s] = $"g{s / covariate.Length}";
            counts[0, s] = new BinomialCount(successes[i], trials[i]);
            y[0, s] = double.NaN;
            for (int t = 0; t < trials[i]; t++) bernoulli.Add((covariate[i], t < successes[i] ? 1 : 0));
        }
        var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, DetectionLink.Identity, residualVariance: new[] { 1e-6 },
            estimator: VarianceEstimator.MaximumLikelihood);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.GroupVariance[0], Is.EqualTo(0), "no variation between groups");
        Assert.That(fit.ResidualVarianceGiven[0], Is.True);
        Assert.That(fit.ResidualVariance[0], Is.EqualTo(1e-6).Within(1e-18));

        var by = new double[1, bernoulli.Count];
        var bd = new double[bernoulli.Count, 2];
        for (int i = 0; i < bernoulli.Count; i++) { by[0, i] = bernoulli[i].outcome; bd[i, 0] = 1; bd[i, 1] = bernoulli[i].x; }
        var glm = LogisticRegression.Fit(by, bd);
        Assert.That(fit.Coefficient(0, 0), Is.EqualTo(glm.Coefficient(0, 0)).Within(1e-4));
        Assert.That(fit.Coefficient(0, 1), Is.EqualTo(glm.Coefficient(0, 1)).Within(1e-4));
        Assert.That(fit.StandardError(0, 1), Is.EqualTo(glm.StandardError(0, 1)).Within(1e-3 * glm.StandardError(0, 1)));
    }

    /// <summary>
    /// 16 animals of 4 runs; a group-level trait (0/1) and an intercept. Latent = −0.5 + 1.0 trait + u + e,
    /// u ~ N(0, 0.6²), e ~ N(0, 0.5²). Every other run is measured; the rest are counted, 5 trials each through the
    /// link (0.4, 0.6), binomial.
    /// </summary>
    private static (double[,] y, BinomialCount?[,] counts, double[,] design, string[] groups) Simulated(int seed, int features = 1)
    {
        int G = 16, k = 4, n = G * k;
        var rng = new Random(seed);
        var design = new double[n, 2];
        var groups = new string[n];
        var y = new double[features, n];
        var counts = new BinomialCount?[features, n];
        for (int s = 0; s < n; s++) { int g = s / k; design[s, 0] = 1; design[s, 1] = g % 2; groups[s] = $"a{g}"; }
        for (int f = 0; f < features; f++)
        {
            var u = Enumerable.Range(0, G).Select(_ => Normal.Sample(rng, 0, 0.6)).ToArray();
            for (int s = 0; s < n; s++)
            {
                double latent = -0.5 + 1.0 * design[s, 1] + u[s / k] + Normal.Sample(rng, 0, 0.5);
                if (s % 2 == 0) { y[f, s] = latent; continue; }
                y[f, s] = double.NaN;
                double g = Expit(0.4 + 0.6 * latent);
                int successes = Enumerable.Range(0, 5).Count(_ => rng.NextDouble() < g);
                counts[f, s] = new BinomialCount(successes, 5);
            }
        }
        return (y, counts, design, groups);
    }

    [Test]
    public void MeasuredAndCountedSamplesRecoverTheTraitEffect()
    {
        var link = new DetectionLink(0.4, 0.6);
        var estimates = new List<(double b, double se)>();
        foreach (int seed in new[] { 1, 2, 3 })
        {
            var (y, counts, design, groups) = Simulated(seed);
            var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, link, coefficientNames: new[] { "Intercept", "trait" });
            Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
            Assert.That(fit.Measured[0], Is.EqualTo(32));
            Assert.That(fit.Counted[0], Is.EqualTo(32));
            Assert.That(fit.DegreesOfFreedom(0, 1), Is.EqualTo(16 - 2), "a group-level trait: G − 1 − 1 df");
            estimates.Add((fit.Coefficient(0, fit.IndexOf("trait")), fit.StandardError(0, 1)));
        }
        foreach (var (b, se) in estimates)
            Assert.That(Math.Abs(b - 1.0), Is.LessThan(3 * se), $"trait {b} ± {se}");

        // The counted samples add information: the trait's SE is smaller than from the measured samples alone.
        var (y1, counts1, design1, groups1) = Simulated(1);
        var both = LogitNormalMixedModel.Fit(y1, counts1, design1, groups1, link);
        var measuredOnly = MixedModel.Fit(y1, design1, groups1);
        Assert.That(both.StandardError(0, 1), Is.LessThan(measuredOnly.StandardError(0, 1)));
    }

    [Test]
    public void SeparatedCountsAreReportedNotEstimated()
    {
        // Animals with trait 1 show every trial succeed, animals with trait 0 none: the trait runs to infinity.
        int G = 6, k = 3, n = G * k;
        var design = new double[n, 2];
        var groups = new string[n];
        var y = new double[1, n];
        var counts = new BinomialCount?[1, n];
        for (int s = 0; s < n; s++)
        {
            int g = s / k;
            design[s, 0] = 1; design[s, 1] = g % 2; groups[s] = $"a{g}";
            y[0, s] = double.NaN;
            counts[0, s] = new BinomialCount(g % 2 == 1 ? 4 : 0, 4);
        }
        var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, DetectionLink.Identity, residualVariance: new[] { 0.25 });
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Separated));
        Assert.That(fit.PValue(0, 1), Is.NaN);
    }

    private static List<LinkCell> LinkCells(DetectionLink truth, int count, int seed)
    {
        var rng = new Random(seed);
        var cells = new List<LinkCell>();
        for (int i = 0; i < count; i++)
        {
            double mean = Normal.Sample(rng, 0, 1.5), sd = 0.5, depth = Normal.Sample(rng, 0, 1);
            double latent = mean + Normal.Sample(rng, 0, sd);
            double g = Expit(truth.Intercept + truth.Slope * latent + truth.DepthSlope * depth);
            int trials = 1 + rng.Next(6);
            // Beta-binomial: one chance per cell drawn from Beta(g s, (1 − g) s), s = (1 − ρ)/ρ; binomial when ρ = 0.
            double s = truth.Overdispersion > 0 ? (1 - truth.Overdispersion) / truth.Overdispersion : double.NaN;
            double chance = truth.Overdispersion > 0 ? Beta.Sample(rng, g * s, (1 - g) * s) : g;
            int successes = Enumerable.Range(0, trials).Count(_ => rng.NextDouble() < chance);
            cells.Add(new LinkCell(mean, sd, successes, trials, depth));
        }
        return cells;
    }

    [Test]
    public void FitLinkRecoversASimulatedLink()
    {
        // On detection (default): whether any trial succeeded.
        var truth = new DetectionLink(-0.5, 0.5, 0.3);
        var link = LogitNormalMixedModel.FitLink(LinkCells(truth, 4000, 5), withDepth: true)!;
        Assert.That(link.Intercept, Is.EqualTo(truth.Intercept).Within(0.15));
        Assert.That(link.Slope, Is.EqualTo(truth.Slope).Within(0.1));
        Assert.That(link.DepthSlope, Is.EqualTo(truth.DepthSlope).Within(0.1));
        Assert.That(link.Overdispersion, Is.EqualTo(0));

        // On the full counts, with correlated trials.
        var overdispersed = truth with { Overdispersion = 0.1 };
        var cells = LinkCells(overdispersed, 4000, 6);
        var full = LogitNormalMixedModel.FitLink(cells, withDepth: true, withOverdispersion: true, detectionOnly: false)!;
        Assert.That(full.Intercept, Is.EqualTo(truth.Intercept).Within(0.15));
        Assert.That(full.Slope, Is.EqualTo(truth.Slope).Within(0.1));
        Assert.That(full.Overdispersion, Is.EqualTo(0.1).Within(0.05));

        var plain = LogitNormalMixedModel.FitLink(cells)!;
        Assert.That(plain.DepthSlope, Is.EqualTo(0));
        Assert.That(plain.Overdispersion, Is.EqualTo(0));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.FitLink(new[] { new LinkCell(0, 1, 3, 2) }));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.FitLink(Array.Empty<LinkCell>()));
    }

    [Test]
    public void ResultsDoNotDependOnTheThreadCount()
    {
        var (y, counts, design, groups) = Simulated(9, features: 4);
        var link = new DetectionLink(0.4, 0.6);
        var one = LogitNormalMixedModel.Fit(y, counts, design, groups, link, maxThreads: 1);
        var many = LogitNormalMixedModel.Fit(y, counts, design, groups, link, maxThreads: 4);
        for (int f = 0; f < 4; f++)
        {
            Assert.That(many.Coefficient(f, 1), Is.EqualTo(one.Coefficient(f, 1)));
            Assert.That(many.StandardError(f, 1), Is.EqualTo(one.StandardError(f, 1)));
            Assert.That(many.GroupVariance[f], Is.EqualTo(one.GroupVariance[f]));
        }
    }

    [Test]
    public void InvalidInputsAreRefused()
    {
        var (y, counts, design, groups) = Simulated(4);
        var link = new DetectionLink(0.4, 0.6);
        var both = (double[,])y.Clone();
        both[0, 1] = 0.5;                                                   // sample 1 is also counted
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(both, counts, design, groups, link));
        var bad = (BinomialCount?[,])counts.Clone();
        bad[0, 1] = new BinomialCount(6, 5);
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(y, bad, design, groups, link));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(y, counts, design, groups, link with { DepthSlope = 0.2 }));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(y, counts, design, groups, link with { Slope = 0 }));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(y, counts, design, groups, link with { Overdispersion = 1 }));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(y, counts, design, groups, link, residualVariance: new[] { 1.0, 2.0 }));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(y, counts, design, groups, link, depth: new[] { 1.0 }));
        Assert.Throws<ArgumentException>(() => LogitNormalMixedModel.Fit(y, new BinomialCount?[1, 3], design, groups, link));
    }
}
