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
/// What <see cref="NestedMixedModel"/> promises beyond the lme4 reference values: which random structure a protein gets,
/// which proteins are not fittable and why, independence from thread count, and refusals.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class NestedMixedModelTests
{
    private static readonly string[] Names = { "intercept", "x" };

    /// <summary>Two peptides in 8 samples, x = 0/1 by sample; deterministic values with sample-level offsets.</summary>
    private static NestedMixedModelProblem Problem(string?[]? individualOfSample = null, int peptides = 2)
    {
        int s = 8, n = s * peptides;
        var y = new double[n];
        var common = new double[n, 2];
        var nuisance = peptides > 1 ? new double[n, peptides - 1] : null;
        var sample = new string[n];
        var individual = individualOfSample is null ? null : new string?[n];
        double[] offset = { 0.3, -0.2, 0.1, -0.35, 0.25, -0.05, 0.15, -0.2 };
        double[] noise = { 0.02, -0.05, 0.04, 0.01, -0.03, 0.05, -0.01, 0.03, -0.04, 0.02, 0.0, -0.02, 0.03, -0.01, 0.04, -0.03 };
        for (int k = 0; k < peptides; k++)
            for (int j = 0; j < s; j++)
            {
                int i = k * s + j;
                common[i, 0] = 1;
                common[i, 1] = j % 2;
                if (k > 0) nuisance![i, k - 1] = 1;
                sample[i] = $"s{j}";
                if (individual is not null) individual[i] = individualOfSample![j];
                y[i] = 20 + 0.5 * (j % 2) + 0.4 * k + offset[j] + noise[i % noise.Length];
            }
        return new NestedMixedModelProblem(y, common, nuisance, sample, individual);
    }

    [Test]
    public void NoIndividualMeansARandomSampleIntercept()
    {
        var fit = NestedMixedModel.Fit(new[] { Problem() }, Names);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Structure[0], Is.EqualTo(RandomStructure.Sample));
        Assert.That(fit.IndividualVariance[0], Is.NaN);
        Assert.That(fit.SampleVariance[0], Is.GreaterThan(0));
    }

    [Test]
    public void IndividualsWithSeveralSamplesGetBothIntercepts()
    {
        var fit = NestedMixedModel.Fit(new[] { Problem(new string?[] { "a", "a", "b", "b", "c", "c", "d", "d" }) }, Names);
        Assert.That(fit.Structure[0], Is.EqualTo(RandomStructure.IndividualAndSample));
        Assert.That(fit.IndividualVariance[0], Is.Not.NaN);
    }

    [Test]
    public void OneSamplePerIndividualIsOneIntercept()
    {
        // Individual and sample are then the same grouping; fitting both would split one variance arbitrarily.
        var fit = NestedMixedModel.Fit(new[] { Problem(new string?[] { "a", "b", "c", "d", "e", "f", "g", "h" }) }, Names);
        Assert.That(fit.Structure[0], Is.EqualTo(RandomStructure.Individual));
        Assert.That(fit.SampleVariance[0], Is.NaN);
    }

    [Test]
    public void ASinglePeptideProteinIsNotFittable()
    {
        // One value per sample: the sample intercept cannot be told from the residual (lme4 refuses it too).
        var fit = NestedMixedModel.Fit(new[] { Problem(peptides: 1) }, Names);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.TooFewGroups));
        Assert.That(fit.Coefficient(0, 1), Is.NaN);
        Assert.That(fit.SatterthwaiteDf(0, new[] { 0.0, 1 }), Is.NaN);
    }

    [Test]
    public void ResultsDoNotDependOnThreadCount()
    {
        var data = LmmData.Read("nested");
        var a = NestedMixedModel.Fit(data.Problems, data.CommonNames, maxThreads: 1);
        var b = NestedMixedModel.Fit(data.Problems, data.CommonNames, maxThreads: 4);
        for (int f = 0; f < a.FeatureCount; f++)
        {
            Assert.That(b.Coefficient(f, 1), Is.EqualTo(a.Coefficient(f, 1)));
            Assert.That(b.RemlCriterion[f], Is.EqualTo(a.RemlCriterion[f]));
            if (a.Status[f] == FeatureFitStatus.Fitted)
                Assert.That(b.SatterthwaiteDf(f, new[] { 0.0, 1, 0 }), Is.EqualTo(a.SatterthwaiteDf(f, new[] { 0.0, 1, 0 })));
        }
    }

    /// <summary>
    /// 24 features, 10 samples of 3 injections each: log2 values around 18-26 on intercept + condition (between samples)
    /// + run order (within samples), a sample SD from 0.15 to 0.6 (every fourth feature none), a residual SD from 0.15 to
    /// 0.4, and a few values missing. Feature 0's noise sums to zero within every sample, so no sample variance is estimated
    /// at all: the boundary. Seven features end on the boundary, the other 17 inside.
    /// </summary>
    private static (double[,] Y, double[,] Design, string[] Sample) SampleOnlyData()
    {
        var rng = new MersenneTwister(20261013);
        int samples = 10, injections = 3, n = samples * injections, features = 24;
        var design = new double[n, 3];
        var sample = new string[n];
        for (int i = 0; i < n; i++)
        {
            int s = i / injections;
            design[i, 0] = 1;
            design[i, 1] = s >= samples / 2 ? 1 : 0;
            design[i, 2] = (i % injections - 1) + 0.1 * s; // run order, centred within a sample, drifting across samples
            sample[i] = $"s{s:D2}";
        }
        var y = new double[features, n];
        for (int f = 0; f < features; f++)
        {
            double baseline = 18 + 8 * rng.NextDouble(), effect = Normal.Sample(rng, 0, 0.5), drift = Normal.Sample(rng, 0, 0.05);
            double tau = f % 4 == 0 ? 0 : 0.15 + 0.45 * rng.NextDouble(), sigma = 0.15 + 0.25 * rng.NextDouble();
            for (int s = 0; s < samples; s++)
            {
                double offset = tau * Normal.Sample(rng, 0, 1);
                var noise = Enumerable.Range(0, injections).Select(_ => sigma * Normal.Sample(rng, 0, 1)).ToArray();
                if (f == 0) { double mean = noise.Average(); for (int k = 0; k < injections; k++) noise[k] -= mean; }
                for (int k = 0; k < injections; k++)
                {
                    int i = s * injections + k;
                    y[f, i] = baseline + effect * design[i, 1] + drift * design[i, 2] + offset + noise[k];
                }
            }
        }
        foreach (var (f, i) in new[] { (3, 4), (7, 0), (7, 1), (11, 29), (15, 12), (15, 13), (15, 14), (22, 8) })
            y[f, i] = double.NaN;
        return (y, design, sample);
    }

    /// <summary>
    /// How closely the two fits can agree is set by <see cref="MixedModel"/>'s criterion, not by either optimiser's stopping
    /// rule (Brent to 1e-12; Newton to rounding). MixedModel profiles through rss = y′V⁻¹y − r′β, which at log2 values near
    /// 22 cancels to about 1e-10 of the criterion (the two criteria differ by up to 2e-10 even at θ = 0, where both are
    /// least squares; NestedMixedModel sums squared residuals). A minimum of a function that noisy is located only to about
    /// the square root, ~1e-5 in log θ. Variances move with θ at first order (observed up to 5e-6 of σ²); standard errors
    /// (2e-6) and the generalised least-squares fixed effects (7e-8) less; the criterion at its optimum only by its noise
    /// (5e-11). Each tolerance is about ten times what is observed.
    /// </summary>
    private const double SameCriterion = 1e-9, SameFixedEffects = 1e-6, SameStandardErrors = 1e-5, SameVariances = 1e-4;

    [Test]
    public void TheSampleOnlyCaseEqualsMixedModel()
    {
        // With one random intercept per sample, no nuisance columns and one design for every feature, NestedMixedModel fits
        // the model MixedModel fits (each sample a group), both by REML. The parametrisations differ (θ = τ/σ here,
        // τ²/σ² there) and so do the optimisers, but the optimum is the same point.
        var (y, design, sample) = SampleOnlyData();
        string[] names = { "intercept", "conditionB", "runOrder" };
        var reference = MixedModel.Fit(y, design, sample, names, VarianceEstimator.Reml);
        var problems = Enumerable.Range(0, y.GetLength(0)).Select(f =>
            new NestedMixedModelProblem(Enumerable.Range(0, y.GetLength(1)).Select(i => y[f, i]).ToArray(), design, null, sample, null))
            .ToList();
        var fit = NestedMixedModel.Fit(problems, names);

        int interior = 0;
        for (int f = 0; f < fit.FeatureCount; f++)
        {
            string at = $"feature {f}";
            Assert.That(reference.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at + " (MixedModel)");
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at);
            Assert.That(fit.Structure[f], Is.EqualTo(RandomStructure.Sample), at);
            Assert.That(fit.Observed[f], Is.EqualTo(reference.Observed[f]), at);
            double sigma2 = reference.ResidualVariance[f];
            Assert.That(fit.ResidualVariance[f], Is.EqualTo(sigma2).Within(SameVariances * sigma2), at + " residual variance");
            Assert.That(fit.SampleVariance[f], Is.EqualTo(reference.GroupVariance[f]).Within(SameVariances * sigma2), at + " sample variance");
            // The boundary is kept by both, exactly: a variance of 0 here is 0 there.
            Assert.That(fit.SampleVariance[f] == 0, Is.EqualTo(reference.GroupVariance[f] == 0), at + " on the boundary");
            if (reference.GroupVariance[f] > 0) interior++;
            double criterion = -2 * reference.LogLikelihood[f];
            Assert.That(fit.RemlCriterion[f], Is.EqualTo(criterion).Within(SameCriterion * Math.Max(1, Math.Abs(criterion))), at + " REML criterion");
            for (int j = 0; j < names.Length; j++)
            {
                double beta = reference.Coefficient(f, j), se = reference.StandardError(f, j);
                Assert.That(fit.Coefficient(f, j), Is.EqualTo(beta).Within(SameFixedEffects * Math.Max(1, Math.Abs(beta))), $"{at} {names[j]}");
                Assert.That(Math.Sqrt(fit.ResidualVariance[f] * fit.UnscaledCovariance(f, j, j)), Is.EqualTo(se).Within(SameStandardErrors * se),
                    $"{at} SE({names[j]})");
            }
        }
        Assert.That(fit.SampleVariance[0], Is.EqualTo(0), "feature 0 is built on the boundary");
        Assert.That(interior, Is.GreaterThanOrEqualTo(fit.FeatureCount / 2), "most features have a sample variance");
    }

    [Test]
    public void MisuseThrows()
    {
        Assert.Throws<ArgumentException>(() => new NestedMixedModelProblem(new double[3], new double[2, 2], null, new[] { "a", "b", "c" }, null),
            "design rows differ from responses");
        Assert.Throws<ArgumentException>(() => new NestedMixedModelProblem(new double[2], new double[2, 2], null, new[] { "a" }, null),
            "sample labels differ from responses");
        Assert.Throws<ArgumentException>(() => new NestedMixedModelProblem(new double[2], new double[,] { { 1, double.NaN }, { 1, 0 } }, null,
            new[] { "a", "b" }, null), "non-finite design");
        Assert.Throws<ArgumentException>(() => NestedMixedModel.Fit(new[] { Problem() }, new[] { "intercept" }), "names vs columns");
        var fit = NestedMixedModel.Fit(new[] { Problem(), Problem() }, Names);
        Assert.Throws<ArgumentException>(() => fit.SatterthwaiteDf(0, new[] { 1.0 }));
        Assert.Throws<ArgumentException>(() => NestedMixedModel.ModerateContrasts(fit,
            new[] { new ContrastWeights("x", new[] { 0.0, 1 }) }, new IReadOnlyList<double>[] { new[] { 1.0 } }), "variance-df row length");
    }
}
