using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
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
