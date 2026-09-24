using MathNet.Numerics.Distributions;
using NUnit.Framework;
using Statistics;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.Statistics;

/// <summary>
/// The maximum-marginal-likelihood variance prior (<see cref="VariancePriorEstimator.MarginalLikelihood"/>):
/// recovery of known hyperparameters with unequal residual df, the d0 = ∞ limit, agreement with the moment
/// estimator where both are consistent, and agreement BOUNDS with limma's current default eBayes on the
/// synthetic fixture (a benchmark, not parity: mzLib does not implement limma's fitFDistUnequalDF1).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MarginalLikelihoodPriorTests
{
    private static string DataDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "Statistics", "ReferenceData");

    private static double Parse(string s) => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static (string[] Header, List<string[]> Rows) ReadTsv(string name)
    {
        var lines = File.ReadAllLines(Path.Combine(DataDirectory, name)).Where(l => l.Length > 0).ToArray();
        return (lines[0].Split('\t'), lines.Skip(1).Select(l => l.Split('\t')).ToList());
    }

    private static double[,] Matrix(string name)
    {
        var (header, rows) = ReadTsv(name);
        var m = new double[rows.Count, header.Length];
        for (int i = 0; i < rows.Count; i++)
            for (int j = 0; j < header.Length; j++)
                m[i, j] = Parse(rows[i][j]);
        return m;
    }

    /// <summary>Draws s²_g from the model: σ²_g = d0 s0² / χ²(d0), s²_g = σ²_g χ²(d_g) / d_g.</summary>
    private static double[] Simulate(Random rng, double[] df, double d0, Func<int, double> s0Squared)
    {
        return df.Select((d, g) =>
        {
            double sigma2 = double.IsPositiveInfinity(d0) ? s0Squared(g) : d0 * s0Squared(g) / ChiSquared.Sample(rng, d0);
            return sigma2 * ChiSquared.Sample(rng, d) / d;
        }).ToArray();
    }

    [Test]
    public void RecoversThePriorWithUnequalResidualDf()
    {
        var rng = new Random(20260924);
        int n = 6000;
        double d0 = 8, s0 = 0.5;
        var df = Enumerable.Range(0, n).Select(g => (double)(1 + g % 17)).ToArray();   // 1..17, as omitting missing values gives
        var s2 = Simulate(rng, df, d0, _ => s0);
        var prior = EmpiricalBayes.FitPrior(s2, df, estimator: VariancePriorEstimator.MarginalLikelihood);
        Assert.That(prior.Estimator, Is.EqualTo(VariancePriorEstimator.MarginalLikelihood));
        Assert.That(prior.Df, Is.EqualTo(d0).Within(0.15 * d0));
        Assert.That(prior.Scale[0], Is.EqualTo(s0).Within(0.03 * s0));
        Assert.That(prior.Scale.Distinct().Count(), Is.EqualTo(1), "no trend: one prior variance for all");
    }

    [Test]
    public void NoPriorVariationGivesInfiniteOrVeryLargeDf()
    {
        // Every feature has the same true variance: the prior is a point mass and d0 = ∞ in the limit.
        var rng = new Random(3);
        var df = Enumerable.Range(0, 4000).Select(g => (double)(2 + g % 10)).ToArray();
        var s2 = Simulate(rng, df, double.PositiveInfinity, _ => 1.3);
        var prior = EmpiricalBayes.FitPrior(s2, df, estimator: VariancePriorEstimator.MarginalLikelihood);
        Assert.That(prior.Df, Is.GreaterThan(200));
        Assert.That(prior.Scale[0], Is.EqualTo(1.3).Within(0.02 * 1.3));
    }

    [Test]
    public void AgreesWithTheMomentEstimatorWhenDfAreEqual()
    {
        // With equal residual df both estimators are consistent for the same (d0, s0²).
        var rng = new Random(5);
        var df = Enumerable.Repeat(6.0, 8000).ToArray();
        var s2 = Simulate(rng, df, 5, _ => 0.2);
        var mom = EmpiricalBayes.FitPrior(s2, df);
        var mml = EmpiricalBayes.FitPrior(s2, df, estimator: VariancePriorEstimator.MarginalLikelihood);
        Assert.That(mom.Estimator, Is.EqualTo(VariancePriorEstimator.MomentsLegacy), "legacy remains the default");
        Assert.That(mml.Df, Is.EqualTo(mom.Df).Within(0.15 * mom.Df));
        Assert.That(mml.Scale[0], Is.EqualTo(mom.Scale[0]).Within(0.03 * mom.Scale[0]));
    }

    [Test]
    public void RecoversAnIntensityTrend()
    {
        var rng = new Random(9);
        int n = 5000;
        var amean = Enumerable.Range(0, n).Select(_ => 18 + 8 * rng.NextDouble()).ToArray();
        var df = Enumerable.Range(0, n).Select(g => (double)(2 + g % 12)).ToArray();
        double Truth(int g) => 0.05 * Math.Exp(-(amean[g] - 22) / 3);
        var s2 = Simulate(rng, df, 10, Truth);
        var prior = EmpiricalBayes.FitPrior(s2, df, amean, estimator: VariancePriorEstimator.MarginalLikelihood);
        Assert.That(prior.Trended, Is.True);
        Assert.That(prior.Df, Is.EqualTo(10).Within(2.0));
        var logRatio = Enumerable.Range(0, n).Select(g => Math.Abs(Math.Log(prior.Scale[g] / Truth(g)))).OrderBy(v => v).ToArray();
        Assert.That(logRatio[n / 2], Is.LessThan(0.05), "median prior variance within 5% of the true trend");
    }

    [TestCase(false, "limma_default_notrend.tsv")]
    [TestCase(true, "limma_default_trend.tsv")]
    public void StaysWithinBoundsOfLimmasCurrentDefault(bool trend, string file)
    {
        // Benchmark, not parity: limma's default (fitFDistUnequalDF1) is a different estimator of the same
        // prior. The synthetic fixture has five distinct residual df (5 to 9).
        var fit = LinearModel.Fit(Matrix("limma_responses.tsv"), Matrix("limma_design.tsv"),
            ReadTsv("limma_design.tsv").Header);
        var test = EmpiricalBayes.Moderate(fit, "age_decades", trend, estimator: VariancePriorEstimator.MarginalLikelihood);
        var (ph, prows) = ReadTsv("limma_default_prior.tsv");
        double d0Limma = Parse(prows.Single(r => r[0] == (trend ? "TRUE" : "FALSE"))[1]);
        var (h, rows) = ReadTsv(file);
        int s0Col = Array.IndexOf(h, "s2_prior_default"), pCol = Array.IndexOf(h, "p_default");
        var fitted = Enumerable.Range(0, fit.FeatureCount).Where(f => fit.Status[f] == FeatureFitStatus.Fitted).ToArray();
        Assert.That(fitted.Length, Is.EqualTo(rows.Count));

        var legacy = EmpiricalBayes.Moderate(fit, "age_decades", trend);
        var scaleRatio = new List<double>();
        var dLogP = new List<double>();
        var dLogPLegacy = new List<double>();
        for (int k = 0; k < fitted.Length; k++)
        {
            int f = fitted[k];
            double pLimma = Parse(rows[k][pCol]);
            scaleRatio.Add(Math.Abs(Math.Log(test.Prior.Scale[f] / Parse(rows[k][s0Col]))));
            dLogP.Add(Math.Abs(Math.Log10(test.PValue[f]) - Math.Log10(pLimma)));
            dLogPLegacy.Add(Math.Abs(Math.Log10(legacy.PValue[f]) - Math.Log10(pLimma)));
        }
        double Median(List<double> v) => v.OrderBy(x => x).ElementAt(v.Count / 2);
        TestContext.WriteLine($"trend={trend}: d0 {test.Prior.Df:G6} vs limma default {d0Limma:G6}; " +
            $"median |log s0 ratio| {Median(scaleRatio):G3}, max {scaleRatio.Max():G3}; median |dlog10 p| {Median(dLogP):G3}, " +
            $"max {dLogP.Max():G3} (legacy: median {Median(dLogPLegacy):G3}, max {dLogPLegacy.Max():G3})");

        Assert.That(test.Prior.Df, Is.EqualTo(d0Limma).Within(0.03 * d0Limma));
        Assert.That(Median(scaleRatio), Is.LessThan(0.02));
        Assert.That(Median(dLogP), Is.LessThan(0.005));
        // With a trend, the prior variance at the ends of the intensity range depends on the curve: a natural
        // spline here, lowess in limma. That, not the prior estimator, sets the worst-case difference, and the
        // legacy estimator (same spline) shows it too. So the extremes are bounded more loosely with a trend.
        Assert.That(scaleRatio.Max(), Is.LessThan(trend ? 0.2 : 0.02));
        Assert.That(dLogP.Max(), Is.LessThan(trend ? 0.5 : 0.05));
        // At the worst feature, never further from limma's default than the legacy estimator is. (At the median
        // on this synthetic set, whose df differ little, legacy is as close or slightly closer: 0.0011 vs 0.0018
        // with a trend. On real label-free data, where df range from 1 to 17, this estimator is 4-6x closer.)
        Assert.That(dLogP.Max(), Is.LessThanOrEqualTo(dLogPLegacy.Max()));
    }

    [Test]
    public void LikelihoodIsConcaveEnoughToFindTheSameOptimumFromAnywhere()
    {
        // The profile over d0 is found by grid + Brent: nudging the data slightly moves d0 slightly.
        var rng = new Random(21);
        var df = Enumerable.Range(0, 3000).Select(g => (double)(1 + g % 9)).ToArray();
        var s2 = Simulate(rng, df, 4, _ => 1.0);
        var a = EmpiricalBayes.FitPrior(s2, df, estimator: VariancePriorEstimator.MarginalLikelihood);
        var b = EmpiricalBayes.FitPrior(s2.Select(v => v * (1 + 1e-9)).ToArray(), df, estimator: VariancePriorEstimator.MarginalLikelihood);
        Assert.That(b.Df, Is.EqualTo(a.Df).Within(1e-5 * a.Df));
        Assert.That(b.Scale[0], Is.EqualTo(a.Scale[0] * (1 + 1e-9)).Within(1e-7 * a.Scale[0]));
    }
}
