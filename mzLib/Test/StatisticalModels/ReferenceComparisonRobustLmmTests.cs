using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using StatisticalModels;

namespace Test.StatisticalModels;

/// <summary>
/// <see cref="NestedMixedModel"/>'s robust option against frozen reference values from
/// <c>ReferenceData/make_robust_lmm_fixtures.R</c>: msqrob2 1.20.0's own reweighting loop (<c>.robust_fitting</c>:
/// <c>psi.huber</c> on <c>resid</c> scaled by <c>mad(res, 0)</c>, <c>refit</c>, stop at a 1e-6 relative change in pwrss),
/// run to convergence (QuantProject GR-22), then lmerTest's Satterthwaite df and GR-21's moderation on the final weighted
/// fit. Two scenarios of 30 proteins with about 6% outlying values. R is never run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReferenceComparisonRobustLmmTests
{
    /// <summary>
    /// Weights, estimates and σ². Each round refits from scratch and the loop stops on a 1e-6 relative change in pwrss, so
    /// the final weights agree to the precision that rule leaves, not to the optimizer's.
    /// </summary>
    private const double Fit = 1e-5;

    /// <summary>Satterthwaite df and everything downstream, as in <see cref="ReferenceComparisonNestedLmmTests"/>.</summary>
    private const double Df = 1e-4;

    [TestCase("nested")]
    [TestCase("sample_only")]
    public void RobustFitAndModerationEqualReference(string scenario)
    {
        var data = LmmData.Read(scenario, "rlmm");
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames, robust: true);
        var moderated = NestedMixedModel.ModerateContrasts(fit, data.Contrasts, data.ConditionRows);
        var expected = LmmData.Rows($"rlmm_{scenario}_results.tsv");
        var unscaled = LmmData.Rows($"rlmm_{scenario}_unscaled.tsv").ToDictionary(r => r["protein"]);
        var weights = LmmData.Rows($"rlmm_{scenario}_weights.tsv")
            .ToDictionary(r => (r["protein"], r["peptide"], r["sample"]), r => LmmData.Parse(r["weight"]));

        int downWeighted = 0;
        foreach (var protein in expected.GroupBy(r => r["protein"]))
        {
            int f = data.Proteins.IndexOf(protein.Key);
            var first = protein.First();
            string at = $"{scenario} {protein.Key}";
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at);
            Assert.That(fit.RobustConverged[f], Is.EqualTo(bool.Parse(first["converged"])), at + " converged");
            Assert.That(Math.Abs(fit.RobustIterations[f] - int.Parse(first["iterations"])), Is.LessThanOrEqualTo(1),
                at + " iterations (the stopping rule can fall either side of 1e-6 by one round)");
            var w = fit.Weights(f)!;
            for (int i = 0; i < data.Keys[f].Count; i++)
            {
                var (peptide, sample) = data.Keys[f][i];
                double expectedWeight = weights[(protein.Key, peptide, sample)];
                Assert.That(w[i], Is.EqualTo(expectedWeight).Within(Fit), $"{at} weight {peptide} {sample}");
                if (expectedWeight < 1) downWeighted++;
            }
            Close(fit.ResidualVariance[f], first["sigma2"], Fit, at + " sigma2");
            Close(fit.RemlCriterion[f], first["reml"], Fit, at + " REML");
            Variance(fit.IndividualVariance[f], first["var_individual"], fit.ResidualVariance[f], at + " var(individual)");
            Variance(fit.SampleVariance[f], first["var_sample"], fit.ResidualVariance[f], at + " var(sample)");
            int p = data.CommonNames.Count;
            var u = unscaled[protein.Key];
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                    Close(fit.UnscaledCovariance(f, i, j), u[$"X{i * p + j + 1}"], Fit, $"{at} unscaled[{i},{j}]");
            foreach (var row in protein)
            {
                var r = moderated.Single(m => m.Name == row["contrast"]);
                var c = data.Contrasts.Single(x => x.Name == row["contrast"]).Weights;
                string where = $"{at} {row["contrast"]}";
                Close(r.Estimate[f], row["estimate"], Fit, where + " estimate");
                Close(fit.SatterthwaiteDf(f, c), row["df_satterthwaite"], Df, where + " Satterthwaite df");
                Close(fit.JointSatterthwaiteDf(f, data.ConditionRows), row["variance_df"], Df, where + " joint DenDF");
                Close(r.PosteriorVariance[f], row["s2_post"], Df, where + " s2.post");
                Close(r.Prior.Df, row["df_prior"], Df, where + " df.prior");
                Close(r.StandardError[f], row["se"], Df, where + " se");
                Close(r.DfTotal[f], row["df"], Df, where + " df");
                Close(r.T[f], row["t"], Df, where + " t");
                Close(r.PValue[f], row["p"], Df, where + " p");
                Close(r.BenjaminiHochbergAdjusted[f], row["adj_p"], Df, where + " adj.p");
                Close(r.ConfidenceLow[f], row["ci_low"], Df, where + " CI low");
                Close(r.ConfidenceHigh[f], row["ci_high"], Df, where + " CI high");
            }
        }
        Assert.That(downWeighted, Is.GreaterThan(20), "the fixture's outliers are down-weighted, so the test exercises reweighting");
    }

    [Test]
    public void WithoutRobustEveryWeightIsOne()
    {
        var data = LmmData.Read("nested", "rlmm");
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames);
        for (int f = 0; f < fit.FeatureCount; f++)
        {
            Assert.That(fit.RobustIterations[f], Is.EqualTo(0));
            Assert.That(fit.RobustConverged[f], Is.Null);
            Assert.That(fit.Weights(f)!.Where(double.IsFinite), Is.All.EqualTo(1.0));
        }
    }

    [Test]
    public void RobustFittingMovesEstimatesTowardTheBulk()
    {
        // One wild value (+8 log2) added to one condition-B observation of each protein: the plain condition estimate
        // follows it, the robust one moves less for every protein (median about an eighth as far; at most 0.6 as far, p009,
        // where the sample's random effect takes part of it) and the wild value is down-weighted. Reweighting that pulled
        // toward outliers would move the robust estimate at least as far as the plain one.
        var data = LmmData.Read("nested", "rlmm");
        var wild = data.Problems.Select(problem =>
        {
            var y = problem.Response.ToArray();
            int k = Enumerable.Range(0, y.Length).First(i => problem.Common[i, 1] == 1);
            y[k] += 8;
            return (Problem: Copy(problem, y), Index: k);
        }).ToList();
        var wildProblems = wild.Select(w => w.Problem).ToList();
        var plain = NestedMixedModel.Fit(data.Problems, data.CommonNames);
        var robust = NestedMixedModel.Fit(data.Problems, data.CommonNames, robust: true);
        var plainWild = NestedMixedModel.Fit(wildProblems, data.CommonNames);
        var robustWild = NestedMixedModel.Fit(wildProblems, data.CommonNames, robust: true);
        var ratios = new List<double>();
        for (int f = 0; f < data.Problems.Count; f++)
        {
            double plainShift = Math.Abs(plainWild.Coefficient(f, 1) - plain.Coefficient(f, 1));
            double robustShift = Math.Abs(robustWild.Coefficient(f, 1) - robust.Coefficient(f, 1));
            Assert.That(plainShift, Is.GreaterThan(0.1), $"{data.Proteins[f]}: the wild value moves the plain fit");
            Assert.That(robustShift, Is.LessThan(0.75 * plainShift), $"{data.Proteins[f]}: robust {robustShift} vs plain {plainShift}");
            Assert.That(robustWild.Weights(f)![wild[f].Index], Is.LessThan(0.5), $"{data.Proteins[f]}: the wild value's weight");
            ratios.Add(robustShift / plainShift);
        }
        Assert.That(ratios.Median(), Is.LessThan(0.2), "typically the robust fit moves a small fraction as far");
        Assert.That(Enumerable.Range(0, robust.FeatureCount).All(f => robust.Weights(f)!.Where(double.IsFinite).All(w => w > 0 && w <= 1)));
    }

    [Test]
    public void ReweightingThatReachesItsRoundLimitIsNotConvergedAndLeftOutOfModeration()
    {
        // The reference needs 7-27 rounds on this scenario; with a limit of 15, proteins that need 17 or more cannot
        // converge (the round count may differ from R's by one) and those that need 13 or fewer still do.
        const int limit = 15;
        var data = LmmData.Read("nested", "rlmm");
        var rounds = LmmData.Rows("rlmm_nested_results.tsv").GroupBy(r => r["protein"])
            .ToDictionary(g => g.Key, g => int.Parse(g.First()["iterations"]));
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames, -1, true, limit);
        var moderated = NestedMixedModel.ModerateContrasts(fit, data.Contrasts, data.ConditionRows).Single(m => m.Name == "conditionB");
        int capped = 0, converged = 0;
        for (int f = 0; f < data.Proteins.Count; f++)
        {
            int needed = rounds[data.Proteins[f]];
            string at = $"{data.Proteins[f]} (needs {needed})";
            if (needed >= limit + 2)
            {
                capped++;
                Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.NotConverged), at);
                Assert.That(fit.RobustConverged[f], Is.False, at);
                Assert.That(fit.RobustIterations[f], Is.EqualTo(limit), at);
                Assert.That(fit.Weights(f)!.Where(double.IsFinite).Any(w => w < 1), at + ": the last round's weights are kept");
                Assert.That(moderated.Estimate[f], Is.NaN, at + ": left out of moderation");
                Assert.That(moderated.PValue[f], Is.NaN, at);
            }
            else if (needed <= limit - 2)
            {
                converged++;
                Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at);
                Assert.That(fit.RobustConverged[f], Is.True, at);
                Assert.That(moderated.Estimate[f], Is.Not.NaN, at);
            }
        }
        Assert.That(capped, Is.GreaterThan(5));
        Assert.That(converged, Is.GreaterThan(5));
    }

    [Test]
    public void ZeroMadGivesNoHuberWeights()
    {
        // More than half the residuals exactly 0: MAD is 0, no residual can be scaled, so the weights are refused and left
        // unchanged (the fit then becomes NotConverged), never all set to 1 with the outlier at full weight.
        var weights = new[] { 0.5, 0.5, 0.5, 0.5, 0.5 };
        Assert.That(NestedMixedModel.HuberWeights(new[] { 0.0, 0.0, 0.0, 1.0, 50.0 }, weights), Is.False);
        Assert.That(weights, Is.All.EqualTo(0.5));

        // A usable scale: MAD₀ = 1.4826 · median(0.5, 1, 1, 2, 10) = 1.4826; |r| / MAD₀ above 1.345 gets 1.345 / scaled.
        Assert.That(NestedMixedModel.HuberWeights(new[] { 0.5, -1.0, 1.0, 2.0, -10.0 }, weights), Is.True);
        Assert.That(weights[0], Is.EqualTo(1));
        Assert.That(weights[1], Is.EqualTo(1));
        Assert.That(weights[3], Is.EqualTo(1.345 * 1.4826 / 2.0).Within(1e-12));
        Assert.That(weights[4], Is.EqualTo(1.345 * 1.4826 / 10.0).Within(1e-12));
    }

    [Test]
    public void AMissingResponseHasNaNWeightInItsOwnPlace()
    {
        // The same protein with one extra, missing observation in the middle fits exactly as without it, and its weight
        // vector has NaN at that observation and the other weights in their own places.
        var data = LmmData.Read("nested", "rlmm");
        var problem = data.Problems[0];
        int n = problem.Response.Count, at = n / 2;
        var rows = Enumerable.Range(0, n + 1).Select(i => i <= at ? i : i - 1).ToArray();
        var y = rows.Select((r, i) => i == at ? double.NaN : problem.Response[r]).ToArray();
        var withMissing = new NestedMixedModelProblem(y, Rows(problem.Common, rows),
            problem.Nuisance is null ? null : Rows(problem.Nuisance, rows),
            rows.Select(r => problem.Sample[r]).ToArray(), problem.Individual is null ? null : rows.Select(r => problem.Individual[r]).ToArray());

        var plain = NestedMixedModel.Fit(new[] { problem, problem }, data.CommonNames, robust: true);
        var fit = NestedMixedModel.Fit(new[] { withMissing, withMissing }, data.CommonNames, robust: true);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Observed[0], Is.EqualTo(n));
        var w = fit.Weights(0)!;
        var expected = plain.Weights(0)!;
        Assert.That(w.Count, Is.EqualTo(n + 1));
        Assert.That(w[at], Is.NaN);
        for (int i = 0; i <= n; i++)
            if (i != at) Assert.That(w[i], Is.EqualTo(expected[rows[i]]), $"observation {i}");
        Assert.That(expected.Distinct().Count(), Is.GreaterThan(1), "the weights differ, so a shifted vector would fail");
        Assert.That(fit.Coefficient(0, 1), Is.EqualTo(plain.Coefficient(0, 1)));
    }

    [Test]
    public void RobustResultsDoNotDependOnTheThreadCount()
    {
        var data = LmmData.Read("sample_only", "rlmm");
        var one = NestedMixedModel.Fit(data.Problems, data.CommonNames, maxThreads: 1, robust: true);
        var many = NestedMixedModel.Fit(data.Problems, data.CommonNames, maxThreads: 4, robust: true);
        for (int f = 0; f < data.Problems.Count; f++)
        {
            Assert.That(many.Status[f], Is.EqualTo(one.Status[f]));
            Assert.That(many.RobustIterations[f], Is.EqualTo(one.RobustIterations[f]));
            Assert.That(many.Weights(f), Is.EqualTo(one.Weights(f)));
            Assert.That(many.ResidualVariance[f], Is.EqualTo(one.ResidualVariance[f]));
            for (int i = 0; i < data.CommonNames.Count; i++)
                Assert.That(many.Coefficient(f, i), Is.EqualTo(one.Coefficient(f, i)));
        }
    }

    private static NestedMixedModelProblem Copy(NestedMixedModelProblem problem, double[] response) =>
        new(response, problem.Common, problem.Nuisance, problem.Sample, problem.Individual);

    private static double[,] Rows(double[,] matrix, int[] rows)
    {
        var copy = new double[rows.Length, matrix.GetLength(1)];
        for (int i = 0; i < rows.Length; i++)
            for (int j = 0; j < matrix.GetLength(1); j++) copy[i, j] = matrix[rows[i], j];
        return copy;
    }

    /// <summary>A variance component, compared absolutely on the residual variance's scale (lme4 can sit at or near 0).</summary>
    private static void Variance(double actual, string expectedText, double sigma2, string where)
    {
        double expected = LmmData.Parse(expectedText);
        if (double.IsNaN(expected)) { Assert.That(actual, Is.NaN, where); return; }
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(Df * sigma2), $"{where}: {actual:R} vs {expected:R}");
    }

    private static void Close(double actual, string expectedText, double tolerance, string where)
    {
        double expected = LmmData.Parse(expectedText);
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Abs(expected))),
            $"{where}: {actual:R} vs {expected:R}");
    }
}
