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
/// fit. Two scenarios of 30 proteins with about 6% outlying values. Reweighting that never settles falls back to msqrob2's
/// default of one round (GR-22 amended 2026-11-03): the "fallback" scenario (2-peptide proteins on 4 samples, some never
/// settling) and the "set" fixture (one SETTEST1 set model) hold R's one-round fits, whose weights the script checks
/// against <c>msqrob2:::.robust_fitting(maxitRob = 1)</c>. R is never run by these tests.
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

    /// <summary>
    /// A one-round fit's weights, estimates and σ²: no stopping rule stands between the two fits, only each optimizer's
    /// precision on the unweighted and the weighted fit.
    /// </summary>
    private const double OneRound = 1e-6;

    [TestCase("nested")]
    [TestCase("sample_only")]
    [TestCase("fallback")]
    public void RobustFitAndModerationEqualReference(string scenario)
    {
        var data = LmmData.Read(scenario, "rlmm");
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames, robust: true);
        var moderated = NestedMixedModel.ModerateContrasts(fit, data.Contrasts, data.ConditionRows);
        var expected = LmmData.Rows($"rlmm_{scenario}_results.tsv");
        var unscaled = LmmData.Rows($"rlmm_{scenario}_unscaled.tsv").ToDictionary(r => r["protein"]);
        var weights = LmmData.Rows($"rlmm_{scenario}_weights.tsv")
            .ToDictionary(r => (r["protein"], r["peptide"], r["sample"]), r => LmmData.Parse(r["weight"]));

        int downWeighted = 0, fellBack = 0;
        foreach (var protein in expected.GroupBy(r => r["protein"]))
        {
            int f = data.Proteins.IndexOf(protein.Key);
            var first = protein.First();
            string at = $"{scenario} {protein.Key}";
            // R's loop either settles or reaches its 100 rounds and falls back to one round; both are ordinary fits.
            bool settled = bool.Parse(first["converged"]);
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at);
            Assert.That(fit.RobustConverged[f], Is.EqualTo(settled), at + " converged");
            Assert.That(fit.RobustFellBack[f], Is.EqualTo(!settled), at + " fell back to one round");
            if (!settled) fellBack++;
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
        Assert.That(downWeighted, Is.GreaterThan(scenario == "fallback" ? 5 : 20),
            "the fixture's outliers are down-weighted, so the test exercises reweighting");
        Assert.That(fellBack, scenario == "fallback" ? Is.GreaterThanOrEqualTo(2) : Is.Zero,
            "the fallback scenario holds proteins that never settle; the others settle");
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
            Assert.That(fit.RobustFellBack[f], Is.Null);
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
    public void ReweightingThatReachesItsRoundLimit_FallsBackToOneRound_AndIsModerated()
    {
        // GR-22 amended 2026-11-03 (was: NotConverged and left out of moderation). The reference needs 7-27 rounds on this
        // scenario; with a limit of 15, proteins that need 17 or more cannot converge (the round count may differ from R's
        // by one) and those that need 13 or fewer still do. A capped protein is fitted exactly as a fit allowed one round.
        const int limit = 15;
        var data = LmmData.Read("nested", "rlmm");
        var rounds = LmmData.Rows("rlmm_nested_results.tsv").GroupBy(r => r["protein"])
            .ToDictionary(g => g.Key, g => int.Parse(g.First()["iterations"]));
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames, -1, true, limit);
        var oneRound = NestedMixedModel.Fit(data.Problems, data.CommonNames, -1, true, 1);
        var moderated = NestedMixedModel.ModerateContrasts(fit, data.Contrasts, data.ConditionRows).Single(m => m.Name == "conditionB");
        var contrast = data.Contrasts.Single(c => c.Name == "conditionB").Weights;
        int capped = 0, converged = 0;
        for (int f = 0; f < data.Proteins.Count; f++)
        {
            int needed = rounds[data.Proteins[f]];
            string at = $"{data.Proteins[f]} (needs {needed})";
            if (needed >= limit + 2)
            {
                capped++;
                Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at);
                Assert.That(fit.RobustConverged[f], Is.False, at);
                Assert.That(fit.RobustFellBack[f], Is.True, at);
                Assert.That(fit.RobustIterations[f], Is.EqualTo(limit), at + ": the rounds run before falling back");
                Assert.That(oneRound.RobustFellBack[f], Is.True, at);
                Assert.That(oneRound.RobustIterations[f], Is.EqualTo(1), at);
                Assert.That(fit.Weights(f), Is.EqualTo(oneRound.Weights(f)), at + ": the first round's weights");
                Assert.That(fit.Weights(f)!.Where(double.IsFinite).Any(w => w < 1), at);
                for (int i = 0; i < data.CommonNames.Count; i++)
                    Assert.That(fit.Coefficient(f, i), Is.EqualTo(oneRound.Coefficient(f, i)), at);
                Assert.That(fit.ResidualVariance[f], Is.EqualTo(oneRound.ResidualVariance[f]), at);
                Assert.That(fit.RemlCriterion[f], Is.EqualTo(oneRound.RemlCriterion[f]), at);
                Assert.That(fit.SatterthwaiteDf(f, contrast), Is.EqualTo(oneRound.SatterthwaiteDf(f, contrast)), at);
                Assert.That(moderated.Estimate[f], Is.EqualTo(fit.Coefficient(f, 1)), at + ": in moderation");
                Assert.That(moderated.PValue[f], Is.Not.NaN, at);
            }
            else if (needed <= limit - 2)
            {
                converged++;
                Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at);
                Assert.That(fit.RobustConverged[f], Is.True, at);
                Assert.That(fit.RobustFellBack[f], Is.False, at);
                Assert.That(moderated.Estimate[f], Is.Not.NaN, at);
            }
        }
        Assert.That(capped, Is.GreaterThan(5));
        Assert.That(converged, Is.GreaterThan(5));
    }

    [Test]
    public void ASetModelWhoseReweightingNeverSettles_FallsBackToOneRound_AsR()
    {
        // SETTEST1 S2's set of three 2-peptide members on 4 samples (rlmm_set_*). Reweighting drives one peptide's weights
        // toward 0 and never settles: R cycles for all 100 rounds; here a later round's refit fails first (round 41). Either
        // way the set is fitted on the first round's weights, which the script checks are msqrob2's maxitRob = 1 weights,
        // and its fit, variance components and Satterthwaite df match lmerTest's fit with those weights in the call.
        var rows = LmmData.Rows("rlmm_set_values.tsv");
        var header = LmmData.Lines("rlmm_set_values.tsv")[0];
        var commonColumns = header.Where(h => h.StartsWith("common_")).ToList();
        var nuisanceColumns = header.Where(h => h.StartsWith("nuisance_")).ToList();
        var names = commonColumns.Select(h => h["common_".Length..]).ToList();
        int n = rows.Count;
        var common = new double[n, commonColumns.Count];
        var nuisance = new double[n, nuisanceColumns.Count];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < commonColumns.Count; j++) common[i, j] = LmmData.Parse(rows[i][commonColumns[j]]);
            for (int j = 0; j < nuisanceColumns.Count; j++) nuisance[i, j] = LmmData.Parse(rows[i][nuisanceColumns[j]]);
        }
        // NestedMixedModel's IndividualAndSample with the sample as the "individual" and the member as the "sample" in it.
        var problem = new NestedMixedModelProblem(rows.Select(r => LmmData.Parse(r["y"])).ToArray(), common, nuisance,
            rows.Select(r => r["member"]).ToArray(), rows.Select(r => (string?)r["sample"]).ToArray());

        var fit = NestedMixedModel.Fit(new[] { problem }, names, maxThreads: 1, robust: true);

        var expected = LmmData.Rows("rlmm_set_results.tsv");
        Assert.That(bool.Parse(expected[0]["converged"]), Is.False, "R never settles");
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Structure[0], Is.EqualTo(RandomStructure.IndividualAndSample));
        Assert.That(fit.RobustConverged[0], Is.False);
        Assert.That(fit.RobustFellBack[0], Is.True);
        Assert.That(fit.RobustIterations[0], Is.LessThan(100),
            "a later round's refit fails before the round limit, so this covers that way of not settling");
        var weights = LmmData.Rows("rlmm_set_weights.tsv");
        var w = fit.Weights(0)!;
        for (int i = 0; i < n; i++)
        {
            Assert.That((weights[i]["peptide"], weights[i]["sample"]), Is.EqualTo((rows[i]["peptide"], rows[i]["sample"])));
            Close(w[i], weights[i]["weight"], OneRound, $"weight {rows[i]["peptide"]} {rows[i]["sample"]}");
        }
        Assert.That(w.Count(x => x < 1), Is.GreaterThan(0), "the one round down-weights");
        var first = expected[0];
        Close(fit.ResidualVariance[0], first["sigma2"], OneRound, "sigma2");
        Close(fit.RemlCriterion[0], first["reml"], OneRound, "REML");
        Variance(fit.IndividualVariance[0], first["var_sample"], fit.ResidualVariance[0], "var (1|sample)");
        Variance(fit.SampleVariance[0], first["var_sample_member"], fit.ResidualVariance[0], "var (1|sample:member)");
        var contrasts = LmmData.Rows("rlmm_set_contrasts.tsv").ToDictionary(c => c["contrast"], c => names.Select(x => LmmData.Parse(c[x])).ToArray());
        foreach (var row in expected)
        {
            var c = contrasts[row["contrast"]];
            double estimate = 0, unscaled = 0;
            for (int i = 0; i < c.Length; i++)
            {
                estimate += c[i] * fit.Coefficient(0, i);
                for (int j = 0; j < c.Length; j++) unscaled += c[i] * c[j] * fit.UnscaledCovariance(0, i, j);
            }
            Close(estimate, row["estimate"], OneRound, row["contrast"] + " estimate");
            Close(Math.Sqrt(unscaled * fit.ResidualVariance[0]), row["se"], OneRound, row["contrast"] + " se");
            Close(fit.SatterthwaiteDf(0, c), row["df_satterthwaite"], Df, row["contrast"] + " Satterthwaite df");
        }
    }

    [Test]
    public void ReweightingAllowedNoRound_HasNoRoundToFallBackTo()
    {
        // A first round that cannot be formed leaves nothing to fall back to: the feature is NotConverged and out of
        // moderation. A limit of no rounds reaches that branch directly (a zero MAD or a failed first refit does too).
        var data = LmmData.Read("nested", "rlmm");
        var fit = NestedMixedModel.Fit(data.Problems.Take(2).ToList(), data.CommonNames, -1, true, 0);
        for (int f = 0; f < 2; f++)
        {
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.NotConverged));
            Assert.That(fit.RobustConverged[f], Is.False);
            Assert.That(fit.RobustFellBack[f], Is.False);
            Assert.That(fit.RobustIterations[f], Is.Zero);
            Assert.That(fit.Weights(f)!.Where(double.IsFinite), Is.All.EqualTo(1.0));
        }
    }

    [Test]
    public void ZeroMadGivesNoHuberWeights()
    {
        // More than half the residuals exactly 0: MAD is 0, no residual can be scaled, so the weights are refused and left
        // unchanged, never all set to 1 with the outlier at full weight. In a first round the fit is then NotConverged; in a
        // later round it falls back to the first round's weights.
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
            Assert.That(many.RobustFellBack[f], Is.EqualTo(one.RobustFellBack[f]));
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
