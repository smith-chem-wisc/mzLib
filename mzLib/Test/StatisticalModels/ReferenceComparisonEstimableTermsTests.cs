using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using StatisticalModels;

namespace Test.StatisticalModels;

/// <summary>
/// <see cref="NestedMixedModel.Fit"/> with <c>estimableRandomTermsOnly</c> (QuantProject STAT1 ruling 16) against R on two
/// real proteins of the MSstatsTMT controlled mixtures whose full model, y ~ ion + condition + (1 | individual / sample),
/// lme4 refuses: each has one ion in every plex run it is in, so a channel within a run (the inner term) holds one value
/// (make_estimable_terms_lmm_fixtures.R; PROVENANCE_estimable_terms_lmm.txt quotes lme4's refusal). A3KMH1 is seen in two
/// or more runs of one mixture, so its individual term stays: lme4 <c>lmer</c> with (1 | individual) and lmerTest's
/// Satterthwaite df. A6NED2 is seen at most once per individual, so no term stays: <c>lm</c>, residual df. Without the
/// option both are refused as lme4 refuses them.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReferenceComparisonEstimableTermsTests
{
    /// <summary>The REML criterion, relative (as the other nested-model reference tests).</summary>
    private const double Reml = 1e-9;
    /// <summary>σ², variance components, coefficients and their standard errors, relative.</summary>
    private const double Fit = 1e-6;
    /// <summary>Satterthwaite df, relative: numerical derivatives, as lmerTest's.</summary>
    private const double Df = 1e-4;

    private static readonly string[] Levels = { "0.5", "0.667", "1" };
    private static readonly string[] Names = { "intercept", "condition0.5", "condition0.667", "condition1" };

    private static NestedMixedModelProblem Problem(string protein)
    {
        var values = LmmData.Rows("lmm_estimable_terms_values.tsv").Where(r => r["protein"] == protein).ToList();
        var common = new double[values.Count, Levels.Length + 1];
        for (int i = 0; i < values.Count; i++)
        {
            common[i, 0] = 1;
            for (int k = 0; k < Levels.Length; k++) common[i, k + 1] = values[i]["condition"] == Levels[k] ? 1 : 0;
        }
        return NestedMixedModelProblem.FromNuisanceLabels(values.Select(v => LmmData.Parse(v["y"])).ToArray(), common,
            values.Select(v => v["ion"]).ToArray(), values.Select(v => v["sample"]).ToArray(),
            values.Select(v => (string?)v["individual"]).ToArray());
    }

    [TestCase("A3KMH1", RandomStructure.Individual)]
    [TestCase("A6NED2", RandomStructure.None)]
    public void OnlyTheTermsTheDataCanEstimateAreFittedAsLme4AndLmDo(string protein, RandomStructure structure)
    {
        var expected = LmmData.Rows("lmm_estimable_terms_results.tsv").Single(r => r["protein"] == protein);
        var problem = Problem(protein);

        var refused = NestedMixedModel.Fit(new[] { problem }, Names);
        Assert.That(refused.Status[0], Is.EqualTo(FeatureFitStatus.TooFewGroups), "the full model, as lme4");

        var fit = NestedMixedModel.Fit(new[] { problem }, Names, estimableRandomTermsOnly: true);

        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Structure[0], Is.EqualTo(structure));
        Assert.That(fit.SampleVariance[0], Is.NaN);
        if (structure == RandomStructure.Individual)
        {
            Close(fit.RemlCriterion[0], expected["reml"], Reml, $"{protein} REML");
            Close(fit.IndividualVariance[0], expected["var_individual"], Fit, $"{protein} var(individual)");
        }
        else
        {
            Assert.That(fit.IndividualVariance[0], Is.NaN);
            Assert.That(fit.RemlCriterion[0], Is.Not.NaN);
        }
        Close(fit.ResidualVariance[0], expected["sigma2"], Fit, $"{protein} sigma2");
        for (int k = 0; k < Levels.Length; k++)
        {
            var unit = Enumerable.Range(0, Names.Length).Select(j => j == k + 1 ? 1.0 : 0).ToArray();
            Close(fit.Coefficient(0, k + 1), expected[$"coef_{Levels[k]}"], Fit, $"{protein} condition {Levels[k]}");
            Close(Math.Sqrt(fit.UnscaledCovariance(0, k + 1, k + 1) * fit.ResidualVariance[0]), expected[$"se_{Levels[k]}"], Fit,
                $"{protein} se {Levels[k]}");
            Close(fit.SatterthwaiteDf(0, unit), expected[$"df_{Levels[k]}"], structure == RandomStructure.None ? 1e-9 : Df,
                $"{protein} df {Levels[k]}");
        }
    }

    /// <summary>
    /// The joint F df of the condition terms (moderation's variance df) is the residual df too for a least-squares fit, and
    /// the fit takes part in moderation like any other.
    /// </summary>
    [Test]
    public void ALeastSquaresFitTakesPartInModeration()
    {
        var problems = new[] { Problem("A3KMH1"), Problem("A6NED2"), Problem("A3KMH1") };
        var fit = NestedMixedModel.Fit(problems, Names, estimableRandomTermsOnly: true);
        var rows = new IReadOnlyList<double>[] { new[] { 0.0, 1, 0, 0 }, new[] { 0.0, 0, 1, 0 }, new[] { 0.0, 0, 0, 1 } };

        Assert.That(fit.JointSatterthwaiteDf(1, rows), Is.EqualTo(25).Within(1e-9));
        var moderated = NestedMixedModel.ModerateContrasts(fit, new[] { new ContrastWeights("1 vs 0.125", new[] { 0.0, 0, 0, 1 }) },
            rows);
        Assert.That(moderated[0].Estimate[1], Is.EqualTo(fit.Coefficient(1, 3)).Within(1e-12));
        Assert.That(moderated[0].PValue.All(double.IsFinite), Is.True);
    }

    /// <summary>
    /// Robust reweighting runs on a least-squares fit as on any other: Huber weights from its residuals, refitted to
    /// convergence; an outlier is down-weighted.
    /// </summary>
    [Test]
    public void RobustReweightingWorksWithoutARandomTerm()
    {
        var values = LmmData.Rows("lmm_estimable_terms_values.tsv").Where(r => r["protein"] == "A6NED2").ToList();
        var problem = Problem("A6NED2");
        var y = problem.Response.ToArray();
        y[0] += 3; // one outlying value
        var shifted = new NestedMixedModelProblem(y, problem.Common, problem.Nuisance, problem.Sample, problem.Individual);

        var fit = NestedMixedModel.Fit(new[] { shifted }, Names, robust: true, estimableRandomTermsOnly: true);

        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Structure[0], Is.EqualTo(RandomStructure.None));
        Assert.That(fit.RobustConverged[0], Is.True);
        Assert.That(fit.Weights(0)![0], Is.LessThan(0.5), "the outlier");
        Assert.That(fit.Weights(0)!.Skip(1).Count(w => w < 1), Is.LessThan(values.Count / 2));
    }

    private static void Close(double actual, string expectedText, double tolerance, string where)
    {
        double expected = LmmData.Parse(expectedText);
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Abs(expected))),
            $"{where}: {actual:R} vs {expected:R}");
    }
}
