using NUnit.Framework;
using StatisticalModels;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// The pooled moderation overloads (several fits, one variance prior) and <see cref="NestedMixedModelFit.NuisanceCoefficients"/>.
/// One part, and features split over parts in their original order, must give exactly a single fit's numbers; with fewer
/// than 3 usable features no prior is fitted (QuantProject STAT1 ruling 5), checked against lmerTest's and limma's
/// unmoderated values; the peptide coefficients are checked against lme4 (make_nuisance_coefficient_fixtures.R). R is never
/// run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PooledModerationTests
{
    /// <summary>Estimates and SEs against lme4, as <c>ReferenceComparisonNestedLmmTests</c> (bobyqa's stopping rule).</summary>
    private const double Fit = 1e-6;

    /// <summary>Satterthwaite df against lmerTest, as <c>ReferenceComparisonNestedLmmTests</c> (numerical derivatives).</summary>
    private const double Df = 1e-4;

    /// <summary>limma's lmFit outputs, as <c>ReferenceComparisonTests</c>.</summary>
    private const double Limma = 1e-8;

    private static void Close(double actual, double expected, double tolerance, string what) =>
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Abs(expected))),
            $"{what}: {actual:R} vs {expected:R}");

    /// <summary>Every per-feature array, the prior and the flag of two results are identical, NaN included.</summary>
    private static void Same(ModeratedContrast actual, ModeratedContrast expected, string what)
    {
        Assert.That(actual.Name, Is.EqualTo(expected.Name), what);
        Assert.That(actual.Prior.Df, Is.EqualTo(expected.Prior.Df), $"{what} prior df");
        Assert.That(actual.Prior.Scale, Is.EqualTo(expected.Prior.Scale), $"{what} prior scale");
        Assert.That(actual.ResidualDfDiffer, Is.EqualTo(expected.ResidualDfDiffer), $"{what} df differ");
        Assert.That(actual.Estimate, Is.EqualTo(expected.Estimate), $"{what} estimate");
        Assert.That(actual.StandardError, Is.EqualTo(expected.StandardError), $"{what} se");
        Assert.That(actual.T, Is.EqualTo(expected.T), $"{what} t");
        Assert.That(actual.DfTotal, Is.EqualTo(expected.DfTotal), $"{what} df");
        Assert.That(actual.PValue, Is.EqualTo(expected.PValue), $"{what} p");
        Assert.That(actual.BenjaminiHochbergAdjusted, Is.EqualTo(expected.BenjaminiHochbergAdjusted), $"{what} adjusted p");
        Assert.That(actual.PosteriorVariance, Is.EqualTo(expected.PosteriorVariance), $"{what} posterior variance");
        Assert.That(actual.ConfidenceLow, Is.EqualTo(expected.ConfidenceLow), $"{what} low");
        Assert.That(actual.ConfidenceHigh, Is.EqualTo(expected.ConfidenceHigh), $"{what} high");
    }

    private static IReadOnlyList<double> Concat(IEnumerable<ModeratedContrast> parts, Func<ModeratedContrast, IReadOnlyList<double>> values) =>
        parts.SelectMany(values).ToArray();

    // ---- NuisanceCoefficients ----

    [TestCase("nested")]
    [TestCase("sample_only")]
    public void NuisanceCoefficientsMatchLme4(string scenario)
    {
        var data = LmmData.Read(scenario);
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames);
        var expected = LmmData.Rows($"lmm_{scenario}_nuisance.tsv");
        Assert.That(expected, Is.Not.Empty);

        foreach (var protein in expected.GroupBy(r => r["protein"]))
        {
            int f = data.Proteins.IndexOf(protein.Key);
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), protein.Key);
            var peptides = data.Keys[f].Select(k => k.Peptide).Distinct().OrderBy(x => x, StringComparer.Ordinal).ToList();
            var mine = fit.NuisanceCoefficients(f)!;
            Assert.That(mine.Count, Is.EqualTo(peptides.Count - 1), $"{protein.Key}: one coefficient per non-reference peptide");
            foreach (var row in protein)
                Close(mine[peptides.IndexOf(row["peptide"]) - 1], LmmData.Parse(row["estimate"]), Fit, $"{scenario} {protein.Key} {row["peptide"]}");
        }
    }

    [Test]
    public void NuisanceCoefficientsAreNullWhenNotFittedAndEmptyWithoutNuisanceColumns()
    {
        var data = LmmData.Read("nested");
        var fit = NestedMixedModel.Fit(data.Problems.Take(2).ToList(), data.CommonNames);
        Assert.That(fit.Status[0], Is.Not.EqualTo(FeatureFitStatus.Fitted), "p001 has one peptide");
        Assert.That(fit.NuisanceCoefficients(0), Is.Null);

        // Two observations per sample and no nuisance column: fitted, with no nuisance coefficient.
        var y = new[] { 1.0, 1.2, 2.1, 1.9, 1.1, 0.8, 2.3, 2.0, 0.9, 1.3, 2.2, 1.7 };
        var common = new double[12, 2];
        var sample = new string[12];
        for (int i = 0; i < 12; i++)
        {
            common[i, 0] = 1;
            common[i, 1] = (i / 2) % 2;
            sample[i] = "s" + (i / 2).ToString(CultureInfo.InvariantCulture);
        }
        var plain = NestedMixedModel.Fit(new[] { new NestedMixedModelProblem(y, common, null, sample, null) }, new[] { "intercept", "b" });
        Assert.That(plain.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(plain.NuisanceCoefficients(0), Is.Empty);
    }

    // ---- pooled mixed-model moderation ----

    [TestCase("nested", "lmm", false)]
    [TestCase("sample_only", "lmm", false)]
    [TestCase("nested", "rlmm", true)]
    public void OneMixedModelPartGivesExactlyTheSingleFitResult(string scenario, string prefix, bool robust)
    {
        var data = LmmData.Read(scenario, prefix);
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames, robust: robust);
        var single = NestedMixedModel.ModerateContrasts(fit, data.Contrasts, data.ConditionRows);
        var pooled = NestedMixedModel.ModerateContrasts(new[] { new MixedModelModerationPart(fit, data.Contrasts, data.ConditionRows) });

        Assert.That(pooled, Has.Count.EqualTo(1));
        Assert.That(pooled[0], Has.Count.EqualTo(single.Count));
        for (int c = 0; c < single.Count; c++) Same(pooled[0][c], single[c], $"{prefix} {scenario} {single[c].Name}");
    }

    [Test]
    public void MixedModelFeaturesSplitOverTwoPartsGiveTheSameNumbersAsOneFit()
    {
        var data = LmmData.Read("nested");
        int half = data.Problems.Count / 2;
        var all = NestedMixedModel.Fit(data.Problems, data.CommonNames);
        var first = NestedMixedModel.Fit(data.Problems.Take(half).ToList(), data.CommonNames);
        var second = NestedMixedModel.Fit(data.Problems.Skip(half).ToList(), data.CommonNames);
        var single = NestedMixedModel.ModerateContrasts(all, data.Contrasts, data.ConditionRows);
        var pooled = NestedMixedModel.ModerateContrasts(new[]
        {
            new MixedModelModerationPart(first, data.Contrasts, data.ConditionRows),
            new MixedModelModerationPart(second, data.Contrasts, data.ConditionRows),
        });

        for (int c = 0; c < single.Count; c++)
        {
            var parts = new[] { pooled[0][c], pooled[1][c] };
            string what = single[c].Name;
            Assert.That(parts.Select(p => p.Prior.Df), Is.All.EqualTo(single[c].Prior.Df), what);
            Assert.That(Concat(parts, p => p.Prior.Scale), Is.EqualTo(single[c].Prior.Scale), what);
            Assert.That(Concat(parts, p => p.StandardError), Is.EqualTo(single[c].StandardError), what);
            Assert.That(Concat(parts, p => p.DfTotal), Is.EqualTo(single[c].DfTotal), $"{what}: the df cap is over both parts");
            Assert.That(Concat(parts, p => p.PValue), Is.EqualTo(single[c].PValue), what);
            Assert.That(Concat(parts, p => p.BenjaminiHochbergAdjusted), Is.EqualTo(single[c].BenjaminiHochbergAdjusted),
                $"{what}: Benjamini-Hochberg is over both parts");
        }
    }

    /// <summary>
    /// p001 (one peptide, not fitted), p002 and p003: two usable features, so no prior (d0 = 0), and each contrast keeps
    /// lmerTest's unmoderated SE and Satterthwaite df (capped at the two features' variance df, which do not bind here).
    /// The single-fit overload, given the same three features, would fit a prior from the two (it counts every feature).
    /// </summary>
    [Test]
    public void FewerThanThreeUsableMixedModelFeaturesFitNoPriorAndKeepLmerTestsUnmoderatedValues()
    {
        var data = LmmData.Read("nested");
        var fit = NestedMixedModel.Fit(data.Problems.Take(3).ToList(), data.CommonNames);
        var pooled = NestedMixedModel.ModerateContrasts(new[] { new MixedModelModerationPart(fit, data.Contrasts, data.ConditionRows) })[0];
        var expected = LmmData.Rows("lmm_nested_results.tsv");
        double cap = expected.Where(r => r["protein"] is "p002" or "p003").GroupBy(r => r["protein"])
            .Sum(g => LmmData.Parse(g.First()["variance_df"]));

        foreach (var c in pooled)
        {
            Assert.That(c.Prior.Df, Is.EqualTo(0), "no prior from fewer than 3 usable features");
            Assert.That(c.Estimate[0], Is.NaN, "p001 is not fitted");
            foreach (string protein in new[] { "p002", "p003" })
            {
                int f = data.Proteins.IndexOf(protein);
                var row = expected.Single(r => r["protein"] == protein && r["contrast"] == c.Name);
                Close(c.Estimate[f], LmmData.Parse(row["estimate"]), Fit, $"{protein} {c.Name} estimate");
                Close(c.StandardError[f], LmmData.Parse(row["se_unmoderated"]), Fit, $"{protein} {c.Name} se");
                Close(c.DfTotal[f], Math.Min(LmmData.Parse(row["df_satterthwaite"]), cap), Df, $"{protein} {c.Name} df");
                Close(c.PosteriorVariance[f], LmmData.Parse(row["sigma2"]), Fit, $"{protein} {c.Name} own variance");
            }
        }
    }

    [Test]
    public void OneUsableMixedModelFeatureKeepsItsOwnVarianceAndNoneGivesNaNWithoutThrowing()
    {
        var data = LmmData.Read("nested");
        var expected = LmmData.Rows("lmm_nested_results.tsv");

        var one = NestedMixedModel.Fit(data.Problems.Take(2).ToList(), data.CommonNames);
        foreach (var c in NestedMixedModel.ModerateContrasts(new[] { new MixedModelModerationPart(one, data.Contrasts, data.ConditionRows) })[0])
        {
            var row = expected.Single(r => r["protein"] == "p002" && r["contrast"] == c.Name);
            Assert.That(c.Prior.Df, Is.EqualTo(0));
            Close(c.StandardError[1], LmmData.Parse(row["se_unmoderated"]), Fit, $"p002 {c.Name} se");
        }

        var none = NestedMixedModel.Fit(data.Problems.Take(1).ToList(), data.CommonNames);
        var results = NestedMixedModel.ModerateContrasts(new[] { new MixedModelModerationPart(none, data.Contrasts, data.ConditionRows) })[0];
        Assert.That(results.SelectMany(c => c.Estimate.Concat(c.StandardError).Concat(c.PValue)), Is.All.NaN);
    }

    [Test]
    public void PartsMayTestDifferentContrastsAndOneWithNoneStillInformsThePrior()
    {
        var data = LmmData.Read("nested");
        int half = data.Problems.Count / 2;
        var first = NestedMixedModel.Fit(data.Problems.Take(half).ToList(), data.CommonNames);
        var second = NestedMixedModel.Fit(data.Problems.Skip(half).ToList(), data.CommonNames);
        var both = NestedMixedModel.ModerateContrasts(new[]
        {
            new MixedModelModerationPart(first, data.Contrasts, data.ConditionRows),
            new MixedModelModerationPart(second, data.Contrasts, data.ConditionRows),
        });
        var silent = NestedMixedModel.ModerateContrasts(new[]
        {
            new MixedModelModerationPart(first, data.Contrasts, data.ConditionRows),
            new MixedModelModerationPart(second, Array.Empty<ContrastWeights>(), data.ConditionRows),
        });

        Assert.That(silent[1], Is.Empty);
        for (int c = 0; c < data.Contrasts.Count; c++)
        {
            Assert.That(silent[0][c].Prior.Df, Is.EqualTo(both[0][c].Prior.Df), "the silent part's variances still enter the prior");
            Assert.That(silent[0][c].StandardError, Is.EqualTo(both[0][c].StandardError));
            Assert.That(silent[0][c].DfTotal, Is.EqualTo(both[0][c].DfTotal), "and the df cap");
        }
    }

    [Test]
    public void APooledMixedModelRequestThatCannotBeReadIsRefusedByPart()
    {
        var data = LmmData.Read("nested");
        var fit = NestedMixedModel.Fit(data.Problems.Take(4).ToList(), data.CommonNames);
        Assert.Throws<ArgumentException>(() => NestedMixedModel.ModerateContrasts(Array.Empty<MixedModelModerationPart>()));
        var wrong = new[] { new ContrastWeights("short", new[] { 0.0, 1.0 }) };
        var ex = Assert.Throws<ArgumentException>(() => NestedMixedModel.ModerateContrasts(new[]
        {
            new MixedModelModerationPart(fit, data.Contrasts, data.ConditionRows),
            new MixedModelModerationPart(fit, wrong, data.ConditionRows),
        }));
        Assert.That(ex!.Message, Does.Contain("part 1"));
        Assert.Throws<ArgumentException>(() => NestedMixedModel.ModerateContrasts(new[]
            { new MixedModelModerationPart(fit, data.Contrasts, Array.Empty<IReadOnlyList<double>>()) }));
        Assert.Throws<ArgumentOutOfRangeException>(() => NestedMixedModel.ModerateContrasts(new[]
            { new MixedModelModerationPart(fit, data.Contrasts, data.ConditionRows) }, confidenceLevel: 1));
    }

    // ---- pooled linear-model moderation ----

    private static (double[,] Responses, double[,] Design, string[] Names) LimmaData()
    {
        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData");
        var responses = File.ReadAllLines(Path.Combine(dir, "limma_responses.tsv")).Skip(1).Where(l => l.Length > 0)
            .Select(l => l.Split('\t').Select(LmmData.Parse).ToArray()).ToArray();
        var designLines = File.ReadAllLines(Path.Combine(dir, "limma_design.tsv")).Where(l => l.Length > 0).ToArray();
        var names = designLines[0].Split('\t');
        var design = designLines.Skip(1).Select(l => l.Split('\t').Select(LmmData.Parse).ToArray()).ToArray();
        var y = new double[responses.Length, responses[0].Length];
        for (int f = 0; f < responses.Length; f++)
            for (int s = 0; s < responses[f].Length; s++) y[f, s] = responses[f][s];
        var x = new double[design.Length, names.Length];
        for (int s = 0; s < design.Length; s++)
            for (int j = 0; j < names.Length; j++) x[s, j] = design[s][j];
        return (y, x, names);
    }

    private static double[,] Rows(double[,] y, int from, int count)
    {
        var part = new double[count, y.GetLength(1)];
        for (int f = 0; f < count; f++)
            for (int s = 0; s < y.GetLength(1); s++) part[f, s] = y[from + f, s];
        return part;
    }

    private static readonly ContrastWeights[] LimmaContrasts =
    {
        new("age", new[] { 0.0, 1.0, 0.0 }),
        new("age_minus_sex", new[] { 0.0, 1.0, -1.0 }),
    };

    [Test]
    public void OneLinearPartGivesExactlyTheSingleFitResult()
    {
        var (y, x, names) = LimmaData();
        var fit = LinearModel.Fit(y, x, names);
        var single = EmpiricalBayes.ModerateContrasts(fit, LimmaContrasts, trend: false);
        var pooled = EmpiricalBayes.ModerateContrasts(new[] { new LinearModelModerationPart(fit, LimmaContrasts) });

        for (int c = 0; c < single.Count; c++) Same(pooled[0][c], single[c], single[c].Name);
    }

    [Test]
    public void LinearFeaturesSplitOverTwoPartsGiveTheSameNumbersAsOneFit()
    {
        var (y, x, names) = LimmaData();
        int n = y.GetLength(0), half = n / 2;
        var single = EmpiricalBayes.ModerateContrasts(LinearModel.Fit(y, x, names), LimmaContrasts, trend: false);
        var pooled = EmpiricalBayes.ModerateContrasts(new[]
        {
            new LinearModelModerationPart(LinearModel.Fit(Rows(y, 0, half), x, names), LimmaContrasts),
            new LinearModelModerationPart(LinearModel.Fit(Rows(y, half, n - half), x, names), LimmaContrasts),
        });

        for (int c = 0; c < single.Count; c++)
        {
            var parts = new[] { pooled[0][c], pooled[1][c] };
            Assert.That(parts.Select(p => p.Prior.Df), Is.All.EqualTo(single[c].Prior.Df));
            Assert.That(Concat(parts, p => p.StandardError), Is.EqualTo(single[c].StandardError));
            Assert.That(Concat(parts, p => p.DfTotal), Is.EqualTo(single[c].DfTotal));
            Assert.That(Concat(parts, p => p.BenjaminiHochbergAdjusted), Is.EqualTo(single[c].BenjaminiHochbergAdjusted));
        }
    }

    /// <summary>
    /// Two fitted features: no prior, so each keeps its own residual variance and df, which is limma's
    /// <c>eBayes</c> given those two features alone: SE = R's sigma × stdev.unscaled, df = df.residual (limma_fit.tsv).
    /// </summary>
    [Test]
    public void TwoFittedLinearFeaturesFitNoPriorAndKeepLimmasOwnVariance()
    {
        var (y, x, names) = LimmaData();
        var age = new[] { LimmaContrasts[0] };
        var result = EmpiricalBayes.ModerateContrasts(new[] { new LinearModelModerationPart(LinearModel.Fit(Rows(y, 0, 2), x, names), age) })[0][0];
        var lines = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData", "limma_fit.tsv"))
            .Where(l => l.Length > 0).Select(l => l.Split('\t')).ToArray();
        int Col(string name) => Array.IndexOf(lines[0], name);

        Assert.That(result.Prior.Df, Is.EqualTo(0));
        for (int f = 0; f < 2; f++)
        {
            var r = lines[f + 1];
            double se = LmmData.Parse(r[Col("sigma")]) * LmmData.Parse(r[Col("unscaled_age")]);
            Close(result.Estimate[f], LmmData.Parse(r[Col("coef_age")]), Limma, $"feature {f} estimate");
            Close(result.StandardError[f], se, Limma, $"feature {f} se");
            Close(result.DfTotal[f], LmmData.Parse(r[Col("df_residual")]), Limma, $"feature {f} df");
        }
    }

    /// <summary>
    /// Every feature's residuals are the same vector, so every variance is equal, the variances spread less than sampling
    /// explains and the prior df is infinite: each feature's df is then the cap, the residual df summed over EVERY part's
    /// features (6 features x 4 df), not over its own part's 3.
    /// </summary>
    [Test]
    public void TheDfCapIsTheSumOverEveryPart()
    {
        var x = new double[6, 2];
        var noise = new[] { 0.3, -0.1, -0.2, 0.25, -0.05, -0.2 };
        for (int s = 0; s < 6; s++) { x[s, 0] = 1; x[s, 1] = s < 3 ? 0 : 1; }
        double[,] Features(double shift)
        {
            var y = new double[3, 6];
            for (int f = 0; f < 3; f++)
                for (int s = 0; s < 6; s++) y[f, s] = 10 + shift + f + x[s, 1] * (f - 1) + noise[s];
            return y;
        }
        var names = new[] { "intercept", "b" };
        var contrast = new[] { new ContrastWeights("b", new[] { 0.0, 1.0 }) };
        var pooled = EmpiricalBayes.ModerateContrasts(new[]
        {
            new LinearModelModerationPart(LinearModel.Fit(Features(0), x, names), contrast),
            new LinearModelModerationPart(LinearModel.Fit(Features(5), x, names), contrast),
        });

        Assert.That(pooled[0][0].Prior.Df, Is.EqualTo(double.PositiveInfinity));
        Assert.That(pooled.SelectMany(p => p[0].DfTotal), Is.All.EqualTo(24.0));
    }
}
