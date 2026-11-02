using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using NUnit.Framework;
using StatisticalModels;

namespace Test.StatisticalModels;

/// <summary>
/// <see cref="NestedMixedModel"/> against frozen reference values from <c>ReferenceData/make_nested_lmm_fixtures.R</c>:
/// per protein, peptide-level log2 values ~ condition terms + peptide + random intercepts, fitted by
/// <c>lme4::lmer</c> (REML); Satterthwaite df from <c>lmerTest::contest</c>; and the MSstatsTMT-style moderation of
/// QuantProject GR-21 (<c>limma::squeezeVar(legacy = TRUE)</c> on σ² with the joint condition F-test's DenDF). Two
/// scenarios of 40 proteins: <c>nested</c> (individual and sample within individual) and <c>sample_only</c>. R is never
/// run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReferenceComparisonNestedLmmTests
{
    /// <summary>
    /// The REML criterion at the optimum: both optimizers reach the same minimum, so it is compared tightly.
    /// </summary>
    private const double Reml = 1e-9;

    /// <summary>
    /// Estimates, unscaled covariance and residual variance; limited by lme4's bobyqa stopping rule (it warned at rhoend
    /// 1e-12 on two proteins).
    /// </summary>
    private const double Fit = 1e-6;

    /// <summary>
    /// Satterthwaite df and everything downstream. lmerTest differentiates numerically (numDeriv: a 0.1 relative step for
    /// the Hessian, 1e-4 for the Jacobian, Richardson-extrapolated); this implementation reproduces those steps, and agrees
    /// far inside this tolerance. Near a variance component close to zero the df is very sensitive to where the optimum
    /// lies, which is why the REML criterion above is evaluated from the residuals, as lme4 does.
    /// </summary>
    private const double Df = 1e-4;

    [TestCase("nested")]
    [TestCase("sample_only")]
    public void FitAndModerationEqualReference(string scenario)
    {
        var data = LmmData.Read(scenario);
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames);
        var moderated = NestedMixedModel.ModerateContrasts(fit, data.Contrasts, data.ConditionRows);
        var expected = LmmData.Rows($"lmm_{scenario}_results.tsv");
        var unscaled = LmmData.Rows($"lmm_{scenario}_unscaled.tsv").ToDictionary(r => r["protein"]);

        var fitted = expected.Select(r => r["protein"]).Distinct().ToHashSet();
        for (int f = 0; f < data.Proteins.Count; f++)
            Assert.That(fit.Status[f] == FeatureFitStatus.Fitted, Is.EqualTo(fitted.Contains(data.Proteins[f])),
                $"{scenario} {data.Proteins[f]}: fitted here iff lme4 fitted it ({fit.Status[f]})");

        foreach (var row in expected)
        {
            int f = data.Proteins.IndexOf(row["protein"]);
            var r = moderated.Single(m => m.Name == row["contrast"]);
            var c = data.Contrasts.Single(x => x.Name == row["contrast"]).Weights;
            string at = $"{scenario} {row["protein"]} {row["contrast"]}";
            Close(fit.ResidualVariance[f], row["sigma2"], Fit, at + " sigma2");
            Close(fit.RemlCriterion[f], row["reml"], Reml, at + " REML");
            Variance(fit.IndividualVariance[f], row["var_individual"], fit.ResidualVariance[f], at + " var(individual)");
            Variance(fit.SampleVariance[f], row["var_sample"], fit.ResidualVariance[f], at + " var(sample)");
            Close(r.Estimate[f], row["estimate"], Fit, at + " estimate");
            Close(fit.SatterthwaiteDf(f, c), row["df_satterthwaite"], Df, at + " Satterthwaite df");
            Close(fit.JointSatterthwaiteDf(f, data.ConditionRows), row["variance_df"], Df, at + " joint DenDF");
            Close(r.PosteriorVariance[f], row["s2_post"], Df, at + " s2.post");
            Close(r.StandardError[f], row["se"], Df, at + " se");
            Close(r.DfTotal[f], row["df"], Df, at + " df");
            Close(r.T[f], row["t"], Df, at + " t");
            Close(r.PValue[f], row["p"], Df, at + " p");
            Close(r.BenjaminiHochbergAdjusted[f], row["adj_p"], Df, at + " adj.p");
            Close(r.ConfidenceLow[f], row["ci_low"], Df, at + " CI low");
            Close(r.ConfidenceHigh[f], row["ci_high"], Df, at + " CI high");
            Close(r.Prior.Df, row["df_prior"], Df, at + " df.prior");
            int p = data.CommonNames.Count;
            var u = unscaled[row["protein"]];
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                    Close(fit.UnscaledCovariance(f, i, j), u[$"X{i * p + j + 1}"], Fit, $"{at} unscaled[{i},{j}]");
        }
    }

    private static void Close(double actual, string expectedText, double tolerance, string where)
    {
        double expected = LmmData.Parse(expectedText);
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Abs(expected))),
            $"{where}: {actual:R} vs {expected:R}");
    }

    /// <summary>A variance component, compared absolutely on the residual variance's scale (lme4 can sit at or near 0).</summary>
    private static void Variance(double actual, string expectedText, double sigma2, string where)
    {
        double expected = LmmData.Parse(expectedText);
        if (double.IsNaN(expected)) { Assert.That(actual, Is.NaN, where); return; }
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(Df * sigma2), $"{where}: {actual:R} vs {expected:R}");
    }
}

/// <summary>Reads the nested-LMM fixtures into <see cref="NestedMixedModelProblem"/>s, building peptide columns as R did.</summary>
[ExcludeFromCodeCoverage]
internal sealed class LmmData
{
    public List<string> Proteins { get; } = new();
    public List<NestedMixedModelProblem> Problems { get; } = new();
    public List<string> CommonNames { get; } = new();
    public List<ContrastWeights> Contrasts { get; } = new();
    public List<IReadOnlyList<double>> ConditionRows { get; } = new();

    public static LmmData Read(string scenario)
    {
        var d = new LmmData();
        var samples = Rows($"lmm_{scenario}_samples.tsv");
        var header = Lines($"lmm_{scenario}_samples.tsv")[0];
        d.CommonNames.AddRange(header.SkipWhile(h => h != "intercept"));
        int p = d.CommonNames.Count;
        var sampleIds = samples.Select(s => s["sample"]).ToList();
        foreach (var c in Rows($"lmm_{scenario}_contrasts.tsv"))
            d.Contrasts.Add(new ContrastWeights(c["contrast"], d.CommonNames.Select(n => Parse(c[n])).ToArray()));
        for (int r = 1; r < p; r++) d.ConditionRows.Add(Enumerable.Range(0, p).Select(j => j == r ? 1.0 : 0.0).ToArray());

        foreach (var protein in Rows($"lmm_{scenario}_values.tsv").GroupBy(v => v["protein"]))
        {
            var observed = protein.Where(v => !double.IsNaN(Parse(v["y"]))).ToList();
            var peptides = observed.Select(v => v["peptide"]).Distinct().OrderBy(x => x, StringComparer.Ordinal).ToList();
            int n = observed.Count;
            var common = new double[n, p];
            var nuisance = new double[n, Math.Max(0, peptides.Count - 1)];
            var y = new double[n];
            var sample = new string[n];
            var individual = new string?[n];
            for (int i = 0; i < n; i++)
            {
                var v = observed[i];
                var s = samples[sampleIds.IndexOf(v["sample"])];
                for (int j = 0; j < p; j++) common[i, j] = Parse(s[d.CommonNames[j]]);
                int k = peptides.IndexOf(v["peptide"]);
                if (k > 0) nuisance[i, k - 1] = 1;
                y[i] = Parse(v["y"]);
                sample[i] = v["sample"];
                individual[i] = s["individual"] is "NA" or "" ? null : s["individual"];
            }
            d.Proteins.Add(protein.Key);
            d.Problems.Add(new NestedMixedModelProblem(y, common, peptides.Count > 1 ? nuisance : null, sample,
                individual.All(x => x is null) ? null : individual));
        }
        return d;
    }

    public static double Parse(string s) => s switch
    {
        "NaN" or "NA" => double.NaN,
        "Inf" => double.PositiveInfinity,
        "-Inf" => double.NegativeInfinity,
        _ => double.Parse(s, CultureInfo.InvariantCulture),
    };

    public static List<Dictionary<string, string>> Rows(string file)
    {
        var lines = Lines(file);
        return lines.Skip(1).Select(r => lines[0].Zip(r).ToDictionary(t => t.First, t => t.Second)).ToList();
    }

    public static List<string[]> Lines(string file) =>
        File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData", file))
            .Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
}
