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
/// Moderated contrasts with confidence intervals (<see cref="EmpiricalBayes.ModerateContrasts"/>) against frozen
/// reference values from <c>ReferenceData/make_contrast_fixtures.R</c>. Three contrasts (C vs B, the age slope,
/// the mean of B and C vs A) on 300 simulated features:
/// <list type="bullet">
/// <item>complete data: limma's <c>lmFit</c> + <c>contrasts.fit</c> + <c>eBayes(legacy = TRUE)</c> +
/// <c>topTable(confint = 0.95)</c>, without and with the intensity trend;</item>
/// <item>about 10% missing values: the exact per-feature unscaled SD of each contrast with limma's
/// <c>squeezeVar(legacy = TRUE)</c> prior. limma's own <c>contrasts.fit</c> approximates that SD when missing values
/// change how a feature's coefficients correlate, so it is not the reference there.</item>
/// </list>
/// R is never run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ContrastReferenceTests
{
    private const double Tolerance = 1e-8;

    [TestCase("complete", false)]
    [TestCase("complete", true)]
    [TestCase("na", false)]
    [TestCase("na", true)]
    public void ModeratedContrastsEqualReference(string data, bool trend)
    {
        var (design, names) = ReadDesign();
        var responses = ReadResponses(data);
        var contrasts = ReadWeights();
        var fit = LinearModel.Fit(responses, design, names);
        var results = EmpiricalBayes.ModerateContrasts(fit, contrasts, trend);
        var expected = Read($"limma_contrast_{data}_{(trend ? "trend" : "none")}.tsv");

        Assert.That(results.Select(r => r.Name), Is.EqualTo(contrasts.Select(c => c.Name)));
        int checkedRows = 0;
        foreach (var row in expected)
        {
            var r = results.Single(x => x.Name == row["contrast"]);
            int f = int.Parse(row["feature"][1..], CultureInfo.InvariantCulture) - 1;
            string where = $"{data}/{(trend ? "trend" : "none")} {row["contrast"]} {row["feature"]}";
            Close(r.Estimate[f], row["estimate"], where + " estimate");
            Close(r.StandardError[f], row["se"], where + " se");
            Close(r.T[f], row["t"], where + " t");
            Close(r.DfTotal[f], row["df_total"], where + " df.total");
            Close(r.PValue[f], row["p"], where + " p");
            Close(r.BenjaminiHochbergAdjusted[f], row["adj_p"], where + " adj.p");
            Close(r.ConfidenceLow[f], row["ci_low"], where + " CI.L");
            Close(r.ConfidenceHigh[f], row["ci_high"], where + " CI.R");
            Close(r.PosteriorVariance[f], row["s2_post"], where + " s2.post");
            Close(r.Prior.Scale[f], row["s2_prior"], where + " s2.prior");
            Close(r.Prior.Df, row["df_prior"], where + " df.prior");
            checkedRows++;
        }
        Assert.That(checkedRows, Is.EqualTo(3 * 300));
        Assert.That(results.All(r => r.ConfidenceLevel == 0.95));
    }

    private static void Close(double actual, string expectedText, string where)
    {
        double expected = Parse(expectedText);
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(Tolerance * Math.Max(1, Math.Abs(expected))),
            $"{where}: {actual:R} vs {expected:R}");
    }

    private static double Parse(string s) => s switch
    {
        "NaN" => double.NaN,
        "Inf" => double.PositiveInfinity,
        "-Inf" => double.NegativeInfinity,
        _ => double.Parse(s, CultureInfo.InvariantCulture),
    };

    private static (double[,] Design, string[] Names) ReadDesign()
    {
        var lines = Lines("limma_contrast_design.tsv");
        var names = lines[0];
        var design = new double[lines.Count - 1, names.Length];
        for (int i = 1; i < lines.Count; i++)
            for (int j = 0; j < names.Length; j++) design[i - 1, j] = Parse(lines[i][j]);
        return (design, names);
    }

    private static double[,] ReadResponses(string data)
    {
        var lines = Lines($"limma_contrast_responses_{data}.tsv");
        int n = lines[0].Length - 1;
        var y = new double[lines.Count - 1, n];
        for (int i = 1; i < lines.Count; i++)
            for (int j = 0; j < n; j++) y[i - 1, j] = Parse(lines[i][j + 1]);
        return y;
    }

    private static List<ContrastWeights> ReadWeights() =>
        Lines("limma_contrast_weights.tsv").Skip(1)
            .Select(r => new ContrastWeights(r[0], r.Skip(1).Select(Parse).ToArray())).ToList();

    private static List<Dictionary<string, string>> Read(string file)
    {
        var lines = Lines(file);
        return lines.Skip(1).Select(r => lines[0].Zip(r).ToDictionary(t => t.First, t => t.Second)).ToList();
    }

    private static List<string[]> Lines(string file) =>
        File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData", file))
            .Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
}
