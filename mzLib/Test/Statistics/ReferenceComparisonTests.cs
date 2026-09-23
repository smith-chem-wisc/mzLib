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
/// Compares the Statistics project against frozen outputs of the published implementations: limma's
/// lmFit and eBayes (legacy = TRUE, the method-of-moments prior this project implements) and metafor's
/// DerSimonian-Laird rma. The outputs, the R script that made them and its versions are in
/// ReferenceData (PROVENANCE.txt). R is never run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReferenceComparisonTests
{
    private const double Tolerance = 1e-8;

    private static string DataDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "Statistics", "ReferenceData");

    private static (string[] Header, List<string[]> Rows) ReadTsv(string name)
    {
        var lines = File.ReadAllLines(Path.Combine(DataDirectory, name)).Where(l => l.Length > 0).ToArray();
        return (lines[0].Split('\t'), lines.Skip(1).Select(l => l.Split('\t')).ToList());
    }

    private static double Parse(string s) => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static double[] Column(string name, string column)
    {
        var (header, rows) = ReadTsv(name);
        int j = Array.IndexOf(header, column);
        Assert.That(j, Is.GreaterThanOrEqualTo(0), $"{name} has no column {column}");
        return rows.Select(r => Parse(r[j])).ToArray();
    }

    private static void AssertClose(IReadOnlyList<double> actual, IReadOnlyList<double> expected, string what)
    {
        Assert.That(actual.Count, Is.EqualTo(expected.Count), what);
        for (int i = 0; i < expected.Count; i++)
        {
            if (double.IsNaN(expected[i]))
            {
                Assert.That(actual[i], Is.NaN, $"{what}[{i}]");
                continue;
            }
            // Relative, so a p-value of 1e-12 is held to the same number of digits as a coefficient.
            double scale = Math.Max(Math.Abs(expected[i]), 1e-300);
            Assert.That(Math.Abs(actual[i] - expected[i]) / scale, Is.LessThan(Tolerance),
                $"{what}[{i}]: {actual[i]:R} vs reference {expected[i]:R}");
        }
    }

    private static LinearModelFit FitReferenceData()
    {
        var (_, responseRows) = ReadTsv("limma_responses.tsv");
        var (designHeader, designRows) = ReadTsv("limma_design.tsv");
        int features = responseRows.Count, samples = designRows.Count, p = designHeader.Length;
        var responses = new double[features, samples];
        for (int f = 0; f < features; f++)
            for (int s = 0; s < samples; s++)
                responses[f, s] = Parse(responseRows[f][s]);
        var design = new double[samples, p];
        for (int s = 0; s < samples; s++)
            for (int j = 0; j < p; j++)
                design[s, j] = Parse(designRows[s][j]);
        return LinearModel.Fit(responses, design, designHeader);
    }

    [Test]
    public void LinearModelMatchesLimmaLmFit()
    {
        var fit = FitReferenceData();
        int n = fit.FeatureCount;
        const string file = "limma_fit.tsv";
        Assert.That(fit.Status.All(s => s == FeatureFitStatus.Fitted), "the fixture keeps every feature fittable");
        AssertClose(Enumerable.Range(0, n).Select(f => fit.Coefficient(f, 0)).ToArray(), Column(file, "coef_intercept"), "intercept");
        AssertClose(Enumerable.Range(0, n).Select(f => fit.Coefficient(f, 1)).ToArray(), Column(file, "coef_age"), "age_decades");
        AssertClose(Enumerable.Range(0, n).Select(f => fit.Coefficient(f, 2)).ToArray(), Column(file, "coef_sex"), "sex");
        AssertClose(Enumerable.Range(0, n).Select(f => fit.StdevUnscaled(f, 1)).ToArray(), Column(file, "unscaled_age"), "stdev.unscaled");
        AssertClose(fit.Sigma, Column(file, "sigma"), "sigma");
        AssertClose(fit.DfResidual.Select(d => (double)d).ToArray(), Column(file, "df_residual"), "df.residual");
        AssertClose(fit.AverageResponse, Column(file, "amean"), "Amean");
    }

    [TestCase(false, "notrend")]
    [TestCase(true, "trend")]
    public void ModerationMatchesLimmaEBayesLegacy(bool trend, string tag)
    {
        var test = EmpiricalBayes.Moderate(FitReferenceData(), "age_decades", trend);
        // The fixture omits missing values, so this is exactly the input on which default limma (>= 3.61)
        // would use its newer estimator; the match below holds for legacy = TRUE only.
        Assert.That(test.ResidualDfDiffer, Is.True);
        string file = $"limma_ebayes_{tag}.tsv";
        AssertClose(new[] { test.Prior.Df }, Column($"limma_ebayes_{tag}_prior.tsv", "df_prior"), "df.prior");
        AssertClose(test.Prior.Scale, Column(file, "s2_prior"), "s2.prior");
        AssertClose(test.PosteriorVariance, Column(file, "s2_post"), "s2.post");
        AssertClose(test.DfTotal, Column(file, "df_total"), "df.total");
        AssertClose(test.T, Column(file, "t_age"), "t");
        AssertClose(test.PValue, Column(file, "p_age"), "p.value");
        AssertClose(test.BenjaminiHochbergAdjusted, Column(file, "bh_age"), "BH");
    }

    [Test]
    public void DerSimonianLairdMatchesMetafor()
    {
        var (inHeader, inRows) = ReadTsv("dl_inputs.tsv");
        var (outHeader, outRows) = ReadTsv("dl_results.tsv");
        int c = Array.IndexOf(inHeader, "case"), yi = Array.IndexOf(inHeader, "yi"), sei = Array.IndexOf(inHeader, "sei");
        double Out(string[] row, string col) => Parse(row[Array.IndexOf(outHeader, col)]);
        foreach (var expected in outRows)
        {
            string name = expected[Array.IndexOf(outHeader, "case")];
            var studies = inRows.Where(r => r[c] == name).ToArray();
            var r = RandomEffectsMeta.Pool(studies.Select(s => Parse(s[yi])).ToArray(),
                                           studies.Select(s => Parse(s[sei])).ToArray());
            AssertClose(new[] { r.Estimate, r.StandardError, r.Q, r.PValue, r.ConfidenceLow, r.ConfidenceHigh },
                new[] { Out(expected, "estimate"), Out(expected, "se"), Out(expected, "q"), Out(expected, "pval"),
                        Out(expected, "ci_lb"), Out(expected, "ci_ub") }, name);
            // tau² and I² are exactly 0 when Q <= k - 1, so compare them absolutely.
            Assert.That(r.Tau2, Is.EqualTo(Out(expected, "tau2")).Within(1e-12), $"{name} tau2");
            Assert.That(r.ISquared, Is.EqualTo(Out(expected, "i2")).Within(1e-12), $"{name} I2");
        }
    }
}
