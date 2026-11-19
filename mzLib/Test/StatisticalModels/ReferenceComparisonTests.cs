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
/// Compares the StatisticalModels project against frozen outputs of the published implementations: limma's
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
        Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData");

    private static (string[] Header, List<string[]> Rows) ReadTsv(string name)
    {
        var lines = File.ReadAllLines(Path.Combine(DataDirectory, name)).Where(l => l.Length > 0).ToArray();
        return (lines[0].Split('\t'), lines.Skip(1).Select(l => l.Split('\t')).ToList());
    }

    private static double Parse(string s) => s switch
    {
        "Inf" => double.PositiveInfinity,                  // how R's sprintf writes an infinite df.prior
        "-Inf" => double.NegativeInfinity,
        _ => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture),
    };

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
            if (double.IsInfinity(expected[i]))
            {
                Assert.That(actual[i], Is.EqualTo(expected[i]), $"{what}[{i}]");
                continue;
            }
            // Relative, so a p-value of 1e-12 is held to the same number of digits as a coefficient.
            double scale = Math.Max(Math.Abs(expected[i]), 1e-300);
            Assert.That(Math.Abs(actual[i] - expected[i]) / scale, Is.LessThan(Tolerance),
                $"{what}[{i}]: {actual[i]:R} vs reference {expected[i]:R}");
        }
    }

    /// <param name="firstFeatures">Fit only the first this many features (all when null).</param>
    /// <param name="responsesFile">The responses fixture, one row per feature.</param>
    /// <param name="designFile">The design fixture, one row per sample.</param>
    private static LinearModelFit FitReferenceData(int? firstFeatures = null, string responsesFile = "limma_responses.tsv",
        string designFile = "limma_design.tsv")
    {
        var (_, responseRows) = ReadTsv(responsesFile);
        var (designHeader, designRows) = ReadTsv(designFile);
        int features = firstFeatures ?? responseRows.Count, samples = designRows.Count, p = designHeader.Length;
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

    /// <summary>
    /// With few features limma's default trend spline is smaller (1 + (n ≥ 3) + (n ≥ 6) + (n ≥ 30) basis
    /// functions), so the first n features of the fixture are moderated at each step of that rule and
    /// compared with <c>eBayes(trend = TRUE, legacy = TRUE)</c> on the same rows (make_small_n_fixtures.R).
    /// At n = 2 limma fits no prior at all (d0 = 0), so both features keep their own variance.
    /// </summary>
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(5)]
    [TestCase(6)]
    [TestCase(20)]
    [TestCase(29)]
    [TestCase(30)]
    public void TrendOnFewFeaturesMatchesLimmaEBayesLegacy(int n)
    {
        var test = EmpiricalBayes.Moderate(FitReferenceData(n), "age_decades", trend: true);
        const string file = "limma_ebayes_trend_small_n.tsv";
        var (header, rows) = ReadTsv(file);
        var mine = rows.Where(r => r[0] == n.ToString(CultureInfo.InvariantCulture)).ToList();
        double[] Expected(string column) => mine.Select(r => Parse(r[Array.IndexOf(header, column)])).ToArray();
        var (priorHeader, priorRows) = ReadTsv("limma_ebayes_trend_small_n_prior.tsv");
        double df0 = Parse(priorRows.Single(r => r[0] == n.ToString(CultureInfo.InvariantCulture))[Array.IndexOf(priorHeader, "df_prior")]);

        Assert.That(test.Prior.SplineBasisCount, Is.EqualTo(EmpiricalBayes.TrendBasisCount(n, n, null)));
        AssertClose(new[] { test.Prior.Df }, new[] { df0 }, "df.prior");
        AssertClose(test.Prior.Scale, Expected("s2_prior"), "s2.prior");
        AssertClose(test.PosteriorVariance, Expected("s2_post"), "s2.post");
        AssertClose(test.DfTotal, Expected("df_total"), "df.total");
        AssertClose(test.T, Expected("t_age"), "t");
        AssertClose(test.PValue, Expected("p_age"), "p.value");
    }

    /// <summary>
    /// When the feature variances spread no more than sampling explains, the prior's degrees of freedom are infinite
    /// and every feature takes the prior variance. Without a trend, limma (since 3.32.0, 2017) sets that variance to
    /// the mean of the features' variances; its earlier versions returned exp(mean(e)), the limit of the finite-df
    /// formula, which is larger (here 0.0419 against 0.0392). make_infinite_prior_fixtures.R; R is not run here.
    /// </summary>
    [Test]
    public void AnInfinitePriorDfTakesTheMeanVarianceAsLimmaDoes()
    {
        var fit = FitReferenceData(null, "limma_infinite_prior_responses.tsv", "limma_infinite_prior_design.tsv");
        var test = EmpiricalBayes.Moderate(fit, "group", trend: false);
        const string file = "limma_infinite_prior_ebayes.tsv";

        Assert.That(test.Prior.Df, Is.EqualTo(double.PositiveInfinity));
        AssertClose(test.Prior.Scale, Column(file, "s2_prior"), "s2.prior");
        AssertClose(test.PosteriorVariance, Column(file, "s2_post"), "s2.post");
        AssertClose(test.DfTotal, Column(file, "df_total"), "df.total");
        AssertClose(test.T, Column(file, "t_group"), "t");
        AssertClose(test.PValue, Column(file, "p_group"), "p.value");
    }

    /// <summary>
    /// The same rule on the prior alone: three variances 0.10, 0.11 and 0.105 on 4 df spread less than sampling
    /// explains, so limma 3.68.5's squeezeVar(legacy = TRUE) gives df.prior = Inf and var.prior = 0.105, their mean;
    /// the pre-2017 limit exp(mean(e)) is 0.1375, 31% higher. With a trend the limit stays exp(fitted), as limma's.
    /// </summary>
    [Test]
    public void AnInfinitePriorDfWithoutATrendIsTheMeanOfTheVariances()
    {
        var prior = EmpiricalBayes.FitPrior(new[] { 0.10, 0.11, 0.105 }, new[] { 4.0, 4.0, 4.0 });

        Assert.That(prior.Df, Is.EqualTo(double.PositiveInfinity));
        AssertClose(prior.Scale, new[] { 0.105, 0.105, 0.105 }, "s2.prior");
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

    /// <summary>metafor <c>rma(method = "DL", test = "knha")</c> at 95% and 90% (make_hartung_knapp_fixtures.R).</summary>
    [Test]
    public void HartungKnappMatchesMetafor()
    {
        var (inHeader, inRows) = ReadTsv("hk_inputs.tsv");
        var (outHeader, outRows) = ReadTsv("hk_results.tsv");
        int c = Array.IndexOf(inHeader, "case"), yi = Array.IndexOf(inHeader, "yi"), sei = Array.IndexOf(inHeader, "sei");
        double Out(string[] row, string col) => Parse(row[Array.IndexOf(outHeader, col)]);
        Assert.That(outRows, Has.Count.EqualTo(14));
        foreach (var expected in outRows)
        {
            string name = expected[Array.IndexOf(outHeader, "case")];
            double level = Out(expected, "level");
            var studies = inRows.Where(r => r[c] == name).ToArray();
            var r = RandomEffectsMeta.Pool(studies.Select(s => Parse(s[yi])).ToArray(),
                studies.Select(s => Parse(s[sei])).ToArray(), level, RandomEffectsInference.HartungKnapp);
            string what = $"{name} at {level}";
            Assert.That(r.Studies, Is.EqualTo((int)Out(expected, "k")), what);
            Assert.That(r.DegreesOfFreedom, Is.EqualTo(Out(expected, "df")), what);
            AssertClose(new[] { r.Estimate, r.StandardError, r.ConfidenceLow, r.ConfidenceHigh, r.PValue },
                new[] { Out(expected, "estimate"), Out(expected, "se"), Out(expected, "ci_lb"), Out(expected, "ci_ub"),
                        Out(expected, "pval") }, what);
            Assert.That(r.Tau2, Is.EqualTo(Out(expected, "tau2")).Within(1e-12), $"{what} tau2");
        }
    }

    /// <summary>
    /// metafor <c>rma(method = "DL", test = "adhoc")</c> at 95% and 90% (make_hartung_knapp_adhoc_fixtures.R): the
    /// Hartung–Knapp standard error never below the DerSimonian–Laird one. The fixture also carries metafor's plain
    /// <c>knha</c> and <c>z</c> standard errors, so each case is checked to be the larger of the two, and the cases where
    /// the correction binds (<c>agreeing_three</c>, <c>equal_estimates</c>, <c>near_identical_two</c>) are checked to
    /// differ from plain Hartung–Knapp.
    /// </summary>
    [Test]
    public void HartungKnappAdHocMatchesMetafor()
    {
        var (inHeader, inRows) = ReadTsv("hk_adhoc_inputs.tsv");
        var (outHeader, outRows) = ReadTsv("hk_adhoc_results.tsv");
        int c = Array.IndexOf(inHeader, "case"), yi = Array.IndexOf(inHeader, "yi"), sei = Array.IndexOf(inHeader, "sei");
        double Out(string[] row, string col) => Parse(row[Array.IndexOf(outHeader, col)]);
        Assert.That(outRows, Has.Count.EqualTo(18));
        int bound = 0;
        foreach (var expected in outRows)
        {
            string name = expected[Array.IndexOf(outHeader, "case")];
            double level = Out(expected, "level");
            var studies = inRows.Where(r => r[c] == name).ToArray();
            double[] b = studies.Select(s => Parse(s[yi])).ToArray(), se = studies.Select(s => Parse(s[sei])).ToArray();
            var r = RandomEffectsMeta.Pool(b, se, level, RandomEffectsInference.HartungKnappAdHoc);
            var plain = RandomEffectsMeta.Pool(b, se, level, RandomEffectsInference.HartungKnapp);
            string what = $"{name} at {level}";
            Assert.That(r.Studies, Is.EqualTo((int)Out(expected, "k")), what);
            Assert.That(r.DegreesOfFreedom, Is.EqualTo(Out(expected, "df")), what);
            AssertClose(new[] { r.Estimate, r.StandardError, r.ConfidenceLow, r.ConfidenceHigh, r.PValue },
                new[] { Out(expected, "estimate"), Out(expected, "se"), Out(expected, "ci_lb"), Out(expected, "ci_ub"),
                        Out(expected, "pval") }, what);
            Assert.That(r.Tau2, Is.EqualTo(Out(expected, "tau2")).Within(1e-12), $"{what} tau2");
            Assert.That(Out(expected, "se"), Is.EqualTo(Math.Max(Out(expected, "se_knha"), Out(expected, "se_dl"))).Within(1e-15),
                $"{what}: metafor's adhoc se is the larger of knha and DL");
            if (Out(expected, "se_dl") > Out(expected, "se_knha") * (1 + 1e-9))
            {
                bound++;
                Assert.That(r.StandardError, Is.GreaterThan(plain.StandardError), $"{what}: the correction binds");
            }
        }
        Assert.That(bound, Is.EqualTo(6), "agreeing_three, equal_estimates and near_identical_two, at two levels");
    }

    /// <summary>The partial conjunction closed form with base R's <c>pchisq</c> (make_partial_conjunction_fixtures.R).</summary>
    [Test]
    public void FisherPartialConjunctionMatchesPchisq()
    {
        var (inHeader, inRows) = ReadTsv("partial_conjunction_inputs.tsv");
        var (outHeader, outRows) = ReadTsv("partial_conjunction_results.tsv");
        double Out(string[] row, string col) => Parse(row[Array.IndexOf(outHeader, col)]);
        foreach (var expected in outRows)
        {
            string name = expected[Array.IndexOf(outHeader, "case")];
            var p = inRows.Where(r => r[Array.IndexOf(inHeader, "case")] == name)
                .Select(r => Parse(r[Array.IndexOf(inHeader, "p")])).ToArray();
            int u = (int)Out(expected, "u");
            var r = PValueCombination.FisherPartialConjunction(p, u);
            Assert.That(r.Studies, Is.EqualTo((int)Out(expected, "n")), name);
            AssertClose(new[] { r.Statistic, r.PValue }, new[] { Out(expected, "statistic"), Out(expected, "pval") }, $"{name} u={u}");
        }
    }

    /// <summary>
    /// Families selected by Benjamini–Hochberg on one p each, then each selected family adjusted within itself and compared
    /// at q·R/m, against base R's <c>p.adjust</c> (make_partial_conjunction_fixtures.R).
    /// </summary>
    [Test]
    public void BenjaminiBogomolovMatchesPAdjust()
    {
        var (inHeader, inRows) = ReadTsv("bogomolov_inputs.tsv");
        var (outHeader, outRows) = ReadTsv("bogomolov_results.tsv");
        int fc = Array.IndexOf(inHeader, "family");
        var families = inRows.GroupBy(r => r[fc]).OrderBy(g => int.Parse(g.Key, CultureInfo.InvariantCulture)).ToList();
        var familyP = families.Select(g => Parse(g.First()[Array.IndexOf(inHeader, "family_p")])).ToArray();
        var within = families
            .Select(g => (IReadOnlyList<double>)g.Select(r => Parse(r[Array.IndexOf(inHeader, "p")])).ToArray()).ToArray();
        var selected = MultipleTesting.BenjaminiHochberg(familyP).Select(q => q <= 0.05).ToArray();

        var adjustment = MultipleTesting.BenjaminiBogomolov(within, selected);

        Assert.That(outRows, Has.Count.EqualTo(within.Sum(f => f.Count)));
        foreach (var expected in outRows)
        {
            int f = int.Parse(expected[Array.IndexOf(outHeader, "family")], CultureInfo.InvariantCulture) - 1;
            int i = int.Parse(expected[Array.IndexOf(outHeader, "index")], CultureInfo.InvariantCulture) - 1;
            double Out(string col) => Parse(expected[Array.IndexOf(outHeader, col)]);
            Assert.That(selected[f], Is.EqualTo(Out("selected") == 1), $"family {f + 1} selected");
            AssertClose(new[] { adjustment.Adjusted[f][i], adjustment.Level }, new[] { Out("adjusted"), Out("level") },
                $"family {f + 1} p {i + 1}");
        }
    }
}
