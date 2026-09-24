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
/// Compares LogisticRegression, MixedModel, SpearmanCorrelation and PValueCombination against frozen
/// outputs of R's glm, nlme::lme, cor.test and pchisq/pnorm. The outputs, the script that made them and
/// its versions are in ReferenceData (make_glm_lmm_fixtures.R, PROVENANCE_glm_lmm.txt). R is never run
/// by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReferenceComparisonGlmMixedTests
{
    /// <summary>Closed forms and IRLS to the MLE agree with R to this relative tolerance.</summary>
    private const double Tight = 1e-8;

    /// <summary>
    /// nlme optimizes the variance parameters numerically, on a log scale and with its own stopping rule, so
    /// agreement is limited by nlme's precision rather than by the method.
    /// </summary>
    private const double Lme = 1e-5;

    private static string DataDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "Statistics", "ReferenceData");

    private static (string[] Header, List<string[]> Rows) ReadTsv(string name)
    {
        var lines = File.ReadAllLines(Path.Combine(DataDirectory, name)).Where(l => l.Length > 0).ToArray();
        return (lines[0].Split('\t'), lines.Skip(1).Select(l => l.Split('\t')).ToList());
    }

    private static double Parse(string s) => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static double[,] Matrix(string name)
    {
        var (header, rows) = ReadTsv(name);
        var m = new double[rows.Count, header.Length];
        for (int i = 0; i < rows.Count; i++)
            for (int j = 0; j < header.Length; j++)
                m[i, j] = Parse(rows[i][j]);
        return m;
    }

    private static void Close(double actual, double expected, double tolerance, string what)
    {
        double scale = Math.Max(Math.Abs(expected), 1e-300);
        Assert.That(Math.Abs(actual - expected) / scale, Is.LessThan(tolerance), $"{what}: {actual:R} vs reference {expected:R}");
    }

    [Test]
    public void LogisticRegressionMatchesGlm()
    {
        var fit = LogisticRegression.Fit(Matrix("glm_responses.tsv"), Matrix("glm_design.tsv"));
        var (h, rows) = ReadTsv("glm_fit.tsv");
        double Col(string[] r, string c) => Parse(r[Array.IndexOf(h, c)]);
        int fitted = 0;
        for (int f = 0; f < rows.Count; f++)
        {
            string status = rows[f][Array.IndexOf(h, "status")];
            if (status == "separated") { Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Separated), $"feature {f}"); continue; }
            if (status == "not_fitted") { Assert.That(fit.Status[f], Is.Not.EqualTo(FeatureFitStatus.Fitted), $"feature {f}"); continue; }
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), $"feature {f}");
            fitted++;
            for (int j = 0; j < 3; j++)
            {
                Close(fit.Coefficient(f, j), Col(rows[f], $"b{j}"), Tight, $"feature {f} coefficient {j}");
                // summary.glm takes SEs from the weights of the LAST IRLS step, i.e. at the previous
                // iterate; mzLib evaluates the information at the estimate itself. With both at the MLE the
                // coefficients agree to 1e-8 and the SEs to about 1e-7.
                Close(fit.StandardError(f, j), Col(rows[f], $"se{j}"), 1e-6, $"feature {f} SE {j}");
            }
            Close(fit.PValue(f, 1), Col(rows[f], "p1"), 1e-6, $"feature {f} p");
            Close(fit.Deviance[f], Col(rows[f], "deviance"), Tight, $"feature {f} deviance");
        }
        Assert.That(fitted, Is.GreaterThan(rows.Count / 2), "most fixture features are fittable");
    }

    [TestCase("REML", VarianceEstimator.Reml)]
    [TestCase("ML", VarianceEstimator.MaximumLikelihood)]
    public void MixedModelMatchesNlme(string method, VarianceEstimator estimator)
    {
        var groups = ReadTsv("lmm_groups.tsv").Rows.Select(r => r[0]).ToArray();
        var fit = MixedModel.Fit(Matrix("lmm_responses.tsv"), Matrix("lmm_design.tsv"), groups, estimator: estimator);
        var (h, all) = ReadTsv("lmm_fit.tsv");
        var rows = all.Where(r => r[Array.IndexOf(h, "method")] == method).ToList();
        double Col(string[] r, string c) => Parse(r[Array.IndexOf(h, c)]);
        int compared = 0;
        foreach (var r in rows)
        {
            int f = (int)Col(r, "feature");
            string status = r[Array.IndexOf(h, "status")];
            if (status == "r_error") continue;
            if (status == "not_fitted") { Assert.That(fit.Status[f], Is.Not.EqualTo(FeatureFitStatus.Fitted), $"feature {f}"); continue; }
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), $"feature {f}");
            compared++;
            double sigma2 = Col(r, "sigma2");
            for (int j = 0; j < 3; j++)
            {
                Close(fit.Coefficient(f, j), Col(r, $"b{j}"), Lme, $"{method} feature {f} coefficient {j}");
                Close(fit.StandardError(f, j), Col(r, $"se{j}"), Lme, $"{method} feature {f} SE {j}");
                Assert.That(fit.DegreesOfFreedom(f, j), Is.EqualTo(Col(r, $"df{j}")), $"{method} feature {f} df {j}");
            }
            // A very small p-value magnifies a tiny difference in t, so p is compared on the log scale.
            Close(Math.Log(fit.PValue(f, 1)), Math.Log(Col(r, "p1")), Lme, $"{method} feature {f} log p(age)");
            Close(Math.Log(fit.PValue(f, 2)), Math.Log(Col(r, "p2")), Lme, $"{method} feature {f} log p(wp)");
            Close(fit.ResidualVariance[f], sigma2, Lme, $"{method} feature {f} sigma2");
            // nlme cannot reach τ² = 0 (it works on log τ), so compare on the scale of σ² rather than relatively.
            Assert.That(Math.Abs(fit.GroupVariance[f] - Col(r, "tau2")), Is.LessThan(Lme * sigma2),
                $"{method} feature {f} tau2: {fit.GroupVariance[f]:R} vs {Col(r, "tau2"):R}");
            Close(fit.LogLikelihood[f], Col(r, "loglik"), Lme, $"{method} feature {f} logLik");
        }
        Assert.That(compared, Is.GreaterThan(rows.Count / 2), "most fixture features are fittable");
    }

    [Test]
    public void SpearmanMatchesCorTest()
    {
        var (ih, inputs) = ReadTsv("spearman_inputs.tsv");
        var (oh, outputs) = ReadTsv("spearman_results.tsv");
        foreach (var o in outputs)
        {
            string name = o[Array.IndexOf(oh, "case")];
            var rows = inputs.Where(r => r[Array.IndexOf(ih, "case")] == name).ToArray();
            var r = SpearmanCorrelation.Correlate(rows.Select(x => Parse(x[Array.IndexOf(ih, "x")])).ToArray(),
                                                  rows.Select(x => Parse(x[Array.IndexOf(ih, "y")])).ToArray());
            string method = o[Array.IndexOf(oh, "method")];
            Assert.That(r.Method, Is.EqualTo(method == "exact" ? SpearmanPValueMethod.Exact : SpearmanPValueMethod.Asymptotic), name);
            Close(r.Rho, Parse(o[Array.IndexOf(oh, "rho")]), Tight, $"{name} rho");
            Close(r.PValue, Parse(o[Array.IndexOf(oh, "p")]), Tight, $"{name} p");
        }
    }

    [Test]
    public void CombinationMatchesPchisqAndPnorm()
    {
        var (ih, inputs) = ReadTsv("combination_inputs.tsv");
        var (oh, outputs) = ReadTsv("combination_results.tsv");
        foreach (var o in outputs)
        {
            string name = o[Array.IndexOf(oh, "case")];
            var rows = inputs.Where(r => r[Array.IndexOf(ih, "case")] == name).ToArray();
            var p = rows.Select(r => Parse(r[Array.IndexOf(ih, "p")])).ToArray();
            var w = rows.Select(r => Parse(r[Array.IndexOf(ih, "weight")])).ToArray();
            double Out(string c) => Parse(o[Array.IndexOf(oh, c)]);
            var fisher = PValueCombination.Fisher(p);
            var stouffer = PValueCombination.Stouffer(p);
            var weighted = PValueCombination.Stouffer(p, w);
            Close(fisher.Statistic, Out("fisher_x"), Tight, $"{name} Fisher X");
            Close(fisher.PValue, Out("fisher_p"), Tight, $"{name} Fisher p");
            Close(stouffer.Statistic, Out("stouffer_z"), Tight, $"{name} Stouffer Z");
            Close(stouffer.PValue, Out("stouffer_p"), Tight, $"{name} Stouffer p");
            Close(weighted.Statistic, Out("weighted_z"), Tight, $"{name} weighted Z");
            Close(weighted.PValue, Out("weighted_p"), Tight, $"{name} weighted p");
        }
    }
}
