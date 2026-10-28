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
/// Compares LogitNormalMixedModel against frozen outputs of lme4 (lmer, glmer) and glmmTMB. The outputs, the script that
/// made them and its versions are in ReferenceData (make_lnmm_fixtures.R, PROVENANCE_lnmm.txt). R is never run by these
/// tests. All three cases share one layout: 8 groups of 2 to 4 samples, a group-level trait and a sample-level covariate.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReferenceComparisonLogitNormalTests
{
    /// <summary>
    /// The σ² given in the counted cases, where R fits no run scatter. The model takes σ² ≥ e^−16; at this value the
    /// scatter moves a log-likelihood by about 1e-6.
    /// </summary>
    private const double TinyResidualVariance = 2e-7;

    private static string DataDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData");

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

    private static double[,] Design => Matrix("lnmm_design.tsv");

    private static string[] Groups => ReadTsv("lnmm_groups.tsv").Rows.Select(r => r[0]).ToArray();

    /// <summary>The counts of one fixture ("NaN" = no count), and an all-missing measured matrix of the same shape.</summary>
    private static (BinomialCount?[,] Counts, double[,] Unmeasured) Counts(string label)
    {
        var k = Matrix($"lnmm_{label}_successes.tsv");
        var n = Matrix($"lnmm_{label}_trials.tsv");
        int features = k.GetLength(0), samples = k.GetLength(1);
        var counts = new BinomialCount?[features, samples];
        var unmeasured = new double[features, samples];
        for (int f = 0; f < features; f++)
            for (int s = 0; s < samples; s++)
            {
                unmeasured[f, s] = double.NaN;
                if (!double.IsNaN(k[f, s])) counts[f, s] = new BinomialCount((int)k[f, s], (int)n[f, s]);
            }
        return (counts, unmeasured);
    }

    /// <summary>One fixture row's columns by name.</summary>
    private sealed class Row(string[] header, string[] values)
    {
        public string Text(string column) => values[Array.IndexOf(header, column)];
        public double this[string column] => Parse(Text(column));
        public int Feature => (int)this["feature"];
    }

    private static List<Row> Rows(string name, Func<Row, bool> where)
    {
        var (h, rows) = ReadTsv(name);
        return rows.Select(r => new Row(h, r)).Where(where).ToList();
    }

    /// <summary>
    /// β in units of R's standard error, standard errors and σ² relatively, τ² on the scale of max(τ², σ²) (or of
    /// max(τ², 1) when σ² is given), the log-likelihood absolutely.
    /// </summary>
    private static void AssertFit(LogitNormalMixedModelFit fit, Row r, double beta, double se, double variance, double logLik, bool sigmaEstimated, string what)
    {
        int f = r.Feature;
        Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), $"{what} feature {f}");
        for (int j = 0; j < 3; j++)
        {
            double rse = r[$"se{j}"];
            Assert.That(Math.Abs(fit.Coefficient(f, j) - r[$"b{j}"]), Is.LessThan(beta * rse),
                $"{what} feature {f} coefficient {j}: {fit.Coefficient(f, j):R} vs {r[$"b{j}"]:R}");
            Assert.That(Math.Abs(fit.StandardError(f, j) - rse) / rse, Is.LessThan(se),
                $"{what} feature {f} SE {j}: {fit.StandardError(f, j):R} vs {rse:R}");
        }
        double scale = 1;
        if (sigmaEstimated)
        {
            scale = r["sigma2"];
            Assert.That(Math.Abs(fit.ResidualVariance[f] - scale) / scale, Is.LessThan(variance),
                $"{what} feature {f} sigma2: {fit.ResidualVariance[f]:R} vs {scale:R}");
        }
        Assert.That(Math.Abs(fit.GroupVariance[f] - r["tau2"]), Is.LessThan(variance * Math.Max(scale, r["tau2"])),
            $"{what} feature {f} tau2: {fit.GroupVariance[f]:R} vs {r["tau2"]:R}");
        Assert.That(Math.Abs(fit.LogLikelihood[f] - r["loglik"]), Is.LessThan(logLik),
            $"{what} feature {f} logLik: {fit.LogLikelihood[f]:R} vs {r["loglik"]:R}");
    }

    /// <summary>
    /// With nothing counted the model is a linear mixed model, and lmer's fit. lme4 stops its variance search on its own
    /// rule, which limits the agreement to about 1e-5.
    /// </summary>
    [TestCase("REML", VarianceEstimator.Reml)]
    [TestCase("ML", VarianceEstimator.MaximumLikelihood)]
    public void WithEverySampleMeasuredItMatchesLmer(string method, VarianceEstimator estimator)
    {
        var responses = Matrix("lnmm_measured_responses.tsv");
        var counts = new BinomialCount?[responses.GetLength(0), responses.GetLength(1)];
        var fit = LogitNormalMixedModel.Fit(responses, counts, Design, Groups, DetectionLink.Identity, estimator: estimator);
        var rows = Rows("lnmm_measured_fit.tsv", r => r.Text("method") == method);
        foreach (var r in rows.Where(r => r.Text("status") == "fitted"))
            AssertFit(fit, r, beta: 1e-5, se: 1e-4, variance: 2e-4, logLik: 1e-5, sigmaEstimated: true, $"{method}");
        Assert.That(rows.Count(r => r.Text("status") == "fitted"), Is.EqualTo(rows.Count), "every fixture feature is fittable");
    }

    /// <summary>
    /// Every sample counted, with no run scatter: glmer's 25-point adaptive quadrature, through the identity link and
    /// through one with an intercept, a slope of 0.6 and a depth term.
    /// <list type="bullet">
    /// <item>ML is glmer's own fit, and exact; it agrees to about 1e-6.</item>
    /// <item>
    /// REML integrates β out of glmer's likelihood by Laplace, with a numerical Hessian. That makes R's criterion noisy at
    /// about 3e-6 and flat to that noise over ±5e-4 in τ², so the agreement is looser there.
    /// </item>
    /// <item>
    /// Fits where glmer's own quadrature is not exact (τ² above about 10; "reference_inexact") are not compared. The
    /// fixture checked each fit's log-likelihood against integrate().
    /// </item>
    /// </list>
    /// </summary>
    [TestCase("identity", "ML", VarianceEstimator.MaximumLikelihood)]
    [TestCase("identity", "REML", VarianceEstimator.Reml)]
    [TestCase("link", "ML", VarianceEstimator.MaximumLikelihood)]
    [TestCase("link", "REML", VarianceEstimator.Reml)]
    public void WithEverySampleCountedItMatchesGlmerQuadrature(string linkName, string method, VarianceEstimator estimator)
    {
        var l = Rows("lnmm_counted_links.tsv", r => r.Text("name") == linkName).Single();
        var link = new DetectionLink(l["a"], l["b"], l["c"]);
        var depthMatrix = Matrix("lnmm_counted_depth.tsv");
        var depth = Enumerable.Range(0, depthMatrix.GetLength(0)).Select(i => depthMatrix[i, 0]).ToArray();
        var (counts, unmeasured) = Counts($"counted_{linkName}");
        var given = Enumerable.Repeat(TinyResidualVariance, counts.GetLength(0)).ToArray();
        var fit = LogitNormalMixedModel.Fit(unmeasured, counts, Design, Groups, link, depth: depth, residualVariance: given, estimator: estimator);

        var rows = Rows("lnmm_counted_fit.tsv", r => r.Text("link") == linkName && r.Text("method") == method);
        var exact = rows.Where(r => r.Text("status") == "fitted").ToList();
        bool ml = estimator == VarianceEstimator.MaximumLikelihood;
        foreach (var r in exact)
            AssertFit(fit, r, beta: ml ? 1e-5 : 2e-3, se: ml ? 5e-5 : 2e-3, variance: ml ? 5e-5 : 2e-3, logLik: ml ? 1e-5 : 5e-5,
                sigmaEstimated: false, $"{linkName} {method}");
        Assert.That(exact.Count, Is.GreaterThanOrEqualTo(rows.Count * 9 / 10), "glmer's quadrature is exact for most fixture features");
    }

    /// <summary>
    /// Every sample counted, with run scatter (σ² estimated): glmmTMB with an observation-level random intercept. glmmTMB
    /// integrates by Laplace, so the two agree only to the Laplace error. Here (20 to 40 trials per sample) that is about
    /// 1% in τ² and in the standard errors, about 5% in σ² (Laplace understates it) and about 0.1 in the log-likelihood.
    /// The fixture's standard errors are glmmTMB's at its estimated variances, as ours are, not its vcov(). Under REML,
    /// glmmTMB reports β as the joint mode of (β, u), not as the maximum of the integrated likelihood, so β is compared
    /// under ML only. The REML and ML variances differ far more than the Laplace error, so the test also requires each of
    /// ours to sit nearer glmmTMB's same-method value than its other-method value.
    /// </summary>
    [Test]
    public void WithRunScatterItMatchesGlmmTmbToTheLaplaceError()
    {
        var (counts, unmeasured) = Counts("scatter");
        var reml = LogitNormalMixedModel.Fit(unmeasured, counts, Design, Groups, DetectionLink.Identity, estimator: VarianceEstimator.Reml);
        var ml = LogitNormalMixedModel.Fit(unmeasured, counts, Design, Groups, DetectionLink.Identity, estimator: VarianceEstimator.MaximumLikelihood);
        var tmbReml = Rows("lnmm_scatter_fit.tsv", r => r.Text("method") == "REML").ToDictionary(r => r.Feature);
        var tmbMl = Rows("lnmm_scatter_fit.tsv", r => r.Text("method") == "ML").ToDictionary(r => r.Feature);
        int compared = 0, discriminating = 0;
        foreach (int f in tmbReml.Keys.Where(f => tmbReml[f].Text("status") == "fitted" && tmbMl[f].Text("status") == "fitted"))
        {
            compared++;
            foreach (var (fit, r, method) in new[] { (reml, tmbReml[f], "REML"), (ml, tmbMl[f], "ML") })
            {
                Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), $"{method} feature {f}");
                Assert.That(Math.Abs(fit.GroupVariance[f] - r["tau2"]), Is.LessThan(0.03 * r["tau2"] + 0.003),
                    $"{method} feature {f} tau2: {fit.GroupVariance[f]:R} vs {r["tau2"]:R}");
                Assert.That(Math.Abs(fit.ResidualVariance[f] - r["sigma2"]), Is.LessThan(0.1 * r["sigma2"] + 0.005),
                    $"{method} feature {f} sigma2: {fit.ResidualVariance[f]:R} vs {r["sigma2"]:R}");
                Assert.That(Math.Abs(fit.LogLikelihood[f] - r["loglik"]), Is.LessThan(0.25),
                    $"{method} feature {f} logLik: {fit.LogLikelihood[f]:R} vs {r["loglik"]:R}");
                for (int j = 0; j < 3; j++)
                {
                    double rse = r[$"se{j}"];
                    Assert.That(Math.Abs(fit.StandardError(f, j) - rse) / rse, Is.LessThan(0.05),
                        $"{method} feature {f} SE {j}: {fit.StandardError(f, j):R} vs {rse:R}");
                    if (method == "ML")
                        Assert.That(Math.Abs(fit.Coefficient(f, j) - r[$"b{j}"]), Is.LessThan(0.03 * rse),
                            $"ML feature {f} coefficient {j}: {fit.Coefficient(f, j):R} vs {r[$"b{j}"]:R}");
                }
            }
            // A feature discriminates when glmmTMB's REML and ML τ² differ by well over the Laplace band above.
            double gap = tmbReml[f]["tau2"] - tmbMl[f]["tau2"];
            if (Math.Abs(gap) <= 0.1 * tmbReml[f]["tau2"] + 0.01) continue;
            discriminating++;
            Assert.That(Math.Abs(reml.GroupVariance[f] - tmbReml[f]["tau2"]), Is.LessThan(Math.Abs(reml.GroupVariance[f] - tmbMl[f]["tau2"])),
                $"feature {f}: our REML tau2 sits nearer glmmTMB's REML than its ML");
            Assert.That(Math.Abs(ml.GroupVariance[f] - tmbMl[f]["tau2"]), Is.LessThan(Math.Abs(ml.GroupVariance[f] - tmbReml[f]["tau2"])),
                $"feature {f}: our ML tau2 sits nearer glmmTMB's ML than its REML");
        }
        Assert.That(compared, Is.GreaterThanOrEqualTo(tmbReml.Count * 9 / 10), "glmmTMB fitted most fixture features");
        Assert.That(discriminating, Is.GreaterThanOrEqualTo(compared / 2), "REML and ML differ in most fixture features");
    }
}
