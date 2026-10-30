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

    /// <summary>
    /// The fixture features the required run fits. Counted: the τ² = 0 boundary feature; the largest τ², where more than
    /// half the groups take the panel rule (identity 4 and 5, link 10); and a feature that fails when the u rule is
    /// centred on the prior instead of the joint mode (identity 4, link 1; of every feature, identity 1 and 4 and link 1
    /// and 7 fail, so the subset catches that mutation in the same cases). Scatter: the one whose σ² is on its boundary, one with a
    /// large σ², and one with a large τ² whose groups take the panel rule (12 and 19 are left out as the slowest, at four
    /// to nine times the median feature). Every feature took about 15 of the
    /// build job's 30 minutes (Debug, with coverage, on 4 cores), so the every-feature cases are explicit: run them
    /// locally (dotnet test --filter TestCategory=EveryFeature) before a push that changes the model.
    /// </summary>
    private static readonly Dictionary<string, int[]> RequiredFeatures = new()
    {
        ["counted_identity"] = new[] { 0, 4, 5 },
        ["counted_link"] = new[] { 0, 1, 10 },
        ["scatter"] = new[] { 2, 4, 14 },
    };

    private const string EveryFeatureCategory = "EveryFeature";
    private const string EveryFeatureReason = "Every fixture feature, about 15 minutes in CI; run locally before a push that changes the model.";

    /// <summary>
    /// The counts of one fixture ("NaN" = no count), and an all-missing measured matrix of the same shape: every feature,
    /// or only <see cref="RequiredFeatures"/>. Features are fitted independently, so fitting a subset changes none of
    /// their fits. Features[i] is the fixture feature on row i.
    /// </summary>
    private static (BinomialCount?[,] Counts, double[,] Unmeasured, int[] Features) Counts(string label, bool everyFeature)
    {
        var k = Matrix($"lnmm_{label}_successes.tsv");
        var n = Matrix($"lnmm_{label}_trials.tsv");
        int samples = k.GetLength(1);
        var features = everyFeature ? Enumerable.Range(0, k.GetLength(0)).ToArray() : RequiredFeatures[label];
        var counts = new BinomialCount?[features.Length, samples];
        var unmeasured = new double[features.Length, samples];
        for (int i = 0; i < features.Length; i++)
            for (int s = 0; s < samples; s++)
            {
                int f = features[i];
                unmeasured[i, s] = double.NaN;
                if (!double.IsNaN(k[f, s])) counts[i, s] = new BinomialCount((int)k[f, s], (int)n[f, s]);
            }
        return (counts, unmeasured, features);
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
    /// max(τ², 1) when σ² is given), the log-likelihood absolutely. <paramref name="f"/> is the fit's row for r's feature.
    /// </summary>
    private static void AssertFit(LogitNormalMixedModelFit fit, int f, Row r, double beta, double se, double variance, double logLik, bool sigmaEstimated, string what)
    {
        int feature = r.Feature;
        Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), $"{what} feature {feature}");
        for (int j = 0; j < 3; j++)
        {
            double rse = r[$"se{j}"];
            Assert.That(Math.Abs(fit.Coefficient(f, j) - r[$"b{j}"]), Is.LessThan(beta * rse),
                $"{what} feature {feature} coefficient {j}: {fit.Coefficient(f, j):R} vs {r[$"b{j}"]:R}");
            Assert.That(Math.Abs(fit.StandardError(f, j) - rse) / rse, Is.LessThan(se),
                $"{what} feature {feature} SE {j}: {fit.StandardError(f, j):R} vs {rse:R}");
        }
        double scale = 1;
        if (sigmaEstimated)
        {
            scale = r["sigma2"];
            Assert.That(Math.Abs(fit.ResidualVariance[f] - scale) / scale, Is.LessThan(variance),
                $"{what} feature {feature} sigma2: {fit.ResidualVariance[f]:R} vs {scale:R}");
        }
        Assert.That(Math.Abs(fit.GroupVariance[f] - r["tau2"]), Is.LessThan(variance * Math.Max(scale, r["tau2"])),
            $"{what} feature {feature} tau2: {fit.GroupVariance[f]:R} vs {r["tau2"]:R}");
        Assert.That(Math.Abs(fit.LogLikelihood[f] - r["loglik"]), Is.LessThan(logLik),
            $"{what} feature {feature} logLik: {fit.LogLikelihood[f]:R} vs {r["loglik"]:R}");
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
            AssertFit(fit, r.Feature, r, beta: 1e-5, se: 1e-4, variance: 2e-4, logLik: 1e-5, sigmaEstimated: true, $"{method}");
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
    /// The required run fits <see cref="RequiredFeatures"/>; the explicit cases fit every feature.
    /// </summary>
    [TestCase("identity", "ML", VarianceEstimator.MaximumLikelihood, false)]
    [TestCase("identity", "REML", VarianceEstimator.Reml, false)]
    [TestCase("link", "ML", VarianceEstimator.MaximumLikelihood, false)]
    [TestCase("link", "REML", VarianceEstimator.Reml, false)]
    [TestCase("identity", "ML", VarianceEstimator.MaximumLikelihood, true, Explicit = true, Reason = EveryFeatureReason, Category = EveryFeatureCategory)]
    [TestCase("identity", "REML", VarianceEstimator.Reml, true, Explicit = true, Reason = EveryFeatureReason, Category = EveryFeatureCategory)]
    [TestCase("link", "ML", VarianceEstimator.MaximumLikelihood, true, Explicit = true, Reason = EveryFeatureReason, Category = EveryFeatureCategory)]
    [TestCase("link", "REML", VarianceEstimator.Reml, true, Explicit = true, Reason = EveryFeatureReason, Category = EveryFeatureCategory)]
    public void WithEverySampleCountedItMatchesGlmerQuadrature(string linkName, string method, VarianceEstimator estimator, bool everyFeature)
    {
        var l = Rows("lnmm_counted_links.tsv", r => r.Text("name") == linkName).Single();
        var link = new DetectionLink(l["a"], l["b"], l["c"]);
        var depthMatrix = Matrix("lnmm_counted_depth.tsv");
        var depth = Enumerable.Range(0, depthMatrix.GetLength(0)).Select(i => depthMatrix[i, 0]).ToArray();
        var (counts, unmeasured, features) = Counts($"counted_{linkName}", everyFeature);
        var given = Enumerable.Repeat(TinyResidualVariance, counts.GetLength(0)).ToArray();
        // The estimates only; the likelihood-ratio tests are WithEverySampleCountedTheLikelihoodRatioMatchesGlmer's.
        var fit = LogitNormalMixedModel.Fit(unmeasured, counts, Design, Groups, link, depth: depth, residualVariance: given, estimator: estimator,
            testedCoefficients: Array.Empty<string>());

        var rows = Rows("lnmm_counted_fit.tsv", r => r.Text("link") == linkName && r.Text("method") == method && features.Contains(r.Feature));
        Assert.That(rows.Count, Is.EqualTo(features.Length));
        var exact = rows.Where(r => r.Text("status") == "fitted").ToList();
        bool ml = estimator == VarianceEstimator.MaximumLikelihood;
        foreach (var r in exact)
            AssertFit(fit, Array.IndexOf(features, r.Feature), r, beta: ml ? 1e-5 : 2e-3, se: ml ? 5e-5 : 2e-3, variance: ml ? 5e-5 : 2e-3,
                logLik: ml ? 1e-5 : 5e-5, sigmaEstimated: false, $"{linkName} {method}");
        Assert.That(exact.Count, Is.GreaterThanOrEqualTo(rows.Count * 9 / 10), "glmer's quadrature is exact for most fixture features");
    }

    /// <summary>
    /// Every sample counted, with run scatter (σ² estimated): glmmTMB with an observation-level random intercept. glmmTMB
    /// integrates by Laplace, so the two agree only to the Laplace error. Here (20 to 40 trials per sample) that is about
    /// 1% in τ² and in the standard errors, about 5% in σ² (Laplace understates it) and about 0.1 in the log-likelihood.
    /// The fixture's standard errors are glmmTMB's at its estimated variances, as ours are, not its vcov(). Under REML,
    /// glmmTMB reports β as the joint mode of (β, u), not as the maximum of the integrated likelihood, so β is compared
    /// under ML only. The REML and ML variances differ far more than the Laplace error, so the test also requires each of
    /// ours to sit nearer glmmTMB's same-method value than its other-method value. The required run fits
    /// <see cref="RequiredFeatures"/>; the explicit case fits every feature.
    /// </summary>
    [TestCase(false)]
    [TestCase(true, Explicit = true, Reason = EveryFeatureReason, Category = EveryFeatureCategory)]
    public void WithRunScatterItMatchesGlmmTmbToTheLaplaceError(bool everyFeature)
    {
        var (counts, unmeasured, features) = Counts("scatter", everyFeature);
        // The estimates only (each likelihood-ratio test is two more fits with σ estimated).
        var none = Array.Empty<string>();
        var reml = LogitNormalMixedModel.Fit(unmeasured, counts, Design, Groups, DetectionLink.Identity, estimator: VarianceEstimator.Reml, testedCoefficients: none);
        var ml = LogitNormalMixedModel.Fit(unmeasured, counts, Design, Groups, DetectionLink.Identity, estimator: VarianceEstimator.MaximumLikelihood, testedCoefficients: none);
        var tmbReml = Rows("lnmm_scatter_fit.tsv", r => r.Text("method") == "REML" && features.Contains(r.Feature)).ToDictionary(r => r.Feature);
        var tmbMl = Rows("lnmm_scatter_fit.tsv", r => r.Text("method") == "ML" && features.Contains(r.Feature)).ToDictionary(r => r.Feature);
        Assert.That(tmbReml.Count, Is.EqualTo(features.Length));
        int compared = 0, discriminating = 0;
        foreach (int feature in tmbReml.Keys.Where(f => tmbReml[f].Text("status") == "fitted" && tmbMl[f].Text("status") == "fitted"))
        {
            compared++;
            int f = Array.IndexOf(features, feature);
            foreach (var (fit, r, method) in new[] { (reml, tmbReml[feature], "REML"), (ml, tmbMl[feature], "ML") })
            {
                Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), $"{method} feature {feature}");
                Assert.That(Math.Abs(fit.GroupVariance[f] - r["tau2"]), Is.LessThan(0.03 * r["tau2"] + 0.003),
                    $"{method} feature {feature} tau2: {fit.GroupVariance[f]:R} vs {r["tau2"]:R}");
                Assert.That(Math.Abs(fit.ResidualVariance[f] - r["sigma2"]), Is.LessThan(0.1 * r["sigma2"] + 0.005),
                    $"{method} feature {feature} sigma2: {fit.ResidualVariance[f]:R} vs {r["sigma2"]:R}");
                Assert.That(Math.Abs(fit.LogLikelihood[f] - r["loglik"]), Is.LessThan(0.25),
                    $"{method} feature {feature} logLik: {fit.LogLikelihood[f]:R} vs {r["loglik"]:R}");
                for (int j = 0; j < 3; j++)
                {
                    double rse = r[$"se{j}"];
                    Assert.That(Math.Abs(fit.StandardError(f, j) - rse) / rse, Is.LessThan(0.05),
                        $"{method} feature {feature} SE {j}: {fit.StandardError(f, j):R} vs {rse:R}");
                    if (method == "ML")
                        Assert.That(Math.Abs(fit.Coefficient(f, j) - r[$"b{j}"]), Is.LessThan(0.03 * rse),
                            $"ML feature {feature} coefficient {j}: {fit.Coefficient(f, j):R} vs {r[$"b{j}"]:R}");
                }
            }
            // A feature discriminates when glmmTMB's REML and ML τ² differ by well over the Laplace band above.
            double gap = tmbReml[feature]["tau2"] - tmbMl[feature]["tau2"];
            if (Math.Abs(gap) <= 0.1 * tmbReml[feature]["tau2"] + 0.01) continue;
            discriminating++;
            Assert.That(Math.Abs(reml.GroupVariance[f] - tmbReml[feature]["tau2"]), Is.LessThan(Math.Abs(reml.GroupVariance[f] - tmbMl[feature]["tau2"])),
                $"feature {feature}: our REML tau2 sits nearer glmmTMB's REML than its ML");
            Assert.That(Math.Abs(ml.GroupVariance[f] - tmbMl[feature]["tau2"]), Is.LessThan(Math.Abs(ml.GroupVariance[f] - tmbReml[feature]["tau2"])),
                $"feature {feature}: our ML tau2 sits nearer glmmTMB's ML than its REML");
        }
        Assert.That(compared, Is.GreaterThanOrEqualTo(tmbReml.Count * 9 / 10), "glmmTMB fitted most fixture features");
        Assert.That(discriminating, Is.GreaterThanOrEqualTo(compared / 2), "REML and ML differ in most fixture features");
    }

    /// <summary>
    /// One group at a time, intercept only, every sample counted: ℓ against the exact integral (make_lnmm_ridge_fixtures.R)
    /// over one-sided and near-one-sided patterns of 1, 2 and 4 trials, β0 from −12 to 12, τ² from 0.0025 to 21904. At a
    /// large τ a one-sided group's integrand is the prior times a smoothed step, which the 40-node Gauss-Hermite rule
    /// missed by up to 0.6 (#1511). σ from 0.1 to 3: a one-sided count's integrand over e is the same shape at a large σ,
    /// which the 12-node e rule missed by 2e-7 at σ = 1, 2e-4 at σ = 2 and 4e-3 at σ = 3.
    /// </summary>
    [Test]
    public void OneGroupLikelihoodsMatchTheExactIntegral()
    {
        int compared = 0;
        var misses = new List<string>();
        double worst = 0;
        foreach (var r in Rows("lnmm_ridge_groups.tsv", r => true))
        {
            var k = r.Text("successes").Split(',').Select(int.Parse).ToArray();
            int n = (int)r["trials"];
            var design = new double[k.Length, 1];
            for (int s = 0; s < k.Length; s++) design[s, 0] = 1;
            double at = LogitNormalMixedModel.LogLikelihoodAt(Enumerable.Repeat(double.NaN, k.Length).ToArray(),
                k.Select(v => (BinomialCount?)new BinomialCount(v, n)).ToArray(), design, Enumerable.Repeat("a", k.Length).ToArray(),
                DetectionLink.Identity, new[] { r["b0"] }, r["tau2"], r["sigma2"]);
            Assert.That(double.IsFinite(at), Is.True);
            worst = Math.Max(worst, Math.Abs(at - r["loglik"]));
            if (!(Math.Abs(at - r["loglik"]) < 1e-8))
                misses.Add($"k {r.Text("successes")} of {n}, b0 {r["b0"]}, tau2 {r["tau2"]}, sigma2 {r["sigma2"]}: {at:R} vs {r["loglik"]:R}");
            compared++;
        }
        foreach (var m in misses.Take(30)) TestContext.Out.WriteLine(m);
        Assert.That(misses, Is.Empty, $"{misses.Count} of {compared} groups miss by 1e-8 or more; worst {worst:G3}");
        Assert.That(compared, Is.EqualTo(9 * 5 * 7 * 5));
        TestContext.Out.WriteLine($"{compared} groups, worst |Δℓ| {worst:G3}");
    }

    /// <summary>
    /// One count's integral over e, against integrate() (make_lnmm_ridge_fixtures.R): counts of 30 or 36 trials, mixed
    /// and one-sided, and small ones, σ from 0.5 to 5. A count of many trials near one end has a log integrand with one
    /// linear tail; the 12-node Gauss-Hermite rule missed 34 of 36 at σ = e by 3.8e-5, so such a rule gets panels.
    /// </summary>
    [Test]
    public void OneCountIntegralsMatchTheExactIntegral()
    {
        var rows = Rows("lnmm_ridge_cells.tsv", r => true);
        foreach (var r in rows)
        {
            double at = LogitNormalMixedModel.LogCountIntegral((int)r["successes"], (int)r["trials"], r["offset"], 1, Math.Sqrt(r["sigma2"]));
            Assert.That(Math.Abs(at - r["loglik"]), Is.LessThan(1e-9),
                $"{r["successes"]} of {r["trials"]}, offset {r["offset"]}, sigma2 {r["sigma2"]}: {at:R} vs {r["loglik"]:R}");
        }
        Assert.That(rows.Count, Is.EqualTo(8 * 3 * 4));
    }

    /// <summary>The exact ML fit of one of #1511's layouts, σ² = 0.25 given (make_lnmm_ridge_fixtures.R).</summary>
    private static (double[,] Unmeasured, BinomialCount?[,] Counts, double[,] Design, string[] Groups, Row Exact) RidgeLayout(string which)
    {
        int G = 8, k = 3, n = G * k;
        var design = new double[n, 2];
        var groups = new string[n];
        var y = new double[1, n];
        var counts = new BinomialCount?[1, n];
        var full = new[] { true, false, false, false, true, true, true, false };
        for (int s = 0; s < n; s++)
        {
            int g = s / k;
            design[s, 0] = 1; design[s, 1] = g < 4 ? 0 : 1; groups[s] = $"a{g}";
            y[0, s] = double.NaN;
            counts[0, s] = new BinomialCount(full[g] ? 2 : 0, 2);
        }
        counts[0, which == "A" ? 0 : k] = new BinomialCount(1, 2);
        return (y, counts, design, groups, Rows("lnmm_ridge_fits.tsv", r => r.Text("layout") == which).Single());
    }

    /// <summary>
    /// #1511's layouts, ML with σ² given: the fit lands on the exact maximum. Before, the u rule's error put layout A's
    /// fit at τ² 233 (exact 142) and called layout B Separated, though its maximum exists. The ridge is flat (ℓ falls by
    /// 0.05 from τ 12 to 16), so β and τ² are compared to 1e-3 relatively and ℓ to 1e-7.
    /// </summary>
    [TestCase("A")]
    [TestCase("B")]
    public void NearOneSidedLayoutsFitToTheExactMaximum(string which)
    {
        var (y, counts, design, groups, r) = RidgeLayout(which);
        // Untested: the estimate itself (the ridge check would blank β; see NearOneSidedLayoutsAreUnboundedWithALikelihoodRatioP).
        var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, DetectionLink.Identity, residualVariance: new[] { r["sigma2"] },
            estimator: VarianceEstimator.MaximumLikelihood, testedCoefficients: Array.Empty<string>());
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.LogLikelihood[0], Is.EqualTo(r["loglik"]).Within(1e-7));
        Assert.That(fit.GroupVariance[0], Is.EqualTo(r["tau2"]).Within(1e-3 * r["tau2"]));
        Assert.That(fit.Coefficient(0, 0), Is.EqualTo(r["b0"]).Within(1e-3 * Math.Abs(r["b0"])));
        Assert.That(fit.Coefficient(0, 1), Is.EqualTo(r["b1"]).Within(1e-3 * Math.Abs(r["b1"])));
    }

    /// <summary>
    /// At #1511's exact estimates, the panel rule's β-gradient vanishes and its H = −∂²ℓ/∂β² is the exact one (numDeriv in
    /// R, good to about 1e-7 relatively): the derivatives the Newton steps and standard errors use are right on the rule.
    /// </summary>
    [TestCase("A")]
    [TestCase("B")]
    public void TheDerivativesAtTheExactEstimateAreExact(string which)
    {
        var (y, counts, design, groups, r) = RidgeLayout(which);
        int n = y.GetLength(1);
        var gradient = new double[2];
        var h = new double[2, 2];
        LogitNormalMixedModel.LogLikelihoodAt(Enumerable.Range(0, n).Select(s => y[0, s]).ToArray(), Enumerable.Range(0, n).Select(s => counts[0, s]).ToArray(),
            design, groups, DetectionLink.Identity, new[] { r["b0"], r["b1"] }, r["tau2"], r["sigma2"], gradient: gradient, negativeHessian: h);
        Assert.That(Math.Abs(gradient[0]) + Math.Abs(gradient[1]), Is.LessThan(1e-5));
        Assert.That(h[0, 0], Is.EqualTo(r["h00"]).Within(1e-5 * r["h00"]));
        Assert.That(h[0, 1], Is.EqualTo(r["h01"]).Within(1e-5 * Math.Abs(r["h01"])));
        Assert.That(h[1, 0], Is.EqualTo(h[0, 1]).Within(1e-12 * Math.Abs(h[0, 1])));
        Assert.That(h[1, 1], Is.EqualTo(r["h11"]).Within(1e-5 * r["h11"]));
    }

    /// <summary>
    /// #1511's layouts, the default tests: the likelihood-ratio statistic is the exact one, 2(ℓ̂ − ℓ̂₀) from the exact ML fits,
    /// read against F(1, 6); and the trait is unbounded (its estimate is past an effect of 10 logits, with the likelihood
    /// still within 1 of its maximum by β1 = 48), so β and its SE are NaN and the p-value stands. Under REML too: REML
    /// has no maximum here (its criterion rises with τ until β passes SeparationLimit), and the ML fits carry the tests.
    /// Before, ML reported p = 0.013 (A) or Separated (B), and REML p = 0.037 and 0.025.
    /// </summary>
    [TestCase("A", VarianceEstimator.MaximumLikelihood)]
    [TestCase("A", VarianceEstimator.Reml)]
    [TestCase("B", VarianceEstimator.MaximumLikelihood)]
    [TestCase("B", VarianceEstimator.Reml)]
    public void NearOneSidedLayoutsAreUnboundedWithALikelihoodRatioP(string which, VarianceEstimator estimator)
    {
        var (y, counts, design, groups, r) = RidgeLayout(which);
        var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, DetectionLink.Identity, residualVariance: new[] { r["sigma2"] },
            estimator: estimator);
        double lr = 2 * (r["loglik"] - r["null_loglik"]);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Unbounded));
        Assert.That(fit.IsUnbounded(0, 1), Is.True);
        Assert.That(fit.Coefficient(0, 1), Is.NaN);
        Assert.That(fit.StandardError(0, 1), Is.NaN);
        Assert.That(fit.DegreesOfFreedom(0, 1), Is.EqualTo(6));
        Assert.That(fit.LikelihoodRatio(0, 1), Is.EqualTo(lr).Within(2e-7));
        double p = 2 * MathNet.Numerics.Distributions.StudentT.CDF(0, 1, 6, -Math.Sqrt(lr));
        Assert.That(double.IsFinite(fit.PValue(0, 1)), Is.True);
        Assert.That(fit.PValue(0, 1), Is.EqualTo(p).Within(1e-6 * p));
        TestContext.Out.WriteLine($"{which} {estimator}: LR {fit.LikelihoodRatio(0, 1):G6}, p {fit.PValue(0, 1):G4}");
    }

    /// <summary>
    /// The ridge check, against exact fits (make_lnmm_ridge_fixtures.R, lnmm_ridge_held.tsv), on two layouts like #1511's
    /// whose trait estimates (5.1 and 3.3) are under the effect bound of 10 logits: held at 10, D's likelihood is still
    /// inside the 95% region (LR below F(1, 6)'s 5.99), so D is unbounded; E's is not, so E is fitted, with its estimate.
    /// </summary>
    [TestCase("D", true)]
    [TestCase("E", false)]
    public void TheRidgeCheckHoldsTheCoefficientAtTheEffectBound(string which, bool unbounded)
    {
        var r = Rows("lnmm_ridge_held.tsv", r => r.Text("layout") == which).Single();
        var pattern = r.Text("pattern").Split(' ').Select(g => g.Select(c => c - '0').ToArray()).ToArray();
        int G = 8, k = 3, n = G * k;
        var design = new double[n, 2];
        var groups = new string[n];
        var y = new double[1, n];
        var counts = new BinomialCount?[1, n];
        for (int s = 0; s < n; s++)
        {
            int g = s / k;
            design[s, 0] = 1; design[s, 1] = g < 4 ? 0 : 1; groups[s] = $"a{g}";
            y[0, s] = double.NaN;
            counts[0, s] = new BinomialCount(pattern[g][s % k], 2);
        }
        double critical = Math.Pow(MathNet.Numerics.Distributions.StudentT.InvCDF(0, 1, 6, 0.975), 2);
        double heldLr = 2 * (r["loglik"] - r["bound_loglik"]);
        Assert.That(heldLr < critical, Is.EqualTo(unbounded), $"the fixture itself: held LR {heldLr:G5} against {critical:G5}");
        Assert.That(Math.Abs(heldLr - critical), Is.GreaterThan(0.05), "the fixture is not on the edge");

        var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, DetectionLink.Identity, residualVariance: new[] { r["sigma2"] },
            estimator: VarianceEstimator.MaximumLikelihood, testedCoefficients: new[] { "c1" });
        Assert.That(fit.LikelihoodRatio(0, 1), Is.EqualTo(2 * (r["loglik"] - r["null_loglik"])).Within(2e-6));
        Assert.That(fit.IsUnbounded(0, 1), Is.EqualTo(unbounded));
        Assert.That(fit.Status[0], Is.EqualTo(unbounded ? FeatureFitStatus.Unbounded : FeatureFitStatus.Fitted));
        if (unbounded) Assert.That(fit.Coefficient(0, 1), Is.NaN);
        else
        {
            Assert.That(fit.Coefficient(0, 1), Is.EqualTo(r["b1"]).Within(1e-3 * Math.Abs(r["b1"])));
            Assert.That(fit.LogLikelihood[0], Is.EqualTo(r["loglik"]).Within(1e-6));
        }
    }

    /// <summary>
    /// The likelihood-ratio statistic against lme4: make_lnmm_fixtures.R's count-only features, glmer at 25-point AGQ with
    /// the trait or the covariate x dropped (lnmm_ridge_glmer_lr.tsv). Only the fits where glmer's quadrature is exact
    /// ("fitted" in lnmm_counted_fit.tsv) are compared; σ² is given tiny here and absent in glmer, which moves ℓ by about
    /// 1e-6. No coefficient there is unbounded: these are ordinary features. The required run fits
    /// <see cref="RequiredFeatures"/>; the explicit cases fit every feature.
    /// </summary>
    [TestCase("identity", false)]
    [TestCase("link", false)]
    [TestCase("identity", true, Explicit = true, Reason = EveryFeatureReason, Category = EveryFeatureCategory)]
    [TestCase("link", true, Explicit = true, Reason = EveryFeatureReason, Category = EveryFeatureCategory)]
    public void WithEverySampleCountedTheLikelihoodRatioMatchesGlmer(string linkName, bool everyFeature)
    {
        var l = Rows("lnmm_counted_links.tsv", r => r.Text("name") == linkName).Single();
        var link = new DetectionLink(l["a"], l["b"], l["c"]);
        var depthMatrix = Matrix("lnmm_counted_depth.tsv");
        var depth = Enumerable.Range(0, depthMatrix.GetLength(0)).Select(i => depthMatrix[i, 0]).ToArray();
        var (counts, unmeasured, features) = Counts($"counted_{linkName}", everyFeature);
        var given = Enumerable.Repeat(TinyResidualVariance, counts.GetLength(0)).ToArray();
        var fit = LogitNormalMixedModel.Fit(unmeasured, counts, Design, Groups, link, depth: depth, residualVariance: given,
            estimator: VarianceEstimator.MaximumLikelihood, coefficientNames: new[] { "intercept", "trait", "x" }, testedCoefficients: new[] { "trait", "x" });
        var exact = Rows("lnmm_counted_fit.tsv", r => r.Text("link") == linkName && r.Text("method") == "ML" && r.Text("status") == "fitted")
            .Select(r => r.Feature).Where(features.Contains).ToHashSet();
        int compared = 0;
        foreach (var r in Rows("lnmm_ridge_glmer_lr.tsv", r => r.Text("link") == linkName && exact.Contains(r.Feature)))
        {
            int f = Array.IndexOf(features, r.Feature), j = fit.IndexOf(r.Text("coefficient"));
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), $"feature {r.Feature}");
            Assert.That(fit.IsUnbounded(f, j), Is.False);
            Assert.That(fit.LikelihoodRatio(f, j), Is.EqualTo(r["lr"]).Within(2e-5),
                $"feature {r.Feature} {r.Text("coefficient")}: {fit.LikelihoodRatio(f, j):R} vs {r["lr"]:R}");
            compared++;
        }
        Assert.That(compared, Is.GreaterThanOrEqualTo(2 * (features.Length - (everyFeature ? 2 : 0))), "both coefficients of (nearly) every feature");
    }
}
