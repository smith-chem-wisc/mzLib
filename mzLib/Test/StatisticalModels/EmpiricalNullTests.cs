using MathNet.Numerics.Distributions;
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
/// <see cref="EmpiricalNull"/>: Efron's maximum-likelihood empirical null, against R's locfdr 1.1.8
/// (<c>locfdr(z, nulltype = 1)</c>, frozen in ReferenceData, PROVENANCE_empirical_null.txt; R is never run here), and its
/// refusals and failure statuses. No z-value is stored: each scenario's are built from its recipe, as the R script built
/// them.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class EmpiricalNullTests
{
    // The exact maximum (mle_* in the fixture) is locmle run to ep = 1e-14 from the window scores' moments; EmpiricalNull
    // solves to 1e-12 on the same windows, so the two agree to rounding.
    private const double ExactTolerance = 1e-9; // absolute, in z units (delta, sigma) and for p0
    // locfdr's own estimate stops when its natural parameters move by less than 1e-5, which leaves delta and sigma within
    // about 1e-5 of the maximum, more as sigma grows (1.3e-5 at sigma 6.5 in clamped_one_side_30).
    private const double LocfdrTolerance = 1e-5; // relative to max(1, sigma0)
    private const double WindowTolerance = 1e-9; // the first window is closed-form in both

    private static List<Dictionary<string, string>> ReadTsv(string name)
    {
        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData");
        var lines = File.ReadAllLines(Path.Combine(dir, name)).Where(l => l.Length > 0).ToArray();
        string[] header = lines[0].Split('\t');
        return lines.Skip(1).Select(l => l.Split('\t')).Select(f => header.Zip(f).ToDictionary(p => p.First, p => p.Second)).ToList();
    }

    private static double D(string s) =>
        s == "NA" ? double.NaN : double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static IEnumerable<TestCaseData> Scenarios() =>
        ReadTsv("empirical_null_results.tsv").Select(r => new TestCaseData(r["scenario"]).SetName("MatchesLocfdr_" + r["scenario"]));

    /// <summary>
    /// A scenario's z-values, built from its recipe exactly as make_empirical_null_fixtures.R builds them: parts separated
    /// by ';', each "n,mean,sd" giving mean + sd * Φ⁻¹((i − 0.5) / n) for i = 1..n, or n copies of mean when sd is 0; an
    /// optional "|k" repeats every value k times.
    /// </summary>
    private static double[] Scores(string scenario)
    {
        string recipe = ReadTsv("empirical_null_results.tsv").Single(r => r["scenario"] == scenario)["recipe"];
        string[] halves = recipe.Split('|');
        int each = halves.Length > 1 ? int.Parse(halves[1], CultureInfo.InvariantCulture) : 1;
        var z = new List<double>();
        foreach (string part in halves[0].Split(';'))
        {
            double[] v = part.Split(',').Select(D).ToArray();
            int n = (int)v[0];
            for (int i = 1; i <= n; i++)
            {
                double value = v[2] == 0 ? v[1] : v[1] + v[2] * Normal.InvCDF(0, 1, (i - 0.5) / n);
                for (int k = 0; k < each; k++) z.Add(value);
            }
        }
        return z.ToArray();
    }

    [TestCaseSource(nameof(Scenarios))]
    public void MatchesLocfdr(string scenario)
    {
        var expected = ReadTsv("empirical_null_results.tsv").Single(r => r["scenario"] == scenario);
        double[] z = Scores(scenario);
        Assert.That(z, Has.Length.EqualTo(int.Parse(expected["n"], CultureInfo.InvariantCulture)));

        var fit = EmpiricalNull.Fit(z);

        Assert.That(fit.Scores, Is.EqualTo(z.Length));
        Assert.That(fit.Status, Is.EqualTo(EmpiricalNullStatus.Estimated));
        Assert.That(fit.FirstDelta, Is.EqualTo(D(expected["mle_first_delta0"])).Within(ExactTolerance), "first pass delta0");
        Assert.That(fit.FirstSigma, Is.EqualTo(D(expected["mle_first_sigma0"])).Within(ExactTolerance), "first pass sigma0");
        Assert.That(fit.Delta, Is.EqualTo(D(expected["mle_delta0"])).Within(ExactTolerance), "delta0");
        Assert.That(fit.Sigma, Is.EqualTo(D(expected["mle_sigma0"])).Within(ExactTolerance), "sigma0");
        Assert.That(fit.P0, Is.EqualTo(D(expected["mle_p0"])).Within(ExactTolerance), "p0");

        // locfdr's own result, where it gave one: equal within its stopping rule. Where it stopped with "ML estimation
        // failed" (clamped_both_25) the maximum exists all the same, and EmpiricalNull finds it.
        if (expected["locfdr_status"] != "estimated") return;
        double tolerance = LocfdrTolerance * Math.Max(1, fit.Sigma);
        Assert.That(fit.Delta, Is.EqualTo(D(expected["locfdr_delta0"])).Within(tolerance), "locfdr delta0");
        Assert.That(fit.Sigma, Is.EqualTo(D(expected["locfdr_sigma0"])).Within(tolerance), "locfdr sigma0");
        Assert.That(fit.P0, Is.EqualTo(D(expected["locfdr_p0"])).Within(tolerance), "locfdr p0");
    }

    [TestCaseSource(nameof(Scenarios))]
    public void FirstWindowIsLocfdrs(string scenario)
    {
        var expected = ReadTsv("empirical_null_results.tsv").Single(r => r["scenario"] == scenario);
        double[] z = Scores(scenario);
        Array.Sort(z);
        int n = z.Length;
        // R's type-7 quantiles, written out, so the window rule is checked apart from the fit.
        double Q(double p) { double h = (n - 1) * p; int lo = (int)Math.Floor(h); return z[lo] + (h - lo) * (z[Math.Min(lo + 1, n - 1)] - z[lo]); }
        double b = 4.3 * Math.Exp(-0.26 * Math.Log10(n));
        Assert.That(Q(0.5), Is.EqualTo(D(expected["first_center"])).Within(WindowTolerance));
        Assert.That(b * (Q(0.75) - Q(0.25)) / (2 * Normal.InvCDF(0, 1, 0.75)), Is.EqualTo(D(expected["first_half_width"])).Within(WindowTolerance));
    }

    [Test]
    public void TheFinalWindowIsTheFirstEstimatePlusOrMinusBTimesItsSpread()
    {
        double[] z = Scores("null_2000");
        var fit = EmpiricalNull.Fit(z);
        double b = 4.3 * Math.Exp(-0.26 * Math.Log10(z.Length));
        Assert.That(fit.WindowLow, Is.EqualTo(fit.FirstDelta - b * fit.FirstSigma).Within(1e-12));
        Assert.That(fit.WindowHigh, Is.EqualTo(fit.FirstDelta + b * fit.FirstSigma).Within(1e-12));
        Assert.That(fit.ScoresInWindow, Is.EqualTo(z.Count(v => v >= fit.WindowLow && v <= fit.WindowHigh)));
    }

    [Test]
    public void AShiftOfEveryScoreShiftsDeltaAndLeavesSigma()
    {
        double[] z = Scores("one_sided_10pct");
        var fit = EmpiricalNull.Fit(z);
        var moved = EmpiricalNull.Fit(z.Select(v => v + 0.25).ToArray());
        Assert.That(moved.Delta, Is.EqualTo(fit.Delta + 0.25).Within(1e-9));
        Assert.That(moved.Sigma, Is.EqualTo(fit.Sigma).Within(1e-9));
        Assert.That(moved.P0, Is.EqualTo(fit.P0).Within(1e-9));
    }

    [Test]
    public void FewerThanThreeScoresAreTooFew()
    {
        var fit = EmpiricalNull.Fit(new[] { 0.1, -0.2 });
        Assert.That(fit.Status, Is.EqualTo(EmpiricalNullStatus.TooFewScores));
        Assert.That(new[] { fit.Delta, fit.Sigma, fit.P0, fit.FirstDelta, fit.FirstSigma }, Has.All.NaN);
        Assert.That(fit.Scores, Is.EqualTo(2));
    }

    [Test]
    public void ScoresThatAllTieHaveNoMaximum()
    {
        // The window then has no width, so no normal can be fitted inside it.
        var fit = EmpiricalNull.Fit(Enumerable.Repeat(0.3, 500).ToArray());
        Assert.That(fit.Status, Is.EqualTo(EmpiricalNullStatus.NoMaximum));
        Assert.That(new[] { fit.Delta, fit.Sigma, fit.P0 }, Has.All.NaN);
    }

    [Test]
    public void ScoresMoreSpreadThanAUniformOverTheWindowHaveNoMaximum()
    {
        // 45% of the scores at -1, 45% at +1 and 10% spread between: the quartiles are -1 and +1, and at 400,000 scores
        // b is about 1, so the window is about +-1.48 and its scores have variance 0.41 in units of the half-width, more
        // than a uniform's 1/3. No normal truncated to the window fits them: the likelihood rises as sigma grows.
        int n = 400_000;
        var z = new double[n];
        for (int i = 0; i < n; i++)
            z[i] = i < 0.45 * n ? -1 : i < 0.9 * n ? 1 : -1 + 2.0 * (i - 0.9 * n) / (0.1 * n);
        var fit = EmpiricalNull.Fit(z);
        Assert.That(fit.Status, Is.EqualTo(EmpiricalNullStatus.NoMaximum));
        Assert.That(new[] { fit.Delta, fit.Sigma, fit.P0 }, Has.All.NaN);
    }

    [Test]
    public void BadInputIsRefused()
    {
        Assert.Throws<ArgumentNullException>(() => EmpiricalNull.Fit(null!));
        Assert.Throws<ArgumentException>(() => EmpiricalNull.Fit(new[] { 0.1, double.NaN, 0.3 }));
        Assert.Throws<ArgumentException>(() => EmpiricalNull.Fit(new[] { 0.1, double.PositiveInfinity, 0.3 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmpiricalNull.Fit(new double[EmpiricalNull.MaximumScores + 1]));
    }
}
