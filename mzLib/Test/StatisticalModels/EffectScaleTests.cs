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
/// <see cref="EffectScale"/>: the scale between two sets of effects with measurement error in both, against mr.raps 0.2's
/// mr.raps.overdispersed (frozen in ReferenceData, PROVENANCE_effect_scale.txt; R is never run here), and its closed
/// forms, refusals and rescaling.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class EffectScaleTests
{
    // mr.raps maximises each profile with R's optimize at tol = sqrt(machine epsilon) (about 1.5e-8 in the argument), so
    // its scale and τ² carry that much error; EffectScale solves to 1e-12. The tolerances below are set from that.
    private const double ScaleTolerance = 1e-6;          // relative
    private const double Tau2Tolerance = 1e-6;           // absolute, in the scaled effects' squared units
    private const double StandardErrorTolerance = 1e-4;  // relative; the sandwich carries τ²'s error
    private const double ClosedFormTolerance = 1e-12;    // relative; same arithmetic as the R script

    private static string DataDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData");

    private static List<Dictionary<string, string>> ReadTsv(string name)
    {
        var lines = File.ReadAllLines(Path.Combine(DataDirectory, name)).Where(l => l.Length > 0).ToArray();
        string[] header = lines[0].Split('\t');
        return lines.Skip(1).Select(l => l.Split('\t')).Select(f => header.Zip(f).ToDictionary(p => p.First, p => p.Second)).ToList();
    }

    private static double D(string s) => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static IEnumerable<TestCaseData> Scenarios() =>
        ReadTsvStatic("effect_scale_results.tsv").Select(r => new TestCaseData(r["scenario"]).SetName("MatchesMrRaps_" + r["scenario"]));

    private static List<Dictionary<string, string>> ReadTsvStatic(string name)
    {
        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData");
        var lines = File.ReadAllLines(Path.Combine(dir, name)).Where(l => l.Length > 0).ToArray();
        string[] header = lines[0].Split('\t');
        return lines.Skip(1).Select(l => l.Split('\t')).Select(f => header.Zip(f).ToDictionary(p => p.First, p => p.Second)).ToList();
    }

    private static (double[] X, double[] Sx, double[] Y, double[] Sy) Inputs(string scenario)
    {
        var rows = ReadTsv("effect_scale_inputs.tsv").Where(r => r["scenario"] == scenario).ToList();
        return (rows.Select(r => D(r["reference"])).ToArray(), rows.Select(r => D(r["reference_se"])).ToArray(),
                rows.Select(r => D(r["scaled"])).ToArray(), rows.Select(r => D(r["scaled_se"])).ToArray());
    }

    [TestCaseSource(nameof(Scenarios))]
    public void MatchesMrRaps(string scenario)
    {
        var expected = ReadTsv("effect_scale_results.tsv").Single(r => r["scenario"] == scenario);
        var (x, sx, y, sy) = Inputs(scenario);

        EffectScaleResult fit = EffectScale.Estimate(x, sx, y, sy);

        Assert.That(fit.Status, Is.EqualTo(EffectScaleStatus.Estimated));
        Assert.That(fit.Features, Is.EqualTo(int.Parse(expected["features"])));
        Assert.That(fit.Scale, Is.EqualTo(D(expected["scale"])).Within(ScaleTolerance * Math.Abs(D(expected["scale"]))), "scale");
        Assert.That(fit.Tau2, Is.EqualTo(D(expected["tau2"])).Within(Tau2Tolerance), "tau2");
        Assert.That(fit.StandardError, Is.EqualTo(D(expected["scale_se"])).Within(StandardErrorTolerance * D(expected["scale_se"])), "scale se");
        Assert.That(fit.Tau2StandardError, Is.EqualTo(D(expected["tau2_se"])).Within(StandardErrorTolerance * D(expected["tau2_se"])), "tau2 se");
        Assert.That(fit.InverseVarianceScale, Is.EqualTo(D(expected["ivw_scale"])).Within(ClosedFormTolerance * Math.Abs(D(expected["ivw_scale"]))));
        Assert.That(fit.BcesScale, Is.EqualTo(D(expected["bces_scale"])).Within(ClosedFormTolerance * Math.Abs(D(expected["bces_scale"]))));
        Assert.That(fit.ISquaredGx, Is.EqualTo(D(expected["i2_gx"])).Within(1e-12));
        Assert.That(fit.ConfidenceLow, Is.EqualTo(fit.Scale - 1.959963984540054 * fit.StandardError).Within(1e-12));
        Assert.That(fit.ConfidenceHigh, Is.EqualTo(fit.Scale + 1.959963984540054 * fit.StandardError).Within(1e-12));
    }

    [Test]
    public void WeakReferenceEffectsShrinkTheInverseVarianceSlopeButNotTheScale()
    {
        // weak_effects: true scale 0.8; the reference effects are mostly their own error (I²_GX 0.61), so a slope that
        // ignores that error is diluted to about half.
        var (x, sx, y, sy) = Inputs("weak_effects");

        EffectScaleResult fit = EffectScale.Estimate(x, sx, y, sy);

        Assert.That(fit.ISquaredGx, Is.LessThan(0.9));
        Assert.That(fit.InverseVarianceScale, Is.LessThan(0.5));
        Assert.That(0.8, Is.InRange(fit.ConfidenceLow, fit.ConfidenceHigh));
    }

    [Test]
    public void ClosedFormsOnAWorkedExample()
    {
        // x = 1, 2, 3 (se 0.1 each), y = 0.5, 1, 1.5 (se 0.2 each): y = x / 2 exactly.
        // Inverse-variance slope: Σ x y / s_y² ÷ Σ x² / s_y² = 7 / 14 = 0.5.
        // BCES: Σ x y / (Σ x² − Σ s_x²) = 7 / (14 − 0.03) = 0.501073729420186...
        // I²_GX: weighted mean of x = 2; Q = Σ (x − 2)² / 0.01 = 200; (Q − 2) / Q = 0.99.
        var fit = EffectScale.Estimate(new[] { 1.0, 2, 3 }, new[] { 0.1, 0.1, 0.1 }, new[] { 0.5, 1, 1.5 }, new[] { 0.2, 0.2, 0.2 });

        Assert.That(fit.InverseVarianceScale, Is.EqualTo(0.5).Within(1e-15));
        Assert.That(fit.BcesScale, Is.EqualTo(7 / 13.97).Within(1e-15));
        Assert.That(fit.ISquaredGx, Is.EqualTo(0.99).Within(1e-15));
        Assert.That(fit.Status, Is.EqualTo(EffectScaleStatus.Estimated));
        Assert.That(fit.Scale, Is.EqualTo(0.5).Within(1e-6), "no scatter beyond the standard errors: the scale is the line's");
    }

    [Test]
    public void BcesIsNaNWhenTheReferenceErrorIsAsLargeAsItsSpread()
    {
        var fit = EffectScale.Estimate(new[] { 0.1, -0.1, 0.05 }, new[] { 1.0, 1, 1 }, new[] { 0.2, 0.1, 0.3 }, new[] { 0.5, 0.5, 0.5 });

        Assert.That(fit.BcesScale, Is.NaN);
        Assert.That(fit.ISquaredGx, Is.EqualTo(0.0));
    }

    [Test]
    public void TooFewFeaturesOrNoReferenceEffectIsNotEstimable()
    {
        foreach (var fit in new[]
                 {
                     EffectScale.Estimate(new[] { 1.0, 2 }, new[] { 0.1, 0.1 }, new[] { 0.5, 1 }, new[] { 0.1, 0.1 }),
                     EffectScale.Estimate(new[] { 0.0, 0, 0 }, new[] { 0.1, 0.1, 0.1 }, new[] { 0.5, 1, 2 }, new[] { 0.1, 0.1, 0.1 }),
                 })
        {
            Assert.That(fit.Status, Is.EqualTo(EffectScaleStatus.TooFewFeatures));
            Assert.That(new[] { fit.Scale, fit.StandardError, fit.ConfidenceLow, fit.ConfidenceHigh, fit.Tau2, fit.Tau2StandardError,
                fit.InverseVarianceScale, fit.BcesScale, fit.ISquaredGx }, Is.All.NaN);
        }
    }

    [Test]
    public void BadInputIsRefused()
    {
        double[] ok = { 1, 2, 3 };
        double[] se = { 0.1, 0.1, 0.1 };
        Assert.Throws<ArgumentNullException>(() => EffectScale.Estimate(null!, se, ok, se));
        Assert.Throws<ArgumentException>(() => EffectScale.Estimate(ok, se, new[] { 1.0, 2 }, se));
        Assert.Throws<ArgumentException>(() => EffectScale.Estimate(new[] { 1, double.NaN, 3 }, se, ok, se));
        Assert.Throws<ArgumentException>(() => EffectScale.Estimate(ok, new[] { 0.1, 0, 0.1 }, ok, se));
        Assert.Throws<ArgumentException>(() => EffectScale.Estimate(ok, se, ok, new[] { 0.1, double.PositiveInfinity, 0.1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => EffectScale.Estimate(ok, se, ok, se, confidenceLevel: 1));
    }

    [Test]
    public void RescaleDividesByTheScaleAndCarriesItsUncertainty()
    {
        var (x, sx, y, sy) = Inputs("typical");
        EffectScaleResult fit = EffectScale.Estimate(x, sx, y, sy);

        var (estimate, standardError) = EffectScale.Rescale(0.35, 0.1, fit);

        double c = fit.Scale, varC = fit.StandardError * fit.StandardError;
        Assert.That(estimate, Is.EqualTo(0.35 / c).Within(1e-15));
        Assert.That(standardError, Is.EqualTo(Math.Sqrt(0.01 / (c * c) + 0.35 * 0.35 * varC / Math.Pow(c, 4))).Within(1e-15));
        Assert.That(standardError, Is.GreaterThan(0.1 / c), "the scale's own uncertainty widens the rescaled effect");
    }

    [Test]
    public void RescaleIsNaNWithoutAnEstimatedScale()
    {
        var notEstimable = EffectScale.Estimate(new[] { 1.0, 2 }, new[] { 0.1, 0.1 }, new[] { 0.5, 1 }, new[] { 0.1, 0.1 });

        var (estimate, standardError) = EffectScale.Rescale(0.35, 0.1, notEstimable);

        Assert.That(estimate, Is.NaN);
        Assert.That(standardError, Is.NaN);
    }
}
