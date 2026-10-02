using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Reflection;
using Chromatography.RetentionTimeCalibration;
using MzLibUtil;
using NUnit.Framework;

namespace Test.RetentionTimeCalibration;

/// <summary>
/// Calibrating one run onto a spectral library's iRT scale: anchors are (observed apex RT in minutes, library iRT) pairs
/// from a first pass's confident identifications, some of them wrong. The model maps run minutes to iRT and back, never
/// the library to minutes.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestIrtCalibration
{
    /// <summary>A nonlinear, increasing run-minutes to iRT curve over 1 to 40 minutes, with negative iRT at the start.</summary>
    private static double TrueIrt(double minutes) => -20 + 3.5 * minutes + 0.05 * minutes * minutes;

    private static List<(RtMinutes, Irt)> Anchors(Func<double, double> truth, int count, double noiseSd, double outlierFraction = 0,
        int seed = 1, double minRt = 1, double maxRt = 40)
    {
        var random = new Random(seed);
        var anchors = new List<(RtMinutes, Irt)>(count);
        for (int i = 0; i < count; i++)
        {
            double rt = minRt + random.NextDouble() * (maxRt - minRt);
            double irt = random.NextDouble() < outlierFraction
                ? -30 + random.NextDouble() * 200                      // a wrong identification: any iRT at all
                : truth(rt) + noiseSd * Gaussian(random);
            anchors.Add((new RtMinutes(rt), new Irt(irt)));
        }
        return anchors;
    }

    private static double Gaussian(Random random) =>
        Math.Sqrt(-2 * Math.Log(1 - random.NextDouble())) * Math.Cos(2 * Math.PI * random.NextDouble());

    private static double MaxInteriorError(IrtCalibrationModel model, Func<double, double> truth) =>
        Enumerable.Range(0, 100).Select(i => 3 + i * 0.34)
            .Max(rt => Math.Abs(model.ToIrt(new RtMinutes(rt)).Value - truth(rt)));

    [Test]
    public void AnExactLineIsRecoveredByEitherKind()
    {
        Func<double, double> line = rt => 6.25 * rt - 107;
        var anchors = Anchors(line, 200, noiseSd: 0);

        var linear = IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Kind: IrtCalibrationKind.Linear));
        var lowess = IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Kind: IrtCalibrationKind.Lowess));

        Assert.That(MaxInteriorError(linear, line), Is.LessThan(1e-6));
        Assert.That(MaxInteriorError(lowess, line), Is.LessThan(0.05));
    }

    /// <summary>A curved run is fitted by LOWESS; a straight line leaves an error the iRT window would have to absorb.</summary>
    [Test]
    public void LowessRecoversACurveThatALineCannot()
    {
        var anchors = Anchors(TrueIrt, 1000, noiseSd: 1.0);

        var lowess = IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Kind: IrtCalibrationKind.Lowess));
        var linear = IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Kind: IrtCalibrationKind.Linear));

        Assert.That(MaxInteriorError(lowess, TrueIrt), Is.LessThan(0.6));
        Assert.That(MaxInteriorError(linear, TrueIrt), Is.GreaterThan(3), "the curve is real; a line misses it");
    }

    /// <summary>A first pass always carries wrong identifications; 30% of them must not bend the curve.</summary>
    [Test]
    public void ThirtyPercentWrongAnchorsDoNotBendTheCurve()
    {
        var anchors = Anchors(TrueIrt, 1500, noiseSd: 1.0, outlierFraction: 0.3);

        var lowess = IrtCalibration.Fit(anchors);

        Assert.That(MaxInteriorError(lowess, TrueIrt), Is.LessThan(1.0));
        Assert.That(lowess.ResidualSd, Is.EqualTo(1.0).Within(0.25), "the robust residual SD reflects the good anchors, not the wrong ones");
    }

    [Test]
    public void ALinearFitIgnoresWrongAnchorsToo()
    {
        Func<double, double> line = rt => 6.25 * rt - 107;
        var model = IrtCalibration.Fit(Anchors(line, 800, noiseSd: 1.0, outlierFraction: 0.3),
            new IrtCalibrationOptions(Kind: IrtCalibrationKind.Linear));

        Assert.That(MaxInteriorError(model, line), Is.LessThan(0.5));
    }

    [Test]
    public void TheMapIsIncreasingAndInvertible()
    {
        var model = IrtCalibration.Fit(Anchors(TrueIrt, 1000, noiseSd: 2.0, outlierFraction: 0.1));

        double previous = double.NegativeInfinity;
        for (double rt = -10; rt <= 60; rt += 0.05)
        {
            var irt = model.ToIrt(new RtMinutes(rt));
            Assert.That(irt.Value, Is.GreaterThan(previous), $"increasing at {rt} min");
            previous = irt.Value;
            Assert.That(model.ToRtMinutes(irt).Value, Is.EqualTo(rt).Within(1e-6), $"round trip at {rt} min");
        }
    }

    /// <summary>Beyond the anchors the map continues as a straight line, and says so.</summary>
    [Test]
    public void OutsideTheAnchorsTheMapIsLinearAndFlagged()
    {
        var model = IrtCalibration.Fit(Anchors(TrueIrt, 1000, noiseSd: 0.5));

        double step1 = model.ToIrt(new RtMinutes(45)).Value - model.ToIrt(new RtMinutes(44)).Value;
        double step2 = model.ToIrt(new RtMinutes(55)).Value - model.ToIrt(new RtMinutes(54)).Value;
        Assert.That(step2, Is.EqualTo(step1).Within(1e-9));
        Assert.That(model.IsExtrapolated(new RtMinutes(45)), Is.True);
        Assert.That(model.IsExtrapolated(new RtMinutes(0.5)), Is.True);
        Assert.That(model.IsExtrapolated(new RtMinutes(20)), Is.False);
    }

    /// <summary>
    /// An isocratic hold: for five minutes nothing elutes at a new iRT, so noisy local fits dip. The map must still be
    /// increasing and invertible there.
    /// </summary>
    [Test]
    public void AFlatStretchStillGivesAnIncreasingInvertibleMap()
    {
        Func<double, double> hold = rt => rt < 18 ? 4 * rt : rt < 23 ? 72 : 72 + 4 * (rt - 23);
        var model = IrtCalibration.Fit(Anchors(hold, 1500, noiseSd: 3.0, seed: 5), new IrtCalibrationOptions(Knots: 200, Bandwidth: 0.05));

        double previous = double.NegativeInfinity;
        for (double rt = 1; rt <= 40; rt += 0.01)
        {
            double irt = model.ToIrt(new RtMinutes(rt)).Value;
            Assert.That(irt, Is.GreaterThan(previous), $"increasing at {rt:F2} min");
            Assert.That(model.ToRtMinutes(new Irt(irt)).Value, Is.EqualTo(rt).Within(1e-6));
            previous = irt;
        }
    }

    /// <summary>
    /// A gradient that bends sharply, with wrong anchors clustered where it bends. Per-bin medians smear the bend, so
    /// starting weights mistrust good anchors there. Robustness iterations re-weight against the fitted curve and recover it.
    /// </summary>
    [Test]
    public void RobustnessIterationsRecoverASharpBendThatMediansSmear()
    {
        Func<double, double> bend = rt => rt < 20 ? 2 * rt : 40 + 9 * (rt - 20);
        var anchors = Anchors(bend, 1500, noiseSd: 0.5, seed: 11);
        var random = new Random(12);
        for (int i = 0; i < 80; i++) // about 30% of the anchors in the bend: a local minority, as wrong IDs are
        {
            double rt = 18 + random.NextDouble() * 5;
            anchors.Add((new RtMinutes(rt), new Irt(bend(rt) + 25 + 10 * random.NextDouble())));
        }

        double Error(int iterations) => Enumerable.Range(0, 60).Select(i => 17 + i * 0.1)
            .Max(rt => Math.Abs(IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Bandwidth: 0.1, RobustnessIterations: iterations))
                .ToIrt(new RtMinutes(rt)).Value - bend(rt)));

        double withoutIterations = Error(0), withIterations = Error(3);
        Assert.That(withIterations, Is.LessThan(withoutIterations), $"iterations must help: {withIterations:F2} vs {withoutIterations:F2}");
        Assert.That(withIterations, Is.LessThan(4.0));
    }

    [Test]
    public void TheFitDoesNotDependOnAnchorOrder()
    {
        var anchors = Anchors(TrueIrt, 600, noiseSd: 1.0, outlierFraction: 0.2);
        var shuffled = anchors.OrderBy(a => a.Item2.Value).ToList();

        var a = IrtCalibration.Fit(anchors);
        var b = IrtCalibration.Fit(shuffled);

        for (double rt = 1; rt <= 40; rt += 1)
            Assert.That(b.ToIrt(new RtMinutes(rt)).Value, Is.EqualTo(a.ToIrt(new RtMinutes(rt)).Value).Within(1e-9));
    }

    [Test]
    public void TooFewAnchorsAreRefused()
    {
        var e = Assert.Throws<ArgumentException>(() => IrtCalibration.Fit(Anchors(TrueIrt, 19, noiseSd: 0)));
        Assert.That(e!.Message, Does.Contain("19"));

        Assert.DoesNotThrow(() => IrtCalibration.Fit(Anchors(TrueIrt, 20, noiseSd: 0)));
    }

    [Test]
    public void AnchorsThatSpanNoTimeAreRefused()
    {
        var sameTime = Enumerable.Range(0, 50).Select(i => (new RtMinutes(10), new Irt(i))).ToList();

        Assert.Throws<ArgumentException>(() => IrtCalibration.Fit(sameTime));
    }

    [Test]
    public void NonFiniteAnchorsAreRefused()
    {
        var anchors = Anchors(TrueIrt, 50, noiseSd: 0);
        anchors[7] = (new RtMinutes(double.NaN), new Irt(3));

        Assert.Throws<ArgumentException>(() => IrtCalibration.Fit(anchors));
    }

    [Test]
    public void OptionsAreChecked()
    {
        var anchors = Anchors(TrueIrt, 100, noiseSd: 0);

        Assert.Throws<ArgumentNullException>(() => IrtCalibration.Fit(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Bandwidth: 0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Bandwidth: 1.5)));
        Assert.Throws<ArgumentOutOfRangeException>(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Knots: 2)));
        Assert.Throws<ArgumentOutOfRangeException>(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(MinimumAnchors: 1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(RobustnessIterations: -1)));
    }

    /// <summary>The smallest legal values of each option are accepted, not only the defaults.</summary>
    [Test]
    public void TheSmallestLegalOptionsAreAccepted()
    {
        var anchors = Anchors(TrueIrt, 100, noiseSd: 0.5);

        Assert.DoesNotThrow(() => IrtCalibration.Fit(anchors.Take(3).ToList(), new IrtCalibrationOptions(Kind: IrtCalibrationKind.Linear, MinimumAnchors: 3)));
        Assert.DoesNotThrow(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Knots: 3)));
        Assert.DoesNotThrow(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(Bandwidth: 1)));
        Assert.DoesNotThrow(() => IrtCalibration.Fit(anchors, new IrtCalibrationOptions(RobustnessIterations: 0)));
    }

    /// <summary>The first and last anchors are inside the fitted range; only beyond them is extrapolated.</summary>
    [Test]
    public void TheAnchorsThemselvesAreNotExtrapolated()
    {
        var model = IrtCalibration.Fit(Anchors(TrueIrt, 200, noiseSd: 0.5));

        Assert.That(model.IsExtrapolated(model.FirstAnchor), Is.False);
        Assert.That(model.IsExtrapolated(model.LastAnchor), Is.False);
        Assert.That(model.IsExtrapolated(new RtMinutes(model.FirstAnchor.Value - 1e-9)), Is.True);
        Assert.That(model.IsExtrapolated(new RtMinutes(model.LastAnchor.Value + 1e-9)), Is.True);
    }

    [Test]
    public void TheModelReportsWhatItWasFittedOn()
    {
        var model = IrtCalibration.Fit(Anchors(TrueIrt, 300, noiseSd: 1.0, minRt: 2, maxRt: 38));

        Assert.That(model.AnchorCount, Is.EqualTo(300));
        Assert.That(model.Kind, Is.EqualTo(IrtCalibrationKind.Lowess));
        Assert.That(model.FirstAnchor.Value, Is.InRange(2, 2.5));
        Assert.That(model.LastAnchor.Value, Is.InRange(37.5, 38));
    }

    /// <summary>
    /// Before any anchors exist, a search needs a provisional map: the straight line through two known points (for example,
    /// the library's iRT range laid over the run's gradient).
    /// </summary>
    [Test]
    public void AProvisionalLineRunsThroughItsTwoPoints()
    {
        var line = IrtCalibration.Line((new RtMinutes(2), new Irt(-10)), (new RtMinutes(42), new Irt(110)));

        Assert.That(line.ToIrt(new RtMinutes(2)).Value, Is.EqualTo(-10).Within(1e-12));
        Assert.That(line.ToIrt(new RtMinutes(22)).Value, Is.EqualTo(50).Within(1e-12));
        Assert.That(line.ToRtMinutes(new Irt(170)).Value, Is.EqualTo(62).Within(1e-12));
        Assert.That(line.Kind, Is.EqualTo(IrtCalibrationKind.Linear));
        Assert.That(line.AnchorCount, Is.Zero);
        Assert.That(line.ResidualSd, Is.NaN, "a line through two stated points has measured nothing");
    }

    [Test]
    public void AProvisionalLineMustIncrease()
    {
        Assert.Throws<ArgumentException>(() => IrtCalibration.Line((new RtMinutes(2), new Irt(50)), (new RtMinutes(42), new Irt(10))));
        Assert.Throws<ArgumentException>(() => IrtCalibration.Line((new RtMinutes(2), new Irt(10)), (new RtMinutes(2), new Irt(50))));
    }

    /// <summary>iRT and run minutes are different quantities; only a calibration model converts between them.</summary>
    [Test]
    public void IrtAndMinutesDoNotConvert()
    {
        foreach (var type in new[] { typeof(Irt), typeof(RtMinutes) })
            Assert.That(type.GetMethods(BindingFlags.Public | BindingFlags.Static).Where(m => m.Name is "op_Implicit" or "op_Explicit"), Is.Empty, type.Name);
        Assert.That(new Irt(12.5).ToString(), Is.EqualTo("12.5 iRT"));
        Assert.That(new RtMinutes(3.25).ToString(), Is.EqualTo("3.25 min"));
    }
}
