using NUnit.Framework;
using StatisticalModels;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// Semi-supervised target-decoy rescoring (mProphet/Percolator-style) with a ridge-regularized Fisher linear discriminant.
/// The requirements that matter are the ones that make a reported FDR honest:
/// <list type="bullet">
/// <item>no candidate is scored by a model trained on it;</item>
/// <item>a group (for example one peptide) never straddles folds;</item>
/// <item>a run with nothing to find reports nothing.</item>
/// </list>
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TargetDecoyRescorerTests
{
    private static double Gaussian(Random random) =>
        Math.Sqrt(-2 * Math.Log(1 - random.NextDouble())) * Math.Cos(2 * Math.PI * random.NextDouble());

    /// <summary>
    /// Candidates with three features. True targets are shifted by <paramref name="shift"/> on each feature; false
    /// targets and decoys are not. No single feature separates well, but together they do.
    /// </summary>
    private static (double[][] Features, bool[] IsDecoy, string[] Groups, bool[] IsTrue) Candidates(
        int trueTargets, int falseTargets, int decoys, double shift, int seed = 1)
    {
        var random = new Random(seed);
        var features = new List<double[]>();
        var isDecoy = new List<bool>();
        var groups = new List<string>();
        var isTrue = new List<bool>();
        void Add(bool decoy, bool real, int index)
        {
            double s = real ? shift : 0;
            features.Add([s + Gaussian(random), s + Gaussian(random), s + Gaussian(random)]);
            isDecoy.Add(decoy);
            isTrue.Add(real);
            groups.Add($"{(decoy ? "D" : "T")}{index / 2}"); // two candidates (e.g. charge states) per peptide
        }
        for (int i = 0; i < trueTargets; i++) Add(false, true, i);
        for (int i = 0; i < falseTargets; i++) Add(false, false, trueTargets + i);
        for (int i = 0; i < decoys; i++) Add(true, false, i);
        return (features.ToArray(), isDecoy.ToArray(), groups.ToArray(), isTrue.ToArray());
    }

    /// <summary>Targets passing q ≤ cutoff under a score: plain target-decoy competition, (D+1)/T, monotone.</summary>
    private static int TargetsAtQ(double[] scores, bool[] isDecoy, double cutoff)
    {
        var order = Enumerable.Range(0, scores.Length).OrderByDescending(i => scores[i]).ToArray();
        var q = new double[order.Length];
        int d = 0, t = 0;
        for (int k = 0; k < order.Length; k++)
        {
            if (isDecoy[order[k]]) d++; else t++;
            q[k] = t == 0 ? 1 : Math.Min(1, (d + 1.0) / t);
        }
        for (int k = order.Length - 2; k >= 0; k--) q[k] = Math.Min(q[k], q[k + 1]);
        return Enumerable.Range(0, order.Length).Count(k => !isDecoy[order[k]] && q[k] <= cutoff);
    }

    #region LinearDiscriminant

    [Test]
    public void TheDiscriminantPointsFromDecoysToTargets()
    {
        var random = new Random(3);
        var features = new List<double[]>();
        var positive = new List<bool>();
        for (int i = 0; i < 500; i++) { features.Add([1 + Gaussian(random), 1 + Gaussian(random)]); positive.Add(true); }
        for (int i = 0; i < 500; i++) { features.Add([Gaussian(random), Gaussian(random)]); positive.Add(false); }

        var fit = LinearDiscriminant.Fit(features, positive);

        Assert.That(fit.Weights[0], Is.GreaterThan(0));
        Assert.That(fit.Weights[1], Is.GreaterThan(0));
        Assert.That(fit.Weights[0] / fit.Weights[1], Is.EqualTo(1).Within(0.25), "equal shifts, equal weights");
        Assert.That(fit.Score([2, 2]), Is.GreaterThan(fit.Score([-1, -1])));
    }

    /// <summary>A constant feature carries no information: weight 0, and never a throw or NaN.</summary>
    [Test]
    public void AConstantFeatureGetsNoWeight()
    {
        var random = new Random(4);
        var features = Enumerable.Range(0, 200).Select(i => new[] { (i < 100 ? 1 : 0) + Gaussian(random), 7.0 }).ToList();
        var positive = Enumerable.Range(0, 200).Select(i => i < 100).ToList();

        var fit = LinearDiscriminant.Fit(features, positive);

        Assert.That(fit.Weights[1], Is.EqualTo(0));
        Assert.That(double.IsFinite(fit.Score([1, 7])));
    }

    /// <summary>Semi-supervised positives are nearly separable by construction. The fit must still return a direction.</summary>
    [Test]
    public void PerfectlySeparableClassesStillGiveAFiniteDirection()
    {
        var features = Enumerable.Range(0, 100).Select(i => new[] { i < 50 ? 10.0 + i * 0.01 : -10.0 - i * 0.01, i * 0.1 }).ToList();
        var positive = Enumerable.Range(0, 100).Select(i => i < 50).ToList();

        var fit = LinearDiscriminant.Fit(features, positive);

        Assert.That(fit.Weights.All(double.IsFinite));
        Assert.That(fit.Weights[0], Is.GreaterThan(0));
    }

    [Test]
    public void DiscriminantArgumentsAreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => LinearDiscriminant.Fit(null!, [true]));
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([[1.0], [2.0]], [true]), "one label per row");
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([[1.0], [2.0]], [true, true]), "both classes needed");
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([[1.0, 2.0], [2.0]], [true, false]), "ragged rows");
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([[double.NaN], [2.0]], [true, false]), "non-finite");
    }

    #endregion

    #region TargetDecoyRescorer

    /// <summary>Combining weak features finds more targets at 1% than the best of them alone.</summary>
    [Test]
    public void RescoringBeatsTheBestSingleFeature()
    {
        var (features, isDecoy, groups, _) = Candidates(1500, 1500, 3000, shift: 1.2);

        var result = TargetDecoyRescorer.Score(features, isDecoy, groups);

        int best = Enumerable.Range(0, 3).Max(f => TargetsAtQ(features.Select(x => x[f]).ToArray(), isDecoy, 0.01));
        int rescored = TargetsAtQ(result.Scores, isDecoy, 0.01);
        Assert.That(result.Status, Is.EqualTo(RescoreStatus.Rescored));
        Assert.That(rescored, Is.GreaterThan(best * 1.5), $"rescored {rescored} vs best single feature {best}");
    }

    /// <summary>
    /// The leakage test: a held-out candidate's own label never reaches the model that scores it, so flipping that
    /// label cannot change its score.
    /// </summary>
    [Test]
    public void ACandidatesOwnLabelNeverReachesTheModelThatScoresIt()
    {
        var (features, isDecoy, groups, _) = Candidates(800, 800, 1600, shift: 1.2);
        var before = TargetDecoyRescorer.Score(features, isDecoy, groups);

        int row = 17;
        var flipped = isDecoy.ToArray();
        flipped[row] = !flipped[row];
        var after = TargetDecoyRescorer.Score(features, flipped, groups);

        Assert.That(after.Folds[row], Is.EqualTo(before.Folds[row]), "fold assignment depends on groups, not labels");
        Assert.That(after.Scores[row], Is.EqualTo(before.Scores[row]).Within(1e-12));
    }

    [Test]
    public void AGroupNeverStraddlesFoldsAndFoldsAreBalanced()
    {
        var (features, isDecoy, groups, _) = Candidates(900, 900, 1800, shift: 1.2);

        var result = TargetDecoyRescorer.Score(features, isDecoy, groups, folds: 3);

        foreach (var group in groups.Select((g, i) => (g, i)).GroupBy(x => x.g))
            Assert.That(group.Select(x => result.Folds[x.i]).Distinct().Count(), Is.EqualTo(1), $"group {group.Key}");
        var decoysPerFold = Enumerable.Range(0, 3).Select(f => Enumerable.Range(0, isDecoy.Length).Count(i => isDecoy[i] && result.Folds[i] == f)).ToArray();
        var targetsPerFold = Enumerable.Range(0, 3).Select(f => Enumerable.Range(0, isDecoy.Length).Count(i => !isDecoy[i] && result.Folds[i] == f)).ToArray();
        Assert.That(decoysPerFold.Max() - decoysPerFold.Min(), Is.LessThanOrEqualTo(2));
        Assert.That(targetsPerFold.Max() - targetsPerFold.Min(), Is.LessThanOrEqualTo(2));
    }

    /// <summary>With nothing real to find, rescoring must not manufacture discoveries.</summary>
    [Test]
    public void ARunWithNothingRealReportsAlmostNothing()
    {
        var (features, isDecoy, groups, _) = Candidates(0, 3000, 3000, shift: 0, seed: 9);

        var result = TargetDecoyRescorer.Score(features, isDecoy, groups);

        Assert.That(TargetsAtQ(result.Scores, isDecoy, 0.01), Is.LessThanOrEqualTo(30), "at most about 1% of 3000 targets");
    }

    [Test]
    public void RescoringIsDeterministicAndIndependentOfRowOrder()
    {
        var (features, isDecoy, groups, _) = Candidates(600, 600, 1200, shift: 1.2);
        var a = TargetDecoyRescorer.Score(features, isDecoy, groups);
        var b = TargetDecoyRescorer.Score(features, isDecoy, groups);

        var order = Enumerable.Range(0, features.Length).Reverse().ToArray();
        var c = TargetDecoyRescorer.Score(order.Select(i => features[i]).ToArray(), order.Select(i => isDecoy[i]).ToArray(), order.Select(i => groups[i]).ToArray());

        Assert.That(b.Scores, Is.EqualTo(a.Scores));
        for (int k = 0; k < order.Length; k++)
            Assert.That(c.Scores[k], Is.EqualTo(a.Scores[order[k]]).Within(1e-9), $"row {order[k]}");
    }

    /// <summary>Without decoys no q-value exists, so nothing is scored; the status says why.</summary>
    [Test]
    public void WithoutDecoysNothingIsScored()
    {
        var (features, _, groups, _) = Candidates(100, 100, 0, shift: 1.2);
        var isDecoy = new bool[features.Length];

        var result = TargetDecoyRescorer.Score(features, isDecoy, groups);

        Assert.That(result.Status, Is.EqualTo(RescoreStatus.NoDecoys));
        Assert.That(result.Scores.All(double.IsNaN));
    }

    /// <summary>The best single feature seeds the first iteration, even when higher values of it are worse.</summary>
    [Test]
    public void AFeatureWhereLowerIsBetterIsUsedTheRightWayRound()
    {
        var (features, isDecoy, groups, _) = Candidates(1500, 1500, 3000, shift: 1.2);
        var flipped = features.Select(x => new[] { -x[0], x[1], x[2] }).ToArray();

        int a = TargetsAtQ(TargetDecoyRescorer.Score(features, isDecoy, groups).Scores, isDecoy, 0.01);
        int b = TargetsAtQ(TargetDecoyRescorer.Score(flipped, isDecoy, groups).Scores, isDecoy, 0.01);

        Assert.That(b, Is.EqualTo(a).Within(a * 0.05 + 5), "negating a feature changes nothing a linear model cannot absorb");
    }

    [Test]
    public void RescorerArgumentsAreChecked()
    {
        var (features, isDecoy, groups, _) = Candidates(50, 50, 100, shift: 1.2);

        Assert.Throws<ArgumentNullException>(() => TargetDecoyRescorer.Score(null!, isDecoy, groups));
        Assert.Throws<ArgumentException>(() => TargetDecoyRescorer.Score(features, isDecoy.Take(10).ToArray(), groups));
        Assert.Throws<ArgumentException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups.Take(10).ToArray()));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, folds: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, iterations: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0));
    }

    #endregion
}
