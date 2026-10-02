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
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([[1.0, double.NaN], [2.0, 1.0]], [true, false]), "one non-finite value in a row");
        Assert.Throws<ArgumentNullException>(() => LinearDiscriminant.Fit([[1.0]], null!));
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([], []), "no observations");
        Assert.Throws<ArgumentOutOfRangeException>(() => LinearDiscriminant.Fit([[1.0], [2.0]], [true, false], ridge: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => LinearDiscriminant.Fit([[1.0], [2.0]], [true, false], ridge: double.NaN));
        Assert.DoesNotThrow(() => LinearDiscriminant.Fit([[1.0], [2.0], [0.0]], [true, false, false], ridge: 0));
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([[1.0], [2.0]], [true, false]).Score([1.0, 2.0]), "score length");
    }

    /// <summary>
    /// Without a ridge the weights are Fisher's, <c>S_w⁻¹ (μ₊ − μ₋)</c> in original units whatever the standardization, and
    /// the bias puts 0 midway between the class means. Worked by hand: both classes deviate from their means by
    /// (−1,−1), (1,0), (0,1), so S_w = [[1, ½], [½, 1]] (scatter over n − 2). With μ₊ = (4,2) and μ₋ = (1,1) this gives
    /// w = (10/3, −2/3) and a bias of −22/3.
    /// </summary>
    [Test]
    public void WithoutARidgeTheWeightsAreFishers()
    {
        double[][] features = [[3, 1], [5, 2], [4, 3], [0, 0], [2, 1], [1, 2]];
        bool[] positive = [true, true, true, false, false, false];

        var fit = LinearDiscriminant.Fit(features, positive, ridge: 0);

        Assert.That(fit.Weights[0], Is.EqualTo(10.0 / 3).Within(1e-9));
        Assert.That(fit.Weights[1], Is.EqualTo(-2.0 / 3).Within(1e-9), "correlation makes the second weight negative");
        Assert.That(fit.Bias, Is.EqualTo(-22.0 / 3).Within(1e-9));
        Assert.That(fit.Score([4, 2]), Is.EqualTo(14.0 / 3).Within(1e-9));
        Assert.That(fit.Score([1, 1]), Is.EqualTo(-14.0 / 3).Within(1e-9));
    }

    /// <summary>
    /// The ridge is added in standardized units, so in one dimension <c>w = Δ / (s_w + ridge · sd²)</c>, with sd the
    /// sample SD over all rows. For 3, 5 against 0, 2: Δ = 3, s_w = 2, sd² = 13/3.
    /// </summary>
    [Test]
    public void TheRidgeShrinksInStandardizedUnits()
    {
        double[][] features = [[3], [5], [0], [2]];
        bool[] positive = [true, true, false, false];

        var plain = LinearDiscriminant.Fit(features, positive, ridge: 0);
        var ridged = LinearDiscriminant.Fit(features, positive, ridge: 1);

        Assert.That(plain.Weights[0], Is.EqualTo(1.5).Within(1e-9));
        Assert.That(ridged.Weights[0], Is.EqualTo(9.0 / 19).Within(1e-9));
        Assert.That(ridged.Score([2.5]), Is.EqualTo(0).Within(1e-9), "0 midway between the class means");
    }

    /// <summary>With unbalanced classes the midpoint of the class means is not the overall mean, and 0 is still at the midpoint.</summary>
    [Test]
    public void ZeroIsMidwayBetweenUnbalancedClassMeans()
    {
        var fit = LinearDiscriminant.Fit([[3], [5], [0], [2], [1]], [true, true, false, false, false]);

        Assert.That(fit.Score([2.5]), Is.EqualTo(0).Within(1e-9));
        Assert.Throws<ArgumentException>(() => LinearDiscriminant.Fit([[1.0], [2.0], [3.0]], [true, false]), "more rows than labels");
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
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 1));
        Assert.Throws<ArgumentNullException>(() => TargetDecoyRescorer.Score(features, null!, groups));
        Assert.Throws<ArgumentNullException>(() => TargetDecoyRescorer.Score(features, isDecoy, null!));
        var ragged = features.ToArray();
        ragged[3] = [1.0, 2.0];
        Assert.Throws<ArgumentException>(() => TargetDecoyRescorer.Score(ragged, isDecoy, groups));
        var partlyNaN = features.ToArray();
        partlyNaN[3] = [1.0, double.NaN, 1.0];
        Assert.Throws<ArgumentException>(() => TargetDecoyRescorer.Score(partlyNaN, isDecoy, groups));
        Assert.DoesNotThrow(() => TargetDecoyRescorer.Score(features, isDecoy, groups, folds: 2, iterations: 1, positiveQValue: 0.5));
    }

    [Test]
    public void EmptyOrAllDecoyInputsReportWhyNothingWasScored()
    {
        var empty = TargetDecoyRescorer.Score([], [], []);
        Assert.That(empty.Status, Is.EqualTo(RescoreStatus.NoDecoys));
        Assert.That(empty.Scores, Is.Empty);

        var allDecoys = TargetDecoyRescorer.Score([[1.0], [2.0]], [true, true], ["a", "b"]);
        Assert.That(allDecoys.Status, Is.EqualTo(RescoreStatus.NoTargets));
        Assert.That(allDecoys.Scores.All(double.IsNaN));
    }

    /// <summary>
    /// With one target per fold, the fold holding it trains on the other one alone, which is too few to fit. That fold falls
    /// back to the best single feature and says so, and the scores stay finite and in the right order.
    /// </summary>
    [Test]
    public void AFoldWithTooFewPositivesFallsBackToTheBestFeature()
    {
        var features = new List<double[]> { new[] { 100.0 }, new[] { 90.0 } };
        var isDecoy = new List<bool> { false, false };
        for (int i = 0; i < 60; i++) { features.Add([i * 0.1]); isDecoy.Add(true); }
        var groups = Enumerable.Range(0, features.Count).Select(i => $"g{i:D3}").ToArray();

        var result = TargetDecoyRescorer.Score(features, isDecoy, groups);

        Assert.That(result.Status, Is.EqualTo(RescoreStatus.FoldStarved));
        Assert.That(result.Scores.All(double.IsFinite));
        double bestDecoy = result.Scores.Skip(2).Max();
        Assert.That(result.Scores[0], Is.GreaterThan(bestDecoy));
        Assert.That(result.Scores[1], Is.GreaterThan(bestDecoy));
    }

    /// <summary>
    /// Each fold's scores are normalized on its training rows: 0 at the q-value cutoff and −1 at the median decoy. Held-out
    /// decoys come from the same distribution, so in every fold their median lands near −1 and few clear 0.
    /// </summary>
    [Test]
    public void FoldScoresAreNormalizedToTheCutoffAndTheMedianDecoy()
    {
        var (features, isDecoy, groups, _) = Candidates(1500, 1500, 3000, shift: 1.5);

        var result = TargetDecoyRescorer.Score(features, isDecoy, groups);

        for (int f = 0; f < 3; f++)
        {
            var decoyScores = Enumerable.Range(0, features.Length).Where(i => isDecoy[i] && result.Folds[i] == f)
                .Select(i => result.Scores[i]).Order().ToArray();
            Assert.That(decoyScores[decoyScores.Length / 2], Is.EqualTo(-1).Within(0.1), $"fold {f} median decoy");
            Assert.That(decoyScores.Count(s => s >= 0), Is.LessThan(decoyScores.Length * 0.02), $"fold {f} decoys above the cutoff");
            int targetsAbove = Enumerable.Range(0, features.Length).Count(i => !isDecoy[i] && result.Folds[i] == f && result.Scores[i] >= 0);
            Assert.That(targetsAbove, Is.GreaterThan(50), $"fold {f} targets above the cutoff");        }
    }

    [Test]
    public void GroupsAreDealtToFoldsInKeyOrder()
    {
        Assert.That(TargetDecoyRescorer.AssignFolds(["d", "b", "a", "c", "a"], 3), Is.EqualTo(new[] { 0, 1, 0, 2, 0 }));
    }

    /// <summary>
    /// (D+1)/T down the ranking, capped at 1 and made monotone, reported in input order. Ranked: T×8, D, T, D. The eighth
    /// target has 1/8, which carries up to the top. The decoy and the ninth target then share min(2/8, 2/9) = 2/9. The
    /// final decoy's 3/9 stands.
    /// </summary>
    [Test]
    public void QValuesAreDecoysPlusOneOverTargetsMonotoneAndInInputOrder()
    {
        double[] ranked = [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1];
        bool[] decoy = [false, false, false, false, false, false, false, false, true, false, true];
        double[] expected = [.125, .125, .125, .125, .125, .125, .125, .125, 2.0 / 9, 2.0 / 9, 3.0 / 9];

        var order = new[] { 5, 0, 10, 3, 8, 1, 9, 2, 7, 4, 6 };
        double[] q = TargetDecoyRescorer.QValues(order.Select(i => ranked[i]).ToArray(), order.Select(i => decoy[i]).ToArray());

        for (int k = 0; k < order.Length; k++)
            Assert.That(q[k], Is.EqualTo(expected[order[k]]).Within(1e-12), $"rank {order[k]}");
    }

    [Test]
    public void QValuesAreCappedAtOne()
    {
        double[] q = TargetDecoyRescorer.QValues([4, 3, 2, 1], [true, false, true, true]);

        Assert.That(q, Is.EqualTo(new[] { 1.0, 1.0, 1.0, 1.0 }));
    }

    /// <summary>
    /// The training cutoff is the requested one when at least 10 targets pass it. Otherwise it relaxes to the first rung
    /// that passes 10, and to the last rung when none does. Decoys never count toward the 10.
    /// </summary>
    [Test]
    public void TheTrainingCutoffRelaxesOnlyUntilTenTargetsPass()
    {
        double[] Q(int count, double value) => Enumerable.Repeat(value, count).Concat(Enumerable.Repeat(0.001, 20)).ToArray();
        bool[] Decoys(int count) => Enumerable.Repeat(false, count).Concat(Enumerable.Repeat(true, 20)).ToArray();

        Assert.That(TargetDecoyRescorer.TrainingCutoff(Q(10, 0.01), Decoys(10), 0.01), Is.EqualTo(0.01), "exactly 10 pass");
        Assert.That(TargetDecoyRescorer.TrainingCutoff(Q(9, 0.01), Decoys(9), 0.01), Is.EqualTo(0.5), "9 pass; decoys don't count");
        Assert.That(TargetDecoyRescorer.TrainingCutoff(Q(12, 0.04), Decoys(12), 0.01), Is.EqualTo(0.05));
        Assert.That(TargetDecoyRescorer.TrainingCutoff(Q(12, 0.3), Decoys(12), 0.01), Is.EqualTo(0.5));
        Assert.That(TargetDecoyRescorer.TrainingCutoff(Q(12, 0.2), Decoys(12), 0.1), Is.EqualTo(0.25), "relaxes only above the request");
    }

    /// <summary>(mean of targets − mean of decoys) / SD of all values: for 3, 5 against 0, 2 that is 3 / √(13/3).</summary>
    [Test]
    public void StandardizedMeanDifferenceIsTargetsMinusDecoysOverTheSd()
    {
        Assert.That(TargetDecoyRescorer.StandardizedMeanDifference([3, 5, 0, 2], [false, false, true, true]),
            Is.EqualTo(3 / Math.Sqrt(13.0 / 3)).Within(1e-12));
        Assert.That(TargetDecoyRescorer.StandardizedMeanDifference([3, 5], [false, false]), Is.EqualTo(0), "no decoys");
        Assert.That(TargetDecoyRescorer.StandardizedMeanDifference([3, 5], [true, true]), Is.EqualTo(0), "no targets");
        Assert.That(TargetDecoyRescorer.StandardizedMeanDifference([2, 2, 2], [false, true, true]), Is.EqualTo(0), "constant");
    }

    /// <summary>
    /// With too few targets for any to pass, the seed is the feature and sign with the larger standardized separation. Here
    /// that is the second feature, taken negatively: −5 / 2.74 beats 2 / 1.26.
    /// </summary>
    [Test]
    public void TheSeedFeatureIsChosenBySeparationWhenNothingPasses()
    {
        double[][] features = [[1, 0], [2, 0], [3, 0], [0, 5], [0, 5], [0, 5]];
        bool[] decoy = [false, false, false, true, true, true];

        var seed = TargetDecoyRescorer.BestSingleFeature(features, decoy, [0, 1, 2, 3, 4, 5], 0.01);

        Assert.That(seed, Is.EqualTo((1, -1)));
    }

    /// <summary>
    /// Passing targets outrank separation. Feature 0 has a clean tail of 25 targets above everything, so 25 pass at 5%.
    /// Feature 1 separates the bulk better (1.50 against 1.34 standardized) but has 3 decoys on top, so nothing passes at
    /// 5%. Nothing passes at 1% on either, so the seed must climb to the 5% rung and choose feature 0. Staying at 1% would
    /// leave a tie that separation breaks the other way.
    /// </summary>
    [Test]
    public void TheSeedClimbsTheCutoffLadderBeforeFallingBackToSeparation()
    {
        var features = new List<double[]>();
        var decoy = new List<bool>();
        for (int i = 0; i < 40; i++) { features.Add([i < 25 ? 100 : i * 0.01, 5 + i * 0.01]); decoy.Add(false); }
        for (int i = 0; i < 40; i++) { features.Add([i * 0.01, i < 3 ? 10 : i * 0.01]); decoy.Add(true); }
        int[] all = Enumerable.Range(0, 80).ToArray();
        Assert.That(TargetDecoyRescorer.StandardizedMeanDifference(features.Select(x => x[1]).ToArray(), decoy.ToArray()),
            Is.GreaterThan(TargetDecoyRescorer.StandardizedMeanDifference(features.Select(x => x[0]).ToArray(), decoy.ToArray())),
            "fixture: feature 1 separates better");

        var seed = TargetDecoyRescorer.BestSingleFeature(features, decoy, all, 0.01);

        Assert.That(seed, Is.EqualTo((0, 1)));
    }

    /// <summary>When targets do pass, the seed is the feature and sign passing the most, here the second one reversed.</summary>
    [Test]
    public void TheSeedFeatureIsTheOnePassingTheMostTargets()
    {
        var features = new List<double[]>();
        var decoy = new List<bool>();
        for (int i = 0; i < 40; i++) { features.Add([i % 7, -100 - i]); decoy.Add(false); }
        for (int i = 0; i < 40; i++) { features.Add([i % 5, i]); decoy.Add(true); }

        var seed = TargetDecoyRescorer.BestSingleFeature(features, decoy, Enumerable.Range(0, 80).ToArray(), 0.01);

        Assert.That(seed, Is.EqualTo((1, -1)));
    }

    #endregion
}
