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
    public void TheDiscriminantIsPinnedToTheBit()
    {
        // Correlated features, one constant: pins weights and bias exactly, so a speed change that must change nothing cannot
        var random = new Random(21);
        var features = new List<double[]>();
        var positive = new List<bool>();
        for (int i = 0; i < 3000; i++)
        {
            bool target = i % 3 == 0;
            double shared = random.NextDouble();
            features.Add(Enumerable.Range(0, 12).Select(j => j == 5 ? 4.0 : shared * (j % 4) + random.NextDouble() * (1 + j) + (target ? 0.3 * j : 0)).ToArray());
            positive.Add(target);
        }

        var fit = LinearDiscriminant.Fit(features, positive);

        long[] bits = fit.Weights.Append(fit.Bias).Select(BitConverter.DoubleToInt64Bits).ToArray();
        TestContext.Out.WriteLine("pinned: " + string.Join(", ", bits));
        Assert.That(bits, Is.EqualTo(PinnedDiscriminantBits));
    }

    private static readonly long[] PinnedDiscriminantBits = [-4629063209439380192, 4590518683595660101, 4586681106985497765, 4593786861658844435,
        4603675463750120478, 0, 4599426323733484923, 4598695496080304840, 4599563267747816079, 4598620038832881131, 4598124823602177322,
        4597038961658343471, -4599363073103104677];

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

    /// <summary>
    /// With candidate groups (for example several candidate peaks of one precursor), only each group's top-scoring row
    /// trains the model, as in pyProphet. Every row is still scored. So extra low-ranked rows added to some groups change
    /// no other row's score.
    /// </summary>
    [Test]
    public void OnlyEachCandidateGroupsTopRowTrainsTheModel()
    {
        var (features, isDecoy, groups, _) = Candidates(800, 800, 1600, shift: 1.2);
        int n = features.Length;
        int[] candidateGroups = Enumerable.Range(0, n).ToArray();
        var before = TargetDecoyRescorer.Score(features, isDecoy, groups, candidateGroups: candidateGroups);

        // A second, far worse candidate for every tenth row, in the same sequence group and candidate group
        var extra = Enumerable.Range(0, n).Where(i => i % 10 == 0).ToArray();
        var moreFeatures = features.Concat(extra.Select(i => features[i].Select(v => v - 20).ToArray())).ToArray();
        var moreDecoy = isDecoy.Concat(extra.Select(i => isDecoy[i])).ToArray();
        var moreGroups = groups.Concat(extra.Select(i => groups[i])).ToArray();
        var moreCandidates = candidateGroups.Concat(extra).ToArray();
        var after = TargetDecoyRescorer.Score(moreFeatures, moreDecoy, moreGroups, candidateGroups: moreCandidates);

        for (int i = 0; i < n; i++)
            Assert.That(after.Scores[i], Is.EqualTo(before.Scores[i]).Within(1e-9), $"row {i}");
        Assert.That(extra.Select((row, k) => after.Scores[n + k] < after.Scores[row]), Is.All.True, "the added rows are scored, and lower");
        Assert.That(after.Status, Is.EqualTo(RescoreStatus.Rescored));
    }

    /// <summary>When every row is its own candidate group, nothing changes.</summary>
    [Test]
    public void SingletonCandidateGroupsChangeNothing()
    {
        var (features, isDecoy, groups, _) = Candidates(600, 600, 1200, shift: 1.2);

        var plain = TargetDecoyRescorer.Score(features, isDecoy, groups);
        var grouped = TargetDecoyRescorer.Score(features, isDecoy, groups, candidateGroups: Enumerable.Range(0, features.Length).ToArray());

        Assert.That(grouped.Scores, Is.EqualTo(plain.Scores).Within(1e-12));
        Assert.Throws<ArgumentException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, candidateGroups: new int[3]));
    }

    /// <summary>
    /// The confident training sample, with a positive q cutoff: only targets that pass it (by the linear ranking) are positives,
    /// at most half the cap, and as many of the top-ranked decoys. Without a cutoff it is the top half-cap of each, as before.
    /// If no target passes, the sample falls back to that.
    /// </summary>
    [Test]
    public void TheConfidentSampleCanKeepOnlyTargetsThatPassACutoff()
    {
        int[] ranked = Enumerable.Range(0, 12).ToArray();
        double[] scores = ranked.Select(i => 12.0 - i).ToArray();
        bool[] decoy = [false, false, false, true, false, true, false, true, true, false, true, false];
        double[] q = TargetDecoyRescorer.QValues(scores, decoy);
        int[] passing = ranked.Where(i => !decoy[i] && q[i] <= 0.4).ToArray();
        Assert.That(passing, Is.Not.Empty.And.Length.LessThan(ranked.Count(i => !decoy[i])), "the cutoff must bite in this example");

        int[] withCutoff = TargetDecoyRescorer.ConfidentTrainingRows(ranked, scores, decoy, cap: 10, positiveQValue: 0.4);
        int[] withoutCutoff = TargetDecoyRescorer.ConfidentTrainingRows(ranked, scores, decoy, cap: 10, positiveQValue: null);
        int[] noneCanPass = TargetDecoyRescorer.ConfidentTrainingRows(ranked, scores, decoy, cap: 10, positiveQValue: 1e-9);

        Assert.That(withCutoff.Where(i => !decoy[i]), Is.EquivalentTo(passing.Take(5)));
        Assert.That(withCutoff.Where(i => decoy[i]), Is.EquivalentTo(ranked.Where(i => decoy[i]).Take(Math.Min(5, passing.Length))));
        Assert.That(withoutCutoff, Is.EquivalentTo(ranked.Where(i => !decoy[i]).Take(5).Concat(ranked.Where(i => decoy[i]).Take(5))));
        Assert.That(withCutoff, Is.Ordered);
        Assert.That(noneCanPass, Is.EquivalentTo(withoutCutoff), "no passing target falls back to the plain confident sample");
    }

    /// <summary>
    /// The pooled confident sample takes the top rows of the ranking whatever their label: one score threshold, so every
    /// score band keeps targets and decoys in their own proportion. Equal halves leave a band above the decoys' cutoff and
    /// below the targets' where only decoys train, and the network learns that faint means decoy.
    /// </summary>
    [Test]
    public void ThePooledConfidentSampleTakesTheTopRowsWhateverTheirLabel()
    {
        int[] ranked = Enumerable.Range(0, 20).ToArray();
        double[] scores = ranked.Select(i => 20.0 - i).ToArray();
        // Targets dominate the top of the ranking, as real ones lift them
        bool[] decoy = ranked.Select(i => i >= 6 && i % 2 == 1).ToArray();

        int[] pooled = TargetDecoyRescorer.PooledTrainingRows(ranked, decoy, cap: 10);
        Assert.That(pooled, Is.EqualTo(ranked.Take(10).Order()));
        Assert.That(pooled.Count(i => !decoy[i]), Is.GreaterThan(5), "more targets than half when targets lead the ranking");

        var random = new Random(5);
        var features = ranked.Select(i => new[] { random.NextDouble() + (decoy[i] ? 0 : 0.5), random.NextDouble() }).ToList();
        var keys = ranked.Select(i => i.ToString()).ToList();
        var result = TargetDecoyRescorer.Score(features, decoy, keys, model: RescoreModel.NeuralNetworkEnsemble, maxNetworkTrainingRows: 10,
            networkTrainingSample: NetworkTrainingSample.ConfidentPooled, networkMembers: 2, networkEpochs: 2);
        Assert.That(result.Scores, Has.All.Matches<double>(double.IsFinite));
    }

    /// <summary>
    /// When the top rows of the ranking are all one label the network cannot train on them, so the pooled sample falls back
    /// to the top targets and the top decoys, half the cap each.
    /// </summary>
    [Test]
    public void ThePooledSampleWithOneLabelFallsBackToEqualHalves()
    {
        int[] ranked = Enumerable.Range(0, 20).ToArray();
        bool[] decoy = ranked.Select(i => i >= 12).ToArray();

        int[] pooled = TargetDecoyRescorer.PooledTrainingRows(ranked, decoy, cap: 10);
        Assert.That(pooled, Is.EqualTo(new[] { 0, 1, 2, 3, 4, 12, 13, 14, 15, 16 }));

        var random = new Random(5);
        var features = ranked.Select(i => new[] { random.NextDouble() + (decoy[i] ? 0 : 2), random.NextDouble() }).ToList();
        var keys = ranked.Select(i => i.ToString()).ToList();
        var result = TargetDecoyRescorer.Score(features, decoy, keys, model: RescoreModel.NeuralNetworkEnsemble, maxNetworkTrainingRows: 4,
            networkTrainingSample: NetworkTrainingSample.ConfidentPooled, networkMembers: 2, networkEpochs: 2);
        Assert.That(result.Scores, Has.All.Matches<double>(double.IsFinite));
    }

    /// <summary>
    /// The confident sample's share of targets can be set: at 0.7 a cap of 10 takes the top 7 targets and the top 3 decoys,
    /// at 0.5 the halves as before, and a share outside (0, 1) is refused.
    /// </summary>
    [Test]
    public void TheConfidentSampleTargetShareCanBeSet()
    {
        int[] ranked = Enumerable.Range(0, 20).ToArray();
        double[] scores = ranked.Select(i => 20.0 - i).ToArray();
        bool[] decoy = ranked.Select(i => i % 2 == 1).ToArray();

        int[] seventy = TargetDecoyRescorer.ConfidentTrainingRows(ranked, scores, decoy, cap: 10, positiveQValue: null, targetFraction: 0.7);
        Assert.That(seventy, Is.EquivalentTo(ranked.Where(i => !decoy[i]).Take(7).Concat(ranked.Where(i => decoy[i]).Take(3))));
        Assert.That(TargetDecoyRescorer.ConfidentTrainingRows(ranked, scores, decoy, cap: 10, positiveQValue: null, targetFraction: 0.5),
            Is.EqualTo(TargetDecoyRescorer.ConfidentTrainingRows(ranked, scores, decoy, cap: 10, positiveQValue: null)));

        var features = ranked.Select(i => new[] { (double)i }).ToList();
        var keys = ranked.Select(i => i.ToString()).ToList();
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, decoy, keys, networkTargetFraction: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, decoy, keys, networkTargetFraction: 1));
    }

    /// <summary>
    /// The network's hidden layers can be chosen (DIA-NN 2020's 25-20-15-10-5 by default): another architecture gives other
    /// scores, a wider one still finds non-linear signal a line misses, and an empty or zero-unit layer list is refused.
    /// </summary>
    [Test]
    public void TheNetworkArchitectureCanBeChosen()
    {
        var random = new Random(21);
        var features = new List<double[]>();
        var isDecoy = new List<bool>();
        var groups = new List<string>();
        void Add(bool decoy, bool real, int i)
        {
            double spread = real ? 3.0 : 1.0;
            features.Add([spread * Gaussian(random), spread * Gaussian(random), Gaussian(random)]);
            isDecoy.Add(decoy);
            groups.Add($"{(decoy ? "D" : "T")}{i}");
        }
        for (int i = 0; i < 2000; i++) Add(false, true, i);
        for (int i = 0; i < 2000; i++) Add(false, false, 2000 + i);
        for (int i = 0; i < 4000; i++) Add(true, false, i);
        bool[] decoys = isDecoy.ToArray();

        var standard = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble, networkMembers: 3);
        var wide = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble, networkMembers: 3,
            networkLayers: [64, 32, 16]);
        var linear = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15);

        Assert.That(wide.Scores, Is.Not.EqualTo(standard.Scores));
        Assert.That(TargetsAtQ(wide.Scores, decoys, 0.01), Is.GreaterThan(TargetsAtQ(linear.Scores, decoys, 0.01) + 200));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, model: RescoreModel.NeuralNetworkEnsemble, networkLayers: []));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, model: RescoreModel.NeuralNetworkEnsemble, networkLayers: [10, 0]));
    }

    /// <summary>
    /// True targets that differ from decoys only non-linearly (a larger spread on two features, the same mean): the network
    /// ensemble finds them, a line cannot. It trains only on each fold's training rows.
    /// </summary>
    [Test]
    public void TheNetworkModelFindsNonLinearSignalTheLineMisses()
    {
        var random = new Random(21);
        var features = new List<double[]>();
        var isDecoy = new List<bool>();
        var groups = new List<string>();
        void Add(bool decoy, bool real, int i)
        {
            double spread = real ? 3.0 : 1.0;
            features.Add([spread * Gaussian(random), spread * Gaussian(random), Gaussian(random)]);
            isDecoy.Add(decoy);
            groups.Add($"{(decoy ? "D" : "T")}{i}");
        }
        for (int i = 0; i < 2000; i++) Add(false, true, i);
        for (int i = 0; i < 2000; i++) Add(false, false, 2000 + i);
        for (int i = 0; i < 4000; i++) Add(true, false, i);
        bool[] decoys = isDecoy.ToArray();

        var linear = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15);
        var network = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble);

        int atOnePercentLinear = TargetsAtQ(linear.Scores, decoys, 0.01);
        int atOnePercentNetwork = TargetsAtQ(network.Scores, decoys, 0.01);
        Assert.That(atOnePercentNetwork, Is.GreaterThan(atOnePercentLinear + 200), $"network {atOnePercentNetwork} vs line {atOnePercentLinear}");
    }

    /// <summary>
    /// The leakage test for the network model, with candidate groups: flipping one held-out row's label, or adding a
    /// row to its candidate group, never changes how that row is scored. Only other folds train the model that scores it.
    /// </summary>
    [Test]
    public void ACandidatesOwnLabelNeverReachesTheNetworkThatScoresIt()
    {
        var (features, isDecoy, groups, _) = Candidates(800, 800, 1600, shift: 1.2);
        int[] candidateGroups = Enumerable.Range(0, features.Length).ToArray();
        var before = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, candidateGroups: candidateGroups,
            model: RescoreModel.NeuralNetworkEnsemble);

        int row = 17;
        var flipped = isDecoy.ToArray();
        flipped[row] = !flipped[row];
        var after = TargetDecoyRescorer.Score(features, flipped, groups, positiveQValue: 0.15, candidateGroups: candidateGroups,
            model: RescoreModel.NeuralNetworkEnsemble);

        Assert.That(after.Folds[row], Is.EqualTo(before.Folds[row]));
        Assert.That(after.Scores[row], Is.EqualTo(before.Scores[row]).Within(1e-9), "the flipped label stayed out of the model scoring it");
        // Every row of the flipped row's fold is scored by the same model, so all of them are unchanged too
        foreach (int i in Enumerable.Range(0, features.Length).Where(i => before.Folds[i] == before.Folds[row]))
            Assert.That(after.Scores[i], Is.EqualTo(before.Scores[i]).Within(1e-9), $"row {i} in the same fold");
    }

    /// <summary>The fixture for the network's training cap: real targets spread wider than decoys (no line separates them).</summary>
    private static (List<double[]> Features, List<bool> IsDecoy, List<string> Groups) Ring()
    {
        var random = new Random(21);
        var features = new List<double[]>();
        var isDecoy = new List<bool>();
        var groups = new List<string>();
        void Add(bool decoy, bool real, int i)
        {
            double spread = real ? 3.0 : 1.0;
            features.Add([spread * Gaussian(random), spread * Gaussian(random), Gaussian(random)]);
            isDecoy.Add(decoy);
            groups.Add($"{(decoy ? "D" : "T")}{i}");
        }
        for (int i = 0; i < 2000; i++) Add(false, true, i);
        for (int i = 0; i < 2000; i++) Add(false, false, 2000 + i);
        for (int i = 0; i < 4000; i++) Add(true, false, i);
        return (features, isDecoy, groups);
    }

    /// <summary>
    /// The network can train on a random subsample of each fold's training rows. DIA-NN trains on 267k of 3M precursors;
    /// on a whole-proteome search, training on all 2.6M rows was two thirds of the time. A cap at least the training size
    /// changes nothing. (Taking the top rows by the linear score instead was tried first and failed ACappedNetworkStillFindsTheSignal:
    /// where the line misses the signal, its top rows are the wrong ones.)
    /// </summary>
    [Test]
    public void ANetworkTrainingCapAtLeastTheTrainingSizeChangesNothing()
    {
        var (features, isDecoy, groups) = Ring();

        var all = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble);
        var capped = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble,
            maxNetworkTrainingRows: 1_000_000);

        Assert.That(capped.Scores, Is.EqualTo(all.Scores));
    }

    /// <summary>A cap well below the training size is applied, and the network still finds what the line misses.</summary>
    [Test]
    public void ACappedNetworkStillFindsTheSignal()
    {
        var (features, isDecoy, groups) = Ring();
        bool[] decoys = isDecoy.ToArray();

        var linear = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15);
        var all = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble);
        var capped = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble,
            maxNetworkTrainingRows: 2000);

        Assert.That(capped.Scores, Is.Not.EqualTo(all.Scores), "about 5,300 training rows per fold, capped at 2,000");
        Assert.That(TargetsAtQ(capped.Scores, decoys, 0.01), Is.GreaterThan(TargetsAtQ(linear.Scores, decoys, 0.01) + 200));
    }

    /// <summary>The leakage guarantee holds with a cap: the cap only ever picks among the fold's own training rows.</summary>
    [Test]
    public void ACandidatesOwnLabelNeverReachesACappedNetwork()
    {
        var (features, isDecoy, groups, _) = Candidates(800, 800, 1600, shift: 1.2);
        var before = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble,
            maxNetworkTrainingRows: 600);

        int row = 17;
        var flipped = isDecoy.ToArray();
        flipped[row] = !flipped[row];
        var after = TargetDecoyRescorer.Score(features, flipped, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble,
            maxNetworkTrainingRows: 600);

        foreach (int i in Enumerable.Range(0, features.Length).Where(i => before.Folds[i] == before.Folds[row]))
            Assert.That(after.Scores[i], Is.EqualTo(before.Scores[i]).Within(1e-9), $"row {i} in the flipped row's fold");
    }

    /// <summary>
    /// The network's random seed is a parameter, so a caller can measure how much a result moves for reasons that are only
    /// random: the same data with another seed. The default seed reproduces earlier results exactly.
    /// </summary>
    [Test]
    public void TheNetworkSeedIsAParameterAndTheDefaultIsUnchanged()
    {
        var (features, isDecoy, groups) = Ring();

        var before = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble);
        var seed0 = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble, randomSeed: 0);
        var seed1 = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble, randomSeed: 1);
        var seed1Again = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble, randomSeed: 1);

        Assert.That(seed0.Scores, Is.EqualTo(before.Scores));
        Assert.That(seed1.Scores, Is.Not.EqualTo(seed0.Scores));
        Assert.That(seed1Again.Scores, Is.EqualTo(seed1.Scores), "a seed is deterministic");
    }

    /// <summary>
    /// The ensemble's size and training length are parameters (DIA-NN uses 12 networks for one epoch; our default is 5 for 10).
    /// The defaults reproduce earlier results exactly; other values change the scores and still find what a line misses.
    /// </summary>
    [Test]
    public void EnsembleSizeAndEpochsAreParameters()
    {
        var (features, isDecoy, groups) = Ring();
        bool[] decoys = isDecoy.ToArray();

        var before = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble);
        var defaults = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble,
            networkMembers: 5, networkEpochs: 10);
        var diann = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble,
            networkMembers: 12, networkEpochs: 1);

        Assert.That(defaults.Scores, Is.EqualTo(before.Scores));
        Assert.That(diann.Scores, Is.Not.EqualTo(before.Scores));
        // One epoch underfits a fixture this small (8,000 rows: 202 targets at 1% against the line's 266); DIA-NN trains on
        // far more rows, so whether one epoch suits a real search is for the real-data comparison to show
        Assert.That(TargetsAtQ(diann.Scores, decoys, 0.01), Is.GreaterThan(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, model: RescoreModel.NeuralNetworkEnsemble, networkMembers: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, model: RescoreModel.NeuralNetworkEnsemble, networkEpochs: 0));
    }

    [Test]
    public void ANetworkTrainingCapBelowTwoIsRefused()
    {
        var (features, isDecoy, groups) = Ring();
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups,
            model: RescoreModel.NeuralNetworkEnsemble, maxNetworkTrainingRows: 1));
    }

    /// <summary>The network model must not manufacture discoveries: with nothing real, about 1% at most, as for the line.</summary>
    [Test]
    public void TheNetworkModelFindsNothingWhereThereIsNothing()
    {
        var (features, isDecoy, groups, _) = Candidates(0, 3000, 3000, shift: 0, seed: 9);

        var result = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, model: RescoreModel.NeuralNetworkEnsemble);

        Assert.That(TargetsAtQ(result.Scores, isDecoy, 0.01), Is.LessThanOrEqualTo(30));
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

    #region second network pass

    /// <summary>
    /// The network trains on each candidate group's top row, picked by the linear model. Where the line cannot tell the real
    /// candidate from noise (here the real one only spreads wider), the network learns from mostly wrong rows. A second pass,
    /// as DIA-NN trains its networks twice, re-picks each group's top row with the first network and trains again.
    /// </summary>
    [Test]
    public void ASecondNetworkPassTrainsOnTheRowsTheFirstNetworkPicked()
    {
        var random = new Random(31);
        var features = new List<double[]>();
        var isDecoy = new List<bool>();
        var groups = new List<string>();
        var candidateGroups = new List<int>();
        int group = 0;
        void AddGroup(bool decoy, bool hasReal, int i)
        {
            for (int c = 0; c < 3; c++)
            {
                double spread = hasReal && c == 0 ? 3.0 : 1.0;
                features.Add([spread * Gaussian(random), spread * Gaussian(random), Gaussian(random)]);
                isDecoy.Add(decoy);
                groups.Add($"{(decoy ? "D" : "T")}{i}");
                candidateGroups.Add(group);
            }
            group++;
        }
        for (int i = 0; i < 1500; i++) AddGroup(false, true, i);
        for (int i = 0; i < 1500; i++) AddGroup(false, false, 1500 + i);
        for (int i = 0; i < 3000; i++) AddGroup(true, false, i);

        int GroupsAtOnePercent(double[] scores)
        {
            var best = Enumerable.Range(0, group).Select(g => Enumerable.Range(3 * g, 3).Max(r => scores[r])).ToArray();
            var decoyGroup = Enumerable.Range(0, group).Select(g => isDecoy[3 * g]).ToArray();
            return TargetsAtQ(best, decoyGroup, 0.01);
        }
        var once = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, candidateGroups: candidateGroups,
            model: RescoreModel.NeuralNetworkEnsemble);
        var twice = TargetDecoyRescorer.Score(features, isDecoy, groups, positiveQValue: 0.15, candidateGroups: candidateGroups,
            model: RescoreModel.NeuralNetworkEnsemble, networkPasses: 2);

        int first = GroupsAtOnePercent(once.Scores), second = GroupsAtOnePercent(twice.Scores);
        Assert.That(second, Is.GreaterThan(first + 100), $"two passes {second} vs one {first}");
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, networkPasses: 0));
    }

    #endregion

    #region confident network training sample

    /// <summary>
    /// A whole-proteome DIA search has about one real target in twenty, so a random training sample teaches the network
    /// mostly noise against decoys. DIA-NN removes low-confidence identifications first and trains on about 156k targets and
    /// 138k decoys. The confident sample takes half its rows from the targets the line ranks highest and half from the decoys
    /// it ranks highest.
    /// </summary>
    [Test]
    public void AConfidentTrainingSampleBeatsARandomOneWhenFewTargetsAreReal()
    {
        var (features, isDecoy, groups, _) = Candidates(500, 9500, 10000, shift: 2.5);

        var random = TargetDecoyRescorer.Score(features, isDecoy, groups, model: RescoreModel.NeuralNetworkEnsemble,
            maxNetworkTrainingRows: 600);
        var confident = TargetDecoyRescorer.Score(features, isDecoy, groups, model: RescoreModel.NeuralNetworkEnsemble,
            maxNetworkTrainingRows: 600, networkTrainingSample: NetworkTrainingSample.Confident);

        int byRandom = TargetsAtQ(random.Scores, isDecoy, 0.01), byConfident = TargetsAtQ(confident.Scores, isDecoy, 0.01);
        Assert.That(byConfident, Is.GreaterThan(byRandom + 30), $"confident {byConfident} vs random {byRandom}");
    }

    #endregion

    #region sampled normalisation

    /// <summary>
    /// Each fold scored every training row only to set its normalisation (the 1% threshold and median decoy), twice the work
    /// of scoring its own held-out rows. A random sample of training groups estimates both. Within a fold the scores remain
    /// an affine map of the full normalisation's, so the fold's ranking is unchanged; only the pooling across folds can move.
    /// </summary>
    [Test]
    public void ASampleOfTrainingGroupsNormalisesEachFold()
    {
        var (features, isDecoy, groups, _) = Candidates(3000, 3000, 6000, shift: 1.5);

        var full = TargetDecoyRescorer.Score(features, isDecoy, groups);
        var sampled = TargetDecoyRescorer.Score(features, isDecoy, groups, normalizationGroups: 4000);

        for (int f = 0; f < 3; f++)
        {
            int[] rows = Enumerable.Range(0, features.Length).Where(i => full.Folds[i] == f).ToArray();
            var a = rows.Select(i => full.Scores[i]).ToArray();
            var b = rows.Select(i => sampled.Scores[i]).ToArray();
            double slope = (b[1] - b[0]) / (a[1] - a[0]);
            for (int k = 0; k < rows.Length; k++)
                Assert.That(b[k], Is.EqualTo(b[0] + slope * (a[k] - a[0])).Within(1e-9 * (1 + Math.Abs(b[k]))), $"fold {f}, row {rows[k]}");
        }
        int byFull = TargetsAtQ(full.Scores, isDecoy, 0.01), bySample = TargetsAtQ(sampled.Scores, isDecoy, 0.01);
        // The threshold is the lowest passing training score, an extreme value: on 4,000 groups its noise moves the pooled count
        // by several percent either way. A real search samples hundreds of thousands, and the real-data A/B judges it.
        Assert.That(bySample, Is.EqualTo(byFull).Within(0.1 * byFull));
        Assert.Throws<ArgumentOutOfRangeException>(() => TargetDecoyRescorer.Score(features, isDecoy, groups, normalizationGroups: 0));
    }

    #endregion
}
