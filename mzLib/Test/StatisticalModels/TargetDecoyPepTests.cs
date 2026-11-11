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
/// Posterior error probabilities read off the target-decoy distribution of a score. What matters is that the number is a
/// probability that a match is wrong: it agrees with an independent isotonic regression, it recovers the true PEP of a
/// known mixture, and it is a function of the score alone (ties share it, row order never changes it).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TargetDecoyPepTests
{
    private static double Gaussian(Random random) =>
        Math.Sqrt(-2 * Math.Log(1 - random.NextDouble())) * Math.Cos(2 * Math.PI * random.NextDouble());

    private static string[][] Fixture(string file) =>
        File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "StatisticalModels", "ReferenceData", file))
            .Skip(1).Select(l => l.Split('\t')).ToArray();

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    [Test]
    public void TheIsotonicFitMatchesScikitLearn()
    {
        // make_target_decoy_pep_fixtures.py, PROVENANCE_target_decoy_pep.txt
        string[][] rows = Fixture("target_decoy_pep_isotonic.tsv");
        double[] fitted = TargetDecoyPep.IsotonicNonDecreasing(rows.Select(r => D(r[0])).ToArray(), rows.Select(r => D(r[1])).ToArray());
        for (int i = 0; i < rows.Length; i++)
            Assert.That(fitted[i], Is.EqualTo(D(rows[i][2])).Within(1e-12), $"position {i}");
    }

    [Test]
    public void TheEstimateMatchesScikitLearnOnASimulatedSearch()
    {
        // Ties (scores to two decimals) and half-decoy candidates included.
        string[][] rows = Fixture("target_decoy_pep_search.tsv");
        var result = TargetDecoyPep.Estimate(rows.Select(r => D(r[0])).ToArray(), rows.Select(r => D(r[1])).ToArray());
        Assert.That(result.Status, Is.EqualTo(TargetDecoyPepStatus.Estimated));
        for (int i = 0; i < rows.Length; i++)
            Assert.That(result.Peps[i], Is.EqualTo(D(rows[i][2])).Within(1e-12), $"row {i}");
    }

    [Test]
    public void OnAKnownMixtureTheEstimateTracksTheTruePosteriorErrorProbability()
    {
        // 80,000 correct targets ~ N(3, 1); 20,000 incorrect targets and 20,000 decoys ~ N(0, 1). For a target at s the true
        // PEP is 20000 phi(s) / (20000 phi(s) + 80000 phi(s - 3)).
        var random = new Random(7);
        var scores = new List<double>();
        var isDecoy = new List<bool>();
        var truePep = new List<double>();
        double Phi(double x) => Math.Exp(-0.5 * x * x);
        double True(double s) => 20000 * Phi(s) / (20000 * Phi(s) + 80000 * Phi(s - 3));
        for (int i = 0; i < 80000; i++) { double s = 3 + Gaussian(random); scores.Add(s); isDecoy.Add(false); truePep.Add(True(s)); }
        for (int i = 0; i < 20000; i++) { double s = Gaussian(random); scores.Add(s); isDecoy.Add(false); truePep.Add(True(s)); }
        for (int i = 0; i < 20000; i++) { scores.Add(Gaussian(random)); isDecoy.Add(true); truePep.Add(double.NaN); }

        double[] pep = TargetDecoyPep.Estimate(scores, isDecoy).Peps;

        int[] targets = Enumerable.Range(0, scores.Count).Where(i => !isDecoy[i]).ToArray();
        double estimatedFalse = targets.Sum(i => pep[i]);
        double expectedFalse = targets.Sum(i => truePep[i]);
        Assert.That(estimatedFalse, Is.EqualTo(expectedFalse).Within(0.03 * expectedFalse), "the PEPs sum to the false targets");
        foreach (var (lo, hi) in new[] { (0.01, 0.1), (0.1, 0.3), (0.3, 0.6) })
        {
            int[] bin = targets.Where(i => truePep[i] >= lo && truePep[i] < hi).ToArray();
            Assert.That(bin.Average(i => pep[i]), Is.EqualTo(bin.Average(i => truePep[i])).Within(0.15 * hi),
                $"mean PEP of targets whose true PEP is in [{lo}, {hi})");
        }
    }

    [Test]
    public void ACandidateAboveEveryDecoyHasPepZeroAndWhereDecoysAreAsFrequentAsTargetsThePepIsOne()
    {
        // Ranked: T T D T D D. The decoy share by score is 0 0 1 0 1 1; pooling the violation gives 0 0 .5 .5 1 1.
        double[] pep = TargetDecoyPep.Estimate(new double[] { 5, 4, 3, 2, 1, 0 }, new[] { false, false, true, false, true, true }).Peps;
        Assert.That(pep, Is.EqualTo(new double[] { 0, 0, 1, 1, 1, 1 }));
    }

    [Test]
    public void AHalfDecoyCountsAsHalfADecoy()
    {
        // Decoy share 0 0 .5 0 pools to 0 0 .25 .25, and .25 / .75 = 1/3.
        double[] pep = TargetDecoyPep.Estimate(new double[] { 4, 3, 2, 1 }, new[] { 0.0, 0.0, 0.5, 0.0 }).Peps;
        Assert.That(pep[0], Is.EqualTo(0));
        Assert.That(pep[1], Is.EqualTo(0));
        Assert.That(pep[2], Is.EqualTo(1.0 / 3).Within(1e-15));
        Assert.That(pep[3], Is.EqualTo(1.0 / 3).Within(1e-15));
    }

    [Test]
    public void ThePepNeverRisesAsTheScoreRises()
    {
        var random = new Random(3);
        double[] scores = Enumerable.Range(0, 5000).Select(_ => Gaussian(random)).ToArray();
        bool[] isDecoy = scores.Select(s => random.NextDouble() < 1 / (1 + Math.Exp(3 * s))).ToArray();
        double[] pep = TargetDecoyPep.Estimate(scores, isDecoy).Peps;
        int[] order = Enumerable.Range(0, scores.Length).OrderByDescending(i => scores[i]).ToArray();
        for (int k = 1; k < order.Length; k++)
            Assert.That(pep[order[k]], Is.GreaterThanOrEqualTo(pep[order[k - 1]]), $"rank {k}");
        Assert.That(pep.All(p => p >= 0 && p <= 1), Is.True, "every PEP is in [0, 1]");
    }

    [Test]
    public void TiedScoresShareOnePepAndRowOrderChangesNothing()
    {
        var random = new Random(11);
        double[] scores = Enumerable.Range(0, 3000).Select(_ => Math.Round(2 * Gaussian(random), 1)).ToArray();
        bool[] isDecoy = scores.Select(s => random.NextDouble() < 1 / (1 + Math.Exp(2 * s))).ToArray();
        double[] pep = TargetDecoyPep.Estimate(scores, isDecoy).Peps;

        foreach (var tie in Enumerable.Range(0, scores.Length).GroupBy(i => scores[i]))
            Assert.That(tie.Select(i => pep[i]).Distinct().Count(), Is.EqualTo(1), $"score {tie.Key}");

        // Reversed, decoys come before tied targets; nothing may change.
        int[] reversed = Enumerable.Range(0, scores.Length).Reverse().ToArray();
        double[] again = TargetDecoyPep.Estimate(reversed.Select(i => scores[i]).ToArray(), reversed.Select(i => isDecoy[i]).ToArray()).Peps;
        for (int k = 0; k < reversed.Length; k++)
            Assert.That(again[k], Is.EqualTo(pep[reversed[k]]), $"row {reversed[k]}");
    }

    [Test]
    public void TheBooleanOverloadEqualsZeroOneDecoyFractions()
    {
        var random = new Random(5);
        double[] scores = Enumerable.Range(0, 1000).Select(_ => Gaussian(random)).ToArray();
        bool[] isDecoy = scores.Select(s => random.NextDouble() < 1 / (1 + Math.Exp(2 * s))).ToArray();
        Assert.That(TargetDecoyPep.Estimate(scores, isDecoy).Peps,
            Is.EqualTo(TargetDecoyPep.Estimate(scores, isDecoy.Select(d => d ? 1.0 : 0.0).ToArray()).Peps));
    }

    [Test]
    public void WithNoDecoysOrNoTargetsEveryPepIsNaNAndTheStatusSaysWhy()
    {
        var noDecoys = TargetDecoyPep.Estimate(new double[] { 3, 2, 1 }, new[] { false, false, false });
        Assert.That(noDecoys.Status, Is.EqualTo(TargetDecoyPepStatus.NoDecoys));
        Assert.That(noDecoys.Peps.All(double.IsNaN), Is.True);

        var noTargets = TargetDecoyPep.Estimate(new double[] { 3, 2 }, new[] { true, true });
        Assert.That(noTargets.Status, Is.EqualTo(TargetDecoyPepStatus.NoTargets));
        Assert.That(noTargets.Peps.All(double.IsNaN), Is.True);

        var empty = TargetDecoyPep.Estimate(Array.Empty<double>(), Array.Empty<bool>());
        Assert.That(empty.Status, Is.EqualTo(TargetDecoyPepStatus.NoDecoys));
        Assert.That(empty.Peps, Is.Empty);
    }

    [Test]
    public void BadArgumentsAreRefused()
    {
        Assert.Throws<ArgumentNullException>(() => TargetDecoyPep.Estimate(null!, new[] { false }), "scores is null");
        Assert.Throws<ArgumentNullException>(() => TargetDecoyPep.Estimate(new double[] { 1 }, (IReadOnlyList<bool>)null!), "isDecoy is null");
        Assert.Throws<ArgumentNullException>(() => TargetDecoyPep.Estimate(new double[] { 1 }, (IReadOnlyList<double>)null!), "decoyFraction is null");
        Assert.Throws<ArgumentException>(() => TargetDecoyPep.Estimate(new double[] { 1, 2 }, new[] { false }), "lengths disagree");
        Assert.Throws<ArgumentException>(() => TargetDecoyPep.Estimate(new[] { 1, double.NaN }, new[] { false, true }), "a score is NaN");
        Assert.Throws<ArgumentException>(() => TargetDecoyPep.Estimate(new[] { 1, double.PositiveInfinity }, new[] { false, true }), "a score is infinite");
        Assert.Throws<ArgumentException>(() => TargetDecoyPep.Estimate(new double[] { 1, 2 }, new[] { 0.0, 1.5 }), "a decoy fraction is above 1");
        Assert.Throws<ArgumentException>(() => TargetDecoyPep.Estimate(new double[] { 1, 2 }, new[] { 0.0, -0.1 }), "a decoy fraction is negative");
        Assert.Throws<ArgumentException>(() => TargetDecoyPep.Estimate(new double[] { 1, 2 }, new[] { 0.0, double.NaN }), "a decoy fraction is NaN");
    }
}
