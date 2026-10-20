using NUnit.Framework;
using StatisticalModels;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// Target-decoy q-values and the picked competition (Savitski et al. 2015, Mol. Cell. Proteomics 14:2394). In the picked
/// competition, a target and its own decoy (sharing a pair key, for example a protein accession with DECOY_ stripped)
/// compete. Only the better of each pair is kept, and q-values are (D+1)/T over the kept entries.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PickedTargetDecoyTests
{
    #region TargetDecoyQValues

    /// <summary>
    /// (D+1)/T down the ranking, capped at 1, made monotone and reported in input order. Ranked: T×8, D, T, D. The eighth
    /// target's 1/8 carries up to the top; the decoy and the ninth target share min(2/8, 2/9) = 2/9; the last decoy has 3/9.
    /// </summary>
    [Test]
    public void QValuesAreDecoysPlusOneOverTargetsMonotoneAndInInputOrder()
    {
        double[] ranked = [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1];
        bool[] decoy = [false, false, false, false, false, false, false, false, true, false, true];
        double[] expected = [.125, .125, .125, .125, .125, .125, .125, .125, 2.0 / 9, 2.0 / 9, 3.0 / 9];
        var order = new[] { 5, 0, 10, 3, 8, 1, 9, 2, 7, 4, 6 };

        double[] q = TargetDecoyQValues.Compute(order.Select(i => ranked[i]).ToArray(), order.Select(i => decoy[i]).ToArray());

        for (int k = 0; k < order.Length; k++)
            Assert.That(q[k], Is.EqualTo(expected[order[k]]).Within(1e-12), $"rank {order[k]}");
    }

    [Test]
    public void QValuesAreCappedAtOne()
    {
        Assert.That(TargetDecoyQValues.Compute([4, 3, 2, 1], [true, false, true, true]), Is.EqualTo(new[] { 1.0, 1.0, 1.0, 1.0 }));
    }

    /// <summary>
    /// A decoy tied with a target ranks first, the conservative choice. Ten targets at 10, then a target and a decoy tied
    /// at 5. With the decoy first, the tied target's q-value is 2/11. With the target first it would be 1/11.
    /// </summary>
    [Test]
    public void ATiedDecoyRanksBeforeTheTarget()
    {
        double[] scores = [.. Enumerable.Repeat(10.0, 10), 5, 5];
        bool[] decoy = [.. Enumerable.Repeat(false, 10), false, true];

        double[] q = TargetDecoyQValues.Compute(scores, decoy);

        Assert.That(q[10], Is.EqualTo(2.0 / 11).Within(1e-12));
        Assert.That(q[11], Is.EqualTo(2.0 / 11).Within(1e-12));
    }

    [Test]
    public void QValueArgumentsAreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => TargetDecoyQValues.Compute(null!, [true]));
        Assert.Throws<ArgumentNullException>(() => TargetDecoyQValues.Compute([1.0], null!));
        Assert.Throws<ArgumentException>(() => TargetDecoyQValues.Compute([1.0, 2.0], [true]));
        Assert.Throws<ArgumentException>(() => TargetDecoyQValues.Compute([double.NaN], [true]));
        Assert.That(TargetDecoyQValues.Compute([], []), Is.Empty);
    }

    #endregion

    #region PickedTargetDecoy

    /// <summary>
    /// The picked competition, worked by hand.
    /// <list type="bullet">
    /// <item>Pairs A–I: the target (20 down to 12) beats its decoy (1–9).</item>
    /// <item>Pair J: the decoy (11) beats its target (3).</item>
    /// <item>K: a target (10) with no decoy competes alone.</item>
    /// </list>
    /// The kept ranking is T×9, D, T. The ninth target's 1/9 carries up to the top, and the decoy and the last target share
    /// 2/10. The losers are not kept and have no q-value.
    /// </summary>
    [Test]
    public void TheBetterOfEachPairCompetesAndQValuesCountOnlyTheWinners()
    {
        var keys = "ABCDEFGHI".Select(c => c.ToString()).ToList();
        var pairKeys = keys.Concat(keys).Concat(["J", "J", "K"]).ToArray();
        double[] scores = [20, 19, 18, 17, 16, 15, 14, 13, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 3, 11, 10];
        bool[] isDecoy = [.. Enumerable.Repeat(false, 9), .. Enumerable.Repeat(true, 9), false, true, false];

        var result = PickedTargetDecoy.Compete(pairKeys, scores, isDecoy);

        for (int i = 0; i < 9; i++)
        {
            Assert.That(result.Kept[i], Is.True, $"target {pairKeys[i]}");
            Assert.That(result.QValues[i], Is.EqualTo(1.0 / 9).Within(1e-12), $"target {pairKeys[i]}");
            Assert.That(result.Kept[9 + i], Is.False, $"decoy {pairKeys[i]}");
            Assert.That(result.QValues[9 + i], Is.NaN);
        }
        Assert.That(result.Kept[18], Is.False, "J's target lost to its decoy");
        Assert.That(result.QValues[18], Is.NaN);
        Assert.That(result.Kept[19], Is.True);
        Assert.That(result.QValues[19], Is.EqualTo(0.2).Within(1e-12));
        Assert.That(result.Kept[20], Is.True, "an unpaired target competes alone");
        Assert.That(result.QValues[20], Is.EqualTo(0.2).Within(1e-12));
    }

    /// <summary>A target tied with its own decoy loses: the conservative choice.</summary>
    [Test]
    public void ATargetTiedWithItsDecoyLoses()
    {
        var result = PickedTargetDecoy.Compete(["X", "X"], [5, 5], [false, true]);

        Assert.That(result.Kept, Is.EqualTo(new[] { false, true }));
    }

    /// <summary>Only the single best entry under a key competes, however many targets or decoys share it.</summary>
    [Test]
    public void OnlyTheBestEntryUnderAKeyCompetes()
    {
        var result = PickedTargetDecoy.Compete(["X", "X", "X"], [4, 9, 6], [false, false, true]);

        Assert.That(result.Kept, Is.EqualTo(new[] { false, true, false }));
    }

    [Test]
    public void TheOutcomeDoesNotDependOnInputOrder()
    {
        var random = new Random(7);
        int n = 400;
        string[] keys = Enumerable.Range(0, n).Select(i => $"P{i / 2}").ToArray();
        double[] scores = Enumerable.Range(0, n).Select(_ => Math.Round(random.NextDouble() * 50, 1)).ToArray();
        bool[] decoy = Enumerable.Range(0, n).Select(i => i % 2 == 1).ToArray();
        var order = Enumerable.Range(0, n).OrderBy(_ => random.Next()).ToArray();

        var a = PickedTargetDecoy.Compete(keys, scores, decoy);
        var b = PickedTargetDecoy.Compete(order.Select(i => keys[i]).ToArray(), order.Select(i => scores[i]).ToArray(), order.Select(i => decoy[i]).ToArray());

        for (int k = 0; k < n; k++)
        {
            Assert.That(b.Kept[k], Is.EqualTo(a.Kept[order[k]]), $"row {order[k]}");
            Assert.That(b.QValues[k], Is.EqualTo(a.QValues[order[k]]).Within(1e-12), $"row {order[k]}");
        }
    }

    [Test]
    public void PickedArgumentsAreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => PickedTargetDecoy.Compete(null!, [1.0], [true]));
        Assert.Throws<ArgumentException>(() => PickedTargetDecoy.Compete(["A", "B"], [1.0], [true]));
        Assert.Throws<ArgumentException>(() => PickedTargetDecoy.Compete(["A"], [1.0], [true, false]));
        Assert.Throws<ArgumentException>(() => PickedTargetDecoy.Compete([null!], [1.0], [true]));
        Assert.Throws<ArgumentException>(() => PickedTargetDecoy.Compete(["A"], [double.PositiveInfinity], [true]));
        Assert.Throws<ArgumentNullException>(() => PickedTargetDecoy.Compete(["A"], null!, [true]));
        Assert.Throws<ArgumentNullException>(() => PickedTargetDecoy.Compete(["A"], [1.0], null!));
        Assert.Throws<ArgumentException>(() => PickedTargetDecoy.Compete(["A", "B"], [1.0], [true, false]), "only the scores are short");
        Assert.Throws<ArgumentException>(() => PickedTargetDecoy.Compete(["A", null!], [1.0, 2.0], [true, false]), "one null key among several");
        Assert.Throws<ArgumentException>(() => PickedTargetDecoy.Compete(["A", "B"], [1.0, double.NaN], [true, false]), "one non-finite score among several");
    }

    /// <summary>Two targets tied under one key: the earlier input competes, so the choice is reproducible.</summary>
    [Test]
    public void TiedEntriesOfOneClassUnderAKeyKeepTheEarlierInput()
    {
        var result = PickedTargetDecoy.Compete(["X", "X", "X"], [5, 5, 1], [false, false, true]);

        Assert.That(result.Kept, Is.EqualTo(new[] { true, false, false }));
    }

    #endregion

    /// <summary>
    /// The q-value sort was rewritten for speed (the rescorer computes it for every feature and sign of every fold); on
    /// heavily tied scores it must reproduce the original stable LINQ order exactly.
    /// </summary>
    [Test]
    public void QValuesMatchTheOriginalOrderingOnTiedScores()
    {
        var random = new Random(9);
        for (int trial = 0; trial < 20; trial++)
        {
            int n = 500 + trial * 37;
            double[] scores = Enumerable.Range(0, n).Select(_ => (double)random.Next(0, 40)).ToArray();
            bool[] decoy = Enumerable.Range(0, n).Select(_ => random.NextDouble() < 0.4).ToArray();
            int[] order = Enumerable.Range(0, n).OrderByDescending(i => scores[i]).ThenByDescending(i => decoy[i]).ThenBy(i => i).ToArray();
            var ranked = new double[n];
            int d = 0, t = 0;
            for (int k = 0; k < n; k++)
            {
                if (decoy[order[k]]) d++; else t++;
                ranked[k] = t == 0 ? 1 : Math.Min(1, (d + 1.0) / t);
            }
            for (int k = n - 2; k >= 0; k--)
                ranked[k] = Math.Min(ranked[k], ranked[k + 1]);
            var expected = new double[n];
            for (int k = 0; k < n; k++)
                expected[order[k]] = ranked[k];

            Assert.That(TargetDecoyQValues.Compute(scores, decoy), Is.EqualTo(expected), $"trial {trial}");
        }
    }
}
