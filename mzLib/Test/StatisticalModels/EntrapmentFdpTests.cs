using NUnit.Framework;
using StatisticalModels;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// Entrapment FDP estimators of Wen et al. 2025 (Nat. Methods 22:1454, eqs. 1, 2 and 4).
/// <para>
/// The sweep is checked against the entrapment project's reference script
/// (<c>results/fdp/fdp_analysis.py</c> at commit 92f3108) on the hand-built search below.
/// Every expected value was also derived by hand. The script is never run by these tests.
/// </para>
/// <para>
/// Two deliberate departures from that script follow the paper and FDRBench instead:
/// <list type="bullet">
/// <item>The paired estimate is reported only for r = 1.</item>
/// <item>It is not capped at 1.</item>
/// </list>
/// </para>
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class EntrapmentFdpTests
{
    private const double Tolerance = 1e-5;

    // Targets T1..T6 and T7 (never discovered), as (score, q).
    private static readonly ScoredIdentification T1 = new(100, 0.001);
    private static readonly ScoredIdentification T2 = new(90, 0.002);
    private static readonly ScoredIdentification T3 = new(80, 0.004);
    private static readonly ScoredIdentification T4 = new(70, 0.008);
    private static readonly ScoredIdentification T5 = new(60, 0.02);
    private static readonly ScoredIdentification T6 = new(50, 0.05);
    private static readonly ScoredIdentification T7 = new(20, 0.5);

    private static List<ScoredIdentification> Targets => [T1, T2, T3, T4, T5, T6, T7];

    // Each entrapment exercises one branch of the paired estimator:
    // E1's partner outscores it, E2's partner scores lower but is discovered,
    // E3 has no partner, and E4's partner is never discovered.
    private static List<ScoredEntrapment> Entrapments =>
    [
        new(85, 0.003, T1),
        new(75, 0.006, T4),
        new(65, 0.01, null),
        new(55, 0.03, T7),
    ];

    /// <summary>
    /// fdp_analysis.py (92f3108) output for this search at r = 1:
    /// qValueThreshold, nTarget, nEntrapment, lowerBound, combined, paired.
    /// </summary>
    private static readonly (double Q, int Nt, int Ne, double Lower, double Combined, double Paired)[] ReferenceR1 =
    [
        (0.003, 2, 1, 0.33333, 0.66667, 0.33333),
        (0.006, 3, 2, 0.40000, 0.80000, 0.60000),
        (0.010, 4, 3, 0.42857, 0.85714, 0.85714),
        (0.020, 5, 3, 0.37500, 0.75000, 0.75000),
        (0.030, 5, 4, 0.44444, 0.88889, 0.88889),
        (0.050, 6, 4, 0.40000, 0.80000, 0.80000),
        (0.100, 6, 4, 0.40000, 0.80000, 0.80000),
    ];

    [Test]
    public void SweepMatchesTheReferenceScriptAtROne()
    {
        var points = EntrapmentFdp.Sweep(Targets, Entrapments, 1.0, ReferenceR1.Select(row => row.Q).ToArray());

        Assert.That(points.Count, Is.EqualTo(ReferenceR1.Length));
        for (int i = 0; i < ReferenceR1.Length; i++)
        {
            var expected = ReferenceR1[i];
            var actual = points[i];
            Assert.That(actual.QValueThreshold, Is.EqualTo(expected.Q), $"row {i}");
            Assert.That(actual.TargetCount, Is.EqualTo(expected.Nt), $"q {expected.Q}");
            Assert.That(actual.EntrapmentCount, Is.EqualTo(expected.Ne), $"q {expected.Q}");
            Assert.That(actual.LowerBound, Is.EqualTo(expected.Lower).Within(Tolerance), $"q {expected.Q}");
            Assert.That(actual.Combined, Is.EqualTo(expected.Combined).Within(Tolerance), $"q {expected.Q}");
            Assert.That(actual.Paired, Is.EqualTo(expected.Paired).Within(Tolerance), $"q {expected.Q}");
        }
    }

    /// <summary>
    /// Combined scales with r (fdp_analysis.py at r = 2 gives 0.64286 at q = 0.01). Paired is not reported,
    /// because eq. 4 needs every target paired with exactly one entrapment.
    /// </summary>
    [Test]
    public void AtRTwoCombinedScalesAndPairedIsNotReported()
    {
        var point = EntrapmentFdp.Sweep(Targets, Entrapments, 2.0, [0.01]).Single();

        Assert.That(point.LowerBound, Is.EqualTo(3.0 / 7).Within(1e-12));
        Assert.That(point.Combined, Is.EqualTo(3 * 1.5 / 7).Within(1e-12));
        Assert.That(point.Paired, Is.Null);
    }

    [Test]
    public void ClosedFormsMatchTheEquations()
    {
        Assert.That(EntrapmentFdp.LowerBound(4, 3), Is.EqualTo(3.0 / 7).Within(1e-12));
        Assert.That(EntrapmentFdp.Combined(4, 3, 1.0), Is.EqualTo(6.0 / 7).Within(1e-12));
        Assert.That(EntrapmentFdp.Combined(4, 3, 4.0), Is.EqualTo(3 * 1.25 / 7).Within(1e-12));
        // N_E + N_{E>=s>T} + 2 N_{E>T>=s} over N_T + N_E
        Assert.That(EntrapmentFdp.Paired(4, 3, 1, 1), Is.EqualTo(6.0 / 7).Within(1e-12));
    }

    /// <summary>
    /// With no discoveries every estimate is 0/0. The threshold is left out of the sweep, as in the reference script.
    /// </summary>
    [Test]
    public void AThresholdWithNoDiscoveriesIsLeftOut()
    {
        var points = EntrapmentFdp.Sweep(Targets, Entrapments, 1.0, [0.0005, 0.01]);

        Assert.That(points.Select(point => point.QValueThreshold), Is.EqualTo(new[] { 0.01 }));
    }

    /// <summary>
    /// Neither the paper nor FDRBench caps the paired estimate. Two entrapments that each beat a discovered
    /// partner give (2 + 0 + 2*2) / (2 + 2) = 1.5.
    /// </summary>
    [Test]
    public void PairedIsNotCappedAtOne()
    {
        var a = new ScoredIdentification(10, 0.001);
        var b = new ScoredIdentification(10, 0.001);
        var point = EntrapmentFdp.Sweep([a, b], [new(20, 0.001, a), new(20, 0.001, b)], 1.0, [0.01]).Single();

        Assert.That(point.Paired, Is.EqualTo(1.5).Within(1e-12));
    }

    /// <summary>
    /// An exact tie (equal q and equal score) counts as E > T, the conservative choice for an upper bound.
    /// FDRBench's ranks are unique; its k-fold path (FDPCalcKFold.getNk, ei &gt;= ti) counts a tie the same way.
    /// N_T = 1, N_E = 1, N_{E>T>=s} = 1: (1 + 0 + 2) / 2.
    /// </summary>
    [Test]
    public void AnExactTieCountsAsTheEntrapmentRankedAhead()
    {
        var target = new ScoredIdentification(50, 0.001);
        var point = EntrapmentFdp.Sweep([target], [new(50, 0.001, target)], 1.0, [0.01]).Single();

        Assert.That(point.Paired, Is.EqualTo(3.0 / 2).Within(1e-12));
    }

    /// <summary>
    /// The pair is ordered as FDRBench ranks discoveries: by q-value, then by score. Here the entrapment has
    /// the lower q but the lower score, so it is ranked ahead and counts as N_{E>T>=s}: (1 + 0 + 2) / 2.
    /// </summary>
    [Test]
    public void AnEntrapmentWithTheLowerQRanksAheadOfAHigherScoringPartner()
    {
        var target = new ScoredIdentification(12, 0.005);
        var point = EntrapmentFdp.Sweep([target], [new(10, 0.001, target)], 1.0, [0.01]).Single();

        Assert.That(point.Paired, Is.EqualTo(3.0 / 2).Within(1e-12));
    }

    /// <summary>
    /// The mirror case: the partner has the lower q but the lower score, so it is ranked ahead and the pair
    /// adds nothing beyond the entrapment: (1 + 0 + 0) / 2.
    /// </summary>
    [Test]
    public void APartnerWithTheLowerQRanksAheadOfAHigherScoringEntrapment()
    {
        var target = new ScoredIdentification(10, 0.001);
        var point = EntrapmentFdp.Sweep([target], [new(12, 0.005, target)], 1.0, [0.01]).Single();

        Assert.That(point.Paired, Is.EqualTo(1.0 / 2).Within(1e-12));
    }

    /// <summary>
    /// With equal q-values the score decides, higher first.
    /// </summary>
    [TestCase(12, 10, 3.0 / 2)]
    [TestCase(10, 12, 1.0 / 2)]
    public void WithEqualQTheHigherScoreRanksAhead(double entrapmentScore, double partnerScore, double expectedPaired)
    {
        var target = new ScoredIdentification(partnerScore, 0.004);
        var point = EntrapmentFdp.Sweep([target], [new(entrapmentScore, 0.004, target)], 1.0, [0.01]).Single();

        Assert.That(point.Paired, Is.EqualTo(expectedPaired).Within(1e-12));
    }

    /// <summary>
    /// A partner that is identified but not discovered at the threshold counts as N_{E>=s>T}, even when its
    /// score beats the entrapment's. Discovery is decided by q, not by score.
    /// </summary>
    [Test]
    public void APartnerAboveTheThresholdCountsAsBelowTheCutoffWhateverItsScore()
    {
        var discovered = new ScoredIdentification(90, 0.001);
        var undiscovered = new ScoredIdentification(95, 0.2);
        var point = EntrapmentFdp.Sweep([discovered, undiscovered], [new(60, 0.005, undiscovered)], 1.0, [0.01]).Single();

        // N_T = 1, N_E = 1, N_{E>=s>T} = 1
        Assert.That(point.Paired, Is.EqualTo((1 + 1) / 2.0).Within(1e-12));
    }

    /// <summary>
    /// Foreign-species entrapment has no partner by construction. With every entrapment foreign, the paired
    /// estimate would ignore r and exceed combined for any r above 1, so it is not reported.
    /// </summary>
    [Test]
    public void PairedIsNotReportedWhenEveryEntrapmentIsForeign()
    {
        var point = EntrapmentFdp.Sweep([T1, T2], [new(85, 0.003, null, IsForeign: true)], 1.0, [0.01]).Single();

        Assert.That(point.Paired, Is.Null);
        Assert.That(point.Combined, Is.EqualTo(2.0 / 3).Within(1e-12));
    }

    [Test]
    public void OneNativeEntrapmentKeepsPairedReported()
    {
        var point = EntrapmentFdp.Sweep([T1, T2], [new(85, 0.003, null, IsForeign: true), new(84, 0.003, T2)], 1.0, [0.01]).Single();

        Assert.That(point.Paired, Is.Not.Null);
    }

    [Test]
    public void DefaultThresholdsRunFromOneInAThousandToOneInTen()
    {
        var thresholds = EntrapmentFdp.DefaultQValueThresholds;

        Assert.That(thresholds.Count, Is.EqualTo(100));
        Assert.That(thresholds[0], Is.EqualTo(0.001));
        Assert.That(thresholds[9], Is.EqualTo(0.01), "exactly 10/1000, so q = 0.01 is discovered at 1%");
        Assert.That(thresholds[^1], Is.EqualTo(0.1));
    }

    [TestCase("unrepairableRunCollision", EntrapmentExclusionKind.NotReallyEntrapment)]
    [TestCase("sharedWithTarget", EntrapmentExclusionKind.NotReallyEntrapment)]
    [TestCase("initiatorMethionineCollision", EntrapmentExclusionKind.NotReallyEntrapment)]
    [TestCase("ambiguous", EntrapmentExclusionKind.Unpairable)]
    public void ExclusionReasonsMapToWhatTheyMeanForTheEstimate(string reason, EntrapmentExclusionKind expected)
    {
        Assert.That(EntrapmentFdp.ClassifyExclusion(reason), Is.EqualTo(expected));
    }

    /// <summary>
    /// The two kinds remove contamination from opposite sides of the estimate, so an unknown reason cannot be guessed.
    /// </summary>
    [Test]
    public void AnUnknownExclusionReasonThrows()
    {
        var e = Assert.Throws<ArgumentException>(() => EntrapmentFdp.ClassifyExclusion("somethingNew"));
        Assert.That(e!.Message, Does.Contain("somethingNew"));
    }

    [TestCase(0.0)]
    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void RMustBeAPositiveFiniteRatio(double r)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => EntrapmentFdp.Combined(4, 3, r));
        Assert.Throws<ArgumentOutOfRangeException>(() => EntrapmentFdp.Sweep(Targets, Entrapments, r, [0.01]));
    }

    /// <summary>
    /// A bad r is refused up front, not only when some threshold has discoveries to estimate.
    /// </summary>
    [Test]
    public void ABadRatioIsRefusedEvenWhenNothingIsDiscovered()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => EntrapmentFdp.Sweep(Targets, Entrapments, 0.0, [0.0001]));
    }

    [Test]
    public void NullInputsThrow()
    {
        Assert.That(Assert.Throws<ArgumentNullException>(() => EntrapmentFdp.Sweep(null!, Entrapments, 1.0, [0.01]))!.ParamName, Is.EqualTo("targets"));
        Assert.That(Assert.Throws<ArgumentNullException>(() => EntrapmentFdp.Sweep(Targets, null!, 1.0, [0.01]))!.ParamName, Is.EqualTo("entrapments"));
        Assert.That(Assert.Throws<ArgumentNullException>(() => EntrapmentFdp.Sweep(Targets, Entrapments, 1.0, null!))!.ParamName, Is.EqualTo("qValueThresholds"));
    }
}
