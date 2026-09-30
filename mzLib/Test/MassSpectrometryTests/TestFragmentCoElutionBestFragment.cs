using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests;

/// <summary>
/// Co-elution against the best fragment. The fragment that correlates best with the others serves as the elution
/// profile, once smoothed, and every fragment is then scored against that profile. One reliable reference is harder to
/// corrupt than an average of all the others, which an interfered fragment drags along.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestFragmentCoElutionBestFragment
{
    private static double[] Gaussian(double apex, double height, int length = 21, double sigma = 2.0) =>
        Enumerable.Range(0, length).Select(s => height * Math.Exp(-0.5 * Math.Pow((s - apex) / sigma, 2))).ToArray();

    /// <summary>Interior points are 0.25/0.5/0.25; an end point renormalizes over the neighbour it has (0.5/0.25 → 2/3, 1/3).</summary>
    [Test]
    public void SmoothingWeightsAQuarterHalfQuarter()
    {
        double[] smoothed = FragmentCoElution.Smooth([0, 4, 0, 8, 0]);

        Assert.That(smoothed, Is.EqualTo(new[] { 4.0 / 3, 2.0, 3.0, 4.0, 8.0 / 3 }).Within(1e-12));
        Assert.That(FragmentCoElution.Smooth([5.0]), Is.EqualTo(new[] { 5.0 }));
        Assert.That(FragmentCoElution.Smooth([]), Is.Empty);
    }

    /// <summary>
    /// The best fragment is the one whose summed correlation with the others is largest. Two fragments elute together and
    /// a third is interference elsewhere, so one of the co-eluting pair wins. On a tie, the lower index wins.
    /// </summary>
    [Test]
    public void TheBestFragmentIsTheOneMostCorrelatedWithTheRest()
    {
        var traces = new List<double[]> { Gaussian(4, 50), Gaussian(10, 100), Gaussian(10, 80), Gaussian(10.5, 60) };

        Assert.That(FragmentCoElution.BestFragment(traces, 0, 20), Is.EqualTo(1));
        Assert.That(FragmentCoElution.BestFragment([Gaussian(10, 5), Gaussian(10, 5)], 0, 20), Is.EqualTo(0), "tie goes to the lower index");
    }

    /// <summary>Correlation with the reference, per fragment; flat, anti-correlated or silent fragments score 0.</summary>
    [Test]
    public void CorrelationsToTheReferenceAreClippedAtZero()
    {
        double[] reference = Gaussian(10, 100);
        var traces = new List<double[]>
        {
            Gaussian(10, 3),                                         // same shape: 1
            new double[21],                                          // silent: 0
            reference.Select(v => 100 - v).ToArray(),                // anti-correlated: 0
            Gaussian(12, 50),                                        // shifted: between 0 and 1
        };

        double[] r = FragmentCoElution.CorrelationsTo(traces, reference, 0, 20);

        Assert.That(r[0], Is.EqualTo(1).Within(1e-12));
        Assert.That(r[1], Is.EqualTo(0));
        Assert.That(r[2], Is.EqualTo(0));
        Assert.That(r[3], Is.GreaterThan(0.3).And.LessThan(0.9));
    }

    /// <summary>Only the range is compared: a fragment matching the reference inside the range scores 1 whatever lies outside.</summary>
    [Test]
    public void OnlyTheRangeIsCompared()
    {
        double[] reference = Gaussian(10, 100);
        double[] trace = Gaussian(10, 7);
        trace[0] = 1000;
        trace[20] = 1000;

        Assert.That(FragmentCoElution.CorrelationsTo([trace], reference, 5, 15)[0], Is.EqualTo(1).Within(1e-12));
        Assert.That(FragmentCoElution.BestFragment([trace, Gaussian(10, 3)], 5, 15), Is.EqualTo(0));
    }

    /// <summary>
    /// Candidate peaks: every local maximum of FindApex's apex score, at least a half-width apart, best first. A
    /// precursor whose true peak scores second to interference is then still in the running.
    /// </summary>
    [Test]
    public void CandidateApexesAreTheDistinctLocalMaximaBestFirst()
    {
        double[] library = [1.0, 0.6, 0.3];
        // Two peak groups in library proportions: a big one at scan 30 and a smaller one at scan 10
        var traces = library.Select(l => Enumerable.Range(0, 41)
            .Select(s => l * (1000 * Math.Exp(-0.5 * Math.Pow((s - 30) / 2.0, 2)) + 200 * Math.Exp(-0.5 * Math.Pow((s - 10) / 2.0, 2)))).ToArray()).ToList();

        int[] apexes = FragmentCoElution.FindApexes(traces, library, halfWidth: 3, maxCount: 5);

        Assert.That(apexes, Is.EqualTo(new[] { 30, 10 }));
        Assert.That(apexes[0], Is.EqualTo(FragmentCoElution.FindApex(traces, library, 3)), "the best candidate is FindApex's apex");
        Assert.That(FragmentCoElution.FindApexes(traces, library, 3, maxCount: 1), Is.EqualTo(new[] { 30 }));
        Assert.That(FragmentCoElution.FindApexes(library.Select(_ => new double[41]).ToList(), library, 3, 5), Is.Empty, "no signal, no candidates");
    }

    /// <summary>Within a half-width of a better scan there is no second candidate, whatever small bumps the score has.</summary>
    [Test]
    public void CandidatesAreAtLeastAHalfWidthApart()
    {
        double[] library = [1.0, 0.5];
        var traces = library.Select(l => Enumerable.Range(0, 41)
            .Select(s => l * (1000 * Math.Exp(-0.5 * Math.Pow((s - 20) / 3.0, 2)) + (s == 22 ? 50 : 0))).ToArray()).ToList();

        int[] apexes = FragmentCoElution.FindApexes(traces, library, halfWidth: 3, maxCount: 5);

        Assert.That(apexes.Zip(apexes.Skip(1)).All(pair => Math.Abs(pair.First - pair.Second) > 3));
        Assert.That(apexes[0], Is.EqualTo(20));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentCoElution.FindApexes(traces, library, 3, maxCount: 0));
    }
    [Test]
    public void ArgumentsAreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => FragmentCoElution.Smooth(null!));
        Assert.Throws<ArgumentNullException>(() => FragmentCoElution.BestFragment(null!, 0, 1));
        Assert.Throws<ArgumentException>(() => FragmentCoElution.BestFragment([], 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentCoElution.BestFragment([Gaussian(10, 1)], 5, 30));
        Assert.Throws<ArgumentNullException>(() => FragmentCoElution.CorrelationsTo([Gaussian(10, 1)], null!, 0, 20));
        Assert.Throws<ArgumentException>(() => FragmentCoElution.CorrelationsTo([Gaussian(10, 1)], new double[5], 0, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentCoElution.CorrelationsTo([Gaussian(10, 1)], Gaussian(10, 1), 10, 5));
    }
}
