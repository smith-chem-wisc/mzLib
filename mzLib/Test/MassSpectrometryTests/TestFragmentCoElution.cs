using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests;

/// <summary>
/// DIA library search reads one precursor's fragments as a dense intensity matrix: fragment traces over the MS2 scans of
/// its isolation window, zero where a fragment was not seen. A real precursor's fragments rise and fall together;
/// coincidental peaks at its fragment m/z do not. These tests pin what DIA needs from mzLib's chromatogram code, on
/// hand-built traces whose answer is known.
/// <para>
/// The first group validates the existing <see cref="ExtractedIonChromatogram.FindPeakBoundaries"/> for this use. It
/// passes, so <see cref="FragmentCoElution.PeakBounds"/> reuses it unchanged. The rest specify what is genuinely new.
/// </para>
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestFragmentCoElution
{
    private const int Scans = 40;

    private static double[] Gaussian(double apex, double height, double sigma = 2.0) =>
        Enumerable.Range(0, Scans).Select(s => height * Math.Exp(-0.5 * Math.Pow((s - apex) / sigma, 2))).ToArray();

    private static double[] Spike(int at, double height) =>
        Enumerable.Range(0, Scans).Select(s => s == at ? height : 0.0).ToArray();

    /// <summary>A dense trace as the existing API takes it: one peak per scan, zeros included.</summary>
    private static List<IIndexedPeak> AsPeaks(double[] trace) =>
        trace.Select((intensity, scan) => (IIndexedPeak)new IndexedMassSpectralPeak(500, intensity, scan, scan * 0.05)).ToList();

    #region Existing FindPeakBoundaries, validated for dense DIA traces

    /// <summary>One clean peak on a zero-padded trace is not cut anywhere inside the peak.</summary>
    [Test]
    public void ExistingBoundariesDoNotCutACleanPeakOnADenseTrace()
    {
        var peaks = AsPeaks(Gaussian(20, 1000));

        var boundaries = ExtractedIonChromatogram.FindPeakBoundaries(peaks, 20);

        Assert.That(boundaries.All(b => b.ZeroBasedScanIndex < 20 - 6 || b.ZeroBasedScanIndex > 20 + 6),
            "no boundary may fall within 3 sigma of the apex");
    }

    /// <summary>Two well-separated peaks on one trace are split at the valley between them.</summary>
    [Test]
    public void ExistingBoundariesSplitTwoSeparatedPeaks()
    {
        var trace = Gaussian(12, 1000).Zip(Gaussian(28, 800), (a, b) => a + b).ToArray();

        var boundaries = ExtractedIonChromatogram.FindPeakBoundaries(AsPeaks(trace), 28);

        Assert.That(boundaries.Any(b => b.ZeroBasedScanIndex is > 14 and < 26), "the valley between the peaks is a boundary");
    }

    /// <summary>
    /// A fragment missed in one scan mid-peak is a dropout, not a valley. DIA traces are dense and zero-filled, and the
    /// existing finder does not end the peak there, because it cuts only after two points below the ratio.
    /// </summary>
    [Test]
    public void ExistingBoundariesDoNotEndAPeakAtASingleScanDropout()
    {
        var trace = Gaussian(20, 1000);
        trace[18] = 0;

        var boundaries = ExtractedIonChromatogram.FindPeakBoundaries(AsPeaks(trace), 20);

        Assert.That(boundaries.All(b => b.ZeroBasedScanIndex < 17 || b.ZeroBasedScanIndex > 23));
    }

    /// <summary>The dense-trace adapter reports exactly the existing finder's peak, as inclusive scan bounds.</summary>
    [Test]
    public void PeakBoundsAreJustInsideTheExistingValleys()
    {
        var trace = Gaussian(12, 1000).Zip(Gaussian(28, 800), (a, b) => a + b).ToArray();

        var (start, end) = FragmentCoElution.PeakBounds(trace, 28);
        var valleys = ExtractedIonChromatogram.FindPeakBoundaries(AsPeaks(trace), 28).Select(b => b.ZeroBasedScanIndex).ToList();

        int expectedStart = valleys.Where(i => i < 28).DefaultIfEmpty(-1).Max() + 1;
        int expectedEnd = valleys.Where(i => i > 28).DefaultIfEmpty(Scans).Min() - 1;
        Assert.That((start, end), Is.EqualTo((expectedStart, expectedEnd)));
        Assert.That(start, Is.GreaterThan(14), "the neighbouring peak is excluded");
    }

    [Test]
    public void PeakBoundsSpanTheWholeTraceWhenNothingCutsIt()
    {
        Assert.That(FragmentCoElution.PeakBounds(Gaussian(20, 1000), 20), Is.EqualTo((0, Scans - 1)));
    }

    [Test]
    public void PeakBoundsRejectAnApexOutsideTheTrace()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentCoElution.PeakBounds(new double[5], 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentCoElution.PeakBounds(new double[5], -1));
    }

    #endregion

    #region FragmentCoElution (new)

    [Test]
    public void FragmentsElutingTogetherScoreNearOne()
    {
        var traces = new[] { 1.0, 0.8, 0.6, 0.5, 0.3, 0.2 }.Select(h => Gaussian(20, 1000 * h)).ToArray();

        Assert.That(FragmentCoElution.Score(traces, 14, 26), Is.GreaterThan(0.99));
    }

    [Test]
    public void FragmentsPeakingAtUnrelatedTimesScoreNearZero()
    {
        var traces = new[] { 5, 12, 19, 26, 33, 38 }.Select(apex => Gaussian(apex, 1000, sigma: 1.0)).ToArray();

        Assert.That(FragmentCoElution.Score(traces, 0, Scans - 1), Is.LessThan(0.2));
    }

    /// <summary>One interfered fragment lowers the score but does not sink a real peak group.</summary>
    [Test]
    public void OneInterferedFragmentLowersButDoesNotSinkTheGroup()
    {
        var clean = new[] { 1.0, 0.8, 0.6, 0.5, 0.3, 0.2 }.Select(h => Gaussian(20, 1000 * h)).ToArray();
        var interfered = clean.Select(t => t.ToArray()).ToArray();
        interfered[2] = Gaussian(20, 600).Zip(Gaussian(24, 3000, sigma: 1.0), (a, b) => a + b).ToArray();

        double cleanScore = FragmentCoElution.Score(clean, 14, 26);
        double interferedScore = FragmentCoElution.Score(interfered, 14, 26);
        Assert.That(interferedScore, Is.LessThan(cleanScore));
        Assert.That(interferedScore, Is.GreaterThan(0.5));
    }

    /// <summary>Silence is not co-elution, and is never NaN.</summary>
    [Test]
    public void AllZeroTracesScoreZero()
    {
        var traces = Enumerable.Range(0, 6).Select(_ => new double[Scans]).ToArray();

        Assert.That(FragmentCoElution.Score(traces, 0, Scans - 1), Is.EqualTo(0));
    }

    /// <summary>A fragment seen nowhere counts against the group.</summary>
    [Test]
    public void AMissingFragmentCountsAgainstTheGroup()
    {
        var all = new[] { 1.0, 0.8, 0.6, 0.5, 0.3, 0.2 }.Select(h => Gaussian(20, 1000 * h)).ToArray();
        var oneMissing = all.Select(t => t.ToArray()).ToArray();
        oneMissing[5] = new double[Scans];

        Assert.That(FragmentCoElution.Score(oneMissing, 14, 26), Is.LessThan(FragmentCoElution.Score(all, 14, 26)));
    }

    [Test]
    public void ScoreRejectsABadRange()
    {
        var traces = new[] { Gaussian(20, 1), Gaussian(20, 1) };

        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentCoElution.Score(traces, 30, 10));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentCoElution.Score(traces, 0, Scans));
        Assert.Throws<ArgumentException>(() => FragmentCoElution.Score([new double[3], new double[4]], 0, 2));
    }

    /// <summary>
    /// The tallest summed-intensity scan is the wrong apex when one fragment carries an interfering spike. The apex is where
    /// the fragments co-elute in the library's proportions.
    /// </summary>
    [Test]
    public void TheApexIsTheCoElutingGroupNotTheTallestSpike()
    {
        double[] library = [1.0, 0.8, 0.6, 0.5, 0.3, 0.2];
        var traces = library.Select(h => Gaussian(25, 1000 * h)).ToArray();
        traces[0] = traces[0].Zip(Spike(8, 50000), (a, b) => a + b).ToArray();

        int apex = FragmentCoElution.FindApex(traces, library, halfWidth: 4);

        Assert.That(apex, Is.InRange(24, 26));
    }

    [Test]
    public void ThereIsNoApexInSilence()
    {
        var traces = Enumerable.Range(0, 3).Select(_ => new double[Scans]).ToArray();

        Assert.That(FragmentCoElution.FindApex(traces, [1.0, 0.5, 0.2], halfWidth: 3), Is.EqualTo(-1));
    }

    [Test]
    public void TopIndicesAreTheLargestInOriginalOrder()
    {
        double[] intensities = [0.1, 0.9, 0.3, 1.0, 0.05, 0.9];

        Assert.That(FragmentCoElution.TopIndices(intensities, 3), Is.EqualTo(new[] { 1, 3, 5 }));
        Assert.That(FragmentCoElution.TopIndices(intensities, 10), Is.EqualTo(new[] { 0, 1, 2, 3, 4, 5 }));
    }

    #endregion
}
