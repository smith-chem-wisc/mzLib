using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests;

/// <summary>
/// A DIA precursor's quantity from its fragment traces: each fragment's area within the peak bounds, DIA-NN's interference
/// correction against the best fragment (Demichev et al. 2020, Nat Methods 17:41), one set of fragments chosen across runs,
/// and the quantity with its kind. A quantity is never 0 for "nothing": no value is NaN, and the kind (QuantProject's
/// DEF-DIA-KIND) says which nothing it is. The traces are hand-built, so every answer is known.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestFragmentQuantity
{
    private const int Scans = 40;

    /// <summary>Peak bounds four sigma either side of the apex at scan 20.</summary>
    private const int Start = 12, End = 28;

    /// <summary>Six fragments at library-like relative heights.</summary>
    private static readonly double[] Heights = [1000, 800, 600, 500, 400, 300];

    private static double[] Gaussian(double apex, double height, double sigma = 2.0) =>
        Enumerable.Range(0, Scans).Select(s => height * Math.Exp(-0.5 * Math.Pow((s - apex) / sigma, 2))).ToArray();

    private static double[] Plus(double[] a, double[] b) => a.Zip(b, (x, y) => x + y).ToArray();

    /// <summary>Six co-eluting fragments at <see cref="Heights"/> times <paramref name="scale"/>, plus seeded uniform noise.</summary>
    private static double[][] Precursor(double scale, Random? noise = null, double noiseHeight = 0) =>
        Heights.Select(h => Gaussian(20, h * scale).Select(v => v + (noise is null ? 0 : noise.NextDouble() * noiseHeight)).ToArray())
            .ToArray();

    private static double SumOfAreas(IEnumerable<double[]> traces) => traces.Sum(t => FragmentQuantity.Area(t, Start, End));

    #region Area

    /// <summary>The area is the trapezoid over the bounds less the baseline under the lower end point, as DIA search integrates.</summary>
    [Test]
    public void AreaIsTheTrapezoidLessTheBaselineUnderTheLowerEndPoint()
    {
        double[] trace = [0, 2, 4, 2, 1];

        Assert.That(FragmentQuantity.Area(trace, 0, 4), Is.EqualTo(8.5).Within(1e-12), "1 + 3 + 3 + 1.5, baseline 0");
        Assert.That(FragmentQuantity.Area(trace, 1, 4), Is.EqualTo(4.5).Within(1e-12), "3 + 3 + 1.5, less 1 x 3 scans");
    }

    /// <summary>A planted Gaussian peak, integrated within the bounds <see cref="FragmentCoElution.PeakBounds"/> finds, recovers its area within 2%.</summary>
    [TestCase(1000, 2.0)]
    [TestCase(250, 1.5)]
    [TestCase(50_000, 3.0)]
    public void APlantedGaussianPeakIsRecoveredWithinTwoPercent(double height, double sigma)
    {
        var trace = Gaussian(20, height, sigma);
        var (start, end) = FragmentCoElution.PeakBounds(trace, 20);

        Assert.That(FragmentQuantity.Area(trace, start, end), Is.EqualTo(height * sigma * Math.Sqrt(2 * Math.PI)).Within(2).Percent);
    }

    /// <summary>
    /// An area that is no measurement is NaN, never 0: bounds of one scan (a one-scan peak whose bounds sit just inside its
    /// valleys), and a trace below its own baseline. A trace silent across valid bounds is a real 0.
    /// </summary>
    [Test]
    public void AnAreaThatIsNoMeasurementIsNaNAndASilentTraceIsARealZero()
    {
        Assert.That(FragmentQuantity.Area(Gaussian(20, 1000), 20, 20), Is.NaN, "one scan");
        Assert.That(FragmentQuantity.Area([5, 0, 0, 0, 5], 0, 4), Is.NaN, "below its baseline");
        Assert.That(FragmentQuantity.Area(new double[Scans], Start, End), Is.EqualTo(0), "silent");
    }

    [Test]
    public void AreaRejectsBoundsOutsideTheTrace()
    {
        var trace = Gaussian(20, 1000);
        Assert.Throws<ArgumentNullException>(() => FragmentQuantity.Area(null!, 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Area(trace, -1, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Area(trace, 5, Scans));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Area(trace, 6, 5));
    }

    #endregion

    #region Interference correction

    /// <summary>Fragments that co-elute cleanly are not touched: each corrected area is its plain area.</summary>
    [Test]
    public void CleanCoElutingFragmentsAreNotCorrected()
    {
        var traces = Precursor(1);

        var corrected = FragmentQuantity.InterferenceCorrectedAreas(traces, Start, End);

        for (int f = 0; f < traces.Length; f++)
            Assert.That(corrected[f], Is.EqualTo(FragmentQuantity.Area(traces[f], Start, End)).Within(1e-9), $"fragment {f}");
    }

    /// <summary>
    /// An interferer beside the peak on one fragment (another ion at its m/z, eluting 2.5 sigma later and twice its
    /// height) puts the plain sum more than 30% high. Its correlation with the reference falls below 0.8, so its ratio to
    /// the best fragment is read near the apex, where the interferer is faint, and the rest is capped at 1.5 times that
    /// ratio times the reference: the corrected sum is within 10%.
    /// </summary>
    [Test]
    public void AnInterfererBesideThePeakIsCappedToTheFragmentsOwnProfile()
    {
        double truth = SumOfAreas(Precursor(1));
        var traces = Precursor(1);
        traces[2] = Plus(traces[2], Gaussian(25, 2000));
        var reference = FragmentCoElution.Smooth(traces[FragmentCoElution.BestFragment(traces, Start, End)]);
        Assert.That(FragmentCoElution.CorrelationsTo(traces, reference, Start, End)[2], Is.LessThan(FragmentQuantity.InterferenceCorrelation));

        double plain = SumOfAreas(traces);
        double corrected = FragmentQuantity.InterferenceCorrectedAreas(traces, Start, End).Sum();

        Assert.That(plain / truth - 1, Is.GreaterThan(0.30), "the interferer is large");
        Assert.That(corrected / truth - 1, Is.EqualTo(0).Within(0.10));
    }

    /// <summary>
    /// A one-scan spike on a fragment that still co-elutes (correlation at least 0.8) keeps the fragment's own ratio, and the
    /// spike is clipped to 1.5 times that ratio times the reference: more than half its excess goes, and the other
    /// fragments are untouched.
    /// </summary>
    [Test]
    public void ASpikeOnACoElutingFragmentIsClipped()
    {
        var clean = Precursor(1);
        var traces = Precursor(1);
        traces[2][23] += 400;
        var reference = FragmentCoElution.Smooth(traces[FragmentCoElution.BestFragment(traces, Start, End)]);
        Assert.That(FragmentCoElution.CorrelationsTo(traces, reference, Start, End)[2], Is.GreaterThanOrEqualTo(FragmentQuantity.InterferenceCorrelation));

        var corrected = FragmentQuantity.InterferenceCorrectedAreas(traces, Start, End);

        double plain = FragmentQuantity.Area(traces[2], Start, End), truth = FragmentQuantity.Area(clean[2], Start, End);
        Assert.That(corrected[2], Is.LessThan(plain));
        Assert.That(corrected[2] - truth, Is.LessThan((plain - truth) / 2));
        foreach (int f in new[] { 0, 1, 3, 4, 5 })
            Assert.That(corrected[f], Is.EqualTo(FragmentQuantity.Area(traces[f], Start, End)).Within(1e-9), $"fragment {f}");
    }

    /// <summary>
    /// A scan beside the apex that every fragment missed (the best fragment is 0 there) gives no ratio to read, so it is
    /// skipped, and the interfered fragment is still capped from the scans that have one.
    /// </summary>
    [Test]
    public void AScanEveryFragmentMissedIsSkippedWhenReadingTheRatio()
    {
        static double[][] Missing21(double[][] traces)
        {
            foreach (var trace in traces)
                trace[21] = 0;
            return traces;
        }
        double truth = SumOfAreas(Missing21(Precursor(1)));
        var traces = Precursor(1);
        traces[2] = Plus(traces[2], Gaussian(25, 2000));
        Missing21(traces);

        double corrected = FragmentQuantity.InterferenceCorrectedAreas(traces, Start, End).Sum();

        Assert.That(SumOfAreas(traces) / truth - 1, Is.GreaterThan(0.30), "the interferer is large");
        Assert.That(corrected / truth - 1, Is.EqualTo(0).Within(0.10));
    }

    /// <summary>Silence is a real 0 on every fragment; bounds of one scan are NaN on every fragment.</summary>
    [Test]
    public void ACorrectedAreaIsZeroForSilenceAndNaNWithoutValidBounds()
    {
        var silent = Enumerable.Range(0, 6).Select(_ => new double[Scans]).ToArray();

        Assert.That(FragmentQuantity.InterferenceCorrectedAreas(silent, Start, End), Is.All.EqualTo(0));
        Assert.That(FragmentQuantity.InterferenceCorrectedAreas(Precursor(1), 20, 20), Is.All.NaN);
    }

    [Test]
    public void InterferenceCorrectionRejectsWhatItCannotRead()
    {
        Assert.Throws<ArgumentNullException>(() => FragmentQuantity.InterferenceCorrectedAreas(null!, Start, End));
        Assert.Throws<ArgumentException>(() => FragmentQuantity.InterferenceCorrectedAreas([], Start, End));
        Assert.Throws<ArgumentException>(() => FragmentQuantity.InterferenceCorrectedAreas([new double[Scans], new double[Scans - 1]], Start, End));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.InterferenceCorrectedAreas(Precursor(1), End, Start));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.InterferenceCorrectedAreas(Precursor(1), Start, End, apexHalfWidth: -1));
    }

    #endregion

    #region Fragments chosen across runs

    /// <summary>
    /// One set of fragments serves every run: those with the highest mean correlation over the runs given. Each run's own
    /// best three differ here, the set does not, and the runs' order does not change it.
    /// </summary>
    [Test]
    public void FragmentsAreChosenOnceByTheirMeanCorrelationAcrossRuns()
    {
        float[][] runs =
        [
            [0.95f, 0.90f, 0.20f, 0.85f, 0.80f, 0.60f],
            [0.90f, 0.40f, 0.95f, 0.88f, 0.70f, 0.86f],
            [0.92f, 0.85f, 0.30f, 0.50f, 0.90f, 0.88f],
        ];

        Assert.That(FragmentQuantity.SelectAcrossRuns(runs, 3), Is.EqualTo(new[] { 0, 4, 5 }));
        Assert.That(FragmentQuantity.SelectAcrossRuns(Enumerable.Reverse(runs).ToArray(), 3), Is.EqualTo(new[] { 0, 4, 5 }));
        Assert.That(FragmentQuantity.SelectAcrossRuns([runs[1], runs[2], runs[0]], 3), Is.EqualTo(new[] { 0, 4, 5 }));
    }

    /// <summary>A fragment not measured in a run (NaN) counts 0 there; it is not averaged over the runs that saw it.</summary>
    [Test]
    public void ANaNCorrelationCountsAsZero()
    {
        float[][] runs = [[0.9f, float.NaN, 0.5f], [0.9f, 0.99f, 0.5f]];

        Assert.That(FragmentQuantity.SelectAcrossRuns(runs, 2), Is.EqualTo(new[] { 0, 2 }));
    }

    /// <summary>Ties go to the lower index; the set is in index order; asking for more than there are gives all.</summary>
    [Test]
    public void TiesGoToTheLowerIndexAndTheSetIsInIndexOrder()
    {
        float[][] runs = [[0.5f, 0.9f, 0.5f, 0.9f]];

        Assert.That(FragmentQuantity.SelectAcrossRuns(runs, 3), Is.EqualTo(new[] { 0, 1, 3 }));
        Assert.That(FragmentQuantity.SelectAcrossRuns(runs, 10), Is.EqualTo(new[] { 0, 1, 2, 3 }));
    }

    [Test]
    public void SelectionRejectsWhatItCannotRead()
    {
        Assert.Throws<ArgumentNullException>(() => FragmentQuantity.SelectAcrossRuns(null!, 3));
        Assert.Throws<ArgumentException>(() => FragmentQuantity.SelectAcrossRuns([], 3), "no run");
        Assert.Throws<ArgumentException>(() => FragmentQuantity.SelectAcrossRuns([[0.5f, 0.5f], [0.5f]], 1), "ragged");
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.SelectAcrossRuns([[0.5f]], -1));
    }

    #endregion

    #region The quantity and its kind

    /// <summary>A value is the sum of the selected fragments' areas, from at least three of them, with the floor recorded.</summary>
    [Test]
    public void AValueIsTheSumOfTheSelectedFragmentsAreas()
    {
        var quantity = FragmentQuantity.Quantify([100f, 200f, 300f, 400f, 500f], [0, 2, 4], floor: 50);

        Assert.That(quantity, Is.EqualTo(new PrecursorQuantity(900, PrecursorQuantityKind.Value, 50)));
    }

    /// <summary>Fewer than three usable (measured) fragments is "not looked": no value, never 0.</summary>
    [Test]
    public void FewerThanThreeUsableFragmentsIsNotLooked()
    {
        var twoOfFour = FragmentQuantity.Quantify([100f, float.NaN, 300f, float.NaN], [0, 1, 2, 3], floor: 10);
        var twoSelected = FragmentQuantity.Quantify([100f, 200f, 300f], [0, 2], floor: 10);

        foreach (var quantity in new[] { twoOfFour, twoSelected })
        {
            Assert.That(quantity.Kind, Is.EqualTo(PrecursorQuantityKind.NotLooked));
            Assert.That(quantity.Value, Is.NaN);
            Assert.That(quantity.Floor, Is.EqualTo(10));
        }
    }

    /// <summary>
    /// Every usable fragment at or below the floor is "below the floor": no value. Fragments seen nowhere in valid bounds are
    /// a real zero, which is below any floor, so it is never a quantity of 0.
    /// </summary>
    [Test]
    public void EveryFragmentAtOrBelowTheFloorIsBelowTheFloor()
    {
        var atFloor = FragmentQuantity.Quantify([10f, 20f, 30f], [0, 1, 2], floor: 30);
        var silent = FragmentQuantity.Quantify([0f, 0f, 0f, 0f], [0, 1, 2, 3], floor: 0);

        Assert.That(atFloor.Kind, Is.EqualTo(PrecursorQuantityKind.BelowFloor));
        Assert.That(atFloor.Value, Is.NaN);
        Assert.That(silent.Kind, Is.EqualTo(PrecursorQuantityKind.BelowFloor));
        Assert.That(silent.Value, Is.NaN);
    }

    /// <summary>A selected fragment that was not measured is left out of the sum, never imputed.</summary>
    [Test]
    public void AnUnmeasuredFragmentIsLeftOutNotImputed()
    {
        var quantity = FragmentQuantity.Quantify([100f, float.NaN, 300f, 400f], [0, 1, 2, 3], floor: 0);

        Assert.That(quantity, Is.EqualTo(new PrecursorQuantity(800, PrecursorQuantityKind.Value, 0)));
    }

    /// <summary>Over seeded random areas, some unmeasured and some silent: a value is finite and above the floor, and anything else is NaN.</summary>
    [Test]
    public void AQuantityIsNeverZeroForNothing()
    {
        var random = new Random(1606);
        for (int i = 0; i < 2000; i++)
        {
            var areas = Enumerable.Range(0, 6).Select(_ => random.Next(4) switch
            {
                0 => float.NaN,
                1 => 0f,
                _ => (float)(random.NextDouble() * 100),
            }).ToArray();
            double floor = random.NextDouble() * 50;

            var quantity = FragmentQuantity.Quantify(areas, [0, 1, 2, 3, 4, 5], floor);

            if (quantity.Kind == PrecursorQuantityKind.Value)
                Assert.That(quantity.Value, Is.GreaterThan(floor).And.LessThan(double.PositiveInfinity), $"case {i}");
            else
                Assert.That(quantity.Value, Is.NaN, $"case {i}");
        }
    }

    /// <summary>One trace's area (an MS1 quantity) has the same three kinds.</summary>
    [Test]
    public void OneTracesAreaHasTheSameKinds()
    {
        Assert.That(FragmentQuantity.Quantify(250, floor: 100), Is.EqualTo(new PrecursorQuantity(250, PrecursorQuantityKind.Value, 100)));
        Assert.That(FragmentQuantity.Quantify(100, floor: 100).Kind, Is.EqualTo(PrecursorQuantityKind.BelowFloor));
        Assert.That(FragmentQuantity.Quantify(0, floor: 0).Kind, Is.EqualTo(PrecursorQuantityKind.BelowFloor));
        Assert.That(FragmentQuantity.Quantify(double.NaN, floor: 100).Kind, Is.EqualTo(PrecursorQuantityKind.NotLooked));
        Assert.That(FragmentQuantity.Quantify(double.NaN, floor: 100).Value, Is.NaN);
    }

    [Test]
    public void QuantifyRejectsWhatItCannotRead()
    {
        Assert.Throws<ArgumentNullException>(() => FragmentQuantity.Quantify(null!, [0], 0));
        Assert.Throws<ArgumentNullException>(() => FragmentQuantity.Quantify([1f], null!, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Quantify([1f, 2f, 3f], [0, 1, 3], 0), "no fragment 3");
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Quantify([1f, 2f, 3f], [-1, 1, 2], 0), "no fragment -1");
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Quantify([1f, 2f, 3f], [0, 1, 2], -1), "negative floor");
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Quantify([1f, 2f, 3f], [0, 1, 2], double.NaN), "NaN floor");
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Quantify(5, double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Quantify([1f, -2f, 3f], [0, 1, 2], 0), "an area is never negative");
        Assert.Throws<ArgumentOutOfRangeException>(() => FragmentQuantity.Quantify(-5, 0), "an area is never negative");
    }

    #endregion

    #region From traces to ratios

    /// <summary>
    /// PLAN M9: planted 1:2 and 1:4 ratios are recovered within 10%. Each run has seeded noise and the same interferer on
    /// one fragment, which compresses the plain sum's ratio by more than 10%. Corrected areas, with three fragments chosen
    /// across both runs, recover the ratio.
    /// </summary>
    [TestCase(2)]
    [TestCase(4)]
    public void PlantedRatiosAreRecoveredWithinTenPercent(double ratio)
    {
        var random = new Random(7);
        var runs = new[] { 1.0, ratio }.Select(scale =>
        {
            var traces = Precursor(scale, random, noiseHeight: 10);
            traces[2] = Plus(traces[2], Gaussian(25, 1500));
            return traces;
        }).ToArray();

        var correlations = runs.Select(traces =>
        {
            var reference = FragmentCoElution.Smooth(traces[FragmentCoElution.BestFragment(traces, Start, End)]);
            return FragmentCoElution.CorrelationsTo(traces, reference, Start, End).Select(r => (float)r).ToArray();
        }).ToArray();
        int[] selected = FragmentQuantity.SelectAcrossRuns(correlations, 3);
        var quantities = runs.Select(traces => FragmentQuantity.Quantify(
            FragmentQuantity.InterferenceCorrectedAreas(traces, Start, End).Select(a => (float)a).ToArray(), selected, floor: 0)).ToArray();

        Assert.That(SumOfAreas(runs[1]) / SumOfAreas(runs[0]), Is.LessThan(0.9 * ratio), "the interferer compresses the plain sum");
        Assert.That(quantities.Select(q => q.Kind), Is.All.EqualTo(PrecursorQuantityKind.Value));
        Assert.That(quantities[1].Value / quantities[0].Value, Is.EqualTo(ratio).Within(10).Percent);
    }

    #endregion
}
