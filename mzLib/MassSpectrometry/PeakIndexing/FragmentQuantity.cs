#nullable enable
using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    /// <summary>
    /// A DIA precursor's quantity from its fragment traces, the dense matrix <see cref="FragmentCoElution"/> scores: each
    /// fragment's area within the peak bounds, DIA-NN's interference correction (Demichev et al. 2020, Nat Methods 17:41),
    /// one set of fragments chosen across runs, and the quantity with its <see cref="PrecursorQuantityKind"/>. These are
    /// functions of traces a search has already extracted, so quantifying needs no second pass over the spectra. A quantity
    /// is never 0 for "nothing": no value is NaN, and its kind says which nothing it is.
    /// </summary>
    public static class FragmentQuantity
    {
        /// <summary>
        /// The fewest measured fragments a quantity is read from, as every published tool that selects fragments requires;
        /// with fewer, the quantity is <see cref="PrecursorQuantityKind.NotLooked"/>.
        /// </summary>
        public const int MinimumFragments = 3;

        /// <summary>Below this correlation with the reference profile, DIA-NN takes a fragment to be interfered.</summary>
        public const double InterferenceCorrelation = 0.8;

        /// <summary>DIA-NN lowers every point of a fragment's trace above this multiple of its ratio times the reference.</summary>
        public const double InterferenceMargin = 1.5;

        /// <summary>
        /// Trapezoid area of the trace over scans <paramref name="start"/> to <paramref name="end"/> inclusive, less the
        /// baseline under the lower end point (that end's intensity across the width), as DIA library search integrates a
        /// fragment. A trace silent across the bounds gives 0, a real zero. NaN when the area is no measurement: bounds of
        /// one scan, which have no width, or a trace that lies below its own baseline.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="trace"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The bounds are not within the trace, or start is after end.</exception>
        public static double Area(IReadOnlyList<double> trace, int start, int end)
        {
            ArgumentNullException.ThrowIfNull(trace);
            CheckBounds(start, end, trace.Count);
            return AreaWithin(trace, start, end);
        }

        /// <summary>
        /// Each fragment's <see cref="Area"/> after DIA-NN's interference correction (Demichev et al. 2020). The best
        /// fragment (<see cref="FragmentCoElution.BestFragment"/>) is the one least affected by interference, and its
        /// smoothed trace (<see cref="FragmentCoElution.Smooth"/>) is the reference profile. A fragment's ratio to the best
        /// is the ratio of their intensities weighted by the squared reference, which stresses the scans near the apex. A
        /// fragment whose correlation with the reference is below <see cref="InterferenceCorrelation"/> is taken to be
        /// interfered, and its ratio is instead the smallest ratio to the best fragment within
        /// <paramref name="apexHalfWidth"/> scans of the best fragment's apex. Every point above
        /// <see cref="InterferenceMargin"/> times the ratio times the reference is lowered to it, and the area of what
        /// remains is the fragment's.
        /// </summary>
        /// <param name="traces">One trace per fragment over the same scans, zero where the fragment was not seen.</param>
        /// <returns>One area per trace, in trace order: 0 for every fragment of a silent precursor, NaN for every fragment
        /// when the bounds hold one scan.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="traces"/> or a trace is null.</exception>
        /// <exception cref="ArgumentException">There are no traces, or they differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The bounds are not within the traces, start is after end, or
        /// <paramref name="apexHalfWidth"/> is negative.</exception>
        public static double[] InterferenceCorrectedAreas(IReadOnlyList<double[]> traces, int start, int end, int apexHalfWidth = 1)
        {
            ArgumentNullException.ThrowIfNull(traces);
            if (traces.Count == 0)
                throw new ArgumentException("There are no fragment traces.", nameof(traces));
            foreach (var trace in traces)
            {
                ArgumentNullException.ThrowIfNull(trace, nameof(traces));
                if (trace.Length != traces[0].Length)
                    throw new ArgumentException("Every fragment trace must cover the same scans.", nameof(traces));
            }
            int length = traces[0].Length;
            CheckBounds(start, end, length);
            ArgumentOutOfRangeException.ThrowIfNegative(apexHalfWidth);

            var areas = new double[traces.Count];
            if (start == end)
            {
                Array.Fill(areas, double.NaN);
                return areas;
            }

            double[] best = traces[FragmentCoElution.BestFragment(traces, start, end)];
            double[] reference = FragmentCoElution.Smooth(best);
            double bestWeight = Weighted(best, reference, start, end);
            if (bestWeight <= 0)
            {
                // The best fragment is silent within the bounds, so there is no profile to correct against
                for (int f = 0; f < traces.Count; f++)
                    areas[f] = AreaWithin(traces[f], start, end);
                return areas;
            }

            double[] correlations = FragmentCoElution.CorrelationsTo(traces, reference, start, end);
            int apex = start;
            for (int s = start + 1; s <= end; s++)
                if (best[s] > best[apex])
                    apex = s;

            var clipped = new double[length];
            for (int f = 0; f < traces.Count; f++)
            {
                double[] trace = traces[f];
                double ratio = Weighted(trace, reference, start, end) / bestWeight;
                if (correlations[f] < InterferenceCorrelation)
                {
                    // The best fragment is positive at its apex, so the smallest ratio near it is always defined
                    ratio = double.PositiveInfinity;
                    for (int s = Math.Max(start, apex - apexHalfWidth); s <= Math.Min(end, apex + apexHalfWidth); s++)
                        if (best[s] > 0)
                            ratio = Math.Min(ratio, trace[s] / best[s]);
                }
                for (int s = start; s <= end; s++)
                    clipped[s] = Math.Min(trace[s], InterferenceMargin * ratio * reference[s]);
                areas[f] = AreaWithin(clipped, start, end);
            }
            return areas;
        }

        /// <summary>
        /// One set of fragments for every run of a precursor: the <paramref name="count"/> fragments with the highest mean
        /// correlation over the runs given, as DIA-NN chooses them (Demichev et al. 2020). The same set in every run keeps a
        /// precursor's quantities comparable across runs. A NaN correlation (a fragment not measured in that run) counts as
        /// 0 there. Ties go to the lower index, and the set is returned in index order.
        /// </summary>
        /// <param name="correlationsByRun">Per run, each fragment's correlation with the precursor's elution profile, in
        /// fragment order; typically only the runs where the precursor passed a q-value threshold.</param>
        /// <exception cref="ArgumentNullException"><paramref name="correlationsByRun"/> or a run is null.</exception>
        /// <exception cref="ArgumentException">There are no runs, or they differ in their number of fragments.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="count"/> is negative.</exception>
        public static int[] SelectAcrossRuns(IReadOnlyList<float[]> correlationsByRun, int count)
        {
            ArgumentNullException.ThrowIfNull(correlationsByRun);
            if (correlationsByRun.Count == 0)
                throw new ArgumentException("There are no runs to choose fragments from.", nameof(correlationsByRun));
            ArgumentNullException.ThrowIfNull(correlationsByRun[0], nameof(correlationsByRun));
            var sums = new double[correlationsByRun[0].Length];
            foreach (var run in correlationsByRun)
            {
                ArgumentNullException.ThrowIfNull(run, nameof(correlationsByRun));
                if (run.Length != sums.Length)
                    throw new ArgumentException("Every run must give one correlation per fragment.", nameof(correlationsByRun));
                for (int f = 0; f < run.Length; f++)
                    if (!float.IsNaN(run[f]))
                        sums[f] += run[f];
            }
            // The sums rank the fragments as their means do
            return FragmentCoElution.TopIndices(sums, count);
        }

        /// <summary>
        /// A precursor's quantity in one run from its fragments' areas, as a search stores them: the summed areas of the
        /// selected fragments that were measured. A NaN area is a fragment not measured, and it is left out, never imputed.
        /// <list type="bullet">
        /// <item><see cref="PrecursorQuantityKind.NotLooked"/>: fewer than <see cref="MinimumFragments"/> selected fragments
        /// were measured.</item>
        /// <item><see cref="PrecursorQuantityKind.BelowFloor"/>: every measured one is at or below <paramref name="floor"/>,
        /// silence included.</item>
        /// <item><see cref="PrecursorQuantityKind.Value"/>: otherwise.</item>
        /// </list>
        /// </summary>
        /// <param name="areas">Each fragment's area (<see cref="Area"/>, <see cref="InterferenceCorrectedAreas"/>), or NaN.</param>
        /// <param name="selected">The distinct fragments to sum (<see cref="SelectAcrossRuns"/>).</param>
        /// <param name="floor">The area at or below which a fragment shows nothing, in the areas' units.</param>
        /// <exception cref="ArgumentNullException"><paramref name="areas"/> or <paramref name="selected"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A selected fragment is not among the areas, an area is negative, or
        /// the floor is negative or not finite.</exception>
        public static PrecursorQuantity Quantify(IReadOnlyList<float> areas, IReadOnlyList<int> selected, double floor)
        {
            ArgumentNullException.ThrowIfNull(areas);
            ArgumentNullException.ThrowIfNull(selected);
            CheckFloor(floor);

            int measured = 0;
            bool aboveFloor = false;
            double sum = 0;
            foreach (int f in selected)
            {
                if (f < 0 || f >= areas.Count)
                    throw new ArgumentOutOfRangeException(nameof(selected), f, $"There is no fragment {f} among {areas.Count}.");
                float area = areas[f];
                if (float.IsNaN(area))
                    continue;
                CheckArea(area);
                measured++;
                sum += area;
                aboveFloor |= area > floor;
            }
            if (measured < MinimumFragments)
                return new PrecursorQuantity(double.NaN, PrecursorQuantityKind.NotLooked, floor);
            return aboveFloor
                ? new PrecursorQuantity(sum, PrecursorQuantityKind.Value, floor)
                : new PrecursorQuantity(double.NaN, PrecursorQuantityKind.BelowFloor, floor);
        }

        /// <summary>
        /// One trace's area as a quantity, such as an MS1 quantity: <see cref="PrecursorQuantityKind.NotLooked"/> when the
        /// area is NaN, <see cref="PrecursorQuantityKind.BelowFloor"/> at or below <paramref name="floor"/>, and a value
        /// above it.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">The area is negative, or the floor is negative or not finite.</exception>
        public static PrecursorQuantity Quantify(double area, double floor)
        {
            CheckFloor(floor);
            if (double.IsNaN(area))
                return new PrecursorQuantity(double.NaN, PrecursorQuantityKind.NotLooked, floor);
            CheckArea(area);
            return area > floor
                ? new PrecursorQuantity(area, PrecursorQuantityKind.Value, floor)
                : new PrecursorQuantity(double.NaN, PrecursorQuantityKind.BelowFloor, floor);
        }

        /// <summary><see cref="Area"/> on bounds already checked.</summary>
        private static double AreaWithin(IReadOnlyList<double> trace, int start, int end)
        {
            if (start == end)
                return double.NaN;
            // Summed relative to the baseline, so a flat trace is exactly 0 rather than a rounding error either side of it
            double baseline = Math.Min(trace[start], trace[end]);
            double area = 0;
            for (int s = start; s < end; s++)
                area += 0.5 * (trace[s] + trace[s + 1]) - baseline;
            return area < 0 ? double.NaN : area;
        }

        /// <summary>The trace's intensities over the bounds, each weighted by the squared reference.</summary>
        private static double Weighted(double[] trace, double[] reference, int start, int end)
        {
            double sum = 0;
            for (int s = start; s <= end; s++)
                sum += trace[s] * reference[s] * reference[s];
            return sum;
        }

        private static void CheckBounds(int start, int end, int length)
        {
            if (start < 0 || end >= length || start > end)
                throw new ArgumentOutOfRangeException(nameof(start), $"Scan range [{start}, {end}] is not within the trace's {length} scans.");
        }

        private static void CheckFloor(double floor)
        {
            if (!double.IsFinite(floor) || floor < 0)
                throw new ArgumentOutOfRangeException(nameof(floor), floor, "The floor must be a finite area of at least 0.");
        }

        private static void CheckArea(double area)
        {
            if (area < 0)
                throw new ArgumentOutOfRangeException(nameof(area), area, "An area is never negative; an area that is no measurement is NaN.");
        }
    }
}
