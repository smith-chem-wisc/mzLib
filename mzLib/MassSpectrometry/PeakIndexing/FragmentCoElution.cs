#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics.Statistics;

namespace MassSpectrometry
{
    /// <summary>
    /// Scores how well a precursor's fragment ions co-elute. The input is a dense intensity matrix: one trace per fragment
    /// over the same consecutive scans (typically the MS2 scans of one DIA isolation window, in retention-time order), with
    /// zero where a fragment was not seen. A real precursor's fragments rise and fall together; coincidental peaks at its
    /// fragment m/z do not.
    /// </summary>
    public static class FragmentCoElution
    {
        /// <summary>
        /// Mean, over fragments, of the Pearson correlation between each fragment's trace and the sum of the other
        /// fragments' traces over scans <paramref name="from"/> to <paramref name="to"/> inclusive. Before summing, each
        /// other trace is scaled to its own maximum, so a single large interfering peak counts as one fragment rather than
        /// swamping the rest.
        /// <para>
        /// Negative correlations count as 0, and so does a fragment with no signal in the range. A missing fragment
        /// therefore lowers the group's score. The result is in [0, 1] and never NaN; silence scores 0.
        /// </para>
        /// </summary>
        /// <param name="traces">One trace per fragment, all the same length.</param>
        /// <exception cref="ArgumentNullException"><paramref name="traces"/> is null.</exception>
        /// <exception cref="ArgumentException">The traces differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The range is empty or outside the traces.</exception>
        public static double Score(IReadOnlyList<double[]> traces, int from, int to)
        {
            int length = ValidateTraces(traces);
            if (from < 0 || to >= length || from > to)
                throw new ArgumentOutOfRangeException(nameof(from), $"Scan range [{from}, {to}] is not within the traces' {length} scans.");

            int points = to - from + 1;
            if (traces.Count < 2 || points < 3)
                return 0;

            // Each trace scaled to its own maximum over the range, and their sum
            var scaled = new double[traces.Count][];
            var total = new double[points];
            for (int f = 0; f < traces.Count; f++)
            {
                double max = 0;
                for (int s = from; s <= to; s++)
                    max = Math.Max(max, traces[f][s]);
                scaled[f] = new double[points];
                if (max <= 0)
                    continue;
                for (int s = 0; s < points; s++)
                {
                    scaled[f][s] = traces[f][from + s] / max;
                    total[s] += scaled[f][s];
                }
            }

            double sum = 0;
            var others = new double[points];
            for (int f = 0; f < traces.Count; f++)
            {
                for (int s = 0; s < points; s++)
                    others[s] = total[s] - scaled[f][s];
                double r = Correlation.Pearson(scaled[f], others);
                // A flat trace (no signal, or no signal elsewhere) has no correlation to report
                if (double.IsFinite(r) && r > 0)
                    sum += r;
            }
            return sum / traces.Count;
        }

        /// <summary>
        /// The scan where the fragments best look like the precursor: observed intensities in library proportions
        /// (cosine), co-eluting over the <paramref name="halfWidth"/> scans on either side (<see cref="Score"/>), weighted
        /// by the logarithm of the summed intensity. This is not simply the scan with the most signal, which a single
        /// interfering spike would win.
        /// </summary>
        /// <param name="libraryIntensities">Expected relative intensity of each fragment, in trace order.</param>
        /// <returns>The apex scan index, or -1 when no scan carries any fragment signal.</returns>
        /// <exception cref="ArgumentException">The traces differ in length, or there is not one library intensity per trace.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="halfWidth"/> is negative.</exception>
        public static int FindApex(IReadOnlyList<double[]> traces, IReadOnlyList<double> libraryIntensities, int halfWidth)
        {
            int length = ValidateTraces(traces);
            ArgumentNullException.ThrowIfNull(libraryIntensities);
            if (libraryIntensities.Count != traces.Count)
                throw new ArgumentException("There must be one library intensity per fragment trace.", nameof(libraryIntensities));
            ArgumentOutOfRangeException.ThrowIfNegative(halfWidth);

            double[] library = libraryIntensities.ToArray();
            var observed = new double[traces.Count];
            int apex = -1;
            double best = 0;
            for (int s = 0; s < length; s++)
            {
                double signal = 0;
                for (int f = 0; f < traces.Count; f++)
                {
                    observed[f] = traces[f][s];
                    signal += observed[f];
                }
                if (signal <= 0)
                    continue;

                double value = SpectralSimilarity.CosineOfAlignedVectors(observed, library)
                    * Score(traces, Math.Max(0, s - halfWidth), Math.Min(length - 1, s + halfWidth))
                    * Math.Log(1 + signal);
                if (value > best)
                {
                    best = value;
                    apex = s;
                }
            }
            return apex;
        }

        /// <summary>
        /// The peak around <paramref name="apex"/> on a dense trace, as inclusive scan bounds. The bounds are found by
        /// <see cref="ExtractedIonChromatogram.FindPeakBoundaries"/> (unchanged; validated for dense, zero-filled traces)
        /// and lie just inside the valleys it reports. A side with no valley extends to the end of the trace.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="apex"/> is outside the trace.</exception>
        public static (int Start, int End) PeakBounds(IReadOnlyList<double> trace, int apex)
        {
            ArgumentNullException.ThrowIfNull(trace);
            if (apex < 0 || apex >= trace.Count)
                throw new ArgumentOutOfRangeException(nameof(apex), apex, $"The apex must be a scan of the {trace.Count}-scan trace.");

            var peaks = new List<IIndexedPeak>(trace.Count);
            for (int s = 0; s < trace.Count; s++)
                peaks.Add(new IndexedMassSpectralPeak(0, trace[s], s, s));

            int start = 0, end = trace.Count - 1;
            foreach (var valley in ExtractedIonChromatogram.FindPeakBoundaries(peaks, apex))
            {
                if (valley.ZeroBasedScanIndex < apex)
                    start = Math.Max(start, valley.ZeroBasedScanIndex + 1);
                else if (valley.ZeroBasedScanIndex > apex)
                    end = Math.Min(end, valley.ZeroBasedScanIndex - 1);
            }
            return (start, end);
        }

        /// <summary>
        /// The indices of the <paramref name="count"/> largest values (earlier index first among ties), returned in
        /// ascending index order. Used to keep a library entry's most intense fragments.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="count"/> is negative.</exception>
        public static int[] TopIndices(IReadOnlyList<double> values, int count)
        {
            ArgumentNullException.ThrowIfNull(values);
            ArgumentOutOfRangeException.ThrowIfNegative(count);
            return Enumerable.Range(0, values.Count)
                .OrderByDescending(i => values[i]).ThenBy(i => i)
                .Take(count)
                .Order()
                .ToArray();
        }

        private static int ValidateTraces(IReadOnlyList<double[]> traces)
        {
            ArgumentNullException.ThrowIfNull(traces);
            if (traces.Count == 0)
                return 0;
            int length = traces[0].Length;
            for (int f = 1; f < traces.Count; f++)
                if (traces[f].Length != length)
                    throw new ArgumentException("Every fragment trace must cover the same scans.", nameof(traces));
            return length;
        }
    }
}
