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
        /// Smooths a trace with weights 0.25/0.5/0.25. An end point uses the one neighbour it has, renormalized (2/3, 1/3),
        /// so a peak at the edge is not pulled toward zero.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="trace"/> is null.</exception>
        public static double[] Smooth(IReadOnlyList<double> trace)
        {
            ArgumentNullException.ThrowIfNull(trace);
            int n = trace.Count;
            var smoothed = new double[n];
            if (n == 1)
                smoothed[0] = trace[0];
            for (int i = 0; i < n && n > 1; i++)
            {
                if (i == 0)
                    smoothed[i] = (0.5 * trace[0] + 0.25 * trace[1]) / 0.75;
                else if (i == n - 1)
                    smoothed[i] = (0.5 * trace[i] + 0.25 * trace[i - 1]) / 0.75;
                else
                    smoothed[i] = 0.25 * trace[i - 1] + 0.5 * trace[i] + 0.25 * trace[i + 1];
            }
            return smoothed;
        }

        private static readonly double[] SavitzkyGolay9Weights = [-21, 14, 39, 54, 59, 54, 39, 14, -21];

        /// <summary>
        /// Smooths a trace with the 9-point quadratic/cubic Savitzky–Golay filter (weights −21, 14, 39, 54, 59, 54, 39, 14, −21
        /// over 231), as EncyclopeDIA and Skyline smooth fragment chromatograms before picking peaks. It keeps peak height and
        /// width better than a moving average. Points within four scans of an end are left unchanged, and negative ringing
        /// beside sharp spikes is clipped to 0.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="trace"/> is null.</exception>
        public static double[] SavitzkyGolay9(IReadOnlyList<double> trace)
        {
            ArgumentNullException.ThrowIfNull(trace);
            var smoothed = trace.ToArray();
            for (int s = 4; s < trace.Count - 4; s++)
            {
                double sum = 0;
                for (int w = 0; w < 9; w++)
                    sum += SavitzkyGolay9Weights[w] * trace[s - 4 + w];
                smoothed[s] = Math.Max(0, sum / 231);
            }
            return smoothed;
        }

        /// <summary>
        /// The fragment whose trace, over scans <paramref name="from"/> to <paramref name="to"/>, has the largest summed
        /// Pearson correlation with the other fragments (negative or undefined correlations count as 0). Its smoothed
        /// trace makes a robust elution profile to score the others against (<see cref="CorrelationsTo"/>): one reliable
        /// fragment is harder to corrupt than an average that an interfered fragment drags along. Ties go to the lower index.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="traces"/> is null.</exception>
        /// <exception cref="ArgumentException">There are no traces, or they differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The range is empty or outside the traces.</exception>
        public static int BestFragment(IReadOnlyList<double[]> traces, int from, int to)
        {
            int length = ValidateTraces(traces);
            if (traces.Count == 0)
                throw new ArgumentException("There are no fragment traces.", nameof(traces));
            CheckRange(from, to, length);

            int best = 0;
            double bestSum = double.NegativeInfinity;
            for (int f = 0; f < traces.Count; f++)
            {
                double sum = 0;
                for (int g = 0; g < traces.Count; g++)
                    if (g != f)
                        sum += ClippedPearson(traces[f], traces[g], from, to);
                if (sum > bestSum)
                {
                    bestSum = sum;
                    best = f;
                }
            }
            return best;
        }

        /// <summary>
        /// Each fragment's Pearson correlation with <paramref name="reference"/> over scans <paramref name="from"/> to
        /// <paramref name="to"/>, in trace order. Negative correlations, and flat or silent traces, count as 0.
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A trace and the reference differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The range is empty or outside the traces.</exception>
        public static double[] CorrelationsTo(IReadOnlyList<double[]> traces, IReadOnlyList<double> reference, int from, int to)
        {
            int length = ValidateTraces(traces);
            ArgumentNullException.ThrowIfNull(reference);
            if (traces.Count > 0 && reference.Count != length)
                throw new ArgumentException($"The reference has {reference.Count} scans but the traces have {length}.", nameof(reference));
            CheckRange(from, to, reference.Count);

            var referenceArray = reference as double[] ?? reference.ToArray();
            return traces.Select(trace => ClippedPearson(trace, referenceArray, from, to)).ToArray();
        }

        private static void CheckRange(int from, int to, int length)
        {
            if (from < 0 || to >= length || from > to)
                throw new ArgumentOutOfRangeException(nameof(from), $"Scan range [{from}, {to}] is not within the traces' {length} scans.");
        }

        /// <summary>Pearson over [from, to]; negative, flat or undefined counts as 0.</summary>
        private static double ClippedPearson(double[] a, double[] b, int from, int to)
        {
            int points = to - from + 1;
            if (points < 2)
                return 0;
            double r = Correlation.Pearson(a.Skip(from).Take(points), b.Skip(from).Take(points));
            return double.IsFinite(r) && r > 0 ? r : 0;
        }

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

            double[] values = ApexScores(traces, libraryIntensities.ToArray(), halfWidth, length);
            int apex = -1;
            double best = 0;
            for (int s = 0; s < length; s++)
            {
                if (values[s] > best)
                {
                    best = values[s];
                    apex = s;
                }
            }
            return apex;
        }

        /// <summary>
        /// Candidate apexes: scans whose apex score (see <see cref="FindApex"/>) is a local maximum, more than
        /// <paramref name="halfWidth"/> scans from any better candidate, best first, at most <paramref name="maxCount"/>.
        /// The first is <see cref="FindApex"/>'s apex. Scoring several candidates lets a caller recover a precursor whose
        /// true peak scores second to interference.
        /// </summary>
        /// <exception cref="ArgumentException">The traces differ in length, or there is not one library intensity per trace.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="halfWidth"/> is negative or <paramref name="maxCount"/> is less than 1.</exception>
        public static int[] FindApexes(IReadOnlyList<double[]> traces, IReadOnlyList<double> libraryIntensities, int halfWidth, int maxCount)
        {
            int length = ValidateTraces(traces);
            ArgumentNullException.ThrowIfNull(libraryIntensities);
            if (libraryIntensities.Count != traces.Count)
                throw new ArgumentException("There must be one library intensity per fragment trace.", nameof(libraryIntensities));
            ArgumentOutOfRangeException.ThrowIfNegative(halfWidth);
            ArgumentOutOfRangeException.ThrowIfLessThan(maxCount, 1);

            return FindApexes(ApexScores(traces, libraryIntensities.ToArray(), halfWidth, length), halfWidth, maxCount);
        }

        /// <summary>
        /// Candidate apexes from precomputed apex scores (<see cref="ApexScores(IReadOnlyList{double[]}, IReadOnlyList{double}, int)"/>),
        /// so a caller that also needs the scores computes them once. Same rules as the overload on traces.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="apexScores"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="halfWidth"/> is negative or <paramref name="maxCount"/> is less than 1.</exception>
        public static int[] FindApexes(IReadOnlyList<double> apexScores, int halfWidth, int maxCount)
        {
            ArgumentNullException.ThrowIfNull(apexScores);
            ArgumentOutOfRangeException.ThrowIfNegative(halfWidth);
            ArgumentOutOfRangeException.ThrowIfLessThan(maxCount, 1);
            double[] values = apexScores as double[] ?? apexScores.ToArray();
            int length = values.Length;
            var chosen = new List<int>();
            foreach (int s in Enumerable.Range(0, length).Where(s => values[s] > 0 && IsLocalMaximum(values, s, halfWidth)).OrderByDescending(s => values[s]).ThenBy(s => s))
            {
                if (chosen.All(c => Math.Abs(c - s) > halfWidth))
                    chosen.Add(s);
                if (chosen.Count == maxCount)
                    break;
            }
            return chosen.ToArray();
        }

        /// <summary>True when no scan within <paramref name="halfWidth"/> scores higher, and none earlier scores the same.</summary>
        private static bool IsLocalMaximum(double[] values, int s, int halfWidth)
        {
            for (int t = Math.Max(0, s - halfWidth); t <= Math.Min(values.Length - 1, s + halfWidth); t++)
                if (values[t] > values[s] || (t < s && values[t] == values[s]))
                    return false;
            return true;
        }

        /// <summary>
        /// The apex score of every scan (see <see cref="FindApex"/>): cosine to the library × co-elution around it ×
        /// log(1 + signal), and 0 without signal. It lets a caller judge how far the best peak stands out from the rest of
        /// the window, as PECAN's deltaSn or a z-score does.
        /// </summary>
        /// <exception cref="ArgumentException">The traces differ in length, or there is not one library intensity per trace.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="halfWidth"/> is negative.</exception>
        public static double[] ApexScores(IReadOnlyList<double[]> traces, IReadOnlyList<double> libraryIntensities, int halfWidth)
        {
            int length = ValidateTraces(traces);
            ArgumentNullException.ThrowIfNull(libraryIntensities);
            if (libraryIntensities.Count != traces.Count)
                throw new ArgumentException("There must be one library intensity per fragment trace.", nameof(libraryIntensities));
            ArgumentOutOfRangeException.ThrowIfNegative(halfWidth);
            return ApexScores(traces, libraryIntensities.ToArray(), halfWidth, length);
        }

        private static double[] ApexScores(IReadOnlyList<double[]> traces, double[] library, int halfWidth, int length)
        {
            var values = new double[length];
            var observed = new double[traces.Count];
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
                values[s] = SpectralSimilarity.CosineOfAlignedVectors(observed, library)
                    * Score(traces, Math.Max(0, s - halfWidth), Math.Min(length - 1, s + halfWidth))
                    * Math.Log(1 + signal);
            }
            return values;
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
