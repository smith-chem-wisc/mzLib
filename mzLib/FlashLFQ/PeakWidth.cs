using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// Why a peak has no usable full width at half maximum.
    /// </summary>
    public enum PeakWidthStatus
    {
        /// <summary>The trace fell below half maximum on both sides of the apex.</summary>
        Measured,

        /// <summary>No apex, so no trace to measure.</summary>
        NoApex,

        /// <summary>
        /// Too few time points on the apex charge state to describe a shape. Matches the floor
        /// FlashLfqEngine.CutPeak uses before it looks for valleys.
        /// </summary>
        TooFewPoints,

        /// <summary>
        /// The trace ends while still above half maximum on one or both sides, so the real width is
        /// larger than anything measurable here. FlashLFQ stops extending an XIC after
        /// MissedScansAllowed consecutive misses and CutPeak truncates at valleys, so this is
        /// expected for wide or co-eluting peaks rather than a failure.
        /// </summary>
        Censored
    }

    /// <summary>
    /// Full width at half maximum of a chromatographic peak, in minutes.
    /// </summary>
    /// <remarks>
    /// Measured on the apex charge state alone, which is the same trace FlashLfqEngine.CutPeak uses
    /// for valley detection. ChromatographicPeak.IsotopicEnvelopes interleaves charge states, so an
    /// RT-ordered walk over all of them alternates between charges and describes a sawtooth rather
    /// than a chromatogram.
    /// </remarks>
    public readonly struct PeakWidth
    {
        /// <summary>Minimum time points on the apex charge state before a width is attempted.</summary>
        public const int MinimumTimePoints = 5;

        public PeakWidthStatus Status { get; init; }

        /// <summary>Width in minutes, or NaN when Status is not Measured.</summary>
        public double FullWidthAtHalfMaximum { get; init; }

        /// <summary>Interpolated retention time of the rising half-maximum crossing, or NaN.</summary>
        public double HalfMaximumStart { get; init; }

        /// <summary>Interpolated retention time of the falling half-maximum crossing, or NaN.</summary>
        public double HalfMaximumEnd { get; init; }

        /// <summary>
        /// Time points on the apex charge state's trace, in full -- not the number falling between
        /// <see cref="HalfMaximumStart"/> and <see cref="HalfMaximumEnd"/>.
        /// </summary>
        /// <remarks>
        /// The distinction matters for the duty-cycle question this type is partly meant to answer.
        /// The usual criterion for a reliable area is a healthy number of MS1 points across the peak,
        /// and "across the peak" conventionally means across the full width at half maximum rather
        /// than across the full base -- so this count is the larger of the two and will read as
        /// comfortably above a threshold that the points within the width would not clear. Count the
        /// envelopes between the two crossings if that is the number you want.
        /// </remarks>
        public int TimePointsOnApexCharge { get; init; }

        public bool IsMeasured => Status == PeakWidthStatus.Measured;

        private static PeakWidth Unmeasured(PeakWidthStatus status, int timePoints) => new PeakWidth
        {
            Status = status,
            FullWidthAtHalfMaximum = double.NaN,
            HalfMaximumStart = double.NaN,
            HalfMaximumEnd = double.NaN,
            TimePointsOnApexCharge = timePoints
        };

        /// <summary>
        /// Measures the full width at half maximum of a peak's apex-charge trace.
        /// </summary>
        public static PeakWidth Measure(ChromatographicPeak peak) => Measure(peak, MinimumTimePoints);

        /// <summary>
        /// Overload taking an explicit point floor, for sensitivity analysis of that floor.
        /// </summary>
        public static PeakWidth Measure(ChromatographicPeak peak, int minimumTimePoints)
        {
            if (peak?.Apex == null || peak.IsotopicEnvelopes == null)
            {
                return Unmeasured(PeakWidthStatus.NoApex, 0);
            }

            // The same filter CutPeak applies, for the same reason.
            List<IsotopicEnvelope> trace = peak.IsotopicEnvelopes
                .Where(p => p.ChargeState == peak.Apex.ChargeState)
                .OrderBy(p => p.IndexedPeak.RetentionTime)
                .ToList();

            if (trace.Count < minimumTimePoints)
            {
                return Unmeasured(PeakWidthStatus.TooFewPoints, trace.Count);
            }

            // Apex within the filtered trace rather than peak.Apex directly: the trace is re-sorted by
            // retention time, and reference equality would break if envelopes were ever copied.
            int apexIndex = 0;
            for (int i = 1; i < trace.Count; i++)
            {
                if (trace[i].Intensity > trace[apexIndex].Intensity)
                {
                    apexIndex = i;
                }
            }

            double halfMaximum = trace[apexIndex].Intensity / 2.0;

            // Walk outwards to the nearest crossing on each side, not the first and last crossing of the
            // whole trace: a shoulder that dips below half maximum and comes back would otherwise be
            // measured as one very wide peak.
            double? risingCrossing = FindCrossing(trace, apexIndex, -1, halfMaximum);
            double? fallingCrossing = FindCrossing(trace, apexIndex, +1, halfMaximum);

            if (risingCrossing == null || fallingCrossing == null)
            {
                return Unmeasured(PeakWidthStatus.Censored, trace.Count);
            }

            return new PeakWidth
            {
                Status = PeakWidthStatus.Measured,
                FullWidthAtHalfMaximum = fallingCrossing.Value - risingCrossing.Value,
                HalfMaximumStart = risingCrossing.Value,
                HalfMaximumEnd = fallingCrossing.Value,
                TimePointsOnApexCharge = trace.Count
            };
        }

        /// <summary>
        /// Retention time at which the trace crosses half maximum, walking from the apex in one
        /// direction. Linearly interpolated between the bracketing points: taking the point itself
        /// would quantise every width to a multiple of the MS1 cycle time. Null if the trace runs out
        /// before crossing.
        /// </summary>
        private static double? FindCrossing(List<IsotopicEnvelope> trace, int apexIndex, int direction, double halfMaximum)
        {
            for (int i = apexIndex + direction; i >= 0 && i < trace.Count; i += direction)
            {
                if (trace[i].Intensity > halfMaximum)
                {
                    continue;
                }

                IsotopicEnvelope below = trace[i];
                IsotopicEnvelope above = trace[i - direction];

                double intensitySpan = above.Intensity - below.Intensity;
                double timeSpan = below.IndexedPeak.RetentionTime - above.IndexedPeak.RetentionTime;

                // Equal intensities would mean above is also at half maximum; take that point.
                if (intensitySpan <= 0)
                {
                    return below.IndexedPeak.RetentionTime;
                }

                double fraction = (above.Intensity - halfMaximum) / intensitySpan;
                return above.IndexedPeak.RetentionTime + (fraction * timeSpan);
            }

            return null;
        }
    }
}
