using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// Median peak width over a group of peaks, which is the form the metric is useful in for quality
    /// control: an individual peak's width says little, the central tendency across a file says whether
    /// the column and gradient are behaving.
    /// </summary>
    /// <remarks>
    /// Binning by retention time is offered because on real data the median is not flat across a
    /// gradient. On a 176 minute mouse run it was 0.372 min early, 0.260 min in the middle and 0.393 min
    /// late, so a single median per file averages over a real trend.
    /// </remarks>
    public readonly struct PeakWidthSummary
    {
        public string FileName { get; init; }

        /// <summary>Apex retention time of the earliest peak in the bin.</summary>
        public double RetentionTimeStart { get; init; }

        /// <summary>Apex retention time of the latest peak in the bin.</summary>
        public double RetentionTimeEnd { get; init; }

        /// <summary>Peaks with a measured width. Peaks without one are not counted.</summary>
        public int MeasuredPeakCount { get; init; }

        /// <summary>Peaks in the bin regardless of whether a width could be measured.</summary>
        public int TotalPeakCount { get; init; }

        public double MedianFullWidthAtHalfMaximum { get; init; }
        public double LowerQuartile { get; init; }
        public double UpperQuartile { get; init; }

        /// <summary>
        /// Median width per file, optionally split into equal-count retention time bins. Peaks whose
        /// width could not be measured are excluded from the median but still counted in
        /// TotalPeakCount, so a caller can see how much of the file the median rests on.
        /// </summary>
        /// <param name="peaks">Peaks to summarise; grouped by spectra file internally.</param>
        /// <param name="retentionTimeBins">
        /// Bins per file. 1 gives one summary per file. Bins hold equal numbers of peaks rather than
        /// spanning equal time, so an early or late stretch with few identifications does not produce a
        /// median resting on a handful of peaks.
        /// </param>
        public static List<PeakWidthSummary> Summarize(IEnumerable<ChromatographicPeak> peaks, int retentionTimeBins = 1)
        {
            if (peaks == null)
            {
                throw new ArgumentNullException(nameof(peaks));
            }

            if (retentionTimeBins < 1)
            {
                throw new ArgumentOutOfRangeException(nameof(retentionTimeBins), "at least one bin is required");
            }

            var summaries = new List<PeakWidthSummary>();

            foreach (var byFile in peaks
                .Where(p => p?.SpectraFileInfo != null)
                // Grouped by the SpectraFileInfo itself, not its base name. Two files in different
                // directories, or two fractions or replicates that happen to share a name, are
                // different runs -- merging them would report one summary covering both, with counts
                // and a retention range belonging to neither.
                .GroupBy(p => p.SpectraFileInfo)
                .OrderBy(g => g.Key.FullFilePathWithExtension, StringComparer.Ordinal))
            {
                // A peak with no apex reports ApexRetentionTime as -1, so it would sort ahead of every
                // real peak, take the first bin's RetentionTimeStart down to -1, and inflate that bin's
                // TotalPeakCount. It cannot contribute a width either, so it is not summarised at all.
                List<ChromatographicPeak> ordered = byFile
                    .Where(p => p.Apex != null)
                    .OrderBy(p => p.ApexRetentionTime)
                    .ToList();

                // Equal-count bins. Integer division leaves a remainder, which goes to the last bin
                // rather than becoming a short extra bin.
                int binCount = Math.Min(retentionTimeBins, ordered.Count);
                if (binCount == 0)
                {
                    continue;
                }

                int binSize = ordered.Count / binCount;

                for (int bin = 0; bin < binCount; bin++)
                {
                    int start = bin * binSize;
                    int count = bin == binCount - 1 ? ordered.Count - start : binSize;
                    List<ChromatographicPeak> inBin = ordered.GetRange(start, count);

                    List<double> widths = inBin
                        .Select(p => p.PeakWidth)
                        .Where(w => w.IsMeasured)
                        .Select(w => w.FullWidthAtHalfMaximum)
                        .OrderBy(x => x)
                        .ToList();

                    summaries.Add(new PeakWidthSummary
                    {
                        FileName = byFile.Key.FilenameWithoutExtension,
                        RetentionTimeStart = inBin.First().ApexRetentionTime,
                        RetentionTimeEnd = inBin.Last().ApexRetentionTime,
                        MeasuredPeakCount = widths.Count,
                        TotalPeakCount = inBin.Count,
                        MedianFullWidthAtHalfMaximum = Quantile(widths, 0.50),
                        LowerQuartile = Quantile(widths, 0.25),
                        UpperQuartile = Quantile(widths, 0.75)
                    });
                }
            }

            return summaries;
        }

        /// <summary>
        /// Linearly interpolated quantile of an already-sorted list. NaN when the list is empty, so an
        /// all-censored bin reports no median rather than a zero that would read as a real measurement.
        /// </summary>
        private static double Quantile(List<double> sorted, double quantile)
        {
            if (sorted.Count == 0)
            {
                return double.NaN;
            }

            double position = quantile * (sorted.Count - 1);
            int lower = (int)Math.Floor(position);
            int upper = (int)Math.Ceiling(position);

            return lower == upper
                ? sorted[lower]
                : sorted[lower] + ((position - lower) * (sorted[upper] - sorted[lower]));
        }

        public override string ToString() =>
            $"{FileName} {RetentionTimeStart:F2}-{RetentionTimeEnd:F2} min " +
            $"n={MeasuredPeakCount}/{TotalPeakCount} median={MedianFullWidthAtHalfMaximum:F3}";
    }
}
