using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace TopDownEngine.Features;

public readonly record struct BinomialMatchScore(int KTrue, double PNull, double PValue);

public static class BinomialScorer
{
    private static readonly MethodInfo GetBinsInRangeMethod =
        typeof(IndexingEngine<IndexedMassSpectralPeak>).GetMethod("GetBinsInRange", BindingFlags.Instance | BindingFlags.NonPublic)
        ?? throw new MzLibException("Unable to locate GetBinsInRange helper on PeakIndexingEngine.");

    public static BinomialMatchScore ScoreMatch(
        FeatureBox donorBox,
        PeakIndexingEngine acceptorIndex,
        Func<double, double> rtWarp,
        IReadOnlyList<double> nullShifts)
    {
        ValidateScoreInputs(donorBox, acceptorIndex, rtWarp, nullShifts);

        int n = donorBox.PeakCount;
        DoubleRange warpedRtRange = WarpRtRange(donorBox.RtRange, rtWarp);
        DoubleRange mzRange = new(donorBox.MzRange.Minimum, donorBox.MzRange.Maximum);

        int kTrue = CountPeaksInBox(acceptorIndex, warpedRtRange, mzRange);
        List<int> nullCounts = new();
        nullCounts.AddRange(GetRtShiftNullCounts(acceptorIndex, warpedRtRange, mzRange, nullShifts));
        nullCounts.AddRange(GetMzShiftNullCounts(acceptorIndex, warpedRtRange, mzRange, nullShifts));
        nullCounts.AddRange(GetCombinedShiftNullCounts(acceptorIndex, warpedRtRange, mzRange, nullShifts));

        return BuildScore(kTrue, n, nullCounts);
    }

    public static BinomialMatchScore ScoreMatchRtOnlyShifts(
        FeatureBox donorBox,
        PeakIndexingEngine acceptorIndex,
        Func<double, double> rtWarp,
        IReadOnlyList<double> nullShifts)
    {
        ValidateScoreInputs(donorBox, acceptorIndex, rtWarp, nullShifts);

        int n = donorBox.PeakCount;
        DoubleRange warpedRtRange = WarpRtRange(donorBox.RtRange, rtWarp);
        DoubleRange mzRange = new(donorBox.MzRange.Minimum, donorBox.MzRange.Maximum);

        int kTrue = CountPeaksInBox(acceptorIndex, warpedRtRange, mzRange);
        List<int> nullCounts = GetRtShiftNullCounts(acceptorIndex, warpedRtRange, mzRange, nullShifts);
        return BuildScore(kTrue, n, nullCounts);
    }

    public static BinomialMatchScore ScoreMatchMzOnlyShifts(
        FeatureBox donorBox,
        PeakIndexingEngine acceptorIndex,
        Func<double, double> rtWarp,
        IReadOnlyList<double> nullShifts)
    {
        ValidateScoreInputs(donorBox, acceptorIndex, rtWarp, nullShifts);

        int n = donorBox.PeakCount;
        DoubleRange warpedRtRange = WarpRtRange(donorBox.RtRange, rtWarp);
        DoubleRange mzRange = new(donorBox.MzRange.Minimum, donorBox.MzRange.Maximum);

        int kTrue = CountPeaksInBox(acceptorIndex, warpedRtRange, mzRange);
        List<int> nullCounts = GetMzShiftNullCounts(acceptorIndex, warpedRtRange, mzRange, nullShifts);
        return BuildScore(kTrue, n, nullCounts);
    }

    public static BinomialMatchScore ScoreMatchCombinedShifts(
        FeatureBox donorBox,
        PeakIndexingEngine acceptorIndex,
        Func<double, double> rtWarp,
        IReadOnlyList<double> nullShifts)
    {
        ValidateScoreInputs(donorBox, acceptorIndex, rtWarp, nullShifts);

        int n = donorBox.PeakCount;
        DoubleRange warpedRtRange = WarpRtRange(donorBox.RtRange, rtWarp);
        DoubleRange mzRange = new(donorBox.MzRange.Minimum, donorBox.MzRange.Maximum);

        int kTrue = CountPeaksInBox(acceptorIndex, warpedRtRange, mzRange);
        List<int> nullCounts = GetCombinedShiftNullCounts(acceptorIndex, warpedRtRange, mzRange, nullShifts);
        return BuildScore(kTrue, n, nullCounts);
    }

    public static int CountPeaksInBox(PeakIndexingEngine index, DoubleRange rtRange, DoubleRange mzRange)
    {
        if (index == null)
        {
            throw new ArgumentNullException(nameof(index));
        }

        if (rtRange == null)
        {
            throw new ArgumentNullException(nameof(rtRange));
        }

        if (mzRange == null)
        {
            throw new ArgumentNullException(nameof(mzRange));
        }

        if (index.ScanInfoArray == null || rtRange.Maximum <= rtRange.Minimum || mzRange.Maximum <= mzRange.Minimum)
        {
            return 0;
        }

        HashSet<int> scanIndicesInWindow = new();
        for (int i = 0; i < index.ScanInfoArray.Length; i++)
        {
            ScanInfo scanInfo = index.ScanInfoArray[i];
            double retentionTime = scanInfo.RetentionTime;
            if (retentionTime >= rtRange.Minimum && retentionTime < rtRange.Maximum)
            {
                scanIndicesInWindow.Add(scanInfo.ZeroBasedScanIndex);
            }
        }

        if (scanIndicesInWindow.Count == 0)
        {
            return 0;
        }

        double centerMz = (mzRange.Minimum + mzRange.Maximum) / 2.0;
        double halfWidth = (mzRange.Maximum - mzRange.Minimum) / 2.0;
        AbsoluteTolerance mzTolerance = new(halfWidth);

        var bins = (List<List<IndexedMassSpectralPeak>>)GetBinsInRangeMethod.Invoke(index, new object[] { centerMz, mzTolerance });

        int count = 0;
        foreach (List<IndexedMassSpectralPeak> bin in bins)
        {
            foreach (IndexedMassSpectralPeak peak in bin)
            {
                if (scanIndicesInWindow.Contains(peak.ZeroBasedScanIndex)
                    && peak.M >= mzRange.Minimum
                    && peak.M < mzRange.Maximum)
                {
                    count++;
                }
            }
        }

        return count;
    }

    private static void ValidateScoreInputs(
        FeatureBox donorBox,
        PeakIndexingEngine acceptorIndex,
        Func<double, double> rtWarp,
        IReadOnlyList<double> nullShifts)
    {
        if (acceptorIndex == null)
        {
            throw new ArgumentNullException(nameof(acceptorIndex));
        }

        if (rtWarp == null)
        {
            throw new ArgumentNullException(nameof(rtWarp));
        }

        if (nullShifts == null)
        {
            throw new ArgumentNullException(nameof(nullShifts));
        }

        if (donorBox.PeakCount <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(donorBox), "donorBox.PeakCount must be positive.");
        }

        if (nullShifts.Count == 0)
        {
            throw new ArgumentException("At least one null shift is required.", nameof(nullShifts));
        }

        if (nullShifts.Any(s => double.IsNaN(s) || double.IsInfinity(s) || s <= 0))
        {
            throw new ArgumentException("All null shifts must be finite and positive.", nameof(nullShifts));
        }
    }

    private static DoubleRange WarpRtRange(DoubleRange donorRtRange, Func<double, double> rtWarp)
    {
        double warpedMin = rtWarp(donorRtRange.Minimum);
        double warpedMax = rtWarp(donorRtRange.Maximum);
        return new DoubleRange(Math.Min(warpedMin, warpedMax), Math.Max(warpedMin, warpedMax));
    }

    private static List<int> GetRtShiftNullCounts(
        PeakIndexingEngine acceptorIndex,
        DoubleRange warpedRtRange,
        DoubleRange mzRange,
        IReadOnlyList<double> nullShifts)
    {
        List<int> counts = new(nullShifts.Count * 2);
        foreach (double shift in nullShifts)
        {
            counts.Add(CountPeaksInBox(acceptorIndex, ShiftRtRange(warpedRtRange, shift), mzRange));
            counts.Add(CountPeaksInBox(acceptorIndex, ShiftRtRange(warpedRtRange, -shift), mzRange));
        }

        return counts;
    }

    private static List<int> GetMzShiftNullCounts(
        PeakIndexingEngine acceptorIndex,
        DoubleRange warpedRtRange,
        DoubleRange mzRange,
        IReadOnlyList<double> nullShifts)
    {
        List<int> counts = new(nullShifts.Count * 2);
        foreach (double shift in nullShifts)
        {
            counts.Add(CountPeaksInBox(acceptorIndex, warpedRtRange, ShiftMzRange(mzRange, shift)));
            counts.Add(CountPeaksInBox(acceptorIndex, warpedRtRange, ShiftMzRange(mzRange, -shift)));
        }

        return counts;
    }

    private static List<int> GetCombinedShiftNullCounts(
        PeakIndexingEngine acceptorIndex,
        DoubleRange warpedRtRange,
        DoubleRange mzRange,
        IReadOnlyList<double> nullShifts)
    {
        List<int> counts = new(nullShifts.Count * 4);
        foreach (double shift in nullShifts)
        {
            counts.Add(CountPeaksInBox(acceptorIndex, ShiftRtRange(warpedRtRange, shift), ShiftMzRange(mzRange, shift)));
            counts.Add(CountPeaksInBox(acceptorIndex, ShiftRtRange(warpedRtRange, shift), ShiftMzRange(mzRange, -shift)));
            counts.Add(CountPeaksInBox(acceptorIndex, ShiftRtRange(warpedRtRange, -shift), ShiftMzRange(mzRange, shift)));
            counts.Add(CountPeaksInBox(acceptorIndex, ShiftRtRange(warpedRtRange, -shift), ShiftMzRange(mzRange, -shift)));
        }

        return counts;
    }

    private static DoubleRange ShiftRtRange(DoubleRange range, double shift)
    {
        return new DoubleRange(range.Minimum + shift, range.Maximum + shift);
    }

    private static DoubleRange ShiftMzRange(DoubleRange range, double shift)
    {
        return new DoubleRange(range.Minimum + shift, range.Maximum + shift);
    }

    private static BinomialMatchScore BuildScore(int kTrue, int n, List<int> nullCounts)
    {
        double pNull = Math.Clamp(nullCounts.Average() / n, 0.0, 1.0);

        int boundedKTrue = Math.Clamp(kTrue, 0, n);
        double pValue = boundedKTrue <= 0
            ? 1.0
            : 1.0 - Binomial.CDF(pNull, n, boundedKTrue - 1);

        return new BinomialMatchScore(
            KTrue: kTrue,
            PNull: pNull,
            PValue: Math.Clamp(pValue, 0.0, 1.0));
    }
}
