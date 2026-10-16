using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using MathNet.Numerics.Statistics;

namespace SpectralAveraging;

/// <summary>
///     Reject outliers from a given set of one dimensional data
///     Methods are adapted from https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#description_003
/// </summary>
public static class OutlierRejection
{
    /// <summary>
    ///     Calls the specific rejection function based upon the current static field RejectionType
    ///     Can be used for any array of values
    /// </summary>
    /// <param name="valuesToCheckForRejection">list of mz values to evaluate</param>
    /// <param name="parameters">how to reject outliers</param>
    /// <returns>double array of unrejected values</returns>
    public static double[] RejectOutliers(double[] valuesToCheckForRejection, SpectralAveragingParameters parameters)
    {
        var trimmedMzValues = valuesToCheckForRejection;
        switch (parameters.OutlierRejectionType)
        {
            case OutlierRejectionType.NoRejection:
                break;

            case OutlierRejectionType.MinMaxClipping:
                trimmedMzValues = MinMaxClipping(valuesToCheckForRejection);
                break;

            case OutlierRejectionType.PercentileClipping:
                trimmedMzValues = PercentileClipping(valuesToCheckForRejection, parameters.Percentile);
                break;

            case OutlierRejectionType.SigmaClipping:
                trimmedMzValues = SigmaClipping(valuesToCheckForRejection, parameters.MinSigmaValue,
                    parameters.MaxSigmaValue);
                break;

            case OutlierRejectionType.WinsorizedSigmaClipping:
                trimmedMzValues = WinsorizedSigmaClipping(valuesToCheckForRejection, parameters.MinSigmaValue,
                    parameters.MaxSigmaValue);
                break;

            case OutlierRejectionType.AveragedSigmaClipping:
                trimmedMzValues = AveragedSigmaClipping(valuesToCheckForRejection, parameters.MinSigmaValue,
                    parameters.MaxSigmaValue);
                break;

            case OutlierRejectionType.BelowThresholdRejection:
                trimmedMzValues = BelowThresholdRejection(valuesToCheckForRejection);
                break;
        }

        return trimmedMzValues;
    }

    /// <summary>
    ///     Overload for internal spectral averaging BinnedPeak structure
    /// </summary>
    /// <param name="peaks">list of peaks to reject outliers</param>
    /// <param name="parameters">how to reject outliers</param>
    /// <returns></returns>
    internal static List<BinnedPeak> RejectOutliers(List<BinnedPeak> peaks, SpectralAveragingParameters parameters)
    {
        int kept = RejectOutliers(CollectionsMarshal.AsSpan(peaks), parameters, new OutlierRejectionBuffers());
        peaks.RemoveRange(kept, peaks.Count - kept);
        return peaks;
    }

    /// <summary>
    ///     Rejects outliers among the intensities of one bin's peaks, without allocating. The peaks that are kept are
    ///     moved, in their original order, to the front of <paramref name="peaks"/>, and their number is returned.
    /// </summary>
    /// <remarks>
    ///     A peak is kept when its intensity is one of the values the rejection method leaves, compared with
    ///     <see cref="IEquatable{T}"/> as <c>Contains</c> on the method's result array compared them (so NaN matches NaN
    ///     and 0 matches -0) - the same rule, applied in the same order, as removing the other peaks from a list one at a
    ///     time. Averaging calls this once per m/z bin, often millions of times per file, which is why it keeps its
    ///     working arrays in <paramref name="buffers"/> instead of building lists.
    /// </remarks>
    internal static int RejectOutliers(Span<BinnedPeak> peaks, SpectralAveragingParameters parameters,
        OutlierRejectionBuffers buffers)
    {
        int count = peaks.Length;
        double[] values = buffers.Values(count);
        for (int i = 0; i < count; i++)
            values[i] = peaks[i].Intensity;

        ReadOnlySpan<double> unrejected = parameters.OutlierRejectionType == OutlierRejectionType.SigmaClipping
            ? values.AsSpan(0, SigmaClippingInPlace(values, count, parameters.MinSigmaValue, parameters.MaxSigmaValue, buffers))
            : RejectOutliers(values.AsSpan(0, count).ToArray(), parameters);

        int kept = 0;
        for (int i = 0; i < count; i++)
        {
            if (unrejected.Contains(peaks[i].Intensity))
                peaks[kept++] = peaks[i];
        }

        return kept;
    }

    #region Rejection Functions

    /// <summary>
    ///     Reject the max and min of the set
    /// </summary>
    /// <param name="initialValues">array of mz values to evaluate</param>
    /// <returns>list of mz values with outliers rejected</returns>
    private static double[] MinMaxClipping(double[] initialValues)
    {
        var max = initialValues.Max();
        var min = initialValues.Min();

        return initialValues.Where(p => p < max && p > min).ToArray();
    }

    /// <summary>
    ///     Removes values that fall outside of the central value by the defined percentile exclusively
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="percentile"></param>
    /// <returns>list of mz values with outliers rejected</returns>
    private static double[] PercentileClipping(double[] initialValues, double percentile)
    {
        var trim = (1 - percentile) / 2;
        var highPercentile = 1 - trim;
        var median = initialValues.Median();
        var highCutoff = median * (1 + highPercentile);
        var lowCutoff = median * (1 - highPercentile);
        // round to 4-6 decimal places
        return initialValues.Where(p => highCutoff > p && p > lowCutoff).ToArray();
    }

    /// <summary>
    ///     Iteratively removes values that fall outside of the central value by the defined StandardDeviation amount
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    private static double[] SigmaClipping(double[] initialValues, double sValueMin, double sValueMax)
    {
        var values = (double[])initialValues.Clone();
        int count = SigmaClippingInPlace(values, values.Length, sValueMin, sValueMax, new OutlierRejectionBuffers());
        return values.AsSpan(0, count).ToArray();
    }

    /// <summary>
    ///     Sigma clipping of the first <paramref name="count"/> entries of <paramref name="values"/>, in place: the values
    ///     that are kept are moved, in their original order, to the front, and their number is returned.
    /// </summary>
    /// <remarks>
    ///     Each pass computes the median and standard deviation of what is left and then drops every value outside the
    ///     limits, until a pass drops nothing - the passes, comparisons and arithmetic of the list-based version this
    ///     replaced, which removed values from a <see cref="List{T}"/> one at a time. The median is MathNet's
    ///     <see cref="ArrayStatistics.MedianInplace(double[])"/> on a copy, which is what <c>Median()</c> on a list does; the
    ///     standard deviation is <see cref="StreamingStandardDeviation"/>, the computation <c>StandardDeviation()</c> on a
    ///     list runs.
    /// </remarks>
    internal static int SigmaClippingInPlace(double[] values, int count, double sValueMin, double sValueMax,
        OutlierRejectionBuffers buffers)
    {
        int removed;
        do
        {
            double[] copy = buffers.ExactLength(count);
            Array.Copy(values, copy, count);
            // the standard deviation reads the values in order, so take it before MedianInplace reorders the copy
            double standardDeviation = StreamingStandardDeviation(copy);
            double median = ArrayStatistics.MedianInplace(copy);

            int kept = 0;
            for (int i = 0; i < count; i++)
            {
                double value = values[i];
                if (!SigmaClipping(value, median, standardDeviation, sValueMin, sValueMax))
                    values[kept++] = value;
            }

            removed = count - kept;
            count = kept;
        } while (removed > 0);

        return count;
    }

    /// <summary>
    ///     The sample standard deviation, computed exactly as MathNet.Numerics 5.0's
    ///     <see cref="StreamingStatistics.StandardDeviation(IEnumerable{double})"/> computes it - the same recurrence,
    ///     in the same order - which is what <c>StandardDeviation()</c> runs for a <see cref="List{T}"/>. Reading a span
    ///     directly avoids an enumerator and an interface call per value. NaN for fewer than two values.
    /// </summary>
    internal static double StreamingStandardDeviation(ReadOnlySpan<double> samples)
    {
        double variance = 0;
        double t = 0;
        ulong j = 0;

        if (samples.Length > 0)
        {
            j++;
            t = samples[0];
        }

        for (int i = 1; i < samples.Length; i++)
        {
            j++;
            double xi = samples[i];
            t += xi;
            double diff = (j * xi) - t;
            variance += (diff * diff) / (j * (j - 1));
        }

        return Math.Sqrt(j > 1 ? variance / (j - 1) : double.NaN);
    }

    /// <summary>
    ///     Iteratively replaces values that fall outside of the central value by the defined StandardDeviation amount with the
    ///     values of the median * that standard deviation amount
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <remarks>Documentation for winsorized sigma clipping can be found here https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#description_003</remarks>
    /// <returns></returns>
    private static double[] WinsorizedSigmaClipping(double[] initialValues, double sValueMin, double sValueMax)
    {
        var values = initialValues.ToList();
        int n;
        const double iterationLimitForHuberLoop = 0.00005;
        do
        {
            if (!values.Any())
                break;
            var median = values.Median();
            var standardDeviation = values.Where(p => p != 0).StandardDeviation();
            double[] toProcess = new double[values.Count];
            values.CopyTo(toProcess, 0);
            double winsorizedStandardDeviation;

            do // calculates a new median and standard deviation based on the values to do sigma clipping with (Huber loop)
            {
                var medianLeftBound = median - 1.5 * standardDeviation; // magic number 1.5: optimized parameter per the documentation this method was adapted from, see link at top of class
                var medianRightBound = median + 1.5 * standardDeviation;
                toProcess.Winsorize(medianLeftBound, medianRightBound);
                median = toProcess.Median();
                winsorizedStandardDeviation = standardDeviation;
                standardDeviation = toProcess.StandardDeviation() * 1.134; // magic number 1.134: optimized parameter per the documentation this method was adapted from, see link at top of class
            } while (Math.Abs(standardDeviation - winsorizedStandardDeviation) / winsorizedStandardDeviation >
                     iterationLimitForHuberLoop);

            n = 0;
            for (var i = 0; i < values.Count; i++)
                if (SigmaClipping(values[i], median, standardDeviation, sValueMin, sValueMax))
                {
                    values.RemoveAt(i);
                    n++;
                    i--;
                }
        } while (n > 0 && values.Count > 1); // break loop if nothing was rejected, or only one value remains

        return values.ToArray();
    }

    /// <summary>
    ///     Iteratively removes values that fall outside of a calculated deviation based upon the median of the values
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <remarks>Documentation for averaged sigma clipping can be found here https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#description_003</remarks>
    /// <returns></returns>
    private static double[] AveragedSigmaClipping(double[] initialValues, double sValueMin, double sValueMax)
    {
        var values = initialValues.ToList();
        var median = initialValues.Median();

        double sum = 0;
        for (int i = 0; i < values.Count; i++)
        {
            sum += Math.Pow((values[i] - median), 2) / median;
        }

        double s = Math.Sqrt(sum / (double)(values.Count - 1));

        int n;
        do
        {
            median = values.Where(p => p != 0).Median();
            double sigma = s * Math.Sqrt(median);

            n = 0;
            for (var i = 0; i < values.Count; i++)
            {
                if (SigmaClipping(values[i], median, sigma, sValueMin, sValueMax))
                {
                    values.RemoveAt(i);
                    n++;
                    i--;
                }
            }
        } while (n > 0);

        return values.ToArray();
    }

    /// <summary>
    ///     Sets the array of mz values to null if they have 20% or fewer values than the number of scans
    /// </summary>
    /// <param name="initialValues">array of mz values to evaluate</param>
    /// <param name="cutoffValue">percent in decimal form of where to cutoff </param>
    /// <returns></returns>
    private static double[] BelowThresholdRejection(double[] initialValues, double cutoffValue = 0.2)
    {
        var scanCount = initialValues.Length;
        if (initialValues.Count(p => p != 0) <= scanCount * cutoffValue) initialValues = new double[scanCount];
        return initialValues;
    }

    #endregion

    #region Private Helpers

    /// <summary>
    ///     Helper method to mutate the array of doubles based upon the median value
    /// </summary>
    /// <param name="initialValues">initial values to process</param>
    /// <param name="medianLeftBound">minimum the element in the data set is allowed to be</param>
    /// <param name="medianRightBound">maximum the element in the data set is allowed to be</param>
    /// <returns></returns>
    private static void Winsorize(this double[] initialValues, double medianLeftBound, double medianRightBound)
    {
        for (var i = 0; i < initialValues.Length; i++)
        {
            if (initialValues[i] < medianLeftBound)
            {
                initialValues[i] = medianLeftBound;
            }

            if (initialValues[i] > medianRightBound)
            {
                initialValues[i] = medianRightBound;
            }
        }
    }

    /// <summary>
    ///     Helper delegate method for sigma clipping
    /// </summary>
    /// <param name="value">the value in question</param>
    /// <param name="median">median of the data set</param>
    /// <param name="standardDeviation">standard dev of the data set</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    private static bool SigmaClipping(double value, double median, double standardDeviation, double sValueMin,
        double sValueMax)
    {
        if ((median - value) / standardDeviation > sValueMin)
            return true;
        if ((value - median) / standardDeviation > sValueMax)
            return true;
        return false;
    }

    #endregion
}