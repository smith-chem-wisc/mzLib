using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MzLibSpectralAveraging;

public static class OutlierRejection
{
    /// <summary>
    /// Perform rejection based on the SpectralAveragingOptions based as an argument.
    /// </summary>
    /// <remarks>This method is optimized by using Parallel.ForEach. This is the correct choice because we are iterating over many
    /// pixel stacks in the spectra, exceeding my general criteria for surpassing the additional overhead (roughly 1000 operations) that parallelization requires. </remarks>
    /// <param name="peakBin">stack which should be checked for outliers</param>
    /// <param name="options">which outlier rejection method to use</param>
    public static void RejectOutliers(PeakBin peakBin, SpectralAveragingOptions options)
    {
        switch (options.OutlierRejectionType)
        {
            case OutlierRejectionType.NoRejection:
                break;

            case OutlierRejectionType.MinMaxClipping:
                MinMaxClipping(peakBin);
                break;

            case OutlierRejectionType.PercentileClipping:
                PercentileClipping(peakBin, options.Percentile);
                break;

            case OutlierRejectionType.SigmaClipping:
                SigmaClipping(peakBin, options.MinSigmaValue, options.MaxSigmaValue);
                break;

            case OutlierRejectionType.WinsorizedSigmaClipping:
                WinsorizedSigmaClipping(peakBin, options.MinSigmaValue, options.MaxSigmaValue);
                break;

            case OutlierRejectionType.AveragedSigmaClipping:
                AveragedSigmaClipping(peakBin, options.MinSigmaValue, options.MaxSigmaValue);
                break;

            case OutlierRejectionType.BelowThresholdRejection:
                BelowThresholdRejection(peakBin);
                break;
        }
    }

    /// <summary>
    /// Perform rejection based on the SpectralAveragingOptions based as an argument.
    /// Overload Method for binned spectra
    /// </summary>
    /// <remarks>This method is optimized by using Parallel.ForEach. This is the correct choice because we are iterating over many
    /// pixel stacks in the spectra, exceeding my general criteria for surpassing the additional overhead (roughly 1000 operations) that parallelization requires. </remarks>
    /// <param name="binnedSpectra">spectra whose outlier should be rejected</param>
    /// <param name="options">which outlier rejection method to use</param>
    public static void RejectOutliers(BinnedSpectra binnedSpectra, SpectralAveragingOptions options)
    {
        Parallel.ForEach(binnedSpectra.PeakBins, pixelStack =>
        {
            RejectOutliers(pixelStack, options);
        });
    }

    /// <summary>
    /// Perform rejection based on the SpectralAveragingOptions based as an argument.
    /// Overload Method for binned spectra
    /// </summary>
    /// <remarks>This method is optimized by using Parallel.ForEach. This is the correct choice because we are iterating over many
    /// pixel stacks in the spectra, exceeding my general criteria for surpassing the additional overhead (roughly 1000 operations) that parallelization requires. </remarks>
    /// <param name="intensityArray">intensity values of spectra whose outliers should be rejected</param>
    /// <param name="options">which outlier rejection method to use</param>
    public static double[] RejectOutliers(double[] intensityArray, SpectralAveragingOptions options)
    {
        var xArray = new double[intensityArray.Length];
        var stack = new PeakBin(xArray, intensityArray);
        RejectOutliers(stack, options);
        return stack.UnrejectedIntensities.ToArray();
    }

    /// <summary>
    /// Reject the max and min of the set
    /// </summary>
    /// <param name="stack">array of mz values to evaluate</param>
    /// <returns>list of mz values with outliers rejected</returns>
    private static void MinMaxClipping(PeakBin stack)
    {
        if (stack.Intensities.Any())
        {
            int max = stack.Intensities.IndexOf(stack.Intensities.Max());
            int min = stack.Intensities.IndexOf(stack.Intensities.Min());
            stack.Reject(max);
            stack.Reject(min);
        }
    }

    /// <summary>
    /// Removes values that fall outside of the central value by the defined percentile exclusively
    /// </summary>
    /// <param name="peakBin">list of mz values to evaluate</param>
    /// <param name="percentile"></param>
    /// <returns>list of mz values with outliers rejected</returns>
    private static void PercentileClipping(PeakBin peakBin, double percentile)
    {
        if (!peakBin.Intensities.Any())
        {
            return; 
        }

        double trim = (1 - percentile) / 2;
        double highPercentile = 1 - trim;
        double median = BasicStatistics.CalculateMedian(peakBin.Intensities);
        double highCutoff = median * (1 + highPercentile);
        double lowCutoff = median * (1 - highPercentile);
        // round to 4-6 decimal places
        for (int i = 0; i < peakBin.Length; i++)
        {
            if (peakBin.Intensities[i] < highCutoff && peakBin.Intensities[i] > lowCutoff)
            {
                continue; 
            }
            peakBin.Reject(i);
        }
    }

    /// <summary>
    /// Iteratively removes values that fall outside of the central value by the defined StandardDeviation amount
    /// </summary>
    /// <param name="peakBin">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    private static void SigmaClipping(PeakBin peakBin, double sValueMin, double sValueMax)
    {
        if (!peakBin.Intensities.Any())
        {
            return;
        }

        int n;
        int iterationN = peakBin.Length; 
        do
        {
            double median = BasicStatistics.CalculateMedian(peakBin.UnrejectedIntensities);
            double standardDeviation = BasicStatistics.CalculateStandardDeviation(peakBin.UnrejectedIntensities);
            n = 0; 
            for (int i = 0; i < peakBin.Length; i++)
            {
                if (peakBin.IsIndexRejected(i)) continue; 
                if (SigmaClipping(peakBin.Intensities[i], median, standardDeviation, sValueMin, sValueMax))
                {
                    peakBin.Reject(i);
                    n++;
                }
            }
            iterationN -= n;

        } while (n > 0 && iterationN > 3);
    }

    /// <summary>
    /// Iteratively removes values that fall outside of a calculated deviation based upon the median of the values
    /// </summary>
    /// <param name="peakBin">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    private static void AveragedSigmaClipping(PeakBin peakBin, double sValueMin, double sValueMax)
    {
        if (!peakBin.Intensities.Any())
        {
            return;
        }
        double median = BasicStatistics.CalculateMedian(peakBin.UnrejectedIntensities);
        
        // calculate s
        double numerator = 0;
        double denominator = peakBin.Length - 1;
        int n = 0;

        for (int i = 0; i < peakBin.Length; i++)
        {
            numerator += (Math.Pow((peakBin.Intensities[i] - median), 2) / median);
        }

        double s = Math.Sqrt(numerator / denominator);

        int iterationUnrejectedLength = peakBin.NonRejectedLength;
        do
        {
            median = BasicStatistics.CalculateMedian(peakBin.UnrejectedIntensities);
            double sigma = s * Math.Sqrt(median);
            n = 0;

            for (int i = 0; i < peakBin.Length; i++)
            {
                if (peakBin.IsIndexRejected(i)) continue;
                if (SigmaClipping(peakBin.Intensities[i], median, sigma, sValueMin, sValueMax))
                {
                    peakBin.Reject(i);
                    n++;
                }
            }
            iterationUnrejectedLength -= n;

        } while (n > 0);
    }

    /// <summary>
    /// Iteratively replaces values that fall outside of the central value by the defined StandardDeviation amount with the values of the median * that standard deviation amount
    /// </summary>
    /// <param name="peakBin">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    private static void WinsorizedSigmaClipping(PeakBin peakBin, double sValueMin, double sValueMax)
    {
        if (!peakBin.Intensities.Any())
        {
            return;
        }

        int n = 0; 
        double iterationLimitforHuberLoop = 0.00005;
        double medianLeftBound;
        double medianRightBound;
        double stddev_previous;
        double stddev_current;
        do
        {
            double median = BasicStatistics.CalculateMedian(peakBin.UnrejectedIntensities);
            stddev_current = BasicStatistics.CalculateStandardDeviation(peakBin.UnrejectedIntensities);
            List<double> tempIntensityValues = peakBin.UnrejectedIntensities.ToList();

            do // calculates a new median and standard deviation based on the values to do sigma clipping with (Huber loop)
            {
                medianLeftBound = median - 1.5 * stddev_current;
                medianRightBound = median + 1.5 * stddev_current;
                Winsorize(tempIntensityValues, medianLeftBound, medianRightBound);
                
                median = BasicStatistics.CalculateMedian(tempIntensityValues);
                
                stddev_previous = stddev_current; 
                stddev_current = BasicStatistics.CalculateStandardDeviation(tempIntensityValues) * 1.134;

            } while (Math.Abs(stddev_current - stddev_previous) / stddev_previous > iterationLimitforHuberLoop);

            n = 0;
            for (int i = 0; i < peakBin.Length; i++)
            {
                if (peakBin.IsIndexRejected(i)) continue;
                if (SigmaClipping(peakBin.Intensities[i], median, 
                        stddev_current, sValueMin, sValueMax))
                {
                    peakBin.Reject(i);
                    n++;
                }
            }
        } while (n > 0 && peakBin.UnrejectedIntensities.Count() > 1); // break loop if nothing was rejected, or only one value remains
    }

    /// <summary>
    /// Sets the array of mz values to null if they have 20% or fewer values than the number of scans
    /// </summary>
    /// <param name="peakBin">array of mz values to evaluate</param>
    /// <param name="cutoffValue">percent in decimal form of where to cutoff </param>
    /// <returns></returns>
    private static void BelowThresholdRejection(PeakBin peakBin, double cutoffValue = 0.2)
    {
        if (!peakBin.Intensities.Any())
        {
            return;
        }
        
        if (peakBin.Intensities.Count(p => p != 0) <= peakBin.Length * cutoffValue)
        {
            peakBin.RejectAll();
        }
    }

    #region Helpers

    /// <summary>
    /// Helper method to mutate the array of doubles based upon the median value
    /// </summary>
    /// <param name="array">initial values to process</param>
    /// <param name="medianLeftBound">minimum the element in the dataset is allowed to be</param>
    /// <param name="medianRightBound">maxamum the element in the dataset is allowed to be</param>
    /// <returns></returns>
    private static void Winsorize(List<double> array, double medianLeftBound, double medianRightBound)
    {
        for (int i = 0; i < array.Count; i++)
        {
            array[i] = WinsorizeElement(array[i], medianLeftBound, medianRightBound); 
        }
    }

    private static double WinsorizeElement(double value, double medianLeftBound, double medianRightBound)
    {
        if (value < medianLeftBound)
        {
            return medianLeftBound; 
        }
        if(value > medianRightBound)
        {
            return medianRightBound; 
        }

        return value; 
    }

    /// <summary>
    /// Helper delegate method for sigma clipping
    /// </summary>
    /// <param name="value">the value in question</param>
    /// <param name="median">median of the dataset</param>
    /// <param name="standardDeviation">standard dev of the dataset</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    private static bool SigmaClipping(double value, double median, double standardDeviation, double sValueMin, double sValueMax)
    {
        if ((median - value) / standardDeviation > sValueMin)
        {
            return true;
        }
        else if ((value - median) / standardDeviation > sValueMax)
        {
            return true;
        }
        else
        {
            return false;
        }
    }


    #endregion
}