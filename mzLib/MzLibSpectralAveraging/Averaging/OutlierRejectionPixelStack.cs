using System;
using System.Collections.Generic;
using System.Linq;

namespace MzLibSpectralAveraging;

public static partial class OutlierRejection
{

    public static void RejectOutliers(PixelStack pixelStack, SpectralAveragingOptions options)
    {
        switch (options.RejectionType)
        {
            case RejectionType.NoRejection:
                break;

            case RejectionType.MinMaxClipping:
                MinMaxClipping(pixelStack);
                break;

            case RejectionType.PercentileClipping:
                PercentileClipping(pixelStack, options.Percentile);
                break;

            case RejectionType.SigmaClipping:
                SigmaClipping(pixelStack, options.MinSigmaValue, options.MaxSigmaValue);
                break;

            case RejectionType.WinsorizedSigmaClipping:
                WinsorizedSigmaClipping(pixelStack, options.MinSigmaValue, options.MaxSigmaValue);
                break;

            case RejectionType.AveragedSigmaClipping:
                AveragedSigmaClipping(pixelStack, options.MinSigmaValue, options.MaxSigmaValue);
                break;

            case RejectionType.BelowThresholdRejection:
                BelowThresholdRejection(pixelStack);
                break;
        }
    }
    /// <summary>
    /// Reject the max and min of the set
    /// </summary>
    /// <param name="initialValues">array of mz values to evaluate</param>
    /// <returns>list of mz values with outliers rejected</returns>
    public static void MinMaxClipping(PixelStack stack)
    {
        if (stack.Intensity.Any())
        {
            int max = stack.Intensity.IndexOf(stack.Intensity.Max());
            int min = stack.Intensity.IndexOf(stack.Intensity.Min());
            stack.Reject(max);
            stack.Reject(min);
        }

    }

    /// <summary>
    /// Removes values that fall outside of the central value by the defined percentile exclusively
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="percentile"></param>
    /// <returns>list of mz values with outliers rejected</returns>
    public static void PercentileClipping(PixelStack pixelStack, double percentile)
    {
        if (!pixelStack.Intensity.Any())
        {
            return; 
        }

        double trim = (1 - percentile) / 2;
        double highPercentile = 1 - trim;
        double median = BasicStatistics.CalculateMedian(pixelStack.Intensity);
        double highCutoff = median * (1 + highPercentile);
        double lowCutoff = median * (1 - highPercentile);
        // round to 4-6 decimal places
        for (int i = 0; i < pixelStack.Length; i++)
        {
            if (pixelStack.Intensity[i] < highCutoff && pixelStack.Intensity[i] > lowCutoff)
            {
                continue; 
            }
            pixelStack.Reject(i);
        }
    }

    /// <summary>
    /// Itteratively removes values that fall outside of the central value by the defined StandardDeviation amount
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    public static void SigmaClipping(PixelStack pixelStack, double sValueMin, double sValueMax)
    {
        if (!pixelStack.Intensity.Any())
        {
            return;
        }

        int n;
        int iterationN = pixelStack.Length; 
        do
        {
            double median = BasicStatistics.CalculateMedian(pixelStack.UnrejectedIntensities);
            double standardDeviation = BasicStatistics.CalculateStandardDeviation(pixelStack.UnrejectedIntensities);
            n = 0; 
            for (int i = 0; i < pixelStack.Length; i++)
            {
                if (pixelStack.IsIndexRejected(i)) continue; 
                if (SigmaClipping(pixelStack.Intensity[i], median, standardDeviation, sValueMin, sValueMax))
                {
                    pixelStack.Reject(i);
                    n++;
                }
            }
            iterationN -= n;

        } while (n > 0 && iterationN > 3);
    }

    /// <summary>
    /// Iteratively removes values that fall outside of a calculated deviation based upon the median of the values
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    public static void AveragedSigmaClipping(PixelStack pixelStack, double sValueMin, double sValueMax)
    {
        if (!pixelStack.Intensity.Any())
        {
            return;
        }
        double median = BasicStatistics.CalculateMedian(pixelStack.UnrejectedIntensities);
        
        // calculate s
        double numerator = 0;
        double denominator = pixelStack.Length - 1;
        int n = 0;

        for (int i = 0; i < pixelStack.Length; i++)
        {
            numerator += (Math.Pow((pixelStack.Intensity[i] - median), 2) / median);
        }

        double s = Math.Sqrt(numerator / denominator);

        int iterationUnrejectedLength = pixelStack.NonRejectedLength;
        do
        {
            median = BasicStatistics.CalculateMedian(pixelStack.UnrejectedIntensities);
            double sigma = s * Math.Sqrt(median);
            n = 0;

            for (int i = 0; i < pixelStack.Length; i++)
            {
                if (pixelStack.IsIndexRejected(i)) continue;
                if (SigmaClipping(pixelStack.Intensity[i], median, sigma, sValueMin, sValueMax))
                {
                    pixelStack.Reject(i);
                    n++;
                }
            }
            iterationUnrejectedLength -= n;

        } while (n > 0 && iterationUnrejectedLength - n > 3);
    }

    /// <summary>
    /// Itteratively replaces values that fall outside of the central value by the defined StandardDeviation amount with the values of the median * that standard deviation amount
    /// </summary>
    /// <param name="initialValues">list of mz values to evaluate</param>
    /// <param name="sValueMin">the lower limit of inclusion in sigma (standard deviation) units</param>
    /// <param name="sValueMax">the higher limit of inclusion in sigma (standard deviation) units</param>
    /// <returns></returns>
    public static void WinsorizedSigmaClipping(PixelStack pixelStack, double sValueMin, double sValueMax)
    {
        if (!pixelStack.Intensity.Any())
        {
            return;
        }

        int n = 0; 
        double iterationLimitforHuberLoop = 0.00005;
        double medianLeftBound;
        double medianRightBound;
        double stddev_previous;
        double stddev_current;
        int iterationN = pixelStack.Length;
        do
        {
            double median = BasicStatistics.CalculateNonZeroMedian(pixelStack.UnrejectedIntensities);
            stddev_current = BasicStatistics.CalculateNonZeroStandardDeviation(pixelStack.UnrejectedIntensities);
            List<double> tempIntensityValues = pixelStack.UnrejectedIntensities.ToList();

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
            for (int i = 0; i < pixelStack.Length; i++)
            {
                if (pixelStack.IsIndexRejected(i)) continue;
                if (SigmaClipping(pixelStack.Intensity[i], median, 
                        stddev_current, sValueMin, sValueMax))
                {
                    pixelStack.Reject(i);
                    n++;
                }
            }
            iterationN -= n;
        } while (n > 0 && iterationN > 3); // break loop if nothing was rejected, or only one value remains
    }
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
    /// Sets the array of mz values to null if they have 20% or fewer values than the number of scans
    /// </summary>
    /// <param name="initialValues">array of mz values to evaluate</param>
    /// <param name="cutoffValue">percent in decimal form of where to cutoff </param>
    /// <returns></returns>
    public static void BelowThresholdRejection(PixelStack pixelStack, double cutoffValue = 0.2)
    {
        if (!pixelStack.Intensity.Any())
        {
            return;
        }
        var cutoffVal = pixelStack.Intensity.Max() * cutoffValue;
        for (int i = 0; i < pixelStack.Length; i++)
        {
            if (pixelStack.Intensity[i] < cutoffValue)
            {
                pixelStack.Reject(i); 
            }
        }
    }
}