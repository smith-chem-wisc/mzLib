using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;

namespace SpectralAveraging
{
    public class BasicStatistics
    {
        public static double CalculateNonZeroMedian(IEnumerable<double> toCalc)
        {
            toCalc = toCalc.Where(p => p != 0).ToList();
            if (!toCalc.Any())
                return 0;

            return toCalc.Median();
        }

        /// <summary>
        ///     Calculates the standard deviation of a list of doubles
        /// </summary>
        /// <param name="toCalc">initial list to calculate from</param>
        /// <param name="average">passable value for the average</param>
        /// <returns>double representation of the standard deviation</returns>
        public static double CalculateStandardDeviation(IEnumerable<double> toCalc, double average = 0)
        {
            double deviation = 0;

            var calcList = toCalc.ToList();
            if (!calcList.Any()) return deviation;
            average = average == 0 ? calcList.Average() : average;
            var sum = calcList.Sum(x => Math.Pow(x - average, 2));
            deviation = Math.Sqrt(sum / (calcList.Count - 1));

            return deviation;
        }

        public static double CalculateStandardDeviation(double[] toCalc, double average = 0)
        {
            return CalculateStandardDeviation((IEnumerable<double>)toCalc, average);
        }

        public static double CalculateNonZeroStandardDeviation(IEnumerable<double> toCalc, double average = 0)
        {
            toCalc = toCalc.Where(p => p != 0).ToList();
            if (!toCalc.Any())
                return 0;
            return CalculateStandardDeviation(toCalc, average);
        }
        /// <summary>
        /// Calculates the median absolute deviation from the median, which is then used in
        /// calculating the biweight midvariance. Original algorithm found here:
        /// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_26__
        /// </summary>
        /// <param name="array">Array of values to be calculated.</param>
        /// <returns>The median absolute deviation from median.</returns>
        public static double MedianAbsoluteDeviationFromMedian(double[] array)
        {
            double arrayMedian = array.Median();
            double[] results = new double[array.Length];
            for (int j = 0; j < array.Length; j++)
            {
                results[j] = Math.Abs(array[j] - arrayMedian);
            }

            return results.Median();
        }
        /// <summary>
        /// Calcultes the biweight midvariance for an array. Algorithm orignally found here:
        /// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_27__
        /// </summary>
        /// <param name="array">Array of doubles.</param>
        /// <returns>The biweight midvariance.</returns>
        public static double BiweightMidvariance(double[] array)
        {
            double[] y_i = new double[array.Length];
            double[] a_i = new double[array.Length];
            double MAD_X = MedianAbsoluteDeviationFromMedian(array);
            double median = array.Median();
            for (int i = 0; i < y_i.Length; i++)
            {
                y_i[i] = (array[i] - median) / (9d * MAD_X);
                if (y_i[i] < 1d)
                {
                    a_i[i] = 1d;
                }
                else
                {
                    a_i[i] = 0;
                }
            }

            // biweight midvariance calculation

            double denomSum = 0;
            double numeratorSum = 0;
            for (int i = 0; i < y_i.Length; i++)
            {
                numeratorSum += a_i[i] * Math.Pow(array[i] - median, 2) * Math.Pow(1 - y_i[i] * y_i[i], 4);
                denomSum += a_i[i] * (1 - 5 * y_i[i] * y_i[i]) * (1 - y_i[i] * y_i[i]);
            }

            return (double)y_i.Length * numeratorSum / Math.Pow(Math.Abs(denomSum), 2);
        }
    }
}