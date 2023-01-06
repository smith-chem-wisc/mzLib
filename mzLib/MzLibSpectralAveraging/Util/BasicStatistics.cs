using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibSpectralAveraging
{
    public class BasicStatistics
    {
        /// <summary>
        /// Calculates the median of a list of doubles
        /// </summary>
        /// <param name="toCalc">initial list to calculate from</param>
        /// <returns>double representation of the median</returns>
        public static double CalculateMedian(IEnumerable<double> toCalc)
        {
            IEnumerable<double> sortedValues = toCalc.OrderByDescending(p => p).ToList();
            double median;
            int count = sortedValues.Count();
            if (count % 2 == 0)
                median = sortedValues.Skip(count / 2 - 1).Take(2).Average();
            else
                median = sortedValues.ElementAt(count / 2);
            return median;
        }

        public static double CalculateNonZeroMedian(IEnumerable<double> toCalc)
        {
            toCalc = toCalc.Where(p => p != 0).ToList();
            if (!toCalc.Any())
                return 0;

            return CalculateMedian(toCalc);
        }

        /// <summary>
        /// Calculates the standard deviation of a list of doubles
        /// </summary>
        /// <param name="toCalc">initial list to calculate from</param>
        /// <param name="average">passable value for the average</param>
        /// <returns>double representation of the standard deviation</returns>
        public static double CalculateStandardDeviation(IEnumerable<double> toCalc, double average = 0)
        {
            double deviation = 0;

            if (toCalc.Any())
            {
                List<double> calcList = toCalc.ToList();
                average = average == 0 ? calcList.Average() : average;
                double sum = calcList.Sum(x => Math.Pow(x - average, 2));
                deviation = Math.Sqrt(sum / (double)(calcList.Count() - 1));
            }
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
    }
}
