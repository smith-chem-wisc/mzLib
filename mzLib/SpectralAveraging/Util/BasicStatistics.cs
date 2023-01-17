using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;

namespace SpectralAveraging;

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
}