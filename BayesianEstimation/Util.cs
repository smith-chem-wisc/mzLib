using MathNet.Numerics.Statistics;
using System;

namespace BayesianEstimation
{
    public class Util
    {
        public static (double hdi_start, double hdi_end) GetHighestDensityInterval(double[] data)
        {
            Array.Sort(data);

            var ci_nbr_of_points = (int)Math.Floor(data.Length * 0.95);
            var min_width_ci = (data.Minimum(), data.Maximum()); // initialize interval

            for (var i = 0; i < data.Length - ci_nbr_of_points; i++)
            {
                var ci_width = data[i + ci_nbr_of_points] - data[i];

                if (ci_width < min_width_ci.Item2 - min_width_ci.Item1)
                {
                    min_width_ci = (data[i], data[i + ci_nbr_of_points]);
                }
            }

            return min_width_ci;
        }
    }
}
