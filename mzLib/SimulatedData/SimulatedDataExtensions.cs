using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Plotly.NET.CSharp;
using SpectralAveraging;
// don't add Plotly.NET as a using. Explicitly refer to it to avoid 
// namespace conflicts.


namespace SimulatedData
{
    public static class SimulatedDataExtensions
    {
        public static Plotly.NET.GenericChart.GenericChart Plot(this SimulatedData simData)
        {
            return Chart.Line<double, double, string>(x: simData.Xarray, y: simData.Yarray)
                .WithXAxisStyle<double, double, string>("m/z")
                .WithYAxisStyle<double, double, string>("Intensity"); 
        }

        public static bool TryMrsNoiseEstimation(this SimulatedData simData, double epsilon, out double noiseEstimate, 
            int maxIterations = 25)
        {
            bool mrsNoiseEstimationSuccess = MRSNoiseEstimator.MRSNoiseEstimation(simData.Yarray, 
                epsilon, out double tempNoiseEstimate, maxIterations);
            if (mrsNoiseEstimationSuccess)
            {
                noiseEstimate = tempNoiseEstimate;
            }
            else
            {
                noiseEstimate = double.NaN; 
            }

            return mrsNoiseEstimationSuccess;
        }

        public static double[] NaiveAverage(this List<SimulatedData> simDataList)
        {
            DoubleArray averagedYArray = new double[simDataList[0].Yarray.Length];
            foreach (SimulatedData simulatedData in simDataList)
            {
                // suppress CS8620 with exclamation point, because simulated data Yarray will never 
                // be null. 
                averagedYArray += simulatedData.Yarray!; 
            }

            return (averagedYArray / simDataList.Count).Array;
        }
        public static void NormalizeByTic(this SimulatedData simData)
        {
            double sumOfYArray = simData.Yarray.Sum();
            simData.ApplyElementwise(i => i / sumOfYArray, simData.Yarray);
        }
    }
    
}
