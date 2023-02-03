using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Plotly.NET.CSharp;
using SpectralAveraging; 

namespace Test.SimulatedData
{
    public static class ConvenienceExtensions
    {
        public static Plotly.NET.GenericChart.GenericChart Plot(this double[] yArray, double[] xArray)
        {
            return Chart.Line<double, double, string>(x: xArray, y: yArray)
                .WithXAxisStyle<double, double, string>("m/z")
                .WithYAxisStyle<double, double, string>("Intensity");
        }

        public static double CalculateEnr(this double[] array, double referenceScale,
            double? noiseEstimateArray, double noiseEstimateReference)
        {
            double scale = Math.Sqrt(BasicStatistics.BiweightMidvariance(array));
            double internalNoiseEstimate = 0d;
            if (!noiseEstimateArray.HasValue)
            {
                bool success = MRSNoiseEstimator.MRSNoiseEstimation(array, 0.0001, out double noiseEstimateMrs);
                internalNoiseEstimate = noiseEstimateMrs;
            }
            else
            {
                internalNoiseEstimate = noiseEstimateArray.Value;
            }

            double k = referenceScale / scale;
            return noiseEstimateReference / (k * internalNoiseEstimate);
        }
    }
}
