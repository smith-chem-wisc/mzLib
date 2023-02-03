using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;


// don't add Plotly.NET as a using. Explicitly refer to it to avoid 
// namespace conflicts.


namespace SimulatedData
{
    public static class SimulatedDataExtensions
    {
        public static double[] NaiveAverage(this List<SimulatedData> simDataList)
        {
            DoubleArray averagedYArray = new double[simDataList[0].Yarray.Length];
            foreach (SimulatedData simulatedData in simDataList)
            {
                // suppress CS8620 with exclamation point, because simulated data Yarray will never 
                // be null. 
                averagedYArray += simulatedData.Yarray!;
            }

            return (averagedYArray / (double)simDataList.Count).Array;
        }

        public static double[] NaiveAverage(this List<GaussianPeakSpectra> simDataList)
        {
            List<SimulatedData> list = simDataList.Cast<SimulatedData>().ToList();
            return NaiveAverage(list); 

        }

        public static double[] NaiveAverage(this List<SimulatedChargeStateEnvelope> simDataList)
        {
            return NaiveAverage(simDataList.Cast<SimulatedData>().ToList()); 
        }

        public static double[][] AverageWithRejection(this List<SimulatedData> dataList, 
            SpectralAveraging.SpectralAveragingParameters specAveragingParams)
        {
            double[][] xArrays = new double[dataList.Count][];
            double[][] yArrays = new double[dataList.Count][];
            // put the x and y arrays into jagged array format 
            for (int i = 0; i < dataList.Count; i++)
            {
                xArrays[i] = new double[dataList[i].Xarray.Length]; 
                yArrays[i] = new double[dataList[i].Yarray.Length]; 

                xArrays[i] = dataList[i].Xarray; 
                yArrays[i] = dataList[i].Yarray;
            }

            return SpectralAveraging.SpectraAveraging.AverageSpectra(xArrays, yArrays, specAveragingParams); 
        }

        public static double[][] AverageWithRejection(this List<GaussianPeakSpectra> dataList,
            SpectralAveraging.SpectralAveragingParameters spectralAveragingParameters)
        {
            return AverageWithRejection(dataList.Cast<SimulatedData>().ToList(), spectralAveragingParameters); 
        }

        public static double[][] AverageWithRejection(this List<SimulatedChargeStateEnvelope> cseList,
            SpectralAveraging.SpectralAveragingParameters spectralAveragingParameters)
        {
            return AverageWithRejection(cseList.Cast<SimulatedData>().ToList(), spectralAveragingParameters);
        }
    }
    
}
