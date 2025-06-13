using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class ExtendedCycleCubicSpline
    {
        public double SplineRtInterval { get; set; } // in minutes
        public int NumberOfPeaksToAdd { get; set; } 
        public double Gap { get; set;  } // in minutes, gap between added peaks

        public ExtendedCycleCubicSpline(double splineRtInterval, int numberOfPeaksToAdd, double gap)
        {
            SplineRtInterval = splineRtInterval;
            NumberOfPeaksToAdd = numberOfPeaksToAdd;
            Gap = gap;
        }

        public (double, double)[] GetSplineXYData(double[] rtArray, double[] intensityArray, double startRT, double endRT)
        {
            AddPeaks(rtArray, intensityArray, out double[] newRtArray, out double[] newIntensityArray);
            var cubicSpline = new XicCubicSpline(SplineRtInterval);
            return cubicSpline.GetSplineXYData(newRtArray, newIntensityArray, startRT, endRT);
        }

        private void AddPeaks(double[] rtArray, double[] intensityArray, out double[] newRtArray, out double[] newIntensityArray)
        {
            newRtArray = new double[rtArray.Length + NumberOfPeaksToAdd * 2];
            newIntensityArray = new double[intensityArray.Length + NumberOfPeaksToAdd * 2];
            for(int i = 0; i < newRtArray.Length; i++)
            {
                if (i < NumberOfPeaksToAdd)
                {
                    newRtArray[i] = rtArray[0] - (NumberOfPeaksToAdd - i) * Gap;
                    newIntensityArray[i] = 0;
                }
                else if (i >= rtArray.Length + NumberOfPeaksToAdd)
                {
                    newRtArray[i] = rtArray[rtArray.Length - 1] + (i - rtArray.Length - NumberOfPeaksToAdd + 1) * Gap;
                    newIntensityArray[i] = 0;
                }
                else
                {
                    newRtArray[i] = rtArray[i];
                    newIntensityArray[i] = intensityArray[i];
                }
            }
        }   
    }
}
