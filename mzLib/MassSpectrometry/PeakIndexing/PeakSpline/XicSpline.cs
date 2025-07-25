using MathNet.Numerics.Interpolation;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;

namespace MassSpectrometry
{
    /// <summary>
    /// Impute/smooth a given XIC based on raw data points.
    /// </summary>
    public abstract class XicSpline
    {
        public double SplineRtInterval { get; set; } //can be in RT or scan cycle
        public int NumberOfPeaksToAdd { get; set; } //number of peaks to add on each side of the chromatogram
        public double Gap { get; set; } // gap between added peaks

        /// <summary>
        /// Get the data points as a list of tuple (rt/scanIndex, intensity) after smoothing/interpolation.
        /// Requires a rt/scanIndex array and an intensity array as input. Start and end of rt/index as optionl parameters.
        /// </summary>
        public abstract (double, double)[] GetXicSplineData(float[] rtArray, float[] intensityArray, double start = -1, double end = -1);

        /// <summary>
        /// Helper method that takes a start and end retention time, a spline interval, and an IInterpolation object to calculate the spline data points.
        /// </summary>
        protected (double, double)[] CalculateSpline(double startRT, double endRT, double splineRtInterval, IInterpolation spline)
        {
            int numPoints = (int)Math.Floor((endRT - startRT) / splineRtInterval + 1e-8) + 1;
            var xyData = new (double, double)[numPoints];
            for (int i = 0; i < numPoints; i++)
            {
                var rt = startRT + i * splineRtInterval;
                var intensity = spline.Interpolate(rt);
                xyData[i] = (rt, intensity);
            }
            return xyData;
        }

        /// <summary>
        /// Set the XYData field of an ExtractedIonChromatogram object as the result of spline.
        /// </summary>
        public void SetXicSplineXYData(ExtractedIonChromatogram xic, bool cycle = false, double start = -1, double end = -1)
        {
            float[] peakRts = null;
            if (cycle)
            {
                peakRts = xic.Peaks.Select(p =>(float)p.ZeroBasedScanIndex).ToArray();
            }
            else
            {
                peakRts = xic.Peaks.Select(p => (float)p.RetentionTime).ToArray();
            }
            var peakIntensities = xic.Peaks.Select(p => (float)p.Intensity).ToArray();
            xic.XYData = GetXicSplineData(peakRts, peakIntensities, start, end);
        }

        /// <summary>
        /// Check if the input arrays meet the requirements of interpolation.
        /// </summary>
        protected void CheckArrays(double[] rtArray, double[] intensityArray)
        {
            if (rtArray.Length != intensityArray.Length)
            {
                throw new MzLibException("Input arrays must have the same length.");
            }
            if (rtArray.Length < 5)
            {
                throw new MzLibException("Input arrays must contain at least 5 points.");
            }
        }

        /// <summary>
        /// Add points with 0 intensity to the beginning and end of the XIC before interpolation.
        /// </summary>
        public void AddPeaks(float[] rtArray, float[] intensityArray, out double[] newRtArray, out double[] newIntensityArray)
        {
            if (NumberOfPeaksToAdd == 0)
            {
                newRtArray = rtArray.Select(x => (double)x).ToArray(); //in order to keep the same precision as in the original data
                newIntensityArray = intensityArray.Select(x => (double)x).ToArray();
                return;
            }
            newRtArray = new double[rtArray.Length + NumberOfPeaksToAdd * 2];
            newIntensityArray = new double[intensityArray.Length + NumberOfPeaksToAdd * 2];
            for (int i = 0; i < newRtArray.Length; i++)
            {
                if (i < NumberOfPeaksToAdd)
                {
                    newRtArray[i] = rtArray[0] - (NumberOfPeaksToAdd - i) * Gap;
                    newIntensityArray[i] = 0;
                }
                else if (i >= rtArray.Length + NumberOfPeaksToAdd)
                {
                    newRtArray[i] = newRtArray[i - 1] + Gap;
                    newIntensityArray[i] = 0;
                }
                else
                {
                    newRtArray[i] = rtArray[i - NumberOfPeaksToAdd];
                    newIntensityArray[i] = intensityArray[i - NumberOfPeaksToAdd];
                }
            }
        }
    }
}
