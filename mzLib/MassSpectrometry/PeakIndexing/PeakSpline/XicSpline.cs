using MathNet.Numerics.Interpolation;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public abstract class XicSpline
    {
        public double SplineRtInterval { get; set; } //can be in RT or scan cycle
        public int NumberOfPeaksToAdd { get; set; } //number of peaks to add on each side of the chromatogram
        public double Gap { get; set; } // gap between added peaks

        public abstract (double, double)[] GetXicSplineData(double[] rtArray, double[] intensityArray, double start = -1, double end = -1);

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

        public void SetXicSplineXYData(ExtractedIonChromatogram xic, bool cycle = false, double start = -1, double end = -1)
        {
            double[] peakRts = null;
            if (cycle)
            {
                peakRts = xic.Peaks.Select(p => (double)p.ZeroBasedScanIndex).ToArray();
            }
            else
            {
                peakRts = xic.Peaks.Select(p => p.RetentionTime).ToArray();
            }
            var peakIntensities = xic.Peaks.Select(p => p.Intensity).ToArray();
            xic.XYData = GetXicSplineData(peakRts, peakIntensities, start, end);
        }

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

        protected void AddPeaks(double[] rtArray, double[] intensityArray, out double[] newRtArray, out double[] newIntensityArray)
        {
            if (NumberOfPeaksToAdd == 0)
            {
                newRtArray = rtArray;
                newIntensityArray = intensityArray;
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
                else if (i >= rtArray.Length + NumberOfPeaksToAdd - 1)
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
