using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using CsvHelper.Configuration.Attributes;
using MathNet.Numerics.Interpolation;
using MzLibUtil;

namespace MassSpectrometry
{
    public class XicCubicSpline : XicSpline
    {
        public XicCubicSpline(double splineRtInterval = 0.05, int numberOfPeaksToAdd = 0, double gap = 1, bool scanIndexBased = false) : base(splineRtInterval, numberOfPeaksToAdd, gap, scanIndexBased)
        {
        }

        public override (double, double)[] GetXicSplineData(float[] rtArray, float[] intensityArray, double start = -1, double end = -1)
        {
            AddPeaks(rtArray, intensityArray, out double[] newRtArray, out double[] newIntensityArray);
            CheckArrays(newRtArray, newIntensityArray);
            var cubicSpline = CubicSpline.InterpolateAkima(newRtArray, newIntensityArray);
            if (start == -1 && end == -1)
            {
                start = newRtArray.Min();
                end = newRtArray.Max();
            }
            var XYData = CalculateSpline(start, end, SplineRtInterval, cubicSpline);

            return XYData;
        }
    }
}
