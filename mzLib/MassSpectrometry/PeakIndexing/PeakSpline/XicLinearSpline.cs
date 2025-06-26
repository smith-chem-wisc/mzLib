using MathNet.Numerics.Interpolation;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class XicLinearSpline : XicSpline
    {
        public XicLinearSpline(double splineRtInterval = 0.05, int numberOfPeaksToAdd = 0, double gap = 1)
        {
            SplineRtInterval = splineRtInterval;
            NumberOfPeaksToAdd = numberOfPeaksToAdd;
            Gap = gap;
        }

        public override (double, double)[] GetXicSplineData(double[] rtArray, double[] intensityArray, double start = -1, double end = -1)
        {
            AddPeaks(rtArray, intensityArray, out double[] newRtArray, out double[] newIntensityArray);
            CheckArrays(newRtArray, newIntensityArray);
            var linearSpline = LinearSpline.InterpolateSorted(newRtArray, newIntensityArray);
            if (start == -1 && end == -1)
            {
                start = newRtArray.Min();
                end = newRtArray.Max();
            }
            var XYData = CalculateSpline(start, end, SplineRtInterval, linearSpline);

            return XYData;
        }
    }
}
