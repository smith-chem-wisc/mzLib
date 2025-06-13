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
        public XicLinearSpline(double splineRtInterval = 0.05)
        {
            SplineRtInterval = splineRtInterval;
        }

        public override (double, double)[] GetSplineXYData(double[] rtArray, double[] intensityArray, double start, double end)
        {
            if (rtArray.Length != intensityArray.Length)
            {
                throw new MzLibException("Input arrays must have the same length");
            }
            if (rtArray.Length < 5)
            {
                throw new MzLibException("Input arrays must contain at least 5 points");
            }
            var linearSpline = LinearSpline.InterpolateSorted(rtArray, intensityArray);
            var XYData = CalculateSpline(start, end, SplineRtInterval, linearSpline);

            return XYData;
        }
    }
}
