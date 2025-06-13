using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using MathNet.Numerics.Interpolation;
using MzLibUtil;

namespace MassSpectrometry
{
    public class XicCubicSpline : XicSpline
    {
        public XicCubicSpline(double splineRtInterval = 0.05)
        {
            SplineRtInterval = splineRtInterval;
        }

        public override (double, double)[] GetSplineXYData(double[] rtArray, double[] intensityArray, double start, double end)
        {
            if (rtArray.Length != intensityArray.Length)
            {
                throw new MzLibException("Input arrays must have the same length.");
            }
            if (rtArray.Length < 5)
            {
                throw new MzLibException("Input arrays must contain at least 5 points.");
            }
            var cubicSpline = CubicSpline.InterpolateAkima(rtArray, intensityArray);
            var XYData = CalculateSpline(start, end, SplineRtInterval, cubicSpline);

            return XYData;
        }
    }
}
