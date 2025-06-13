using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public abstract class XicSpline
    {
        public double SplineRtInterval { get; set; } //in minutes

        public abstract (double, double)[] GetSplineXYData(double[] rtArray, double[] intensityArray, double startRT, double endRT);

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

        public void SetXicSplineXYData(ExtractedIonChromatogram xic, double start, double end)
        {
            var peakRts = xic.Peaks.Select(p => p.RetentionTime).ToArray();
            var peakIntensities = xic.Peaks.Select(p => p.Intensity).ToArray();
            xic.XYData = GetSplineXYData(peakRts, peakIntensities, start, end);
        }
    }
}
