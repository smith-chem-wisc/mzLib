using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.PeakIndexing.Interfaces
{
    public class PeakTrace<T> where T : IIndexedPeak
    {
        public List<T> Peaks { get; set; }

        public double ApexRT => Peaks.OrderByDescending(p => p.Intensity).First().RetentionTime;
        public double ApexIntensity => Peaks.OrderByDescending(p => p.Intensity).First().Intensity;
        public int ApexCycle => Peaks.OrderByDescending(p => p.Intensity).First().ZeroBasedScanIndex;
        public int MsLevel { get; }
        public double AveragedM => AverageM();
        public (double, double)[] XYData { get; set; }

        public double StartRT => Peaks.Min(p => p.RetentionTime);
        public double EndRT => Peaks.Max(p => p.RetentionTime);
        public double AverageM()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            double averagedM = 0;
            foreach (var peak in Peaks)
            {
                double weight = peak.Intensity / sumIntensity;
                averagedM += weight * peak.M;
            }
            return averagedM;
        }

        public PeakTrace(List<T> peaks)
        {
            Peaks = peaks;
        }

        public void GetRawXYData()
        {
            XYData = new (double, double)[Peaks.Count];
            for (int i = 0; i < Peaks.Count; i++)
            {
                XYData[i] = (Peaks[i].ZeroBasedScanIndex, Peaks[i].Intensity);
            }
        }

        public void GetLinearSplineXYData(double splineRtInterval, Dictionary<int, double> rtIndexMap = null)
        {
            if (Peaks.Count < 5)
            {
                return;
            }
            var sortedPeaks = Peaks.OrderBy(p => p.RetentionTime).ToList();
            var rtArray = sortedPeaks.Select(p => p.RetentionTime).ToArray();
            var intensityArray = sortedPeaks.Select(p => p.Intensity).ToArray();
            var linearSpline = LinearSpline.InterpolateSorted(rtArray, intensityArray);
            XYData = CalculateSpline(StartRT, EndRT, splineRtInterval, linearSpline);
        }

        public void GetCubicSplineXYData(double splineRtInterval, Dictionary<int, double> rtIndexMap = null)
        {
            if (Peaks.Count < 5)
            {
                return;
            }
            var sortedPeaks = Peaks.OrderBy(p => p.RetentionTime).ToList();
            var rtArray = sortedPeaks.Select(p => p.RetentionTime).ToArray();
            var intensityArray = sortedPeaks.Select(p => p.Intensity).ToArray();
            var cubicSpline = CubicSpline.InterpolateAkima(rtArray, intensityArray);
            XYData = CalculateSpline(StartRT, EndRT, splineRtInterval, cubicSpline);
        }

        private (double, double)[] CalculateSpline(double startRT, double endRT, double splineRtInterval, IInterpolation spline)
        {
            int numPoints = (int)Math.Floor((endRT - startRT) / splineRtInterval) + 1;
            var xyData = new (double, double)[numPoints];
            for (int i = 0; i < numPoints; i++)
            {
                var rt = startRT + i * splineRtInterval;
                var intensity = spline.Interpolate(rt);
                xyData[i] = (rt, intensity);
            }
            return xyData;
        }

        //compare Umpire Spline
    }
}
