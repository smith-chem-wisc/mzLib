using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    /// <summary>
    /// A generic XIC class for all IIndexedPeak objects (mz peak, isotopic envelope, etc.) that can be traced across retention time.
    /// </summary>
    public class ExtractedIonChromatogram
    {
        public List<IIndexedPeak> Peaks { get; set; }

        public double ApexRT;
        public int ApexScanIndex;
        public double AveragedM => AverageM();
        public (double, double)[] XYData { get; set; }
        public double[] NormalizedPeakIntensities { get; set; }
        public double StartRT => Peaks.Min(p => p.RetentionTime);
        public double EndRT => Peaks.Max(p => p.RetentionTime);
        public int StartCycle => Peaks.Min(p => p.ZeroBasedScanIndex);
        public int EndCycle => Peaks.Max(p => p.ZeroBasedScanIndex);
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

        public ExtractedIonChromatogram(List<IIndexedPeak> peaks)
        {
            Peaks = peaks;
            ApexRT = Peaks.OrderByDescending(p => p.Intensity).First().RetentionTime;
            ApexScanIndex = Peaks.OrderByDescending(p => p.Intensity).First().ZeroBasedScanIndex;
        }

        public void SetNormalizedPeakIntensities()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            NormalizedPeakIntensities = Peaks.Select(p => p.Intensity / sumIntensity * 100).ToArray();
        }

    }
}
