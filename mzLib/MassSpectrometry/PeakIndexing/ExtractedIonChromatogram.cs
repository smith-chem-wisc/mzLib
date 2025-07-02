using Easy.Common.Extensions;
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
        public double StartRT { get; set; }
        public double EndRT { get; set; }
        public int StartScanIndex { get; set; }
        public int EndScanIndex { get; set; }
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
            StartRT = Peaks.Min(p => p.RetentionTime);
            EndRT = Peaks.Max(p => p.RetentionTime);
            ApexScanIndex = Peaks.OrderByDescending(p => p.Intensity).First().ZeroBasedScanIndex;
            StartScanIndex = Peaks.Min(p => p.ZeroBasedScanIndex);
            EndScanIndex = Peaks.Max(p => p.ZeroBasedScanIndex);
        }

        public void SetNormalizedPeakIntensities()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            NormalizedPeakIntensities = Peaks.Select(p => p.Intensity / sumIntensity * 100).ToArray();
        }

        /// <summary>
        /// Determines whether a peak should be cut based on the intensity of the surrounding time points.
        /// </summary>
        /// <typeparam name="T">The type of the time points, which must implement ISingleScanDatum.</typeparam>
        /// <param name="timePoints">The list of time points</param>
        /// <param name="apexTimepointIndex">The index of the apex (most intense, best, whatever) in the list of time points.</param>
        /// <param name="discriminationFactorToCutPeak">The discrimination factor to determine if the peak should be cut. Default is 0.6.</param>
        /// <returns>True if the peak should be cut, otherwise false.</returns>
        public static List<IIndexedPeak> FindPeakBoundaries(List<IIndexedPeak> timePoints, int apexTimepointIndex, double discriminationFactorToCutPeak = 0.6)
        {
            List<IIndexedPeak> peakBoundaries = new List<IIndexedPeak>();
            HashSet<int> zeroIndexedScanNumbers = timePoints.Select(p => p.ZeroBasedScanIndex).ToHashSet();

            // -1 checks the left side, +1 checks the right side
            int[] directions = { -1, 1 };
            foreach (int direction in directions)
            {
                IIndexedPeak valley = null;
                int indexOfValley = 0;

                for (int i = apexTimepointIndex + direction; i < timePoints.Count && i >= 0; i += direction)
                {
                    IIndexedPeak timepoint = timePoints[i];

                    // Valley envelope is the lowest intensity point that has been encountered thus far
                    if (valley == null || timepoint.Intensity < valley.Intensity)
                    {
                        valley = timepoint;
                        indexOfValley = i;
                    }

                    double discriminationFactor = (timepoint.Intensity - valley.Intensity) / timepoint.Intensity;

                    // If the time point is at least discriminationFactor times more intense than the valley
                    // We perform an additional check to see if the time point is more intense than the point next to the valley
                    if (discriminationFactor > discriminationFactorToCutPeak &&
                        (indexOfValley + direction < timePoints.Count && indexOfValley + direction >= 0))
                    {
                        IIndexedPeak secondValleyTimepoint = timePoints[indexOfValley + direction];

                        discriminationFactor =
                            (timepoint.Intensity - secondValleyTimepoint.Intensity) / timepoint.Intensity;

                        // If the current timepoint is more intense than the second valley, we cut the peak
                        // If the scan following the valley isn't in the timePointsForApexZ list (i.e., no isotopic envelope is observed in the scan immediately after the valley), we also cut the peak
                        if (discriminationFactor > discriminationFactorToCutPeak || !zeroIndexedScanNumbers.Contains(valley.ZeroBasedScanIndex + direction))
                        {
                            peakBoundaries.Add(valley);
                            break;
                        }
                    }
                }
            }
            return peakBoundaries;
        }

        /// <summary>
        /// Cuts ExtractedIonChromatogram peak list, starting at the ApexRT and removing all IIndexedPeak
        /// that occur before or after the peakBoundaries
        /// </summary>
        /// <param name="peakBoundaries"></param>
        /// <param name="separationValueAtPeakCenter"></param>
        public void RemovePeaks(List<IIndexedPeak> peakBoundaries, double separationValueAtPeakCenter)
        {
            if (peakBoundaries.IsNotNullOrEmpty())
            {
                foreach (var boundary in peakBoundaries)
                {
                    if (boundary.RetentionTime > separationValueAtPeakCenter)
                    {
                        Peaks.RemoveAll(d => d.RetentionTime >= boundary.RetentionTime);
                    }
                    else
                    {
                        Peaks.RemoveAll(d => d.RetentionTime <= boundary.RetentionTime);
                    }
                }
            }
        }

        public void CutPeak(double discriminationFactorToCutPeak = 0.6)
        {
            var apexPeak = Peaks.FirstOrDefault(p => p.ZeroBasedScanIndex == ApexScanIndex);
            var peakBoundaries = FindPeakBoundaries(Peaks, Peaks.IndexOf(apexPeak), discriminationFactorToCutPeak);
            CutPeak(peakBoundaries, apexPeak.RetentionTime);
        }

        public List<IIndexedPeak> CutPeak(List<IIndexedPeak> peaks, double discriminationFactorToCutPeak = 0.6)
        {
            var apexPeak = peaks.OrderByDescending(p => p.Intensity).FirstOrDefault();
            var sortedPeaks = peaks.OrderBy(p => p.RetentionTime).ToList();
            var peakBoundaries = FindPeakBoundaries(sortedPeaks, sortedPeaks.IndexOf(apexPeak), discriminationFactorToCutPeak);
            RemovePeaks(peakBoundaries, apexPeak.RetentionTime);
        }
    }
}
