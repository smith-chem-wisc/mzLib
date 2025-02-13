using Easy.Common.Extensions;
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Interfaces
{
    /// <summary>
    /// Represent an object that can be traced over a separation domain
    /// E.g., A chromatographic peak trace composed of multiple IsotopicEnvelopes
    /// </summary>
    public abstract class TraceablePeak<T> where T : ISingleScanDatum
    {
        /// <summary>
        /// The most intense point in the trace
        /// </summary>
        public T Apex { get; }
        /// <summary>
        /// A list of data points that compose the ITraceable object
        /// This list must be ordered by the separation domain in ascending order! (e.g., retention time)
        /// </summary>
        public List<T> ScanOrderedPoints { get; }

        /// <summary>
        /// Determines whether a peak should be cut based on the intensity of the surrounding time points.
        /// </summary>
        /// <typeparam name="T">The type of the time points, which must implement ISingleScanDatum.</typeparam>
        /// <param name="timePoints">The list of time points</param>
        /// <param name="indexOfPeakCenterInTimePoints">The index of the apex (most intense, best, whatever) in the list of time points.</param>
        /// <param name="valley">The valley point, which is the lowest intensity point encountered.</param>
        /// <param name="discriminationFactorToCutPeak">The discrimination factor to determine if the peak should be cut. Default is 0.6.</param>
        /// <returns>True if the peak should be cut, otherwise false.</returns>
        public static List<T> FindPeakBoundaries<T>(List<T> timePoints, int indexOfPeakCenterInTimePoints, double discriminationFactorToCutPeak = 0.6) where T : ISingleScanDatum
        {
            List<T> peakBoundaries = new List<T>();
            HashSet<int> zeroIndexedScanNumbers = timePoints.Select(p => p.ZeroBasedScanIndex).ToHashSet();

            // -1 checks the left side, +1 checks the right side
            int[] directions = { -1, 1 };
            foreach (int direction in directions)
            {
                T valley = default(T);
                int indexOfValley = 0;

                for (int i = indexOfPeakCenterInTimePoints + direction; i < timePoints.Count && i >= 0; i += direction)
                {
                    ISingleScanDatum timepoint = timePoints[i];

                    // Valley envelope is the lowest intensity point that has been encountered thus far
                    if (valley == null || timepoint.Intensity < valley.Intensity)
                    {
                        valley = (T)timepoint;
                        indexOfValley = i;
                    }

                    double discriminationFactor = (timepoint.Intensity - valley.Intensity) / timepoint.Intensity;

                    // If the time point is at least discriminationFactor times more intense than the valley
                    // We perform an additional check to see if the time point is more intense than the point next to the valley
                    if (discriminationFactor > discriminationFactorToCutPeak &&
                        (indexOfValley + direction < timePoints.Count && indexOfValley + direction >= 0))
                    {
                        ISingleScanDatum secondValleyTimepoint = timePoints[indexOfValley + direction];

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

        public List<T> FindPeakBoundaries(double separationValueAtPeakCenter, double discriminationFactorToCutPeak = 0.6)
        {
            if(ScanOrderedPoints.Count < 5) return null;
            T peakCenter = ScanOrderedPoints.MinBy(d => Math.Abs(d.RelativeSeparationValue -  separationValueAtPeakCenter));
            return FindPeakBoundaries(ScanOrderedPoints, ScanOrderedPoints.IndexOf(peakCenter), discriminationFactorToCutPeak);
        }

        /// <summary>
        /// Recursively cuts ITraceable objects, removing all datapoints
        /// that occur before or after potential "valleys" surrounding the "separationValueAtPeakCenter",
        /// which in the case of a Chromatographic peak will be the associated identification's
        /// MS2 retention time.
        /// </summary>
        /// <param name="separationValueAtPeakCenter"> Time representing the center of the peak </param>
        /// <param name="discriminationFactorToCutPeak"> The discrimination factor to determine if the peak should be cut. Default is 0.6. </param>
        public virtual void CutPeak(double separationValueAtPeakCenter, double discriminationFactorToCutPeak = 0.6)
        {
            var peakBoundaries = FindPeakBoundaries(separationValueAtPeakCenter, discriminationFactorToCutPeak);
            CutPeak(peakBoundaries, separationValueAtPeakCenter);
        }

        /// <summary>
        /// Recursively cuts ITraceable objects, removing all datapoints outside of the peak boundaries, inclusive
        /// </summary>
        /// <param name="peakBoundaries"></param>
        /// <param name="separationValueAtPeakCenter"></param>
        public void CutPeak(List<T> peakBoundaries, double separationValueAtPeakCenter)
        {
            if (peakBoundaries.IsNotNullOrEmpty())
            {
                foreach (var boundary in peakBoundaries)
                {
                    if (boundary.RelativeSeparationValue > separationValueAtPeakCenter)
                    {
                        ScanOrderedPoints.RemoveAll(d => d.RelativeSeparationValue >= boundary.RelativeSeparationValue);
                    }
                    else
                    {
                        ScanOrderedPoints.RemoveAll(d => d.RelativeSeparationValue <= boundary.RelativeSeparationValue);
                    }
                }
            }
        }


    }
}
