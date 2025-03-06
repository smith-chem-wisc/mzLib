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
        public abstract T Apex { get; }
        /// <summary>
        /// A list of data points that compose the ITraceable object
        /// This list must be ordered by the separation domain in ascending order! (e.g., retention time)
        /// </summary>
        public abstract List<T> ScanOrderedPoints { get; }

        /// <summary>
        /// Determines whether a peak should be cut based on the intensity of the surrounding time points.
        /// </summary>
        /// <typeparam name="T">The type of the time points, which must implement ISingleScanDatum.</typeparam>
        /// <param name="timePoints">The list of time points</param>
        /// <param name="apexTimepointIndex">The index of the apex (most intense, best, whatever) in the list of time points.</param>
        /// <param name="discriminationFactorToCutPeak">The discrimination factor to determine if the peak should be cut. Default is 0.6.</param>
        /// <returns>True if the peak should be cut, otherwise false.</returns>
        public static List<T> FindPeakSplits<T>(List<T> timePoints, int apexTimepointIndex, double discriminationFactorToCutPeak = 0.6) where T : ISingleScanDatum
        {
            List<T> peakBoundaries = new List<T>();
            foreach(int idx in FindPeakSplitIndices(timePoints, apexTimepointIndex, discriminationFactorToCutPeak))
            {
                peakBoundaries.Add(timePoints[idx]);
            }
            return peakBoundaries;
        }

        public static List<int> FindPeakSplitIndices<T>(List<T> timePoints, int apexTimepointIndex, double discriminationFactorToCutPeak = 0.6, int allowedMissedScans = 2) where T : ISingleScanDatum
        {
            List<int> peakSplitIndices = new List<int>();
            HashSet<int> zeroIndexedScanNumbers = timePoints.Select(p => p.ZeroBasedScanIndex).ToHashSet();

            // -1 checks the left side, +1 checks the right side
            int[] directions = { -1, 1 };
            foreach (int direction in directions)
            {
                if(apexTimepointIndex < 0)
                {
                    Console.WriteLine("xys");
                    break;
                }
                T valley = default(T);
                int indexOfValley = 0;
                int previousZeroBasedScanIndex = timePoints[apexTimepointIndex].ZeroBasedScanIndex;

                for (int i = apexTimepointIndex + direction; i < timePoints.Count && i >= 0; i += direction)
                {
                    ISingleScanDatum timepoint = timePoints[i];

                    //if (Math.Abs(previousZeroBasedScanIndex - timepoint.ZeroBasedScanIndex) > allowedMissedScans)
                    //{
                    //    peakSplitIndices.Add(i - direction); // Split at previous point
                    //    break;
                    //}
                    //else 
                    previousZeroBasedScanIndex = timepoint.ZeroBasedScanIndex;

                    // Valley envelope is the lowest intensity point that has been encountered thus far
                    if (EqualityComparer<T>.Default.Equals(valley, default(T)) || timepoint.Intensity < valley.Intensity)
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
                            peakSplitIndices.Add(indexOfValley);
                            break;
                        }
                    }
                }
            }

            return peakSplitIndices;
        }

        /// <summary>
        /// Recursively cuts ITraceable objects, removing all datapoints
        /// that occur before or after potential "valleys" surrounding the Apex.
        /// </summary>
        /// <param name="discriminationFactorToCutPeak"> The discrimination factor to determine if the peak should be cut. Default is 0.6. </param>
        public virtual void CutPeak(double discriminationFactorToCutPeak = 0.6)
        {
            var peakBoundaries = FindPeakSplits(ScanOrderedPoints, ScanOrderedPoints.IndexOf(Apex), discriminationFactorToCutPeak);
            CutPeak(peakBoundaries, Apex.RelativeSeparationValue);
        }

        /// <summary>
        /// Cuts a TraceablePeak objects, starting at the separationValueAtPeakCenter and removing all datapoints
        /// that occur before or after the peakBoundaries
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
