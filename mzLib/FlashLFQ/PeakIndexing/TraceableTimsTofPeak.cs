using Easy.Common.Extensions;
using FlashLFQ.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.PeakIndexing
{
    internal class TraceableTimsTofPeak : TraceablePeak<IonMobilityPeak>
    {
        /// <summary>
        /// List 
        /// </summary>
        public List<IonMobilityPeak> IonMobilityPeaks { get; init; }
        public override List<IonMobilityPeak> ScanOrderedPoints => IonMobilityPeaks;
        public override IonMobilityPeak Apex => _apex;
        private IonMobilityPeak _apex;
        public double RetentionTime { get; init; }
        public int ZeroBasedMs1FrameIndex { get; init; }
        //public int SummedIntensity { get; private set; }
        public double Mz => _apex.Mz;

        public TraceableTimsTofPeak(int zeroBasedMs1FrameIndex, double retentionTime, IonMobilityPeak ionMobilityPeak) : this(zeroBasedMs1FrameIndex, retentionTime)
        {
            AddIonMobilityPeak(ionMobilityPeak);
        }

        public TraceableTimsTofPeak(int zeroBasedMs1FrameIndex, double retentionTime)
        {
            IonMobilityPeaks = new List<IonMobilityPeak>(16);
            RetentionTime = retentionTime;
            ZeroBasedMs1FrameIndex = zeroBasedMs1FrameIndex;
            _apex = default(IonMobilityPeak);
        }

        /// <summary>
        /// Very strong assumption that the peaks are added in order of increasing scan number!!!
        /// </summary>
        /// <param name="ionMobilityPeak"></param>
        /// <exception cref="ArgumentException"></exception>
        internal void AddIonMobilityPeak(IonMobilityPeak ionMobilityPeak)
        {
            //if(ionMobilityPeak.OneBasedTimsScanNumber <= IonMobilityPeaks[^1].OneBasedTimsScanNumber)
            //{
            //    throw new ArgumentException("Ion Mobility Peaks must be added sequentially by scan number when building a TraceablePeak");
            //}
            IonMobilityPeaks.Add(ionMobilityPeak);
            if (_apex.Equals(default(IonMobilityPeak)) || ionMobilityPeak.IntegerIntensity > Apex.IntegerIntensity)
            {
                _apex = ionMobilityPeak;
            }
        }

        public List<IndexedMassSpectralPeak> GetIndexedPeaks()
        {
            if (ScanOrderedPoints.Count < 3)
            {
                return null;
            }

            //List<int> allPeakBoundaries = GetPeakBoundariesRecursive();
            List<int> allPeakBoundaries = GetPeakBoundariesLinear();
            int previousPeakIdx = 0;
            List<IndexedMassSpectralPeak> indexedTimsTofPeaks = new List<IndexedMassSpectralPeak>(allPeakBoundaries.Count-1);
            for (int i = 1; i < allPeakBoundaries.Count; i++)
            {
                int index = allPeakBoundaries[i];
                IndexedTimsTofPeak indexedTimsTofPeak;
                if (index - previousPeakIdx > 7) // If we have more than 7 peaks, there is a chance we can split them further
                {
                    var splitIndices = FindPeakSplitIndices(
                        timePoints: IonMobilityPeaks[previousPeakIdx..index],
                        apexTimepointIndex: IonMobilityPeaks[previousPeakIdx..index].IndexOf(IonMobilityPeaks[previousPeakIdx..index].MaxBy(p => p.Intensity)),
                        allowedMissedScans: _allowedMissedScans);
                    if (splitIndices.Any())
                    {
                        foreach (int splitIndex in splitIndices)
                        {
                            indexedTimsTofPeak = GetSinglePeak(previousPeakIdx, splitIndex);
                            if (indexedTimsTofPeak != null)
                                indexedTimsTofPeaks.Add(indexedTimsTofPeak);
                            previousPeakIdx = splitIndex;
                        }
                    }
                    // While it is technically possible that we could encounter a case where a peak needs to be split recursively,
                    // There's a large performance penalty for doing a recursive split for every peak. So, we only attempt to split once and hope that is good
                }
                indexedTimsTofPeak = GetSinglePeak(previousPeakIdx, index);
                if (indexedTimsTofPeak != null)
                    indexedTimsTofPeaks.Add(indexedTimsTofPeak);
                previousPeakIdx = index;
                
            }
            indexedTimsTofPeaks.Sort((x, y) => x.Mz.CompareTo(y.Mz));
            return indexedTimsTofPeaks;
        }

        private IndexedTimsTofPeak GetSinglePeak(int startIndex, int endIndexExclusive)
        {
            if (endIndexExclusive - startIndex >= 5) // Need at least three scans for a peak
            {
                int intensity = IonMobilityPeaks[startIndex..endIndexExclusive].Sum(p => p.IntegerIntensity);
                return new IndexedTimsTofPeak(
                    IonMobilityPeaks[startIndex..endIndexExclusive].MaxBy(p => p.IntegerIntensity).Mz,
                    intensity,
                    ZeroBasedMs1FrameIndex,
                    RetentionTime,
                    GetWeightedAverageTimsIndex(intensity, IonMobilityPeaks[startIndex..endIndexExclusive]));
            }
            return null;
        }

        private int _allowedMissedScans = 20;

        public List<int> GetPeakBoundariesLinear()
        {
            List<int> peakBoundaries = new List<int> { 0 };
            int previousTimsScanNumber = IonMobilityPeaks[0].OneBasedTimsScanNumber;
            for(int i = 1; i < IonMobilityPeaks.Count; i++)
            {
                if (IonMobilityPeaks[i].OneBasedTimsScanNumber - previousTimsScanNumber > _allowedMissedScans)
                    peakBoundaries.Add(i);
                previousTimsScanNumber = IonMobilityPeaks[i].OneBasedTimsScanNumber;
            }
            peakBoundaries.Add(IonMobilityPeaks.Count);
            return peakBoundaries;
        }

        private double GetWeightedAverageMz(int totalIntensity, List<IonMobilityPeak> peaks)
        {
            double mz = 0;
            foreach(var peak in peaks)
            {
                mz += peak.Mz * peak.Intensity / totalIntensity;
            }
            return mz;
        }

        private double GetWeightedAverageTimsIndex(int totalIntensity, List<IonMobilityPeak> peaks)
        {
            double index = 0;
            foreach (var peak in peaks)
            {
                index += (double)peak.OneBasedTimsScanNumber * peak.Intensity / totalIntensity;
            }
            return index;
        }

        private List<int> GetPeakBoundariesRecursive()
        {
            var allBoundarySet = GetPeakBoundariesRecursiveHelper(
                new List<int> { 0, IonMobilityPeaks.Count } ).ToHashSet();
            return allBoundarySet.OrderBy(i => i).ToList();
        }
        // Would be way mor efficient just to do all this by index

        private List<int> GetPeakBoundariesRecursiveHelper(List<int> specificPeakBoundaries)
        {
            int startIndex = specificPeakBoundaries.First();
            int endIndex = specificPeakBoundaries.Last();
            if(endIndex <= startIndex) // Base case
            {
                return specificPeakBoundaries;
            }
            if(endIndex - startIndex <= 3) //minimum length 
            {
                // Split if there are more than _allowedMissedScans between two peaks
                int previousIndex = startIndex;
                for (int i = startIndex + 1; i < endIndex; i++)
                {
                    if (IonMobilityPeaks[i].OneBasedTimsScanNumber - IonMobilityPeaks[previousIndex].OneBasedTimsScanNumber > _allowedMissedScans)
                    {
                        specificPeakBoundaries.Add(i);
                    }
                    previousIndex = i;
                }
                return specificPeakBoundaries;
            }

            var peakSplits = FindPeakSplitIndices(
                timePoints: IonMobilityPeaks[startIndex..endIndex],
                apexTimepointIndex: IonMobilityPeaks[startIndex..endIndex].IndexOf(IonMobilityPeaks[startIndex..endIndex].MaxBy(p => p.Intensity)),
                allowedMissedScans: _allowedMissedScans);

            // base case, can't split any further
            if (!peakSplits.Any())
                return specificPeakBoundaries;

            // recursive part
            int previousIdx = 0;
            foreach(int split in peakSplits) 
            { 
                specificPeakBoundaries.AddRange(GetPeakBoundariesRecursiveHelper(new List<int> { previousIdx, split }));
            }
            if (previousIdx < specificPeakBoundaries.Last())
                specificPeakBoundaries.AddRange(GetPeakBoundariesRecursiveHelper(new List<int> { previousIdx, specificPeakBoundaries.Last() }));

            return specificPeakBoundaries;
        }
    }
}
