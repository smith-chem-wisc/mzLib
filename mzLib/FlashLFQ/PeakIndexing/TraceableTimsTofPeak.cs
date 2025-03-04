using Easy.Common.Extensions;
using FlashLFQ.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
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
        public override IonMobilityPeak Apex => IonMobilityPeaks.MaxBy(p => p.Intensity);
        public double RetentionTime { get; init; }
        public int ZeroBasedMs1FrameIndex { get; init; }
        public int SummedIntensity { get; private set; }
        public double Mz { get; private set; }

        public TraceableTimsTofPeak(int zeroBasedMs1FrameIndex, double retentionTime, IonMobilityPeak ionMobilityPeak) : this(zeroBasedMs1FrameIndex, retentionTime)
        {
            AddIonMobilityPeak(ionMobilityPeak);
        }

        public TraceableTimsTofPeak(int zeroBasedMs1FrameIndex, double retentionTime)
        {
            IonMobilityPeaks = new List<IonMobilityPeak>(32);
            RetentionTime = retentionTime;
            ZeroBasedMs1FrameIndex = zeroBasedMs1FrameIndex;
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
            int newIntensity = SummedIntensity + ionMobilityPeak.IntegerIntensity;
            // Weighted average for new mz
            Mz = (ionMobilityPeak.Mz * ionMobilityPeak.IntegerIntensity / newIntensity) +(Mz * SummedIntensity / newIntensity);
            SummedIntensity = newIntensity;
        }

        public IEnumerable<IndexedMassSpectralPeak> GetIndexedPeaks()
        {
            List<int> allPeakBoundaries = GetPeakBoundariesRecursive();
            int previousPeakIdx = 0;
            for (int i = 1; i < allPeakBoundaries.Count; i++)
            {
                int index = allPeakBoundaries[i];
                if (index - previousPeakIdx >= 3) // Need at least three scans for a peak
                {
                    int intensity = IonMobilityPeaks[previousPeakIdx..index].Sum(p => p.IntegerIntensity);
                    yield return new IndexedTimsTofPeak(
                        GetWeightedAverageMz(intensity, IonMobilityPeaks[previousPeakIdx..index]),
                        intensity, 
                        ZeroBasedMs1FrameIndex, 
                        RetentionTime,
                        GetWeightedAverageTimsIndex(intensity, IonMobilityPeaks[previousPeakIdx..index]));
                }
                previousPeakIdx = index;
            }
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

        private int _allowedMissedScans = 20;

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
