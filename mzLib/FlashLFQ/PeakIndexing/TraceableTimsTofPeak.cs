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
        public List<IonMobilityPeak> IonMobilityPeaks { get; init; }
        public override List<IonMobilityPeak> ScanOrderedPoints => IonMobilityPeaks;
        public override IonMobilityPeak Apex => IonMobilityPeaks.MaxBy(p => p.Intensity);
        public double RetentionTime { get; init; }
        public int ZeroBasedMs1FrameIndex { get; init; }


        public TraceableTimsTofPeak(int zeroBasedMs1FrameIndex, double retentionTime, IonMobilityPeak ionMobilityPeak) : this(zeroBasedMs1FrameIndex, retentionTime)
        {
            IonMobilityPeaks.Add(ionMobilityPeak);
        }

        public TraceableTimsTofPeak(int zeroBasedMs1FrameIndex, double retentionTime)
        {
            IonMobilityPeaks = new List<IonMobilityPeak>(32);
            RetentionTime = retentionTime;
            ZeroBasedMs1FrameIndex = zeroBasedMs1FrameIndex;
        }


        public IEnumerable<IndexedTimsTofPeak> GetIndexedPeaks()
        {
            List<int> allPeakBoundaries = GetPeakBoundariesRecursive();
            int previousPeakIdx = 0;
            foreach (var index in allPeakBoundaries.Skip(1))
            {
                int intensity = IonMobilityPeaks[previousPeakIdx..index].Sum(p => p.IntegerIntensity);
                yield return new IndexedTimsTofPeak(GetWeightedAverageMz(intensity, IonMobilityPeaks[previousPeakIdx..index]),
                    intensity, ZeroBasedMs1FrameIndex, RetentionTime);
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

        private List<int> GetPeakBoundariesRecursive()
        {
            var allBoundarySet = GetPeakBoundariesRecursiveHelper(
                new List<int> { 0, IonMobilityPeaks.Count() } ).ToHashSet();
            return allBoundarySet.OrderBy(i => i).ToList();
        }
        // Would be way mor efficient just to do all this by index

        private List<int> GetPeakBoundariesRecursiveHelper(List<int> specificPeakBoundaries)
        {
            int startIndex = specificPeakBoundaries.First();
            int endIndex = specificPeakBoundaries.Last();
            var peakSplits = FindPeakSplitIndices(IonMobilityPeaks[startIndex..endIndex],
                IonMobilityPeaks[startIndex..endIndex].IndexOf(IonMobilityPeaks[startIndex..endIndex].MaxBy(p => p.Intensity)));

            // base case, can't split any further
            if (!peakSplits.Any())
                return specificPeakBoundaries;

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
