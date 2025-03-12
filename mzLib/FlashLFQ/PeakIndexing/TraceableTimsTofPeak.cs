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
    internal class TraceableTimsTofPeak 
    {
        /// <summary>
        /// List 
        /// </summary>
        public List<IonMobilityPeak> IonMobilityPeaks { get; init; }
        //public override List<IonMobilityPeak> ScanOrderedPoints => IonMobilityPeaks;
        public  IonMobilityPeak Apex => _apex;
        private IonMobilityPeak _apex;
        public double RetentionTime { get; init; }
        public int ZeroBasedMs1FrameIndex { get; init; }
        //public int SummedIntensity { get; private set; }
        //public double Mz => _apex.Mz;
        public uint TofIndex => Apex.TofIndex;

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
            if (_apex.Equals(default(IonMobilityPeak)) || ionMobilityPeak.Intensity > Apex.Intensity)
            {
                _apex = ionMobilityPeak;
            }
        }

        public List<IndexedTimsTofPeak> GetIndexedPeaks()
        {
            if (IonMobilityPeaks.Count < 2)
            {
                return null;
            }
            List<int> allPeakBoundaries = GetPeakBoundariesLinear();
            int previousPeakIdx = 0;
            List<IndexedTimsTofPeak> indexedTimsTofPeaks = new List<IndexedTimsTofPeak>(allPeakBoundaries.Count-1);
            for (int i = 1; i < allPeakBoundaries.Count; i++)
            {
                int index = allPeakBoundaries[i];
                if (TryGetSinglePeak(previousPeakIdx, index, out IndexedTimsTofPeak ittPeak))
                    indexedTimsTofPeaks.Add(ittPeak);

                previousPeakIdx = index; 
            }

            return indexedTimsTofPeaks;
        }

        private bool TryGetSinglePeak(int startIndex, int endIndexExclusive, out IndexedTimsTofPeak indexedTimsTofPeak)
        {
            indexedTimsTofPeak = default(IndexedTimsTofPeak);
            if (endIndexExclusive - startIndex < 2) // Need at least two scans for a peak
                return false;

            int intensity = IonMobilityPeaks[startIndex..endIndexExclusive].Sum(p => p.Intensity);
            if (intensity < 200)
                return false;

            var localApex = IonMobilityPeaks[startIndex..endIndexExclusive].MaxBy(p => p.Intensity);
            indexedTimsTofPeak = new IndexedTimsTofPeak(
                localApex.TofIndex,
                localApex.TimsIndex,
                localApex.Intensity,
                ZeroBasedMs1FrameIndex);
            return true;

        }

        private int _allowedMissedScans = 4;

        public List<int> GetPeakBoundariesLinear()
        {
            List<int> peakBoundaries = new List<int> { 0 };
            int previousTimsScanNumber = IonMobilityPeaks[0].TimsIndex;
            for(int i = 1; i < IonMobilityPeaks.Count; i++)
            {
                if (IonMobilityPeaks[i].TimsIndex - previousTimsScanNumber > _allowedMissedScans)
                    peakBoundaries.Add(i);
                previousTimsScanNumber = IonMobilityPeaks[i].TimsIndex;
            }
            peakBoundaries.Add(IonMobilityPeaks.Count);
            return peakBoundaries;
        }
    }
}
