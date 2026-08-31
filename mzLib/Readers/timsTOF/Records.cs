using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    internal readonly struct Ms1Record
    {
        internal int PrecursorId { get; }
        internal int ScanStart { get; }
        internal int ScanEnd { get; }
        internal double ScanMedian { get; }

        public Ms1Record(int precursorId, int scanStart, int scanEnd, double scanMedian)
        {
            PrecursorId = precursorId;
            ScanStart = scanStart;
            ScanEnd = scanEnd;
            ScanMedian = scanMedian;
        }
    }

    internal readonly struct MrmRecord
    {
        internal long FrameId { get; }
        internal int? ScanStart { get; }
        internal int? ScanEnd { get; }
        internal float IsolationMz { get; }
        internal float IsolationWidth { get; }
        internal float CollisionEnergy { get; }
        
        public MrmRecord(long frame, int? scanStart, int? scanEnd, float isolationMz, float isolationWidth, float collisionEnergy)
        {
            FrameId = frame;
            ScanStart = scanStart;
            ScanEnd = scanEnd;
            IsolationMz = isolationMz;
            IsolationWidth = isolationWidth;
            CollisionEnergy = collisionEnergy;
        }
    }

    internal readonly struct PasefRecord
    {
        internal IEnumerable<long> FrameList { get; }
        internal int PrecursorId { get; }
        internal int ScanStart { get; }
        internal int ScanEnd { get; }
        internal double ScanMedian { get; }
        internal float IsolationMz { get; }
        internal float IsolationWidth { get; }
        internal float CollisionEnergy { get; }
        internal float MostAbundantPrecursorMz { get; }
        internal float PrecursorMonoisotopicMz { get; }
        internal int Charge { get; }
        internal float PrecursorIntensity { get; }

        public PasefRecord(
            IEnumerable<long> frameList,
            int precursorId,
            int scanStart,
            int scanEnd,
            double scanMedian,
            float isolationMz,
            float isolationWidth,
            float collisionEnergy,
            float mostAbundantPrecursorMz,
            float precursorMonoisotopicMz,
            int charge, 
            float precursorIntensity)
        {
            FrameList = frameList ?? throw new ArgumentNullException(nameof(frameList));
            PrecursorId = precursorId;
            ScanStart = scanStart;
            ScanEnd = scanEnd;
            ScanMedian = scanMedian;
            IsolationMz = isolationMz;
            IsolationWidth = isolationWidth;
            CollisionEnergy = collisionEnergy;
            MostAbundantPrecursorMz = mostAbundantPrecursorMz;
            PrecursorMonoisotopicMz = precursorMonoisotopicMz;
            Charge = charge;
            PrecursorIntensity = precursorIntensity;
        }
    }
}
