using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    /// <summary>
    /// For TimsTofPeaks, the ZeroBasedMs1ScanIndex refers to the MS1 Frame index, not the ion mobility scan index
    /// ITTPs store indices instead of mz values. They must be converted to MZ values using the TimsTofFileReader
    /// </summary>
    [Serializable]
    public struct IndexedTimsTofPeak : IComparable<IndexedTimsTofPeak>
    {
        public uint TofIndex { get; init; }
        public int TimsIndex { get; init; }
        public int Intensity { get; init; }
        public int ZeroBasedMs1FrameIndex { get; init; }

        /// <summary>
        /// Stores the information associated with a specific m/z value in one timsTOF frame
        /// The given m/z can be observed in multiple ion mobility scans, each with a different intensity
        /// This information is stored in the IonMobilityPeaks list
        /// </summary>
        public IndexedTimsTofPeak(uint tofIndex, int timsIndex, int intensity, int zeroBasedMs1FrameIndex) 
        {
            TofIndex = tofIndex;
            TimsIndex = timsIndex;
            Intensity = intensity;
            ZeroBasedMs1FrameIndex = zeroBasedMs1FrameIndex;
        }

        public override bool Equals(object obj)
        {
            var otherPeak = (IndexedTimsTofPeak)obj;

            return otherPeak.TimsIndex == TimsIndex
                && otherPeak.TofIndex == TofIndex 
                && otherPeak.Intensity == Intensity
                && otherPeak.ZeroBasedMs1FrameIndex == ZeroBasedMs1FrameIndex;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(TofIndex, ZeroBasedMs1FrameIndex, TimsIndex);
        }

        public override string ToString()
        {
            return TofIndex + "; " + TimsIndex;
        }

        public int CompareTo(IndexedTimsTofPeak other)
        {
            return ZeroBasedMs1FrameIndex.CompareTo(other.ZeroBasedMs1FrameIndex);
        }
    }

    /// <summary>
    /// Stores the scan number and intensity of an ion mobility peak
    /// in the TIMS-TOF data, there is a one-to-one correspondence between
    /// the oneBasedTimsScanNumber and 1/K0.
    /// </summary>
    public readonly record struct IonMobilityPeak(uint tofIndex, int timsIndex, int intensity) : IComparable<IonMobilityPeak>, IEquatable<IonMobilityPeak>
    {
        public uint TofIndex { get; } = tofIndex;
        public int TimsIndex { get; } = timsIndex;
        public int Intensity { get; } = intensity;

        public int CompareTo(IonMobilityPeak other)
        {
            return TofIndex.CompareTo(other.TofIndex);
        }

        public bool Equals(IonMobilityPeak other)
        {
            return TofIndex == other.TofIndex 
                && TimsIndex == other.TimsIndex
                && Intensity == other.Intensity;
        }
    }

}
