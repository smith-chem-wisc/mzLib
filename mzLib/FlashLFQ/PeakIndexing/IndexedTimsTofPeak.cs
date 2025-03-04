using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    /// <summary>
    /// For TimsTofPeaks, the ZeroBasedMs1ScanIndex refers to the MS1 Frame index, not the ion mobility scan index
    /// </summary>
    [Serializable]
    public class IndexedTimsTofPeak : IndexedMassSpectralPeak
    {
        public double TimsIndex { get; init; }

        /// <summary>
        /// Stores the information associated with a specific m/z value in one timsTOF frame
        /// The given m/z can be observed in multiple ion mobility scans, each with a different intensity
        /// This information is stored in the IonMobilityPeaks list
        /// </summary>
        /// <param name="mz"></param>
        /// <param name="zeroBasedMs1FrameIndex">One</param>
        /// <param name="retentionTime"></param>
        /// <param name="ionMobilityPeak"></param>
        public IndexedTimsTofPeak(double mz, double intensity, int zeroBasedMs1FrameIndex, double retentionTime, double timsIndex) : 
            base(mz, intensity, zeroBasedMs1FrameIndex, retentionTime)
        {
            TimsIndex = timsIndex;
        }

        public IndexedTimsTofPeak(double mz, double intensity, int zeroBasedMs1FrameIndex, double retentionTime):
            base(mz, intensity, zeroBasedMs1FrameIndex, retentionTime) { }

        public override bool Equals(object obj)
        {
            var otherPeak = (IndexedTimsTofPeak)obj;

            return otherPeak != null
                && otherPeak.Mz == Mz
                && otherPeak.ZeroBasedMs1ScanIndex == ZeroBasedMs1ScanIndex
                && otherPeak.TimsIndex == TimsIndex;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(Mz, ZeroBasedMs1ScanIndex);
        }

        public override string ToString()
        {
            return Mz.ToString("F3") + "; " + ZeroBasedMs1ScanIndex;
        }
    }

    /// <summary>
    /// Stores the scan number and intensity of an ion mobility peak
    /// in the TIMS-TOF data, there is a one-to-one correspondence between
    /// the oneBasedTimsScanNumber and 1/K0.
    /// This struct does not store information about the m/z of the peak!
    /// </summary>
    [Serializable]
    public readonly struct IonMobilityPeak(double mz, int intensity, int oneBasedTimsScanNumber) : ISingleScanDatum, IComparable<IonMobilityPeak>, IEquatable<IonMobilityPeak>
    {
        public int OneBasedTimsScanNumber { get; } = oneBasedTimsScanNumber;
        public int IntegerIntensity { get; } = intensity;
        public double Mz { get; } = mz;
        public double Intensity => (double)IntegerIntensity;
        public double RelativeSeparationValue => OneBasedTimsScanNumber;
        public int ZeroBasedScanIndex => OneBasedTimsScanNumber - 1;

        public int CompareTo(IonMobilityPeak other)
        {
            return OneBasedTimsScanNumber.CompareTo(other.OneBasedTimsScanNumber);
        }

        public bool Equals(IonMobilityPeak other)
        {
            return OneBasedTimsScanNumber == other.OneBasedTimsScanNumber
                && IntegerIntensity == other.IntegerIntensity
                && Math.Abs(Mz - other.Mz) < 0.0001;
        }
    }

}
