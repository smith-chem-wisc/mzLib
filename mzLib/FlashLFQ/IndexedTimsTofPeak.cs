using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    [Serializable]
    internal class IndexedTimsTofPeak : IndexedMassSpectralPeak
    {
        public List<IonMobilityPeak> IonMobilityPeaks { get; init; }
        public int ZeroBasedMs1FrameIndex => ZeroBasedMs1ScanIndex;

        /// <summary>
        /// Stores the information associated with a specific m/z value in one timsTOF frame
        /// The given m/z can be observed in multiple ion mobility scans, each with a different intensity
        /// This information is stored in the IonMobilityPeaks list
        /// </summary>
        /// <param name="mz"></param>
        /// <param name="zeroBasedMs1FrameIndex">One</param>
        /// <param name="retentionTime"></param>
        /// <param name="ionMobilityPeak"></param>
        public IndexedTimsTofPeak(double mz, int zeroBasedMs1FrameIndex, double retentionTime, IonMobilityPeak ionMobilityPeak) : 
            base(mz, intensity: ionMobilityPeak.Intensity, zeroBasedMs1FrameIndex, retentionTime)
        {
            // Initialize the list of ion mobility peaks with the first peak
            // Set initial size to 32 to minimize resizing
            IonMobilityPeaks = new List<IonMobilityPeak>(32);
            IonMobilityPeaks.Add(ionMobilityPeak);
        }

        /// <summary>
        /// Adds an additional ion mobility peak to the list of ion mobility peaks
        /// and updated the intensity 
        /// </summary>
        public void AddIonMobilityPeak(IonMobilityPeak ionMobilityPeak)
        {
            IonMobilityPeaks.Add(ionMobilityPeak);
            Intensity += ionMobilityPeak.Intensity;
        }

        public override bool Equals(object obj)
        {
            var otherPeak = (IndexedTimsTofPeak)obj;

            return otherPeak != null
                && otherPeak.Mz == this.Mz
                && otherPeak.ZeroBasedMs1ScanIndex == this.ZeroBasedMs1ScanIndex;
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
    internal readonly struct IonMobilityPeak(int oneBasedTimsScanNumber, double intensity)
    {
        public readonly int OneBasedTimsScanNumber { get; } = oneBasedTimsScanNumber;
        public readonly double Intensity { get; } = intensity;
    }

}
