using MassSpectrometry;

namespace FlashLFQ
{
    public class IndexedMassSpectralPeak
    {
        public readonly MassSpectralPeak mainPeak;
        public IMsDataScan<IMzSpectrum<IMzPeak>> scan { get; private set; }
        public readonly int zeroBasedIndexOfPeakInScan;
        public readonly double massSpectralPeakIntensity;
        public readonly double retentionTime;
        public readonly int oneBasedScanNumber;

        public IndexedMassSpectralPeak(MassSpectralPeak peak, IMsDataScan<IMzSpectrum<IMzPeak>> scan, int index)
        {
            mainPeak = peak;
            this.scan = scan;
            zeroBasedIndexOfPeakInScan = index;
            massSpectralPeakIntensity = peak.Intensity;
            oneBasedScanNumber = scan.OneBasedScanNumber;
            retentionTime = scan.RetentionTime;
        }

        public void Compress()
        {
            scan = null;
        }

        public override string ToString()
        {
            if (mainPeak != null)
                return System.Math.Round(mainPeak.Mz, 5) + "; " + System.Math.Round(retentionTime, 2) + "; " + oneBasedScanNumber;
            else
                return "--";
        }
    }
}
