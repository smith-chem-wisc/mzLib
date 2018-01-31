using MassSpectrometry;

namespace FlashLFQ
{
    public class IndexedMassSpectralPeak
    {
        #region Public Fields

        public readonly MassSpectralPeak mainPeak;
        public readonly int zeroBasedIndexOfPeakInScan;
        public readonly double massSpectralPeakIntensity;
        public readonly double retentionTime;
        public readonly int oneBasedScanNumber;

        #endregion Public Fields

        #region Public Constructors

        public IndexedMassSpectralPeak(MassSpectralPeak peak, IMsDataScan<IMzSpectrum<IMzPeak>> scan, int index)
        {
            mainPeak = peak;
            this.scan = scan;
            zeroBasedIndexOfPeakInScan = index;
            massSpectralPeakIntensity = peak.Intensity;
            oneBasedScanNumber = scan.OneBasedScanNumber;
            retentionTime = scan.RetentionTime;
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataScan<IMzSpectrum<IMzPeak>> scan { get; private set; }

        #endregion Public Properties

        #region Public Methods

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

        #endregion Public Methods
    }
}