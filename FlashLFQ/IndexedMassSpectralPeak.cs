namespace FlashLFQ
{
    public class IndexedMassSpectralPeak
    {
        #region Public Fields

        public readonly int zeroBasedIndexOfPeakInScan;
        public readonly double intensity;
        public readonly int oneBasedScanNumber;
        public readonly double mz;

        #endregion Public Fields

        #region Public Constructors

        public IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedIndexOfPeakInScan, int oneBasedScanNumber)
        {
            this.mz = mz;
            this.intensity = intensity;
            this.zeroBasedIndexOfPeakInScan = zeroBasedIndexOfPeakInScan;
            this.oneBasedScanNumber = oneBasedScanNumber;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return mz.ToString("F3") + "; " + oneBasedScanNumber + "; " + zeroBasedIndexOfPeakInScan;
        }

        #endregion Public Methods
    }
}