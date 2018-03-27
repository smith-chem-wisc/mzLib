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

        public IndexedMassSpectralPeak(double mz, double intensity, int indexInScan, int oneBasedScanIndex)
        {
            this.mz = mz;
            zeroBasedIndexOfPeakInScan = indexInScan;
            this.intensity = intensity;
            oneBasedScanNumber = oneBasedScanIndex;
        }

        #endregion Public Constructors
    }
}