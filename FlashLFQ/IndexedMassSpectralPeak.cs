namespace FlashLFQ
{
    public class IndexedMassSpectralPeak
    {
        public readonly int zeroBasedIndexOfPeakInScan;
        public readonly double intensity;
        public readonly int oneBasedScanNumber;
        public readonly double mz;

        public IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedIndexOfPeakInScan, int oneBasedScanNumber)
        {
            this.mz = mz;
            this.intensity = intensity;
            this.zeroBasedIndexOfPeakInScan = zeroBasedIndexOfPeakInScan;
            this.oneBasedScanNumber = oneBasedScanNumber;
        }

        public override string ToString()
        {
            return mz.ToString("F3") + "; " + oneBasedScanNumber + "; " + zeroBasedIndexOfPeakInScan;
        }
    }
}