namespace FlashLFQ
{
    public class IndexedMassSpectralPeak
    {
        public readonly int ZeroBasedIndexOfPeakInScan;
        public readonly double Intensity;
        public readonly int OneBasedScanNumber;
        public readonly double Mz;

        public IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedIndexOfPeakInScan, int oneBasedScanNumber)
        {
            this.Mz = mz;
            this.Intensity = intensity;
            this.ZeroBasedIndexOfPeakInScan = zeroBasedIndexOfPeakInScan;
            this.OneBasedScanNumber = oneBasedScanNumber;
        }

        public override string ToString()
        {
            return Mz.ToString("F3") + "; " + OneBasedScanNumber + "; " + ZeroBasedIndexOfPeakInScan;
        }
    }
}