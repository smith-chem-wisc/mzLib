namespace MassSpectrometry
{
    public class ScanInfo
    {
        public readonly int OneBasedScanNumber;
        public readonly int ZeroBasedScanIndex;
        public readonly double RetentionTime;

        public ScanInfo(int oneBasedScanNumber, int zeroBasedScanIndex, double retentionTime)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            ZeroBasedScanIndex = zeroBasedScanIndex;
            RetentionTime = retentionTime;
        }

        public override string ToString()
        {
            return ZeroBasedScanIndex + "; " + OneBasedScanNumber + "; " + RetentionTime;
        }
    }
}
