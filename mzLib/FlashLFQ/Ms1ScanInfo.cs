namespace FlashLFQ
{
    public class Ms1ScanInfo
    {
        public readonly int OneBasedScanNumber;
        public readonly int ZeroBasedMs1ScanIndex;
        public readonly double RetentionTime;

        public Ms1ScanInfo(int oneBasedScanNumber, int zeroBasedMs1ScanIndex, double retentionTime)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            ZeroBasedMs1ScanIndex = zeroBasedMs1ScanIndex;
            RetentionTime = retentionTime;
        }

        public override string ToString()
        {
            return ZeroBasedMs1ScanIndex + "; " + OneBasedScanNumber + "; " + RetentionTime;
        }
    }
}
