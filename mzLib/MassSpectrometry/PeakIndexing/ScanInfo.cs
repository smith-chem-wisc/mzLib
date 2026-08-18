namespace MassSpectrometry
{
    public class ScanInfo
    {
        public readonly int OneBasedScanNumber;
        public readonly int ZeroBasedScanIndex;
        public readonly double RetentionTime;
        public readonly int MsnOrder;

        public ScanInfo(int oneBasedScanNumber, int zeroBasedScanIndex, double retentionTime, int msnOrder)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            ZeroBasedScanIndex = zeroBasedScanIndex;
            RetentionTime = retentionTime;
            MsnOrder = msnOrder;
        }

        public override string ToString()
        {
            return ZeroBasedScanIndex + "; " + OneBasedScanNumber + "; " + RetentionTime + "; " + MsnOrder;
        }
    }
}
