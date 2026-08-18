namespace MassSpectrometry
{
    public class ScanInfo
    {
        public readonly int OneBasedScanNumber;
        public readonly int ZeroBasedScanIndex;
        public readonly double RetentionTime;
        public readonly int MsnOrder;

        /// <summary>
        /// Total ion current of the scan, or <see cref="double.NaN"/> when the indexing engine that built
        /// this instance did not record one. Callers that write it out are expected to check for NaN.
        /// </summary>
        public readonly double TotalIonCurrent;

        public ScanInfo(int oneBasedScanNumber, int zeroBasedScanIndex, double retentionTime, int msnOrder)
            : this(oneBasedScanNumber, zeroBasedScanIndex, retentionTime, msnOrder, double.NaN)
        {
        }

        public ScanInfo(int oneBasedScanNumber, int zeroBasedScanIndex, double retentionTime, int msnOrder, double totalIonCurrent)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            ZeroBasedScanIndex = zeroBasedScanIndex;
            RetentionTime = retentionTime;
            MsnOrder = msnOrder;
            TotalIonCurrent = totalIonCurrent;
        }

        public override string ToString()
        {
            return ZeroBasedScanIndex + "; " + OneBasedScanNumber + "; " + RetentionTime + "; " + MsnOrder;
        }
    }
}
