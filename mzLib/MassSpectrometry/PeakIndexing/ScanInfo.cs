namespace MassSpectrometry
{
    public class ScanInfo
    {
        public readonly int OneBasedScanNumber;
        public readonly int ZeroBasedScanIndex;
        public readonly double RetentionTime;
        public readonly int MsnOrder;

        /// <summary>
        /// The scan's total ion current, or NaN when it was recorded by a caller that did not have the scan in hand.
        /// Carried here so that summarizing a run's scans does not mean re-reading the file that was just indexed.
        /// </summary>
        public readonly double TotalIonCurrent;

        public ScanInfo(int oneBasedScanNumber, int zeroBasedScanIndex, double retentionTime, int msnOrder,
            double totalIonCurrent = double.NaN)
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
