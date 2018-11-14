namespace FlashLFQ
{
    public class RetentionTimeCalibDataPoint
    {
        public readonly ChromatographicPeak DonorFilePeak;
        public readonly ChromatographicPeak AcceptorFilePeak;
        public readonly double RtDiff;

        public RetentionTimeCalibDataPoint(ChromatographicPeak donorFilePeak, ChromatographicPeak acceptorFilePeak)
        {
            DonorFilePeak = donorFilePeak;
            AcceptorFilePeak = acceptorFilePeak;
            RtDiff = acceptorFilePeak.Apex.IndexedPeak.RetentionTime - donorFilePeak.Apex.IndexedPeak.RetentionTime;
        }

        // for debugging
        public override string ToString()
        {
            return "DonorRT: " + DonorFilePeak.Apex.IndexedPeak.RetentionTime.ToString("F3") 
                 + " AcceptorRT: " + AcceptorFilePeak.Apex.IndexedPeak.RetentionTime.ToString("F3") 
                 + " Diff: " + RtDiff.ToString("F3");
        }
    }
}
