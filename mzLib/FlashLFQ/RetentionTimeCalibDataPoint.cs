using System;

namespace FlashLFQ
{
    public class RetentionTimeCalibDataPoint : IComparable
    {
        public readonly ChromatographicPeak DonorFilePeak;
        public readonly ChromatographicPeak AcceptorFilePeak;
        public readonly double RtDiff;

        public RetentionTimeCalibDataPoint(ChromatographicPeak donorFilePeak, ChromatographicPeak acceptorFilePeak)
        {
            DonorFilePeak = donorFilePeak;
            AcceptorFilePeak = acceptorFilePeak;

            if (donorFilePeak != null && acceptorFilePeak != null)
            {
                RtDiff = acceptorFilePeak.Apex.IndexedMzPeak.RetentionTime - donorFilePeak.Apex.IndexedMzPeak.RetentionTime;
            }
            else
            {
                RtDiff = double.NaN;
            }
        }

        public int CompareTo(object obj)
        {
            var otherPoint = (RetentionTimeCalibDataPoint)obj;

            return this.DonorFilePeak.Apex.IndexedMzPeak.RetentionTime.CompareTo(otherPoint.DonorFilePeak.Apex.IndexedMzPeak.RetentionTime);
        }

        // for debugging
        public override string ToString()
        {
            return "DonorRT: " + DonorFilePeak.Apex.IndexedMzPeak.RetentionTime.ToString("F3")
                 + " AcceptorRT: " + AcceptorFilePeak.Apex.IndexedMzPeak.RetentionTime.ToString("F3")
                 + " Diff: " + RtDiff.ToString("F3");
        }
    }
}
