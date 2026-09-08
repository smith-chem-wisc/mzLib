using System;
using System.Linq;

namespace FlashLFQ
{
    public class RetentionTimeCalibDataPoint : IComparable
    {
        public readonly ChromatographicPeak DonorFilePeak;
        public readonly ChromatographicPeak AcceptorFilePeak;
        public double RtDiff { get; protected set; }

        public RetentionTimeCalibDataPoint(ChromatographicPeak donorFilePeak, ChromatographicPeak acceptorFilePeak)
        {
            DonorFilePeak = donorFilePeak;
            AcceptorFilePeak = acceptorFilePeak;

            if (donorFilePeak != null && acceptorFilePeak != null)
            {
                RtDiff = donorFilePeak.Apex.IndexedPeak.RetentionTime - acceptorFilePeak.Apex.IndexedPeak.RetentionTime;
            }
            else
            {
                RtDiff = double.NaN;
            }
        }

        public int CompareTo(object obj)
        {
            var otherPoint = (RetentionTimeCalibDataPoint)obj;

            return this.DonorFilePeak.Apex.IndexedPeak.RetentionTime.CompareTo(otherPoint.DonorFilePeak.Apex.IndexedPeak.RetentionTime);
        }

        // for debugging
        public override string ToString()
        {
            return "DonorRT: " + DonorFilePeak.Apex.IndexedPeak.RetentionTime.ToString("F3")
                 + " AcceptorRT: " + AcceptorFilePeak.Apex.IndexedPeak.RetentionTime.ToString("F3")
                 + " Diff: " + RtDiff.ToString("F3");
        }
    }

    public class RetentionTimeCalibrationCurve
    {
        public readonly RetentionTimeCalibDataPoint[] DataPoints;
        public RetentionTimeCalibrationCurve(RetentionTimeCalibDataPoint[] dataPoints)
        {
            dataPoints.Order();
            DataPoints = dataPoints;
        }
    }
}
