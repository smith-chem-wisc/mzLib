using System;
using System.Collections.Generic;
using System.Linq;

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
                RtDiff = donorFilePeak.Apex.IndexedPeak.RetentionTime - acceptorFilePeak.Apex.IndexedPeak.RetentionTime;
            }
            else
            {
                RtDiff = double.NaN;
            }
        }

        /// <summary>
        /// Test seam: builds a data point carrying only an <see cref="RtDiff"/>, with no underlying peaks.
        /// Lets tests exercise the prediction-error math in <see cref="MbrScorer.AddRtPredErrorDistribution"/>
        /// without constructing full chromatographic peaks. Points built this way have a null donor peak, so
        /// they sort before points that have one (see <see cref="CompareTo"/>).
        /// </summary>
        internal RetentionTimeCalibDataPoint(double rtDiff)
        {
            RtDiff = rtDiff;
        }

        /// <summary>
        /// Orders data points by the donor peak's apex retention time. Points without a donor apex
        /// retention time (e.g. the probe point used for a binary search) sort before those that have one;
        /// two such points compare equal, so a stable sort leaves them in their original order.
        /// </summary>
        public int CompareTo(object obj)
        {
            var otherPoint = (RetentionTimeCalibDataPoint)obj;

            double? thisRt = DonorFilePeak?.Apex?.IndexedPeak.RetentionTime;
            double? otherRt = otherPoint.DonorFilePeak?.Apex?.IndexedPeak.RetentionTime;

            return Nullable.Compare(thisRt, otherRt);
        }

        // for debugging
        public override string ToString()
        {
            return "DonorRT: " + DonorFilePeak.Apex.IndexedPeak.RetentionTime.ToString("F3")
                 + " AcceptorRT: " + AcceptorFilePeak.Apex.IndexedPeak.RetentionTime.ToString("F3")
                 + " Diff: " + RtDiff.ToString("F3");
        }
    }

    /// <summary>
    /// The set of anchor peptides shared between a donor and an acceptor file, kept ordered by the donor
    /// peak's apex retention time. Match-between-runs relies on that ordering both to build the local
    /// alignment (a binary search over <see cref="DataPoints"/>) and to estimate its prediction error, so
    /// the ordering is established once here rather than at each use site.
    /// </summary>
    public class RetentionTimeCalibrationCurve
    {
        public readonly RetentionTimeCalibDataPoint[] DataPoints;

        public RetentionTimeCalibrationCurve(IEnumerable<RetentionTimeCalibDataPoint> dataPoints)
        {
            DataPoints = dataPoints.Order().ToArray();
        }

        public int Count => DataPoints.Length;
    }
}
