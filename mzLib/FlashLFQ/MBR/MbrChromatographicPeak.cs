using FlashLFQ.PEP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class MbrChromatographicPeak : ChromatographicPeak
    {
        public double PredictedRetentionTime { get; init; }
        public double MbrScore { get; set; }
        public double PpmScore { get; set; }
        public double IntensityScore { get; set; }
        public double RtScore { get; set; }
        public int ScanCount { get; set; }
        public double ScanCountScore { get; set; }
        public double IsotopicDistributionScore { get; set; }
        /// <summary>
        /// Stores the pearson correlation between the apex isotopic envelope and the theoretical isotopic distribution
        /// </summary>
        public double IsotopicPearsonCorrelation { get; set; }
        public double RtPredictionError { get; set; }
        internal double MbrQValue { get; set; }
        public ChromatographicPeakData PepPeakData { get; set; }
        public double? MbrPep { get; set; }
        /// <summary>
        /// Bool that describes whether the retention time of this peak was randomized
        /// If true, implies that this peak is a decoy peak identified by the MBR algorithm
        /// </summary>
        public bool RandomRt { get; }

        public MbrChromatographicPeak(Identification id, SpectraFileInfo fileInfo, double predictedRetentionTime, bool randomRt)
            : base(id, fileInfo, detectionType: DetectionType.MBR)
        {
            PredictedRetentionTime = predictedRetentionTime;
            RandomRt = randomRt;
        }
    }
}
