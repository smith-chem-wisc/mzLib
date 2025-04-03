using FlashLFQ.PEP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace FlashLFQ
{
    public class MbrChromatographicPeak : ChromatographicPeak
    {
        public double PredictedRetentionTime { get; init; }
        public double MbrScore { get; set; }
        public double PpmScore { get; set; }
        public double IntensityScore { get; set; }
        public double RtScore { get; set; }
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

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(SpectraFileInfo.FilenameWithoutExtension + "\t");
            sb.Append(string.Join("|", Identifications.Select(p => p.BaseSequence).Distinct()) + '\t');
            sb.Append(string.Join("|", Identifications.Select(p => p.ModifiedSequence).Distinct()) + '\t');
            var t = Identifications.SelectMany(p => p.ProteinGroups.Select(v => v.ProteinGroupName)).Distinct().OrderBy(p => p);
            if (t.Any())
            {
                sb.Append(string.Join(";", t) + '\t');
                sb.Append(string.Join(";", Identifications.SelectMany(id => id.ProteinGroups).Select(p => p.Organism).Distinct()) + '\t');
            }
            else
            {
                sb.Append('\t');
                sb.Append('\t');
            }

            sb.Append(Identifications.First().MonoisotopicMass + '\t');
            sb.Append('\t'); // No Ms2 ID means no MS2 ID Retention time
            sb.Append(Identifications.First().PrecursorChargeState + '\t');
            sb.Append(ClassExtensions.ToMz(Identifications.First().MonoisotopicMass, Identifications.First().PrecursorChargeState) + '\t');
            sb.Append(Intensity + "\t");

            if (Apex != null)
            {
                sb.Append(IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime) + "\t");
                sb.Append(Apex.IndexedPeak.RetentionTime + "\t");
                sb.Append(IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime) + "\t");
                sb.Append(Apex.IndexedPeak.Mz + "\t");
                sb.Append(Apex.ChargeState + "\t");
            }
            else
            {
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
            }

            sb.Append(NumChargeStatesObserved + "\t");
            sb.Append("MBR" + "\t");

            // MBR Exclusive fields
            sb.Append(MbrQValue + "\t");
            sb.Append(MbrPep + "\t");

            sb.Append(Identifications.Count + "\t");
            sb.Append(NumIdentificationsByBaseSeq + "\t");
            sb.Append(NumIdentificationsByFullSeq + "\t");
            sb.Append(SplitRT + "\t");
            sb.Append(MassError + "\t");
            sb.Append(DecoyPeptide + "\t");
            sb.Append(RandomRt); // Because this isn't an MBR peak, the Random RT Field will always be false

            return sb.ToString();
        }
    }
}
