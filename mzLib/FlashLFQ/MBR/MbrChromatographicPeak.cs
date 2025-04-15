using FlashLFQ.PEP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using System.Globalization;

namespace FlashLFQ
{
    public class MbrChromatographicPeak : ChromatographicPeak
    {
        public double PredictedRetentionTime { get; init; }
        public double RtPredictionError { get; set; }
        internal double MbrQValue { get; set; }
        public ChromatographicPeakData PepPeakData { get; set; }
        public double? MbrPep { get; set; }
        /// <summary>
        /// Bool that describes whether the retention time of this peak was randomized
        /// If true, implies that this peak is a decoy peak identified by the MBR algorithm
        /// </summary>
        public bool RandomRt { get; }
        //These scores are populated by the MBR Scorer class. More details can be found there
        /// <summary>
        /// A score ranging between 0 and 100. Higher scores are better.
        /// </summary>
        public double MbrScore { get; set; }
        public double PpmScore { get; set; }
        public double IntensityScore { get; set; }
        public double RtScore { get; set; }
        public double ScanCountScore { get; set; }
        public double IsotopicDistributionScore { get; set; }

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

            sb.Append(Identifications.First().MonoisotopicMass.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append('\t'); // No Ms2 ID means no MS2 ID Retention time
            sb.Append(Identifications.First().PrecursorChargeState.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(ClassExtensions.ToMz(Identifications.First().MonoisotopicMass, Identifications.First().PrecursorChargeState).ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(Intensity.ToString(CultureInfo.InvariantCulture) + "\t");

            if (Apex != null)
            {
                sb.Append(IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime).ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(Apex.IndexedPeak.RetentionTime.ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime).ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(Apex.IndexedPeak.Mz.ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(Apex.ChargeState.ToString(CultureInfo.InvariantCulture) + "\t");
            }
            else
            {
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
            }

            sb.Append(NumChargeStatesObserved.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append("MBR" + "\t");

            // MBR Exclusive fields
            sb.Append(MbrQValue.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(MbrPep is double pep ? pep.ToString(CultureInfo.InvariantCulture) + "\t" : '\t');

            sb.Append(Identifications.Count.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(NumIdentificationsByBaseSeq.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(NumIdentificationsByFullSeq.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(SplitRT.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(MassError.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(DecoyPeptide.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(RandomRt); // Because this isn't an MBR peak, the Random RT Field will always be false

            return sb.ToString();
        }
    }
}
