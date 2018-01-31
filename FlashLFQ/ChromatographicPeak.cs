using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public class ChromatographicPeak
    {
        public double intensity;
        public IsotopeCluster apexPeak;
        public bool isMbrFeature;
        public string fileName = "";
        public List<Identification> identifyingScans;
        public List<IsotopeCluster> isotopeClusters;
        public int numIdentificationsByBaseSeq { get; private set; }
        public int numIdentificationsByFullSeq { get; private set; }
        public double splitRT;
        public double massError { get; private set; }
        public int numChargeStatesObserved;
        
        public ChromatographicPeak()
        {
            splitRT = 0;
            numChargeStatesObserved = 0;
            massError = double.NaN;
            numIdentificationsByBaseSeq = 1;
            numIdentificationsByFullSeq = 1;
            identifyingScans = new List<Identification>();
            isotopeClusters = new List<IsotopeCluster>();
        }

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("File Name" + "\t");
                sb.Append("Base Sequence" + "\t");
                sb.Append("Full Sequence" + "\t");
                sb.Append("Protein Group" + "\t");
                sb.Append("Peptide Monoisotopic Mass" + "\t");
                sb.Append("MS2 Retention Time" + "\t");
                sb.Append("Precursor Charge" + "\t");
                sb.Append("Theoretical MZ" + "\t");
                sb.Append("Peak intensity" + "\t");
                sb.Append("Peak RT Start" + "\t");
                sb.Append("Peak RT Apex" + "\t");
                sb.Append("Peak RT End" + "\t");
                sb.Append("Peak MZ" + "\t");
                sb.Append("Peak Charge" + "\t");
                sb.Append("Num Charge States Observed" + "\t");
                sb.Append("Peak Detection Type" + "\t");
                sb.Append("PSMs Mapped" + "\t");
                sb.Append("Base Sequences Mapped" + "\t");
                sb.Append("Full Sequences Mapped" + "\t");
                sb.Append("Peak Split Valley RT" + "\t");
                sb.Append("Peak Apex Mass Error (ppm)");
                return sb.ToString();
            }
        }

        public void CalculateIntensityForThisFeature(bool integrate)
        {
            if (isotopeClusters.Any())
            {
                apexPeak = isotopeClusters.Where(p => p.isotopeClusterIntensity == isotopeClusters.Max(v => v.isotopeClusterIntensity)).FirstOrDefault();

                if (integrate)
                    intensity = isotopeClusters.Select(p => p.isotopeClusterIntensity).Sum();
                else
                    intensity = apexPeak.isotopeClusterIntensity;

                massError = ((ClassExtensions.ToMass(apexPeak.peakWithScan.mainPeak.Mz, apexPeak.chargeState) - identifyingScans.First().massToLookFor) / identifyingScans.First().massToLookFor) * 1e6;
                numChargeStatesObserved = isotopeClusters.Select(p => p.chargeState).Distinct().Count();
            }
            else
                apexPeak = null;
        }

        public void MergeFeatureWith(IEnumerable<ChromatographicPeak> otherFeatures, bool integrate)
        {
            var thisFeaturesPeaks = this.isotopeClusters.Select(p => p.peakWithScan);

            foreach (var otherFeature in otherFeatures)
            {
                if (otherFeature != this)
                {
                    this.identifyingScans = this.identifyingScans.Union(otherFeature.identifyingScans).Distinct().ToList();
                    ResolveIdentifications();
                    this.isotopeClusters.AddRange(otherFeature.isotopeClusters.Where(p => !thisFeaturesPeaks.Contains(p.peakWithScan)));
                    otherFeature.intensity = -1;
                }
            }
            this.CalculateIntensityForThisFeature(integrate);
        }

        public void ResolveIdentifications()
        {
            this.numIdentificationsByBaseSeq = identifyingScans.Select(v => v.BaseSequence).Distinct().Count();
            this.numIdentificationsByFullSeq = identifyingScans.Select(v => v.FullSequence).Distinct().Count();
        }

        override public string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(fileName + "\t");
            sb.Append(string.Join("|", identifyingScans.Select(p => p.BaseSequence).Distinct()) + '\t');
            sb.Append(string.Join("|", identifyingScans.Select(p => p.FullSequence).Distinct()) + '\t');

            var t = identifyingScans.SelectMany(p => p.proteinGroups).Distinct().Select(p => p.proteinGroupName).OrderBy(p => p);
            if (t.Any())
                sb.Append(string.Join(";", t) + '\t');
            else
                sb.Append("" + '\t');

            sb.Append("" + identifyingScans.First().monoisotopicMass + '\t');
            if (!isMbrFeature)
                sb.Append("" + identifyingScans.First().ms2RetentionTime + '\t');
            else
                sb.Append("" + '\t');
            sb.Append("" + identifyingScans.First().chargeState + '\t');
            sb.Append("" + ClassExtensions.ToMz(identifyingScans.First().monoisotopicMass, identifyingScans.First().chargeState) + '\t');
            sb.Append("" + intensity + "\t");

            if (apexPeak != null)
            {
                sb.Append("" + isotopeClusters.Select(p => p.peakWithScan.retentionTime).Min() + "\t");
                sb.Append("" + apexPeak.peakWithScan.retentionTime + "\t");
                sb.Append("" + isotopeClusters.Select(p => p.peakWithScan.retentionTime).Max() + "\t");

                sb.Append("" + apexPeak.peakWithScan.mainPeak.Mz + "\t");
                sb.Append("" + apexPeak.chargeState + "\t");
            }
            else
            {
                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");

                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");
            }

            sb.Append("" + numChargeStatesObserved + "\t");

            if (isMbrFeature)
                sb.Append("" + "MBR" + "\t");
            else
                sb.Append("" + "MSMS" + "\t");

            sb.Append("" + identifyingScans.Count + "\t");
            sb.Append("" + numIdentificationsByBaseSeq + "\t");
            sb.Append("" + numIdentificationsByFullSeq + "\t");
            sb.Append("" + splitRT + "\t");
            sb.Append("" + massError);

            return sb.ToString();
        }
    }
}
