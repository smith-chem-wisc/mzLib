using Chemistry;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public class ChromatographicPeak
    {
        #region Public Fields

        public double intensity;
        public IsotopeCluster apexPeak;
        public bool isMbrFeature;
        public RawFileInfo rawFileInfo;
        public List<Identification> identifyingScans;
        public List<IsotopeCluster> isotopeClusters;
        public double splitRT;
        public int numChargeStatesObserved;

        #endregion Public Fields

        #region Public Constructors

        public ChromatographicPeak()
        {
            splitRT = 0;
            numChargeStatesObserved = 0;
            MassError = double.NaN;
            NumIdentificationsByBaseSeq = 1;
            NumIdentificationsByFullSeq = 1;
            identifyingScans = new List<Identification>();
            isotopeClusters = new List<IsotopeCluster>();
        }

        #endregion Public Constructors

        #region Public Properties

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

        public int NumIdentificationsByBaseSeq { get; private set; }
        public int NumIdentificationsByFullSeq { get; private set; }
        public double MassError { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public void CalculateIntensityForThisFeature(bool integrate)
        {
            if (isotopeClusters.Any())
            {
                apexPeak = isotopeClusters.Where(p => p.isotopeClusterIntensity == isotopeClusters.Max(v => v.isotopeClusterIntensity)).FirstOrDefault();

                if (integrate)
                    intensity = isotopeClusters.Select(p => p.isotopeClusterIntensity).Sum();
                else
                    intensity = apexPeak.isotopeClusterIntensity;

                MassError = ((ClassExtensions.ToMass(apexPeak.indexedPeak.mz, apexPeak.chargeState) - identifyingScans.First().massToLookFor) / identifyingScans.First().massToLookFor) * 1e6;
                numChargeStatesObserved = isotopeClusters.Select(p => p.chargeState).Distinct().Count();
            }
            else
                apexPeak = null;
        }

        public void MergeFeatureWith(IEnumerable<ChromatographicPeak> otherFeatures, bool integrate)
        {
            var thisFeaturesPeaks = this.isotopeClusters.Select(p => p.indexedPeak);

            foreach (var otherFeature in otherFeatures)
            {
                if (otherFeature != this)
                {
                    this.identifyingScans = this.identifyingScans.Union(otherFeature.identifyingScans).Distinct().ToList();
                    ResolveIdentifications();
                    this.isotopeClusters.AddRange(otherFeature.isotopeClusters.Where(p => !thisFeaturesPeaks.Contains(p.indexedPeak)));
                    otherFeature.intensity = -1;
                }
            }
            this.CalculateIntensityForThisFeature(integrate);
        }

        public void ResolveIdentifications()
        {
            this.NumIdentificationsByBaseSeq = identifyingScans.Select(v => v.BaseSequence).Distinct().Count();
            this.NumIdentificationsByFullSeq = identifyingScans.Select(v => v.ModifiedSequence).Distinct().Count();
        }

        override public string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(rawFileInfo.filenameWithoutExtension + "\t");
            sb.Append(string.Join("|", identifyingScans.Select(p => p.BaseSequence).Distinct()) + '\t');
            sb.Append(string.Join("|", identifyingScans.Select(p => p.ModifiedSequence).Distinct()) + '\t');

            var t = identifyingScans.SelectMany(p => p.proteinGroupNames).Distinct().OrderBy(p => p);
            if (t.Any())
                sb.Append(string.Join(";", t) + '\t');
            else
                sb.Append("" + '\t');

            sb.Append("" + identifyingScans.First().monoisotopicMass + '\t');
            if (!isMbrFeature)
                sb.Append("" + identifyingScans.First().ms2RetentionTimeInMinutes + '\t');
            else
                sb.Append("" + '\t');
            sb.Append("" + identifyingScans.First().precursorChargeState + '\t');
            sb.Append("" + ClassExtensions.ToMz(identifyingScans.First().monoisotopicMass, identifyingScans.First().precursorChargeState) + '\t');
            sb.Append("" + intensity + "\t");

            if (apexPeak != null)
            {
                sb.Append("" + isotopeClusters.Select(p => p.retentionTime).Min() + "\t");
                sb.Append("" + apexPeak.retentionTime + "\t");
                sb.Append("" + isotopeClusters.Select(p => p.retentionTime).Max() + "\t");

                sb.Append("" + apexPeak.indexedPeak.mz + "\t");
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
            sb.Append("" + NumIdentificationsByBaseSeq + "\t");
            sb.Append("" + NumIdentificationsByFullSeq + "\t");
            sb.Append("" + splitRT + "\t");
            sb.Append("" + MassError);

            return sb.ToString();
        }

        #endregion Public Methods
    }
}