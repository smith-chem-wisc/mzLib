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
        public IsotopicEnvelope apex { get; private set; }
        public readonly bool isMbrFeature;
        public readonly SpectraFileInfo rawFileInfo;
        public List<Identification> identifications { get; private set; }
        public List<IsotopicEnvelope> isotopicEnvelopes;
        public double splitRT;
        public int numChargeStatesObserved { get; private set; }

        #endregion Public Fields

        #region Public Constructors

        public ChromatographicPeak(Identification id, bool isMbrFeature, SpectraFileInfo fileInfo)
        {
            splitRT = 0;
            numChargeStatesObserved = 0;
            MassError = double.NaN;
            NumIdentificationsByBaseSeq = 1;
            NumIdentificationsByFullSeq = 1;
            identifications = new List<Identification>() { id };
            isotopicEnvelopes = new List<IsotopicEnvelope>();
            this.isMbrFeature = isMbrFeature;
            this.rawFileInfo = fileInfo;
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
            if (isotopicEnvelopes.Any())
            {
                double maxIntensity = isotopicEnvelopes.Max(p => p.intensity);
                apex = isotopicEnvelopes.Where(p => p.intensity == maxIntensity).First();

                if (integrate)
                {
                    intensity = isotopicEnvelopes.Sum(p => p.intensity);
                }
                else
                {
                    intensity = apex.intensity;
                }

                MassError = identifications.Min(p => ((ClassExtensions.ToMass(apex.indexedPeak.mz, apex.chargeState) - p.monoisotopicMass) / p.monoisotopicMass) * 1e6);
                numChargeStatesObserved = isotopicEnvelopes.Select(p => p.chargeState).Distinct().Count();
            }
            else
            {
                intensity = 0;
                MassError = double.NaN;
                numChargeStatesObserved = 0;
                apex = null;
            }
        }

        public void MergeFeatureWith(IEnumerable<ChromatographicPeak> otherFeatures, bool integrate)
        {
            var thisFeaturesPeaks = this.isotopicEnvelopes.Select(p => p.indexedPeak);

            foreach (var otherFeature in otherFeatures)
            {
                if (otherFeature != this)
                {
                    this.identifications = this.identifications.Union(otherFeature.identifications).Distinct().ToList();
                    ResolveIdentifications();
                    this.isotopicEnvelopes.AddRange(otherFeature.isotopicEnvelopes.Where(p => !thisFeaturesPeaks.Contains(p.indexedPeak)));
                    otherFeature.intensity = -1;
                }
            }
            this.CalculateIntensityForThisFeature(integrate);
        }

        public void ResolveIdentifications()
        {
            this.NumIdentificationsByBaseSeq = identifications.Select(v => v.BaseSequence).Distinct().Count();
            this.NumIdentificationsByFullSeq = identifications.Select(v => v.ModifiedSequence).Distinct().Count();
        }

        override public string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(rawFileInfo.filenameWithoutExtension + "\t");
            sb.Append(string.Join("|", identifications.Select(p => p.BaseSequence).Distinct()) + '\t');
            sb.Append(string.Join("|", identifications.Select(p => p.ModifiedSequence).Distinct()) + '\t');

            var t = identifications.SelectMany(p => p.proteinGroups.Select(v => v.ProteinGroupName)).Distinct().OrderBy(p => p);
            if (t.Any())
                sb.Append(string.Join(";", t) + '\t');
            else
                sb.Append("" + '\t');

            sb.Append("" + identifications.First().monoisotopicMass + '\t');
            if (!isMbrFeature)
                sb.Append("" + identifications.First().ms2RetentionTimeInMinutes + '\t');
            else
                sb.Append("" + '\t');
            sb.Append("" + identifications.First().precursorChargeState + '\t');
            sb.Append("" + ClassExtensions.ToMz(identifications.First().monoisotopicMass, identifications.First().precursorChargeState) + '\t');
            sb.Append("" + intensity + "\t");

            if (apex != null)
            {
                sb.Append("" + isotopicEnvelopes.Select(p => p.retentionTime).Min() + "\t");
                sb.Append("" + apex.retentionTime + "\t");
                sb.Append("" + isotopicEnvelopes.Select(p => p.retentionTime).Max() + "\t");

                sb.Append("" + apex.indexedPeak.mz + "\t");
                sb.Append("" + apex.chargeState + "\t");
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

            sb.Append("" + identifications.Count + "\t");
            sb.Append("" + NumIdentificationsByBaseSeq + "\t");
            sb.Append("" + NumIdentificationsByFullSeq + "\t");
            sb.Append("" + splitRT + "\t");
            sb.Append("" + MassError);

            return sb.ToString();
        }

        #endregion Public Methods
    }
}