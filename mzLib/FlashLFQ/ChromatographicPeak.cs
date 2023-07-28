using Chemistry;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public class ChromatographicPeak
    {
        public double Intensity;
        public readonly SpectraFileInfo SpectraFileInfo;
        public List<IsotopicEnvelope> IsotopicEnvelopes;
        public double SplitRT;
        public readonly bool IsMbrPeak;
        public double MbrScore;

        public ChromatographicPeak(Identification id, bool isMbrPeak, SpectraFileInfo fileInfo)
        {
            SplitRT = 0;
            NumChargeStatesObserved = 0;
            MassError = double.NaN;
            NumIdentificationsByBaseSeq = 1;
            NumIdentificationsByFullSeq = 1;
            Identifications = new List<Identification>() { id };
            IsotopicEnvelopes = new List<IsotopicEnvelope>();
            IsMbrPeak = isMbrPeak;
            SpectraFileInfo = fileInfo;
        }

        public IsotopicEnvelope Apex { get; private set; }
        public List<Identification> Identifications { get; private set; }
        public int NumChargeStatesObserved { get; private set; }
        public int NumIdentificationsByBaseSeq { get; private set; }
        public int NumIdentificationsByFullSeq { get; private set; }
        public double MassError { get; private set; }
        /// <summary>
        /// Expected retention time for MBR acceptor peaks (mean)
        /// </summary>
        public double? RtHypothesis { get; private set; }
        /// <summary>
        /// Std. Dev of retention time differences between MBR acceptor file and donor file, used if # calibration points < 6
        /// </summary>
        public double? RtStdDev { get; private set;  }
        /// <summary>
        /// Interquartile range of retention time differences between MBR acceptor file and donor file, used if # calibration points >= 6
        /// </summary>
        public double? RtInterquartileRange { get; private set; }

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
                sb.Append("MBR Score" + "\t");
                sb.Append("PSMs Mapped" + "\t");
                sb.Append("Base Sequences Mapped" + "\t");
                sb.Append("Full Sequences Mapped" + "\t");
                sb.Append("Peak Split Valley RT" + "\t");
                sb.Append("Peak Apex Mass Error (ppm)" + "\t");
                //sb.Append("Timepoints");
                return sb.ToString();
            }
        }

        /// <summary>
        /// Sets retention time information for a given peak. Used for MBR peaks
        /// </summary>
        /// <param name="rtHypothesis"> Expected retention time for peak, based on alignment between a donor and acceptor file </param>
        /// <param name="rtStdDev"> Standard deviation in the retention time differences between aligned peaks </param>
        /// <param name="rtInterquartileRange"> Interquartile range og the retention time differences between aligned peaks</param>
        internal void SetRtWindow(double rtHypothesis, double? rtStdDev, double? rtInterquartileRange)
        {
            RtHypothesis = rtHypothesis;
            RtStdDev = rtStdDev;
            RtInterquartileRange = rtInterquartileRange;
        }

        public void CalculateIntensityForThisFeature(bool integrate)
        {
            if (IsotopicEnvelopes.Any())
            {
                double maxIntensity = IsotopicEnvelopes.Max(p => p.Intensity);
                Apex = IsotopicEnvelopes.First(p => p.Intensity == maxIntensity);

                if (integrate)
                {
                    Intensity = IsotopicEnvelopes.Sum(p => p.Intensity);
                }
                else
                {
                    Intensity = Apex.Intensity;
                }

                MassError = double.NaN;

                foreach (Identification id in Identifications)
                {
                    double massErrorForId = ((ClassExtensions.ToMass(Apex.IndexedPeak.Mz, Apex.ChargeState) - id.PeakfindingMass) / id.PeakfindingMass) * 1e6;

                    if (double.IsNaN(MassError) || Math.Abs(massErrorForId) < Math.Abs(MassError))
                    {
                        MassError = massErrorForId;
                    }
                }

                NumChargeStatesObserved = IsotopicEnvelopes.Select(p => p.ChargeState).Distinct().Count();
            }
            else
            {
                Intensity = 0;
                MassError = double.NaN;
                NumChargeStatesObserved = 0;
                Apex = null;
            }
        }

        /// <summary>
        /// Merges ChromatographicPeaks by combining Identifications and IsotopicEnvelopes,
        /// then recalculates feature intensity.
        /// </summary>
        /// <param name="otherFeature"> Peak to be merged in. This peak is not modified</param>
        public void MergeFeatureWith(ChromatographicPeak otherFeature, bool integrate)
        {
            if (otherFeature != this)
            {
                var thisFeaturesPeaks = new HashSet<IndexedMassSpectralPeak>(IsotopicEnvelopes.Select(p => p.IndexedPeak));
                this.Identifications = this.Identifications
                    .Union(otherFeature.Identifications)
                    .Distinct()
                    .OrderBy(p => p.PosteriorErrorProbability).ToList();
                ResolveIdentifications();
                this.IsotopicEnvelopes.AddRange(otherFeature.IsotopicEnvelopes
                    .Where(p => !thisFeaturesPeaks.Contains(p.IndexedPeak)));
                this.CalculateIntensityForThisFeature(integrate);
            }
        }

        /// <summary>
        /// Sets two ChromatographicPeak properties: NumIdentificationsByBaseSeq and NumIdentificationsByFullSeq
        /// </summary>
        public void ResolveIdentifications()
        {
            this.NumIdentificationsByBaseSeq = Identifications.Select(v => v.BaseSequence).Distinct().Count();
            this.NumIdentificationsByFullSeq = Identifications.Select(v => v.ModifiedSequence).Distinct().Count();
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
            }
            else
            {
                sb.Append("" + '\t');
            }

            sb.Append("" + Identifications.First().MonoisotopicMass + '\t');
            if (!IsMbrPeak)
            {
                sb.Append("" + Identifications.First().Ms2RetentionTimeInMinutes + '\t');
            }
            else
            {
                sb.Append("" + '\t');
            }

            sb.Append("" + Identifications.First().PrecursorChargeState + '\t');
            sb.Append("" + ClassExtensions.ToMz(Identifications.First().MonoisotopicMass, Identifications.First().PrecursorChargeState) + '\t');
            sb.Append("" + Intensity + "\t");

            if (Apex != null)
            {
                sb.Append("" + IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime) + "\t");
                sb.Append("" + Apex.IndexedPeak.RetentionTime + "\t");
                sb.Append("" + IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime) + "\t");

                sb.Append("" + Apex.IndexedPeak.Mz + "\t");
                sb.Append("" + Apex.ChargeState + "\t");
            }
            else
            {
                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");

                sb.Append("" + "-" + "\t");
                sb.Append("" + "-" + "\t");
            }

            sb.Append("" + NumChargeStatesObserved + "\t");

            if (IsMbrPeak)
            {
                sb.Append("" + "MBR" + "\t");
            }
            else
            {
                sb.Append("" + "MSMS" + "\t");
            }

            sb.Append("" + (IsMbrPeak ? MbrScore.ToString() : "") + "\t");

            sb.Append("" + Identifications.Count + "\t");
            sb.Append("" + NumIdentificationsByBaseSeq + "\t");
            sb.Append("" + NumIdentificationsByFullSeq + "\t");
            sb.Append("" + SplitRT + "\t");
            sb.Append("" + MassError);
            
            return sb.ToString();
        }
    }
}