using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using ClassExtensions = Chemistry.ClassExtensions;
using MassSpectrometry;

namespace FlashLFQ
{
    public class ChromatographicPeak : IEquatable<ChromatographicPeak>
    {
        public double Intensity { get; private set; }
        public IsotopicEnvelope Apex { get; private set; }
        public double ApexRetentionTime => Apex?.IndexedPeak.RetentionTime ?? -1;
        public double IsotopicPearsonCorrelation => Apex?.PearsonCorrelation ?? -1;
        public SpectraFileInfo SpectraFileInfo { get; init; }
        public List<IsotopicEnvelope> IsotopicEnvelopes { get; set; }
        public int ScanCount => IsotopicEnvelopes.Count;
        /// <summary>
        /// This is a legacy field that is used to store the RT of the split valley between two peaks.
        /// It is written to the QuantifiedPeaks output, but unused otherwise and should be removed at some point
        /// </summary>
        public double SplitRT { get; set; }
        public List<int> ChargeList { get; set; }
        public DetectionType DetectionType { get; set; }
        public List<Identification> Identifications { get; private set; }
        public int NumChargeStatesObserved { get; private set; }
        public int NumIdentificationsByBaseSeq { get; private set; }
        public int NumIdentificationsByFullSeq { get; private set; }
        public double MassError { get; private set; }
        public bool DecoyPeptide => Identifications.First().IsDecoy;

        public ChromatographicPeak(Identification id, SpectraFileInfo fileInfo, DetectionType detectionType = DetectionType.MSMS) :
            this(new List<Identification>() { id }, fileInfo, detectionType) { }

        /// <summary>
        /// overloaded constructor for Isobaric_ambiguity peaks. In this case, the peak is identified by multiple identifications
        /// </summary>
        /// <param name="ids"></param>
        /// <param name="isMbrPeak"></param>
        /// <param name="fileInfo"></param>
        /// <param name="randomRt"></param>
        public ChromatographicPeak(List<Identification> ids, SpectraFileInfo fileInfo, DetectionType detectionType)
        { 
            SplitRT = 0;
            NumChargeStatesObserved = 0;
            MassError = double.NaN;

            DetectionType = detectionType; // default to imputed
            Identifications = ids;
            ResolveIdentifications();
            IsotopicEnvelopes = new List<IsotopicEnvelope>();
            SpectraFileInfo = fileInfo;
        }

        public bool Equals(ChromatographicPeak peak)
        {
            return SpectraFileInfo.Equals(peak.SpectraFileInfo) 
                && Identifications.First().ModifiedSequence.Equals(peak.Identifications.First().ModifiedSequence)
                && ApexRetentionTime == peak.ApexRetentionTime;
        }

        public void CalculateIntensityForThisFeature(bool integrate)
        {
            if (IsotopicEnvelopes.Any())
            {
                Apex = IsotopicEnvelopes.MaxBy(p => p.Intensity);

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
                    double massErrorForId = ((ClassExtensions.ToMass(Apex.IndexedPeak.M, Apex.ChargeState) - id.PeakfindingMass) / id.PeakfindingMass) * 1e6;

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
                var thisFeaturesPeaks = new HashSet<IIndexedPeak>(IsotopicEnvelopes.Select(p => p.IndexedPeak));
                this.Identifications = this.Identifications
                    .Union(otherFeature.Identifications)
                    .Distinct()
                    .ToList();
                ResolveIdentifications();
                this.IsotopicEnvelopes.AddRange(otherFeature.IsotopicEnvelopes
                    .Where(p => !thisFeaturesPeaks.Contains(p.IndexedPeak)));
                this.CalculateIntensityForThisFeature(integrate);
            }

            // If any merge happens on Isobaric peak, the detection type should be set to MSMSAmbiguousPeakfinding
            if (DetectionType == DetectionType.IsoTrack_MBR || DetectionType == DetectionType.IsoTrack_MSMS)
            {
                DetectionType = DetectionType.MSMSAmbiguousPeakfinding;
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
        /// <summary>
        /// The ordinary header, plus the three columns describing the isotope peaks each envelope was summed from.
        /// Written instead of <see cref="TabSeparatedHeader"/> when FlashLFQ was run with verbose peaks on.
        /// </summary>
        public static string VerboseTabSeparatedHeader =>
            TabSeparatedHeader + "\tIsotope Peak Intensity\tIsotope Peak m/z\tIsotope Peak RTs";

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("File Name" + "\t");
                sb.Append("Base Sequence" + "\t");
                sb.Append("Full Sequence" + "\t");
                sb.Append("Protein Group" + "\t");
                sb.Append("Organism" + '\t');
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
                sb.Append("PIP Q-Value" + "\t");
                sb.Append("PIP PEP" + "\t");
                sb.Append("PSMs Mapped" + "\t");
                sb.Append("Base Sequences Mapped" + "\t");
                sb.Append("Full Sequences Mapped" + "\t");
                sb.Append("Peak Split Valley RT" + "\t");
                sb.Append("Peak Apex Mass Error (ppm)" + "\t");
                sb.Append("Decoy Peptide" + "\t");
                sb.Append("Random RT");
                return sb.ToString();
            }
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(SpectraFileInfo.FilenameWithoutExtension + "\t");
            sb.Append(string.Join("|", Identifications.Select(p => p.BaseSequence).Distinct()) + '\t');
            sb.Append(string.Join("|", Identifications.Select(p => p.ModifiedSequence).Distinct()) + '\t');

            //The semi-colon here splitting the protein groups requires some explanation
            //During protein parsimony, you can get situations where all peptides are shared between two or more proteins. In other words, there is no unique peptide that could resolve the parsimony.
            //In this case you would see something like P00001 | P00002.

            //That’s the easy part and you already understand that.

            //    Now imagine another scenario where you have some other peptides(that are not in either P00001 or P00002) that give you a second group, like the one above.Let’s call it P00003 | P00004.
            // Everything is still fine her.

            //    Now you have two protein groups each with two proteins. 

            //    Here is where the semi - colon comes in.
            //Imagine you now find a new peptide(totally different from any of the peptides used to create the two original protein groups) that is shared across all four proteins.The original peptides
            //require that two different protein groups exist, but this new peptide could come from either or both.We don’t know. So, the quantification of that peptide must be allowed to be
            //either/ both groups. For this peptide, the protein accession in the output will be P00001 | P00002; P00003 | P00004.

            //    You could see an output that looks like P0000A; P0000B.Here there is only one protein in each protein group(as decided by parsimony).And you have a peptide that is shared.This would
            // not ever be reported as P0000A | P0000B because each protein has a unique peptide that confirms its existence.

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
            sb.Append(DetectionType == DetectionType.MSMS ? Identifications.First().Ms2RetentionTimeInMinutes.ToString(CultureInfo.InvariantCulture) + '\t' : '\t');
            sb.Append(Identifications.First().PrecursorChargeState.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(ClassExtensions.ToMz(Identifications.First().MonoisotopicMass, Identifications.First().PrecursorChargeState).ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(Intensity.ToString(CultureInfo.InvariantCulture) + "\t");

            if (Apex != null)
            {
                sb.Append(IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime).ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(Apex.IndexedPeak.RetentionTime.ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime).ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(Apex.IndexedPeak.M.ToString(CultureInfo.InvariantCulture) + "\t");
                sb.Append(Apex.ChargeState.ToString(CultureInfo.InvariantCulture) + "\t");
            }
            else
            {
                sb.Append("-" +"\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
                sb.Append("-" + "\t");
            }

            sb.Append(NumChargeStatesObserved.ToString(CultureInfo.InvariantCulture) + "\t");

            // temporary way to distinguish between MBR, MBR_IsoTrack, IsoTrack_Ambiguous, MSMSAmbiguousPeakfinding and MSMS peaks
            switch (this.DetectionType)
            {
                case DetectionType.IsoTrack_MBR:
                    sb.Append("MBR_IsoTrack" + "\t");
                    break;
                case DetectionType.IsoTrack_Ambiguous:
                    sb.Append("IsoTrack_Ambiguous" + "\t");
                    break;
                case DetectionType.MSMSAmbiguousPeakfinding:
                    sb.Append("MSMSAmbiguousPeakfinding" + "\t");
                    break;
                default:
                    sb.Append("MSMS" + "\t");
                    break;
            }

            // MBR Exclusive fields
            sb.Append('\t');
            sb.Append('\t');

            sb.Append(Identifications.Count.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(NumIdentificationsByBaseSeq.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(NumIdentificationsByFullSeq.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(SplitRT.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append(MassError.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append( DecoyPeptide.ToString(CultureInfo.InvariantCulture) + "\t");
            sb.Append("False"); // Because this isn't an MBR peak, the Random RT Field will always be false

            return sb.ToString();
        }

        /// <summary>
        /// The peak as a row of QuantifiedPeaks.tsv, optionally carrying the three isotope peak columns.
        /// </summary>
        /// <param name="verbose"> append the columns described by <see cref="VerboseTabSeparatedHeader"/> </param>
        public string ToString(bool verbose)
        {
            return verbose ? ToString() + "\t" + GetIsotopeInformation() : ToString();
        }

        /// <summary>
        /// The three isotope peak columns as one tab separated fragment: the intensity of every isotope peak, its
        /// m/z, and the retention times the columns of both are indexed by.
        ///
        /// The intensity and m/z columns hold one bracketed list per isotope peak of per charge state, labelled
        /// i&lt;isotope&gt;z&lt;charge&gt;, and each list holds one entry per retention time in the third column - a
        /// dash where that isotope was not observed at that time. So the three columns together are a matrix of the
        /// envelope over its elution, which the summed intensity in the ordinary columns collapses away.
        ///
        /// Empty columns for a peak with no verbose envelopes, which is every peak when FlashLFQ was not run with
        /// verbose peaks on.
        /// </summary>
        internal string GetIsotopeInformation()
        {
            List<VerboseIsotopicEnvelope> envelopes = IsotopicEnvelopes
                .OfType<VerboseIsotopicEnvelope>()
                .ToList();

            if (envelopes.Count == 0)
            {
                return "\t\t";
            }

            // one column per retention time the peak was observed at, ascending, so that the same index means the
            // same scan in every isotope's list
            List<double> retentionTimes = envelopes
                .Select(e => e.RetentionTime)
                .Distinct()
                .OrderBy(rt => rt)
                .ToList();
            var columnByRetentionTime = retentionTimes
                .Select((rt, index) => (rt, index))
                .ToDictionary(pair => pair.rt, pair => pair.index);

            // and one row per isotope peak of per charge state, ordered so the columns read the same way every run
            var rows = new Dictionary<(int IsotopeNumber, int ChargeState), IIndexedPeak[]>();
            foreach (VerboseIsotopicEnvelope envelope in envelopes)
            {
                int column = columnByRetentionTime[envelope.RetentionTime];
                foreach (KeyValuePair<int, IIndexedPeak> isotopePeak in envelope.PeakDictionary)
                {
                    var key = (isotopePeak.Key, envelope.ChargeState);
                    if (!rows.TryGetValue(key, out IIndexedPeak[] row))
                    {
                        row = new IIndexedPeak[retentionTimes.Count];
                        rows[key] = row;
                    }
                    row[column] = isotopePeak.Value;
                }
            }

            var orderedRows = rows
                .OrderBy(row => row.Key.ChargeState)
                .ThenBy(row => row.Key.IsotopeNumber)
                .ToList();

            var intensities = new StringBuilder();
            var mzValues = new StringBuilder();
            foreach (var row in orderedRows)
            {
                string label = VerboseIsotopicEnvelope.GetIsotopePeakName(row.Key);
                intensities.Append(FormatRow(label, row.Value, peak => peak.Intensity));
                mzValues.Append(FormatRow(label, row.Value, peak => peak.M));
            }

            return "\"" + intensities + "\"\t\"" + mzValues + "\"\t["
                + string.Join(", ", retentionTimes.Select(rt => rt.ToString(CultureInfo.InvariantCulture))) + "]";
        }

        /// <summary>
        /// One isotope's row as [label: v1, v2, -, v4];, with a dash standing in for a retention time at which that
        /// isotope was not observed.
        /// </summary>
        private static string FormatRow(string label, IIndexedPeak[] row, Func<IIndexedPeak, float> value)
        {
            return "[" + label + ": "
                + string.Join(", ", row.Select(peak =>
                    peak == null ? "-" : value(peak).ToString(CultureInfo.InvariantCulture)))
                + "];";
        }
    }
}