using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ClassExtensions = Chemistry.ClassExtensions;
using FlashLFQ.PEP;
using FlashLFQ.Interfaces;
using MathNet.Numerics;
using Easy.Common.Extensions;

namespace FlashLFQ
{
    public class ChromatographicPeak : TraceablePeak<IsotopicEnvelope>
    {
        public double Intensity;
        public double ApexRetentionTime => Apex?.IndexedPeak.RetentionTime ?? -1;
        public readonly SpectraFileInfo SpectraFileInfo;
        public List<IsotopicEnvelope> IsotopicEnvelopes { get; set; }
        public override List<IsotopicEnvelope> ScanOrderedPoints => IsotopicEnvelopes;
        public int ScanCount => IsotopicEnvelopes.Count;
        public double SplitRT { get; private set; }
        public readonly bool IsMbrPeak;
        public double MbrScore;
        public double PpmScore { get; set; }
        public double IntensityScore { get; set; }
        public double RtScore { get; set; }
        public double ScanCountScore { get; set; }
        public double IsotopicDistributionScore { get; set; }
        /// <summary>
        /// Stores the pearson correlation between the apex isotopic envelope and the theoretical isotopic distribution
        /// </summary>
        public double IsotopicPearsonCorrelation => Apex?.PearsonCorrelation ?? -1;
        public double RtPredictionError { get; set; }
        public List<int> ChargeList { get; set; }
        internal double MbrQValue { get; set; }
        public ChromatographicPeakData PepPeakData { get; set; }
        public double? MbrPep { get; set; }
        public override IsotopicEnvelope Apex => _apex;
        private IsotopicEnvelope _apex;
        public List<Identification> Identifications { get; private set; }
        public int NumChargeStatesObserved { get; private set; }
        public int NumIdentificationsByBaseSeq { get; private set; }
        public int NumIdentificationsByFullSeq { get; private set; }
        public double MassError { get; private set; }
        /// <summary>
        /// Bool that describes whether the retention time of this peak was randomized
        /// If true, implies that this peak is a decoy peak identified by the MBR algorithm
        /// </summary>
        public bool RandomRt { get; }
        public bool DecoyPeptide => Identifications.First().IsDecoy;

        public ChromatographicPeak(Identification id, bool isMbrPeak, SpectraFileInfo fileInfo, bool randomRt = false)
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
            RandomRt = randomRt;
        }

        public bool Equals(ChromatographicPeak peak)
        {
            return SpectraFileInfo.Equals(peak.SpectraFileInfo) 
                && Identifications.First().ModifiedSequence.Equals(peak.Identifications.First().ModifiedSequence)
                && ApexRetentionTime == peak.ApexRetentionTime;
        }

        /// <summary>
        /// Cuts the peak in such a way as to ensure that the final peak contains the ms2IdRetention time
        /// </summary>
        /// <param name="ms2IdRetentionTime"> The time at which the MS2 scan matched to a peptide was collected</param>
        /// <param name="integrate"> Passed to CalculateIntensityForThisFeature </param>
        public void CutPeak(double ms2IdRetentionTime, double discriminationFactorToCutPeak = 0.6, bool integrate = false)
        {
            // Find all envelopes associated with the apex peak charge state
            List<IsotopicEnvelope> envelopesForApexCharge = IsotopicEnvelopes.Where(e => e.ChargeState == Apex.ChargeState).ToList();

            if (envelopesForApexCharge.Count == 0)
            {
                Console.WriteLine("Something has gone wrong during the peak-cutting procedure");
                return;
            }
            if(envelopesForApexCharge.Count < 5)
            {
                return;
            }

            // Find the boundaries of the peak, based on the apex charge state envelopes
            var peakSplits = FindPeakSplits(envelopesForApexCharge, envelopesForApexCharge.IndexOf(Apex), discriminationFactorToCutPeak);

            // This is done in a while loop, as it's possible that there are multiple distinct peaks in the trace
            while(peakSplits.IsNotNullOrEmpty())
            {
                // Remove all points outside the peak boundaries (inclusive)
                CutPeak(peakSplits, ms2IdRetentionTime);
                if(IsotopicEnvelopes.Count == 0)
                {
                    Intensity = 0;
                    _apex = null;
                    return;
                }
                
                // Recalculate the intensity of the peak and update the split RT
                CalculateIntensityForThisFeature(integrate);
                // technically, there could be two "SplitRT" values if the peak is split into two separate peaks
                // however, this edge case wasn't handled in the legacy code and won't be handled here
                SplitRT = peakSplits.Count > 0 ? peakSplits.Last().RelativeSeparationValue : 0;

                // Update the list of envelopes for the apex charge state (as some were removed during the peak cutting)
                envelopesForApexCharge = IsotopicEnvelopes.Where(e => e.ChargeState == Apex.ChargeState).ToList();
                peakSplits = FindPeakSplits(envelopesForApexCharge, envelopesForApexCharge.IndexOf(Apex), discriminationFactorToCutPeak);
            }
        }

        /// <summary>
        /// Calculated the intensity of the peak. If Integrate is false, this is equal to the summed intensity of all
        /// m/z peak in the the most intense isotopic envelope of the most intense charge state. If Integrate is true,
        /// intensity is  set as the summed intensity of every m/z peak in every isotopic envelope in every charge state.
        /// </summary>
        /// <param name="integrate"></param>
        public void CalculateIntensityForThisFeature(bool integrate)
        {
            if (IsotopicEnvelopes.Any())
            {
                _apex = IsotopicEnvelopes.MaxBy(p => p.Intensity);

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
                _apex = null;
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
        }

        /// <summary>
        /// Sets two ChromatographicPeak properties: NumIdentificationsByBaseSeq and NumIdentificationsByFullSeq
        /// </summary>
        public void ResolveIdentifications()
        {
            this.NumIdentificationsByBaseSeq = Identifications.Select(v => v.BaseSequence).Distinct().Count();
            this.NumIdentificationsByFullSeq = Identifications.Select(v => v.ModifiedSequence).Distinct().Count();
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
                sb.Append("" + '\t');
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

            sb.Append("" + (IsMbrPeak ? MbrQValue.ToString() : "") + "\t");
            sb.Append("" + (IsMbrPeak ? MbrPep.ToString() : "") + "\t");

            sb.Append("" + Identifications.Count + "\t");
            sb.Append("" + NumIdentificationsByBaseSeq + "\t");
            sb.Append("" + NumIdentificationsByFullSeq + "\t");
            sb.Append("" + SplitRT + "\t");
            sb.Append("" + MassError);
            sb.Append("\t" + DecoyPeptide);
            sb.Append("\t" + RandomRt);

            return sb.ToString();
        }
    }
}