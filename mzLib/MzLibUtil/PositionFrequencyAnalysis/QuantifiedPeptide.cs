using Easy.Common.Extensions;
using System;
using System.Collections.Generic;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    /// <summary>
    /// A class to store information about a quantified peptides sharing the same base sequence.
    /// </summary>
    public class QuantifiedPeptide
    {
        public HashSet<string> FullSequences { get; set; }
        public string BaseSequence { get; set; }
        public QuantifiedProtein ParentProtein { get; set; }
        public int OneBasedStartIndexInProtein { get; set; }
        public Dictionary<int, Dictionary<string, QuantifiedModification>> ModifiedAminoAcidPositions { get; set; }
        public double Intensity { get; set; }

        public QuantifiedPeptide(string fullSequence, int oneBasedStartIndexInProtein = -1, double intensity = 0, string modPattern = null)
        {
            ModifiedAminoAcidPositions = new Dictionary<int, Dictionary<string, QuantifiedModification>>();
            OneBasedStartIndexInProtein = oneBasedStartIndexInProtein; // -1 means that the position in the protein is unknown
            Intensity = intensity;
            FullSequences = new HashSet<string> { fullSequence };
            _SetBaseSequence(fullSequence, modPattern);
            _SetModifications(fullSequence, intensity);
        }

        public void AddFullSequence(string fullSeq, double intensity = 0, string modPattern = null)
        {
            if (BaseSequence.Equals(fullSeq.GetBaseSequenceFromFullSequence()))
            {
                FullSequences.Add(fullSeq);
                Intensity += intensity;
                _SetModifications(fullSeq, intensity); // updating the intensity is done here
            }
            else
            {
                throw new Exception("The base sequence of the peptide being added does not match the base sequence of this peptide.");
            }
        }

        public void MergePeptide(QuantifiedPeptide peptideToMerge)
        {
            if (peptideToMerge == null || peptideToMerge.BaseSequence != BaseSequence)
            {
                throw new Exception("The base sequence of the peptide being added does not match the base sequence of this peptide.");
            }
            foreach (var fullSeq in peptideToMerge.FullSequences)
            {
                FullSequences.Add(fullSeq);
                _SetModifications(fullSeq, peptideToMerge.Intensity); // updating the intensity is done here
            }
            Intensity += peptideToMerge.Intensity;
        }

        private void _SetModifications(string fullSeq, double intensity = 0)
        {
            var mods = fullSeq.ParseModifications();

            if (mods.IsNotNullOrEmpty())
            {
                foreach (var modpos in mods.Keys)
                {
                    var mod = mods[modpos];
                    if (!ModifiedAminoAcidPositions.ContainsKey(modpos))
                    {
                        ModifiedAminoAcidPositions[modpos] = new Dictionary<string, QuantifiedModification>();
                    }

                    if (!ModifiedAminoAcidPositions[modpos].ContainsKey(mod))
                    {
                        var modLocalization = modpos == 0 ? "N-terminus" : modpos == BaseSequence.Length + 1 ? "C-terminus" : BaseSequence[modpos - 1].ToString();
                        ModifiedAminoAcidPositions[modpos][mod] = new QuantifiedModification(mod, modpos, modLocalization: modLocalization, intensity: 0);
                    }
                    ModifiedAminoAcidPositions[modpos][mod].Intensity += intensity;

                    // Maybe should update/pass position in protein from here, too.
                }
            }
        }

        private void _SetBaseSequence(string fullSeq, string modPattern)
        {
            BaseSequence = fullSeq.GetBaseSequenceFromFullSequence(modPattern: modPattern);
        }

        public Dictionary<int, Dictionary<string, QuantifiedModification>> GetModStoichiometryForPeptide()
        {
            var aaModsStoichiometry = ModifiedAminoAcidPositions;

            foreach (var modpos in aaModsStoichiometry)
            {
                foreach (var mod in modpos.Value.Values)
                {
                    mod.Intensity /= Intensity;
                }
            }
            return aaModsStoichiometry;
        }
    }
}
