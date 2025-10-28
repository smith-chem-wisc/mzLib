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
        public int ZeroBasedStartIndexInProtein { get; set; }

        /// <summary>
        /// Dictionary mapping zero-based amino acid positions in the peptide to dictionaries of
        /// modification IDs and their corresponding QuantifiedModification objects. This property 
        /// stores ALL of the modifications observed for this peptide across all full sequences.
        /// </summary>
        public Dictionary<int, Dictionary<string, QuantifiedModification>> ModifiedAminoAcidPositions { get; set; }
        public double Intensity { get; set; }

        /// <summary>
        /// Constructor for a QuantifiedPeptide object. The base sequence and modifications are parsed from the full sequence.
        /// </summary>
        /// <param name="fullSequence"></param>
        /// <param name="zeroBasedStartIndexInProtein"></param>
        /// <param name="intensity"></param>
        /// <param name="modPattern"></param>
        public QuantifiedPeptide(string fullSequence, int zeroBasedStartIndexInProtein = -1, double intensity = 0, string modPattern = null)
        {
            ModifiedAminoAcidPositions = new Dictionary<int, Dictionary<string, QuantifiedModification>>();
            ZeroBasedStartIndexInProtein = zeroBasedStartIndexInProtein; // -1 means that the position in the protein is unknown
            Intensity = intensity;
            FullSequences = new HashSet<string> { fullSequence };
            _SetBaseSequence(fullSequence, modPattern);
            _SetModifications(fullSequence, intensity);
        }

        /// <summary>
        /// Adds a new full sequence to the peptide, updating modifications and intensity accordingly.
        /// </summary>
        /// <param name="fullSeq"></param>
        /// <param name="intensity"></param>
        /// <param name="modPattern"></param>
        /// <exception cref="Exception"></exception>
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

        /// <summary>
        /// Merges another QuantifiedPeptide object into this one, combining their full sequences and intensities.
        /// </summary>
        /// <param name="peptideToMerge"></param>
        /// <exception cref="Exception"></exception>
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
                        ModifiedAminoAcidPositions[modpos][mod] = new QuantifiedModification(mod, modpos, intensity: 0);
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

        /// <summary>
        /// Returns the modification stoichiometry for this peptide as a dictionary mapping
        /// zero-based amino acid positions in the peptide to dictionaries of modification IDs and their corresponding
        /// QuantifiedModification objects with normalized intensities (i.e., divided by the total peptide intensity).
        /// </summary>
        /// <returns></returns>
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
