using Easy.Common.Extensions;
using System.Collections.Generic;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    public class PositionFrequencyAnalysis
    {
        public Dictionary<string, QuantifiedProteinGroup> ProteinGroups { get; private set; }

        //public Dictionary<string, (QuantifiedPeptide QuantifiedPeptide, string ProteinGroups)> Peptides { get; private set; }

        /// <summary>
        /// Calculates the occupancy of post-translational modifications at the peptide level. 
        /// </summary>
        /// <param name="peptides"> A List of Tuples whose entries are ordered as (string FullSequence, string BaseSequence, List<string> ProteinGroups, Intensity) for each peptide.</param>
        /// <returns> A nested dictionary whose key mappings are as follows: string ProteinGroup-> string Protein-> string BaseSequence-> int ModifiedAminoAcidIndex-> string ModificationName-> double Intensity
        /// Note: Each BaseSequence dictionary contains a ModifiedAminoAcidIndex key of -1 that then contains a ModificationName key called "Total" that is used to track the total intensity observed for 
        /// all of the amino acids in that peptide.</returns>
        ///
        public void SetUpQuantificationObjectsFromFullSequences(List<(string fullSeq, List<string> proteinGroups, double intensity)> peptides, Dictionary<string, string> proteinSequences=null)
        {
            ProteinGroups = new Dictionary<string, QuantifiedProteinGroup>();

            // Go through the peptides given
            foreach (var pep in peptides)
            {
                string baseSeq = pep.fullSeq.GetBaseSequenceFromFullSequence();

                // Go through the peptide's protein groups
                foreach (var pg in pep.proteinGroups)
                {
                    // If have not seen that protein group, store it
                    if (!ProteinGroups.ContainsKey(pg))
                    {
                        ProteinGroups[pg] = new QuantifiedProteinGroup(pg);
                    }
                    var proteinGroup = ProteinGroups[pg];

                    // Go through the proteins in each protein group
                    foreach (var proteinName in pg.Split('|'))
                    {
                        // Add the protein to the protein group's dictionary if it has not been added
                        if (!proteinGroup.Proteins.ContainsKey(proteinName))
                        {
                            proteinGroup.Proteins[proteinName] = new QuantifiedProtein(proteinName);
                            if (proteinSequences.IsNotNullOrEmpty() && proteinSequences.ContainsKey(proteinName))
                            {
                                proteinGroup.Proteins[proteinName].Sequence = proteinSequences[proteinName];
                            }
                        }
                        var protein = proteinGroup.Proteins[proteinName];

                        // If the peptide's base sequence has not been seen, add it to the protein's dictionary
                        if (!protein.Peptides.ContainsKey(baseSeq))
                        {
                            protein.Peptides[baseSeq] = new QuantifiedPeptide(pep.fullSeq, intensity: pep.intensity);
                        }
                        else
                        {
                            // If the peptide's base sequence has been seen, add the new full sequence to the existing peptide
                            protein.Peptides[baseSeq].AddFullSequence(pep.fullSeq, intensity: pep.intensity);
                        }
                    }
                }
            }
        }
    }
}
