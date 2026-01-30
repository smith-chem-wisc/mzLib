using Easy.Common.Extensions;
using System.Collections.Generic;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    /// <summary>
    /// Handles analysis and organization of protein group quantification from peptide records.
    /// </summary>
    public class PositionFrequencyAnalysis
    {
        /// <summary>
        /// Dictionary mapping protein group names to their quantification data.
        /// </summary>
        public Dictionary<string, QuantifiedProteinGroup> ProteinGroups { get; private set; }

        /// <summary>
        /// Populates protein groups with their respective proteins and peptides from a list of quantifide peptide records. 
        /// The resulting protein groups are stored in the ProteinGroups property with the protein group name strings as keys.
        /// </summary>
        /// <param name="peptides"> A list of QuantifiedPeptideRecord, which store a peptide's full sequence, mapped protein groupsm and intensity.</param>
        /// <param name="proteinSequences"> An optional dictionary of protein sequences to use for mapping peptides to proteins. 
        /// If not provided, the protein sequences will be left null in the QuantifiedProtein objects. However, this parameter should not be null if what we want
        /// is a protein stoichiometry, since it is needed to align the peptides to the parent protein.</param>"
        public void SetUpQuantificationFromQuantifiedPeptideRecords(List<QuantifiedPeptideRecord> peptides, Dictionary<string, string> proteinSequences=null)
        {
            ProteinGroups = new Dictionary<string, QuantifiedProteinGroup>();
            foreach (var peptide in peptides)
            {
                // Iterate through the peptide's protein groups in case it is a shared peptide protein groups.
                // We want to map the peptide separately to each protein group it belongs to, primarily due to 
                // each protein group is reported  separately in MetaMorpheus.
                foreach (var pg in peptide.ProteinGroups)
                {
                    // If have not seen that protein group, store it
                    if (!ProteinGroups.ContainsKey(pg))
                    {
                        ProteinGroups[pg] = new QuantifiedProteinGroup(pg);
                    }
                    var proteinGroup = ProteinGroups[pg];

                    foreach (var proteinName in pg.SplitProteinAccessions())
                    {
                        // Add the protein to the protein group's dictionary if it has not been added
                        if (!proteinGroup.Proteins.ContainsKey(proteinName))
                        {
                            proteinGroup.Proteins[proteinName] = new QuantifiedProtein(proteinName);
                            if (proteinSequences.IsNotNullOrEmpty() && proteinSequences.TryGetValue(proteinName, out var sequence))
                            {
                                proteinGroup.Proteins[proteinName].Sequence = sequence;
                            }
                        }
                        var protein = proteinGroup.Proteins[proteinName];

                        // If the peptide's base sequence has not been seen, add it to the protein's dictionary; otherwise, update the existing entry
                        if (!protein.Peptides.TryGetValue(peptide.BaseSequence, out var quantifiedPeptide))
                        {
                            quantifiedPeptide = new QuantifiedPeptide(peptide.FullSequence, intensity: peptide.Intensity);
                            protein.Peptides[peptide.BaseSequence] = quantifiedPeptide;
                        }
                        else
                        {
                            // If the peptide's base sequence has been seen, add the new full sequence to the existing peptide
                            quantifiedPeptide.AddFullSequence(peptide.FullSequence, intensity: peptide.Intensity);
                        }
                    }
                }
            }
        }
    }
}
