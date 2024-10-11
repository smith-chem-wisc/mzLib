using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;

namespace MzLibUtil
{
    public class PositionFrequencyAnalysis
    {
        // TODO: Write PeptideToProteinPTMOccupancy which takes a dictionary of ProteinAccession:ProteinSequence for determining index offsets to the peptide modification indices. 

        /// <summary>
        /// Calculates the occupancy of post-translational modifications at the peptide level. 
        /// </summary>
        /// <param name="peptides"> A List of Tuples whose entries are ordered as (string FullSequence, string BaseSequence, HashSet<string> ProteinGroups, Intensity) for each peptide.</param>
        /// <param name="IncludeNTerminus"> If true, the index of modifications at the N-terminus will be 0 (zero-based indexing). Otherwise, it is the index of the first amino acid (one-based indexing).</param>
        /// <param name="IncludeCTerminus"> If true, the index of modifications at the C-terminus will be one more than the index of the last amino acid. Otherwise, it is the index of the last amino acid.</param>
        /// <returns> A nested dictionary whose key mappings are as follows: string ProteinGroup-> string Protein-> string BaseSequence-> int ModifiedAminoAcidIndex-> string ModificationName-> double Intensity
        /// Note: Each BaseSequence dictionary contains a ModifiedAminoAcidIndex key of -1 that then contains a ModificationName key called "Total" that is used to track the total intensity observed for 
        /// all of the amino acids in that peptide.</returns>
        public static Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<int, Dictionary<string, double>>>>> PeptidePTMOccupancy(List<Tuple<string, string, List<string>, double>> peptides, bool IncludeNTerminus=true, bool IncludeCTerminus=true)
        {
            var occupancy = new Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<int, Dictionary<string, double>>>>>();
            
            // Go through the peptides given
            foreach (var peptide in peptides)
            {
                double peptideIntensity = peptide.Item4;
                
                string baseSeq;
                if (!peptide.Item2.IsNotNullOrEmpty())
                {
                    baseSeq = new string(peptide.Item1.ToCharArray());
                    ClassExtensions.RemoveSpecialCharacters(ref baseSeq, @"", @"[.*?]");
                }
                else { baseSeq = peptide.Item2; }

                // Go through the peptide's protein groups
                foreach(var pg in peptide.Item3)
                {
                    // If have not seen that protein group, store it
                    if (!occupancy.ContainsKey(pg))
                    {
                        occupancy.Add(pg, new Dictionary<string, Dictionary<string, Dictionary<int, Dictionary<string, double>>>>());
                    }

                    // Go through the proteins in each protein group
                    foreach (var protein in pg.Split('|'))
                    {
                        // Add the protein to the protein group's dictionary if it has not been added
                        if (!occupancy[pg].ContainsKey(protein))
                        {
                            occupancy[pg].Add(protein, new Dictionary<string, Dictionary<int, Dictionary<string, double>>>());
                        }
                        // If the peptide's base sequence has not been seen, add it to the protein's dictionary
                        if (!occupancy[pg][protein].ContainsKey(baseSeq))
                        {
                            occupancy[pg][protein].Add(baseSeq, new Dictionary<int, Dictionary<string, double>>());
                            occupancy[pg][protein][baseSeq].Add(-1, new Dictionary<string, double>());
                            occupancy[pg][protein][baseSeq][-1].Add("Total", 0);
                        }
                        // Increase the total intensity of the peptide base sequence to track the total intensity of all amino acids in that sequence
                        occupancy[pg][protein][baseSeq][-1]["Total"] += peptideIntensity;

                        // Want both arguments passed here to be true if need to later filter out peptide terminal mods that are not protein terminal mods 
                        Dictionary<int, List<string>> peptideMods = peptide.Item1.ParseModifications(IncludeNTerminus, IncludeCTerminus); 
                        // Go through the modified positions found froum the full sequence
                        foreach (var modpos in peptideMods)
                        {
                            // If that position has not been recorded as containing a modification, add it to the base sequence's dictonary
                            if (!occupancy[pg][protein][baseSeq].ContainsKey(modpos.Key))
                            {
                                occupancy[pg][protein][baseSeq].Add(modpos.Key, new Dictionary<string, double>());
                            }
                            // Go through the modifications found at a modified amino acid index
                            foreach (var mod in modpos.Value)
                            {
                                //If the name of that modification has not been seen, record that modification in the index's dictionary with an intensity of 0
                                if (!occupancy[pg][protein][baseSeq][modpos.Key].ContainsKey(mod))
                                {
                                    occupancy[pg][protein][baseSeq][modpos.Key].Add(mod, 0);
                                }
                                // Increase the intensity of the modification by the intensity of the peptide
                                occupancy[pg][protein][baseSeq][modpos.Key][mod] += peptideIntensity;
                            }
                        }
                    }
                }
            }
            return occupancy;
        }
    }
}
