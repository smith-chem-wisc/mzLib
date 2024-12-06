using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using Easy.Common.Extensions;

namespace MzLibUtil
{
    // Should this have all of the parent data (i.e. protein group, protein, peptide, peptide position)? Unnecessary for now, but probably useful later.
    public class UtilModification
    {
        public string IdWithMotif { get; set; }
        public int PeptidePositionZeroIsNTerminus { get; set; } //NEED TO ENFORCE THIS EVERYWHERE OR CHECK IF ZERO OR ONE


        public double Intensity { get; set; }

        public UtilModification(string name, int position, double intensity)
        {
            IdWithMotif = name;
            PeptidePositionZeroIsNTerminus = position;
            Intensity = intensity;
        }

    }
    public class UtilPeptide
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; set; }
        public UtilProtein ParentProtein { get; set; }
        public int IndexInProtein { get; set; }
        public Dictionary<int, Dictionary<string, UtilModification>> ModifiedAminoAcidPositions { get; set; }
        public double Intensity { get; set; } 

        public UtilPeptide(string fullSequence, Dictionary<int, Dictionary<string, UtilModification>> mods = null) 
        {
            FullSequence = fullSequence;
            ModifiedAminoAcidPositions = mods.IsNotNullOrEmpty() ? mods : new Dictionary<int, Dictionary<string, UtilModification>>();
            SetBaseSequence();
        }
        public void SetBaseSequence(string modPattern = @"\[(.+?)\](?<!\[I+\])")
        {
            Regex regexSpecialChar = new(modPattern);
            BaseSequence = regexSpecialChar.Replace(FullSequence, @"");
        }
        public void AddModifications(Dictionary<int, string> mods)
        {
            throw new NotImplementedException();
        }
        public void PeptideToProteinPositions(int offset=0)
        {
            offset = offset != 0 ? offset : ParentProtein.Sequence.IndexOf(BaseSequence);
            foreach (var modpos in ModifiedAminoAcidPositions.Keys)
            {
                int positionInProtein = modpos + offset;
                Dictionary<string, UtilModification> mods = ModifiedAminoAcidPositions[modpos];
                foreach (var mod in mods.Values)
                {
                    mod.PeptidePositionZeroIsNTerminus = positionInProtein;
                }
                ModifiedAminoAcidPositions.Add(positionInProtein, mods);
                ModifiedAminoAcidPositions.Remove(modpos);
            }
        }
    }
    
    public class UtilProtein
    {
        public string Name { get; set; }
        public string Sequence { get; set; }
        public Dictionary<string, UtilPeptide> Peptides { get; set; }

        public UtilProtein(string name, Dictionary<string, UtilPeptide> peptides=null)
        {
            Name = name;
            if (peptides != null) Peptides = peptides;
            else Peptides= new Dictionary<string, UtilPeptide>();
        }
    }

    public class UtilProteinGroup
    {
        public string Name { get; set;}
        public Dictionary<string, UtilProtein> Proteins {  get; set; }

        public UtilProteinGroup(string name, Dictionary<string, UtilProtein> proteins = null)
        {
            Name = name;
            if (proteins != null) Proteins = proteins;
            else Proteins= new Dictionary<string, UtilProtein>();
        }
    }
    public class PositionFrequencyAnalysis
    {
        // TODO: Write PeptideToProteinPTMOccupancy which takes a dictionary of ProteinAccession:ProteinSequence for determining index offsets to the peptide modification indices. 

        /// <summary>
        /// Calculates the occupancy of post-translational modifications at the peptide level. 
        /// </summary>
        /// <param name="peptides"> A List of Tuples whose entries are ordered as (string FullSequence, string BaseSequence, List<string> ProteinGroups, Intensity) for each peptide.</param>
        /// <param name="modOnNTerminus"> If true, the index of modifications at the N-terminus will be 0 (zero-based indexing). Otherwise, it is the index of the first amino acid (one-based indexing).</param>
        /// <param name="modOnCTerminus"> If true, the index of modifications at the C-terminus will be one more than the index of the last amino acid. Otherwise, it is the index of the last amino acid.</param>
        /// <returns> A nested dictionary whose key mappings are as follows: string ProteinGroup-> string Protein-> string BaseSequence-> int ModifiedAminoAcidIndex-> string ModificationName-> double Intensity
        /// Note: Each BaseSequence dictionary contains a ModifiedAminoAcidIndex key of -1 that then contains a ModificationName key called "Total" that is used to track the total intensity observed for 
        /// all of the amino acids in that peptide.</returns>
        /// 

        public Dictionary<string, UtilProteinGroup> Occupancy { get; private set; }
        public string OccupancyLevel { get; private set; }

        
        public void PeptidePTMOccupancy(List<Tuple<string, string, List<string>, double>> peptides, bool modOnNTerminus = true, bool modOnCTerminus = true)
        {
            var proteinGroups = new Dictionary<string, UtilProteinGroup>();
            
            // Go through the peptides given
            foreach (var pep in peptides)
            {
                string fullSeq = pep.Item1;
                string baseSeq = pep.Item2.IsNotNullOrEmpty() ? pep.Item2 : new string(pep.Item1.ToCharArray());
                ClassExtensions.RemoveSpecialCharacters(ref baseSeq, @"", @"\[(.+?)\](?<!\[I+\])"); // in case it is null or empty and we need to get the base sequence from the full sequence
                List<string> pgs = pep.Item3;
                double peptideIntensity = pep.Item4;

                // Go through the peptide's protein groups
                foreach (var pgName in pgs)
                {
                    // If have not seen that protein group, store it
                    if (!proteinGroups.ContainsKey(pgName))
                    {
                        proteinGroups[pgName] = new UtilProteinGroup(pgName);
                    }
                    var proteinGroup = proteinGroups[pgName];

                    // Go through the proteins in each protein group
                    foreach (var proteinName in pgName.Split('|'))
                    {
                        // Add the protein to the protein group's dictionary if it has not been added
                        if (!proteinGroup.Proteins.ContainsKey(proteinName))
                        {
                            proteinGroup.Proteins[proteinName] = new UtilProtein(proteinName);
                        }
                        var protein = proteinGroup.Proteins[proteinName];

                        // If the peptide's base sequence has not been seen, add it to the protein's dictionary
                        if (!protein.Peptides.ContainsKey(baseSeq))
    {
                            protein.Peptides[baseSeq] = new UtilPeptide(fullSeq);
                            protein.Peptides[baseSeq].Intensity = 0;
                        }

                        // Increase the total intensity of the peptide base sequence to track the total intensity of all amino acids in that sequence
                        protein.Peptides[baseSeq].Intensity += peptideIntensity;
                        var peptide = protein.Peptides[baseSeq];

                        // Want both arguments passed here to be true if need to later filter out peptide terminal mods that are not protein terminal mods 
                        Dictionary<int, List<string>> peptideMods = fullSeq.ParseModifications(modOnNTerminus, modOnCTerminus);
                        // Go through the modified positions found froum the full sequence
                        foreach (var modpos in peptideMods)
        {
                            // If that position has not been recorded as containing a modification, add it to the base sequence's dictonary
                            if (!peptide.ModifiedAminoAcidPositions.ContainsKey(modpos.Key))
                            {
                                peptide.ModifiedAminoAcidPositions[modpos.Key] = new Dictionary<string, UtilModification>();
                            }
                            var modifiedPosition = peptide.ModifiedAminoAcidPositions[modpos.Key];
            
                            // Go through the modifications found at a modified amino acid index
                            foreach (var mod in modpos.Value)
                            {
                                //If the name of that modification has not been seen, record that modification in the index's dictionary with an intensity of 0
                                if (!modifiedPosition.ContainsKey(mod))
            {
                                    modifiedPosition[mod] = new UtilModification(mod, modpos.Key, 0);
                                }
                                // Increase the intensity of the modification by the intensity of the peptide
                                modifiedPosition[mod].Intensity += peptideIntensity;
                            }
                        }
                    }
                }
            }
            Occupancy = proteinGroups;
            OccupancyLevel = "peptide";
        }

        public void PeptideToProteinPTMOccupancy(Dictionary<string, string> proteinSequences) // combine this to previous method.
        {
            foreach (var pg in Occupancy.Keys) 
            {
                UtilProteinGroup proteinGroup = Occupancy[pg];
                foreach (var prot in proteinGroup.Proteins.Keys)
                {
                    UtilProtein protein = proteinGroup.Proteins[prot];
                    foreach (var pep in protein.Peptides.Keys)
                    {
                        UtilPeptide peptide = protein.Peptides[pep];
                        peptide.ParentProtein = protein;
                        peptide.PeptideToProteinPositions();
                    }
                }
            }
            OccupancyLevel = "protein";
        }
    }
}
