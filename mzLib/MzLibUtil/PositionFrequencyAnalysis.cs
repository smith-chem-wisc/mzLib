using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using Easy.Common.Extensions;

namespace MzLibUtil
{
    public class UtilModification
    {
        public string IdWithMotif { get; set; }
        public int PeptidePositionZeroIsNTerminus { get; set; }
        public int ProteinPositionZeroIsNTerminus { get; set; }
        public double Intensity { get; set; }

        public UtilModification(string name, int positionInPeptide, int? positionInProtein=null, double intensity=0)
        {
            IdWithMotif = name;
            PeptidePositionZeroIsNTerminus = positionInPeptide;
            ProteinPositionZeroIsNTerminus = positionInProtein ?? -1;
            Intensity = intensity;
        }
    }

    public class UtilPeptide
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; set; }
        public UtilProtein ParentProtein { get; set; }
        public int OneBasedStartIndexInProtein { get; set; }
        public Dictionary<int, Dictionary<string, UtilModification>> ModifiedAminoAcidPositions { get; set; }
        public double Intensity { get; set; }

        public UtilPeptide(string fullSequence, Dictionary<int, Dictionary<string, UtilModification>> mods = null, int oneBasedStartIndexInProtein = 1, double intensity = 0) 
        {
            FullSequence = fullSequence;
            ModifiedAminoAcidPositions = mods.IsNotNullOrEmpty() ? mods : new Dictionary<int, Dictionary<string, UtilModification>>();
            OneBasedStartIndexInProtein = oneBasedStartIndexInProtein;
            Intensity = intensity;
            SetBaseSequence();
        }
        public void SetBaseSequence(string modPattern = @"\[(.+?)\](?<!\[I+\])")
        {
            Regex regexSpecialChar = new(modPattern);
            BaseSequence = regexSpecialChar.Replace(FullSequence, @"");
        }

        public Dictionary<int, Dictionary<string, UtilModification>> GetModStoichiometryFromPeptideMods()
        {
            var aaModsStoichiometry = ModifiedAminoAcidPositions;
            foreach (var modpos in aaModsStoichiometry)
            {
                foreach (var mod in modpos.Value.Values)
                {
                    mod.Intensity = mod.Intensity / Intensity;
                }
            }
            return aaModsStoichiometry;
        }
    }
    
    public class UtilProtein
    {
        public string Accession { get; set; }
        public string Sequence { get; set; }
        public Dictionary<string, UtilPeptide> Peptides { get; set; }
        public Dictionary<int, Dictionary<string, UtilModification>> ModifiedAminoAcidPositionsInProtein { get; set; }
        public Dictionary<int, List<UtilPeptide>> PeptidesByProteinPosition { get; set; }

        public UtilProtein(string accession, Dictionary<string, UtilPeptide> peptides=null)
        {
            Accession = accession;
            if (peptides != null) Peptides = peptides;
            else Peptides= new Dictionary<string, UtilPeptide>();
        }

        public void SetProteinModsFromPeptides()
        {
            ModifiedAminoAcidPositionsInProtein = new Dictionary<int, Dictionary<string, UtilModification>>();
            PeptidesByProteinPosition = new Dictionary<int, List<UtilPeptide>>();

            foreach (var peptide in Peptides.Values)
            {
                foreach (var modpos in peptide.ModifiedAminoAcidPositions)
                {
                    var modPositionInProtein = modpos.Key + peptide.OneBasedStartIndexInProtein - 1;
                    if (Sequence.IsNotNullOrEmpty()) // if the protein sequence is known, ignore terminal modifications that are not protein terminal modifications
                    {
                        if ((modPositionInProtein != 0 && modpos.Key == 0) // if the mod is at the N-terminus of the peptide, but not the protein
                            || (modPositionInProtein != Sequence.Length + 1 && modpos.Key == peptide.BaseSequence.Length + 1)) // if the mod is at the C-terminus of the peptide, but not the protein
                        {
                            continue;
                        }
                    }
                    else // if the protein sequence is not known, ignore peptide terminal modifications
                    {
                        if (modpos.Key == 0 || modpos.Key == peptide.BaseSequence.Length + 1)
                        {
                            continue;
                        }
                    }

                    if (!ModifiedAminoAcidPositionsInProtein.ContainsKey(modPositionInProtein))
                    {
                        ModifiedAminoAcidPositionsInProtein[modPositionInProtein] = new Dictionary<string, UtilModification>();
                        PeptidesByProteinPosition[modPositionInProtein] = new List<UtilPeptide>();
                    }

                    PeptidesByProteinPosition[modPositionInProtein].Add(peptide);

                    foreach (var mod in modpos.Value.Values)
                    {
                        mod.ProteinPositionZeroIsNTerminus = modPositionInProtein;
                        if (!ModifiedAminoAcidPositionsInProtein[modPositionInProtein].ContainsKey(mod.IdWithMotif))
                        {
                            ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.IdWithMotif] = new UtilModification(mod.IdWithMotif, mod.PeptidePositionZeroIsNTerminus, modPositionInProtein, 0);
                        }
                        ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.IdWithMotif].Intensity += mod.Intensity;
                    }
                }
            }
        }

        public Dictionary<int, Dictionary<string, UtilModification>> GetModStoichiometryFromProteinMods()
        {
            if (ModifiedAminoAcidPositionsInProtein == null)
            {
                SetProteinModsFromPeptides();
            }

            var aaModsStoichiometry = ModifiedAminoAcidPositionsInProtein;
            foreach (var modpos in aaModsStoichiometry)
            {
                foreach (var mod in modpos.Value.Values)
                {
                    mod.Intensity = mod.Intensity / PeptidesByProteinPosition[modpos.Key].Select(x => x.Intensity).Sum();
                }
            }
            return aaModsStoichiometry;
        }
    }

    public class UtilProteinGroup
    {
        public string Name { get; set;}
        public Dictionary<string, UtilProtein> Proteins {  get; set; }
        public string OccupancyLevel { get; set; }

        public UtilProteinGroup(string name, Dictionary<string, UtilProtein> proteins = null)
        {
            Name = name;
            if (proteins != null) Proteins = proteins;
            else Proteins= new Dictionary<string, UtilProtein>();
        }
    }
    public class PositionFrequencyAnalysis
    { 

        public Dictionary<string, UtilProteinGroup> Occupancy { get; private set; }

        /// <summary>
        /// Calculates the occupancy of post-translational modifications at the peptide level. 
        /// </summary>
        /// <param name="peptides"> A List of Tuples whose entries are ordered as (string FullSequence, string BaseSequence, List<string> ProteinGroups, Intensity) for each peptide.</param>
        /// <param name="ignoreTerminusMod"> If true, terminal modifications will be ignored.</param>
        /// <returns> A nested dictionary whose key mappings are as follows: string ProteinGroup-> string Protein-> string BaseSequence-> int ModifiedAminoAcidIndex-> string ModificationName-> double Intensity
        /// Note: Each BaseSequence dictionary contains a ModifiedAminoAcidIndex key of -1 that then contains a ModificationName key called "Total" that is used to track the total intensity observed for 
        /// all of the amino acids in that peptide.</returns>
        ///
        public void ProteinGroupsOccupancyByPeptide(List<(string fullSeq, string baseSeq, List<string> proteinGroup, double intensity)> peptides, bool ignoreTerminusMod=false)
        {
            // ToDo: change first argument to Dictionary<IPeptide, intensity>
            var proteinGroups = new Dictionary<string, UtilProteinGroup>();
            
            // Go through the peptides given
            foreach (var pep in peptides)
            {
                string baseSeq = pep.Item2.IsNotNullOrEmpty() ? pep.Item2 : new string(pep.Item1.ToCharArray()); // in case it is null or empty and we need to get the base sequence from the full sequence
                ClassExtensions.RemoveSpecialCharacters(ref baseSeq, @"", @"\[(.+?)\](?<!\[I+\])"); 

                // Go through the peptide's protein groups
                foreach (var pg in pep.proteinGroup)
                {
                    // If have not seen that protein group, store it
                    if (!proteinGroups.ContainsKey(pg))
                    {
                        proteinGroups[pg] = new UtilProteinGroup(pg);
                        proteinGroups[pg].OccupancyLevel = "peptide";
                    }
                    var proteinGroup = proteinGroups[pg];

                    // Go through the proteins in each protein group
                    foreach (var proteinName in pg.Split('|'))
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
                            protein.Peptides[baseSeq] = new UtilPeptide(pep.fullSeq);
                            protein.Peptides[baseSeq].Intensity = 0;
                        }

                        // Increase the total intensity of the peptide base sequence to track the total intensity of all amino acids in that sequence
                        protein.Peptides[baseSeq].Intensity += pep.intensity;
                        var peptide = protein.Peptides[baseSeq];

                        // Want both arguments passed here to be true if need to later filter out peptide terminal mods that are not protein terminal mods 
                        Dictionary<int, List<string>> peptideMods = pep.fullSeq.ParseModifications(ignoreTerminusMod);
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
                                modifiedPosition[mod].Intensity += pep.intensity;
                            }
                        }
                    }
                }
            }
            Occupancy = proteinGroups;
        }

        //public void ProteinGroupsOccupancyByProtein(Dictionary<string, string> proteinSequences) // Dictionary<accession, sequence>

    }
}
