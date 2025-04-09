using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text.RegularExpressions;
using Easy.Common.Extensions;

namespace MzLibUtil
{
    public class QuantifiedModification
    {
        public string IdWithMotif { get; set; }
        public int PeptidePositionZeroIsNTerminus { get; set; }
        public int ProteinPositionZeroIsNTerminus { get; set; }
        public double Intensity { get; set; }

        public QuantifiedModification(string name, int positionInPeptide, int? positionInProtein=null, double intensity=0)
        {
            IdWithMotif = name;
            PeptidePositionZeroIsNTerminus = positionInPeptide;
            ProteinPositionZeroIsNTerminus = positionInProtein ?? -1; // -1 means that the position in the protein is unknown
            Intensity = intensity;
        }
    }
    /// <summary>
    /// A class to store information about a quantified peptide full sequence
    /// UtilPeptide object with modification information from all full sequences stored in the same object 
    /// </summary>
    public class QuantifiedPeptide
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; set; }
        public QuantifiedProtein ParentProtein { get; set; }
        public int OneBasedStartIndexInProtein { get; set; }
        public Dictionary<int, QuantifiedModification> ModifiedAminoAcidPositions { get; set; } 
        public double Intensity { get; set; }

        public QuantifiedPeptide(string fullSequence, Dictionary<int, QuantifiedModification> mods = null, int oneBasedStartIndexInProtein = -1, double intensity = 0, string modPattern=null) 
        {
            FullSequence = fullSequence;
            ModifiedAminoAcidPositions = mods.IsNotNullOrEmpty() ? mods : new Dictionary<int, QuantifiedModification>();
            OneBasedStartIndexInProtein = oneBasedStartIndexInProtein; // -1 means that the position in the protein is unknown
            Intensity = intensity;
            SetBaseSequence(modPattern);
        }
        public void SetBaseSequence(string modPattern)
        {
            Regex mods = modPattern != null ? new(modPattern) : new(ClassExtensions.modificationPattern);
            BaseSequence = mods.Replace(FullSequence, @"");
        }

        public Dictionary<int, QuantifiedModification> GetModStoichiometryForPeptide()
        {
            var aaModsStoichiometry = ModifiedAminoAcidPositions;
            aaModsStoichiometry.ForEach(x => x.Value.Intensity = x.Value.Intensity / Intensity);
            return aaModsStoichiometry;
        }
    }
    
    public class QuantifiedProtein
    {
        public string Accession { get; set; }
        public string Sequence { get; set; }
        public Dictionary<string, QuantifiedPeptide> Peptides { get; set; } //Unique by full sequence
        public Dictionary<int, Dictionary<string, QuantifiedModification>> ModifiedAminoAcidPositionsInProtein { get; set; }
        public Dictionary<int, List<QuantifiedPeptide>> PeptidesByProteinPosition { get; set; }

        public QuantifiedProtein(string accession, string sequence=null, Dictionary<string, QuantifiedPeptide> peptides=null)
        {
            Accession = accession;
            Sequence = sequence;
            Peptides = peptides ?? new Dictionary<string, QuantifiedPeptide>();
        }   

        public void SetProteinModsFromPeptides()
        {
            if (!Sequence.IsNotNullOrEmpty() || !Peptides.IsNotNullOrEmpty())
            {
                throw new Exception("The protein sequence is unknown, or there're no peptides.");
            }   

            ModifiedAminoAcidPositionsInProtein = new Dictionary<int, Dictionary<string, QuantifiedModification>>();
            PeptidesByProteinPosition = new Dictionary<int, List<QuantifiedPeptide>>();

            foreach (var peptide in Peptides.Values)
            {
                foreach (var modpos in peptide.ModifiedAminoAcidPositions)
                {
                    var modPositionInProtein = modpos.Key + peptide.OneBasedStartIndexInProtein - 1;

                    if ((modPositionInProtein != 0 && modpos.Key == 0) // if the mod is at the N-terminus of the peptide, but not the protein.
                        || (modPositionInProtein != Sequence.Length + 1 && modpos.Key == peptide.BaseSequence.Length + 1)) // if the mod is at the C-terminus of the peptide, but not the protein.
                    {
                        continue;
                    }

                    if (!ModifiedAminoAcidPositionsInProtein.ContainsKey(modPositionInProtein))
                    {
                        ModifiedAminoAcidPositionsInProtein[modPositionInProtein] = new Dictionary<string, QuantifiedModification>();
                        PeptidesByProteinPosition[modPositionInProtein] = new List<QuantifiedPeptide>();
                    }

                    PeptidesByProteinPosition[modPositionInProtein].Add(peptide);

                    foreach (var mod in modpos.Value.Values)
                    {
                        mod.ProteinPositionZeroIsNTerminus = modPositionInProtein;
                        if (!ModifiedAminoAcidPositionsInProtein[modPositionInProtein].ContainsKey(mod.IdWithMotif))
                        {
                            ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.IdWithMotif] = new QuantifiedModification(mod.IdWithMotif, mod.PeptidePositionZeroIsNTerminus, modPositionInProtein, 0);
                        }
                        ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.IdWithMotif].Intensity += mod.Intensity;
                    }
                }
            }
        }

        public Dictionary<int, Dictionary<string, QuantifiedModification>> GetModStoichiometryFromProteinMods()
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

    public class QuantifiedProteinGroup
    {
        public string Name { get; set;}
        public Dictionary<string, QuantifiedProtein> Proteins {  get; set; }
        public string OccupancyLevel { get; set; }

        public QuantifiedProteinGroup(string name, Dictionary<string, QuantifiedProtein> proteins = null)
        {
            Name = name;
            if (proteins != null) Proteins = proteins;
            else Proteins= new Dictionary<string, QuantifiedProtein>();
        }
    }
    public class PositionFrequencyAnalysis
    { 

        public Dictionary<string, QuantifiedProteinGroup> Occupancy { get; private set; }

        /// <summary>
        /// Calculates the occupancy of post-translational modifications at the peptide level. 
        /// </summary>
        /// <param name="peptides"> A List of Tuples whose entries are ordered as (string FullSequence, string BaseSequence, List<string> ProteinGroups, Intensity) for each peptide.</param>
        /// <param name="ignoreTerminusMod"> If true, terminal modifications will be ignored.</param>
        /// <returns> A nested dictionary whose key mappings are as follows: string ProteinGroup-> string Protein-> string BaseSequence-> int ModifiedAminoAcidIndex-> string ModificationName-> double Intensity
        /// Note: Each BaseSequence dictionary contains a ModifiedAminoAcidIndex key of -1 that then contains a ModificationName key called "Total" that is used to track the total intensity observed for 
        /// all of the amino acids in that peptide.</returns>
        ///
        public void ProteinGroupsOccupancyByPeptide(List<(string fullSeq, List<string> proteinGroup, double intensity)> peptides, bool ignoreTerminusMod=false)
        {
            // ToDo: change first argument to Dictionary<IPeptide, intensity>
            Occupancy = new Dictionary<string, QuantifiedProteinGroup>();
            
            // Go through the peptides given
            foreach (var pep in peptides)
            {
                //string baseSeq = pep.Item2.IsNotNullOrEmpty() ? pep.Item2 : new string(pep.Item1.ToCharArray()); // in case it is null or empty and we need to get the base sequence from the full sequence
                //ClassExtensions.RemoveSpecialCharacters(ref baseSeq, @"", ClassExtensions.modificationPattern); 
                string baseSeq = pep.fullSeq.ParseBaseSequence();

                // Go through the peptide's protein groups
                foreach (var pg in pep.proteinGroup)
                {
                    // If have not seen that protein group, store it
                    if (!Occupancy.ContainsKey(pg))
                    {
                        Occupancy[pg] = new QuantifiedProteinGroup(pg);
                        Occupancy[pg].OccupancyLevel = "peptide";
                    }
                    var proteinGroup = Occupancy[pg];

                    // Go through the proteins in each protein group
                    foreach (var proteinName in pg.Split('|'))
                    {
                        // Add the protein to the protein group's dictionary if it has not been added
                        if (!proteinGroup.Proteins.ContainsKey(proteinName))
                        {
                            proteinGroup.Proteins[proteinName] = new QuantifiedProtein(proteinName);
                        }
                        var protein = proteinGroup.Proteins[proteinName];

                        // If the peptide's base sequence has not been seen, add it to the protein's dictionary
                        if (!protein.Peptides.ContainsKey(baseSeq))
                        {
                            protein.Peptides[baseSeq] = new QuantifiedPeptide(baseSeq);
                            protein.Peptides[baseSeq].Intensity = 0;
                        }

                        // Increase the total intensity of the peptide base sequence to track the total intensity of all amino acids in that sequence
                        protein.Peptides[baseSeq].Intensity += pep.intensity;
                        var peptide = protein.Peptides[baseSeq];

                        // Want both arguments passed here to be true if need to later filter out peptide terminal mods that are not protein terminal mods 
                        Dictionary<int, string> peptideMods = pep.fullSeq.ParseModifications(ignoreTerminusMod);
                        // Go through the modified positions found froum the full sequence
                        foreach (var modpos in peptideMods)
                        {
                            // If that position has not been recorded as containing a modification, add it to the base sequence's dictonary
                            if (!peptide.ModifiedAminoAcidPositions.ContainsKey(modpos.Key))
                            {
                                peptide.ModifiedAminoAcidPositions[modpos.Key] = new Dictionary<string, QuantifiedModification>();
                            }
                            var modifiedPosition = peptide.ModifiedAminoAcidPositions[modpos.Key];

                            if (!modifiedPosition.ContainsKey(modpos.Value))
                            {
                                modifiedPosition[modpos.Value] = new QuantifiedModification(modpos.Value, modpos.Key, 0);
                            }
                            
                            // Increase the intensity of the modification by the intensity of the peptide
                            modifiedPosition[modpos.Value].Intensity += pep.intensity;
                        }
                    }
                }
            }
        }
    }
}
