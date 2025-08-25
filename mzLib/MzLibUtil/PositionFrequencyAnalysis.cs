using Easy.Common.Extensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace MzLibUtil
{
    public class QuantifiedModification
    {
        public string IdWithMotif { get; set; }
        public string ModificationLocalization { get; set; } // e.g. "N-terminus", "C-terminus", or amino acid name
        public int PeptidePositionZeroIsNTerminus { get; set; }
        public int ProteinPositionZeroIsNTerminus { get; set; }
        public double Intensity { get; set; }

        public QuantifiedModification(string idWithMotif, int positionInPeptide, int? positionInProtein = null, string modLocalization = null, double intensity = 0)
        {
            IdWithMotif = idWithMotif;
            PeptidePositionZeroIsNTerminus = positionInPeptide;
            ProteinPositionZeroIsNTerminus = positionInProtein ?? -1; // -1 means that the position in the protein is unknown
            ModificationLocalization = modLocalization ?? "Unknown"; 
            Intensity = intensity;
        }
    }
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
                        var modLocalization = modpos == 0 ? "N-terminus" : (modpos == BaseSequence.Length + 1 ? "C-terminus" : BaseSequence[modpos - 1].ToString());
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

    public class QuantifiedProtein
    {
        public string Accession { get; set; }
        public string Sequence { get; set; }
        public Dictionary<string, QuantifiedPeptide> Peptides { get; set; }
        public Dictionary<int, Dictionary<string, QuantifiedModification>> ModifiedAminoAcidPositionsInProtein { get; set; }
        public Dictionary<int, HashSet<string>> PeptidesByProteinPosition { get; set; }

        public QuantifiedProtein(string accession, string sequence = null, Dictionary<string, QuantifiedPeptide> peptides = null)
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
            PeptidesByProteinPosition = new Dictionary<int, HashSet<string>>();

            foreach (var peptide in Peptides.Values)
            {
                // if peptide position in protein is unknown, set it using the protein sequence
                if (peptide.OneBasedStartIndexInProtein == -1)
                {
                    peptide.OneBasedStartIndexInProtein = Sequence.IndexOf(peptide.BaseSequence) + 1;
                }
                // if peptide has no modifications, add to all its positions
                if (!peptide.ModifiedAminoAcidPositions.IsNotNullOrEmpty())
                {
                    for (int i = 0; i < peptide.BaseSequence.Length; i++)
                    {
                        var pos = peptide.OneBasedStartIndexInProtein + i;
                        if (!ModifiedAminoAcidPositionsInProtein.ContainsKey(pos))
                        {
                            ModifiedAminoAcidPositionsInProtein[pos] = new Dictionary<string, QuantifiedModification>();
                            PeptidesByProteinPosition[pos] = new HashSet<string>();
                        }
                        PeptidesByProteinPosition[pos].Add(peptide.BaseSequence);
                    }
                    continue;
                }

                else // if peptide has modifications, add to modified positions
                {
                    foreach (var modpos in peptide.ModifiedAminoAcidPositions.Keys)
                    {
                        var modPositionInProtein = modpos + peptide.OneBasedStartIndexInProtein - 1;

                        // Ignore peptide terminal modifications that are not at the protein terminal
                        if ((modPositionInProtein != 0 && modpos == 0) // if the mod is at the N-terminus of the peptide, but not the protein.
                            || (modPositionInProtein != Sequence.Length + 1 && modpos == peptide.BaseSequence.Length + 1)) // if the mod is at the C-terminus of the peptide, but not the protein.
                        {
                            continue;
                        }

                        if (!ModifiedAminoAcidPositionsInProtein.ContainsKey(modPositionInProtein))
                        {
                            ModifiedAminoAcidPositionsInProtein[modPositionInProtein] = new Dictionary<string, QuantifiedModification>();
                            PeptidesByProteinPosition[modPositionInProtein] = new HashSet<string>();
                        }
                        PeptidesByProteinPosition[modPositionInProtein].Add(peptide.BaseSequence);

                        foreach (var mod in peptide.ModifiedAminoAcidPositions[modpos].Values)
                        {
                            mod.ProteinPositionZeroIsNTerminus = modPositionInProtein;

                            if (!ModifiedAminoAcidPositionsInProtein[modPositionInProtein].ContainsKey(mod.IdWithMotif))
                            {
                                ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.IdWithMotif] = new QuantifiedModification(mod.IdWithMotif, mod.PeptidePositionZeroIsNTerminus, modPositionInProtein, null, 0);
                            }
                            ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.IdWithMotif].Intensity += mod.Intensity;
                        }
                    }
                }
            }

            // clean up the dictionary to remove any empty modifications
            var noModPositions = ModifiedAminoAcidPositionsInProtein.Where(x => !x.Value.IsNotNullOrEmpty()).ToDictionary().Keys;
            foreach (var pos in noModPositions)
            {
                ModifiedAminoAcidPositionsInProtein.Remove(pos);
                PeptidesByProteinPosition.Remove(pos);
            }

        }

        public Dictionary<int, Dictionary<string, double>> GetModStoichiometryFromProteinMods()
        {
            SetProteinModsFromPeptides();
            var aaModsStoichiometry = new Dictionary<int, Dictionary<string, double>>();
            foreach (var modpos in ModifiedAminoAcidPositionsInProtein.Keys)
            {
                aaModsStoichiometry[modpos] = new Dictionary<string, double>();

                double totalPositionIntensity = Peptides.Where(pep => PeptidesByProteinPosition[modpos].Contains(pep.Key)).Sum(x => x.Value.Intensity);
                foreach (var mod in ModifiedAminoAcidPositionsInProtein[modpos].Values)
                {
                    double modFraction = mod.Intensity / totalPositionIntensity;
                    aaModsStoichiometry[modpos].Add(mod.IdWithMotif, modFraction);
                }
            }
            return aaModsStoichiometry;
        }
    }

    public class QuantifiedProteinGroup
    {
        public string Name { get; set; }
        public Dictionary<string, QuantifiedProtein> Proteins { get; set; }
        public string OccupancyLevel { get; set; }

        public QuantifiedProteinGroup(string name, Dictionary<string, QuantifiedProtein> proteins = null)
        {
            Name = name;
            if (proteins != null) Proteins = proteins;
            else Proteins = new Dictionary<string, QuantifiedProtein>();
        }
    }
    public class PositionFrequencyAnalysis
    {

        public Dictionary<string, QuantifiedProteinGroup> ProteinGroupOccupancies { get; private set; }
        public Dictionary<string, (QuantifiedPeptide QuantifiedPeptide, string ProteinGroups)> PeptideOccupancies { get; private set; }

        /// <summary>
        /// Calculates the occupancy of post-translational modifications at the peptide level. 
        /// </summary>
        /// <param name="peptides"> A List of Tuples whose entries are ordered as (string FullSequence, string BaseSequence, List<string> ProteinGroups, Intensity) for each peptide.</param>
        /// <param name="ignoreTerminusMod"> If true, terminal modifications will be ignored.</param>
        /// <returns> A nested dictionary whose key mappings are as follows: string ProteinGroup-> string Protein-> string BaseSequence-> int ModifiedAminoAcidIndex-> string ModificationName-> double Intensity
        /// Note: Each BaseSequence dictionary contains a ModifiedAminoAcidIndex key of -1 that then contains a ModificationName key called "Total" that is used to track the total intensity observed for 
        /// all of the amino acids in that peptide.</returns>
        ///
        public void CalculateOccupancies(List<(string fullSeq, List<string> proteinGroups, double intensity)> peptides, bool ignoreTerminusMod = false)
        {
            // ToDo: change first argument to Dictionary<IPeptide, intensity>
            ProteinGroupOccupancies = new Dictionary<string, QuantifiedProteinGroup>();
            PeptideOccupancies = new Dictionary<string, (QuantifiedPeptide, string)>();

            // Go through the peptides given
            foreach (var pep in peptides)
            {
                //string baseSeq = pep.Item2.IsNotNullOrEmpty() ? pep.Item2 : new string(pep.Item1.ToCharArray()); // in case it is null or empty and we need to get the base sequence from the full sequence
                //ClassExtensions.RemoveSpecialCharacters(ref baseSeq, @"", ClassExtensions.modificationPattern); 
                string baseSeq = pep.fullSeq.GetBaseSequenceFromFullSequence();

                if (!PeptideOccupancies.ContainsKey(pep.fullSeq))
                {
                    // Need to make sure clustering of proteingroups is correct
                    string proteinGroupsJoined = string.Join(";", pep.proteinGroups);
                    PeptideOccupancies[pep.fullSeq] = (new QuantifiedPeptide(pep.fullSeq, intensity: pep.intensity), proteinGroupsJoined);
                }
                else
                {
                    PeptideOccupancies[pep.fullSeq].QuantifiedPeptide.AddFullSequence(pep.fullSeq, intensity: pep.intensity);
                }

                // Go through the peptide's protein groups
                foreach (var pg in pep.proteinGroups)
                {
                    // If have not seen that protein group, store it
                    if (!ProteinGroupOccupancies.ContainsKey(pg))
                    {
                        ProteinGroupOccupancies[pg] = new QuantifiedProteinGroup(pg);
                        ProteinGroupOccupancies[pg].OccupancyLevel = "peptide";
                    }
                    var proteinGroup = ProteinGroupOccupancies[pg];

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
