using System;
using System.Collections.Generic;
using System.Linq;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    /// <summary>
    /// A class to store information about a quantified protein. The protein contains peptides 
    /// clustered by their base sequence, rather than by their full sequence. Full sequences are stored
    /// in the QuantifiedPeptide objects.
    /// </summary>
    public class QuantifiedProtein
    {
        public string Accession { get; set; }
        public string Sequence { get; set; }

        /// <summary>
        /// Dictionary mapping peptide base sequences to their corresponding QuantifiedPeptide objects.
        /// </summary>
        public Dictionary<string, QuantifiedPeptide> Peptides { get; set; }

        /// <summary>
        /// Dictionary mapping zero-based amino acid positions in the protein to dictionaries of 
        /// modification IDs and their corresponding QuantifiedModification objects.
        /// Note: the modification positions are 0-based with the N-terminus of the protein being position 0.
        /// </summary>
        public Dictionary<int, Dictionary<string, QuantifiedModification>> ModifiedAminoAcidPositionsInProtein { get; set; }

        /// <summary>
        /// Dictionary mapping zero-based amino acid positions in the protein to sets of peptide base sequences
        /// This is useful to know which peptides contribute to the modification and total intensity at a given position.
        /// </summary>
        public Dictionary<int, HashSet<string>> PeptidesByProteinPosition { get; set; }

        public QuantifiedProtein(string accession, string sequence = null, Dictionary<string, QuantifiedPeptide> peptides = null)
        {
            Accession = accession;
            Sequence = sequence;
            Peptides = peptides ?? new Dictionary<string, QuantifiedPeptide>();
        }

        /// <summary>
        /// Parses and aggregates modifications from the protein's peptides to set the ModifiedAminoAcidPositionsInProtein property.
        /// </summary>
        /// <exception cref="Exception"></exception>
        public void SetProteinModsFromPeptides()
        {
            if (Sequence.IsNullOrEmpty())
            {
                throw new Exception("The protein sequence is unknown.");
            }

            ModifiedAminoAcidPositionsInProtein = new Dictionary<int, Dictionary<string, QuantifiedModification>>();
            PeptidesByProteinPosition = new Dictionary<int, HashSet<string>>();

            if (Peptides.IsNullOrEmpty())
            {
                return;
            }

            foreach (var peptide in Peptides.Values)
            {
                // always recompute start position from this protein's sequence — the peptide instance
                // may be shared across multiple proteins, so a cached value could be stale.    
                int idx = Sequence.IndexOf(peptide.BaseSequence);
                if (idx == -1)
                    throw new InvalidOperationException(
                        $"Peptide '{peptide.BaseSequence}' was not found in protein '{Accession}' sequence.");
                peptide.ZeroBasedStartIndexInProtein = idx + 1;

                // update protein prosition total observations with observed aminoacids from this peptide
                int startIndex = peptide.ZeroBasedStartIndexInProtein == 1 ? 0 : peptide.ZeroBasedStartIndexInProtein; // if the peptide is at the N-terminus of the protein, the start position should include the N-terminus (0).
                int endIndex = peptide.ZeroBasedStartIndexInProtein + peptide.BaseSequence.Length - 1 == Sequence.Length
                    ? peptide.ZeroBasedStartIndexInProtein + peptide.BaseSequence.Length       // C-terminal peptide: extend to include the protein C-terminus position (Sequence.Length + 1)
                    : peptide.ZeroBasedStartIndexInProtein + peptide.BaseSequence.Length - 1;  // non-C-terminal: last amino acid position only
                for (int pos = startIndex; pos <= endIndex; pos++)
                {
                    if (!PeptidesByProteinPosition.ContainsKey(pos))
                    {
                        PeptidesByProteinPosition[pos] = new HashSet<string>();
                    }
                    PeptidesByProteinPosition[pos].Add(peptide.BaseSequence);
                }

                // if no mods in peptide, no need to update the ModifiedAminoAcidPositionsInProtein
                if (peptide.ModifiedAminoAcidPositions.IsNullOrEmpty())
                {
                    continue;
                }

                else
                {
                    foreach (var modpos in peptide.ModifiedAminoAcidPositions.Keys)
                    {
                        var modPositionInProtein = modpos + peptide.ZeroBasedStartIndexInProtein - 1;

                        // Ignore peptide terminal modifications that are not at the protein terminal
                        if ((modPositionInProtein != 0 && modpos == 0) // if the mod is at the N-terminus of the peptide, but not the protein.
                            || (modPositionInProtein != Sequence.Length + 1 && modpos == peptide.BaseSequence.Length + 1)) // if the mod is at the C-terminus of the peptide, but not the protein.
                        {
                            continue;
                        }

                        if (!ModifiedAminoAcidPositionsInProtein.TryGetValue(modPositionInProtein, out _))
                        {
                            ModifiedAminoAcidPositionsInProtein[modPositionInProtein] = new();
                        }

                        foreach (var mod in peptide.ModifiedAminoAcidPositions[modpos].Values)
                        {
                            if (!ModifiedAminoAcidPositionsInProtein[modPositionInProtein].ContainsKey(mod.Name))
                            {
                                ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.Name] = new QuantifiedModification(mod.Name, mod.PeptidePositionZeroIsNTerminus, modPositionInProtein, 0);
                            }
                            ModifiedAminoAcidPositionsInProtein[modPositionInProtein][mod.Name].Intensity += mod.Intensity;
                        }
                    }
                }
            }

            // clean up the dictionary to remove any empty modifications
            var noModPositions = ModifiedAminoAcidPositionsInProtein.Where(x => x.Value.IsNullOrEmpty()).Select(kvp => kvp.Key);
            var alwaysUnmodifiedPositions = PeptidesByProteinPosition.Where(x => !ModifiedAminoAcidPositionsInProtein.ContainsKey(x.Key)).Select(x => x.Key);
            var removablePositions = noModPositions.Concat(alwaysUnmodifiedPositions).Distinct().ToList();
            foreach (var pos in removablePositions)
            {
                ModifiedAminoAcidPositionsInProtein.Remove(pos);
                PeptidesByProteinPosition.Remove(pos);
            }

        }

        /// <summary>
        /// Calculates the stoichiometry of modifications at each amino acid position in the protein.
        /// The output is a dictionary keyed by zero-based amino acid positions in the protein and
        /// and the modification names with their corresponding stoichiometry values (fractions).
        /// </summary>
        /// <returns> 
        /// A dictionary where keys are zero-based amino acid positions in the protein and values are dictionaries
        /// mapping modification names to their stoichiometry (fraction of the total intensity at that position). For example:
        /// { 0: {"Acetyl": 0.5, "Unmodified": 0.5}, 15: {"Phospho": 1.0} } 
        /// indicates that at position 0, 50% of the intensity is from acetylated peptides and 50% from unmodified peptides, 
        /// while at position 15, all of the intensity is from phosphorylated peptides.
        /// </returns>
        public Dictionary<int, Dictionary<string, double>> GetModStoichiometryFromProteinMods()
        {
            SetProteinModsFromPeptides();
            var aaModsStoichiometry = new Dictionary<int, Dictionary<string, double>>();
            foreach (var modpos in ModifiedAminoAcidPositionsInProtein.Keys)
            {
                aaModsStoichiometry[modpos] = new Dictionary<string, double>();

                double totalPositionIntensity = Peptides
                    .Where(pep => PeptidesByProteinPosition[modpos].Contains(pep.Key))
                    .Sum(x => x.Value.Intensity);

                foreach (var mod in ModifiedAminoAcidPositionsInProtein[modpos].Values)
                {
                    double modFraction = totalPositionIntensity > 0
                        ? mod.Intensity / totalPositionIntensity
                        : 0.0;
                    aaModsStoichiometry[modpos].Add(mod.Name, modFraction);
                }
            }
            return aaModsStoichiometry;
        }
    }
}
