using Easy.Common.Extensions;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MzLibUtil.PositionFrequencyAnalysis
{
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
                        if (modPositionInProtein != 0 && modpos == 0 // if the mod is at the N-terminus of the peptide, but not the protein.
                            || modPositionInProtein != Sequence.Length + 1 && modpos == peptide.BaseSequence.Length + 1) // if the mod is at the C-terminus of the peptide, but not the protein.
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
}
