using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;

namespace Omics
{
    /// <summary>
    /// Interface for modified and unmodified precursor ions
    /// </summary>
    /// <remarks>
    /// Proteins -> PeptideWithSetModifications : ProteolyticPeptide
    /// Nucleic Acids -> OligoWithSetMods : NucleolyticOligo
    /// </remarks>
    public interface IBioPolymerWithSetMods : IHasChemicalFormula
    {
        string BaseSequence { get; }
        string FullSequence { get; }
        double MostAbundantMonoisotopicMass { get; }
        string SequenceWithChemicalFormulas { get; }

        int OneBasedStartResidue { get; }
        int OneBasedEndResidue { get; }
        int MissedCleavages { get; }
        CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; }
        char PreviousResidue { get; }
        char NextResidue { get; }
        IDigestionParams DigestionParams { get; }



        Dictionary<int, Modification> AllModsOneIsNterminus { get; }
        int NumMods { get; }
        int NumFixedMods { get; }
        int NumVariableMods { get; }
        int Length { get; }
        char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];


        IBioPolymer Parent { get; }

        public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
            List<Product> products);

        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
            List<Product> products);

        public string DetermineFullSequence()
        {
            var subSequence = new StringBuilder();

            // modification on peptide N-terminus
            if (AllModsOneIsNterminus.TryGetValue(1, out Modification mod))
            {
                subSequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            for (int r = 0; r < Length; r++)
            {
                subSequence.Append(this[r]);

                // modification on this residue
                if (AllModsOneIsNterminus.TryGetValue(r + 2, out mod))
                {
                    subSequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
                }
            }

            // modification on peptide C-terminus
            if (AllModsOneIsNterminus.TryGetValue(Length + 2, out mod))
            {
                subSequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            return subSequence.ToString();
        }

        /// <summary>
        /// This method returns the full sequence with mass shifts INSTEAD OF PTMs in brackets []
        /// Some external tools cannot parse PTMs, instead requiring a numerical input indicating the mass of a PTM in brackets
        /// after the position of that modification
        /// N-terminal mas shifts are in brackets prior to the first amino acid and apparently missing the + sign
        /// </summary>
        /// <returns></returns>
        public string DetermineFullSequenceWithMassShifts()
        {
            var subsequence = new StringBuilder();

            // modification on peptide N-terminus
            if (AllModsOneIsNterminus.TryGetValue(1, out Modification mod))
            {
                subsequence.Append('[' + mod.MonoisotopicMass.RoundedDouble(6).ToString() + ']');
            }

            for (int r = 0; r < Length; r++)
            {
                subsequence.Append(this[r]);

                // modification on this residue
                if (AllModsOneIsNterminus.TryGetValue(r + 2, out mod))
                {
                    if (mod.MonoisotopicMass > 0)
                    {
                        subsequence.Append("[+" + mod.MonoisotopicMass.RoundedDouble(6).ToString() + ']');
                    }
                    else
                    {
                        subsequence.Append("[" + mod.MonoisotopicMass.RoundedDouble(6).ToString() + ']');
                    }
                }
            }

            // modification on peptide C-terminus
            if (AllModsOneIsNterminus.TryGetValue(Length + 2, out mod))
            {
                if (mod.MonoisotopicMass > 0)
                {
                    subsequence.Append("[+" + mod.MonoisotopicMass.RoundedDouble(6).ToString() + ']');
                }
                else
                {
                    subsequence.Append("[" + mod.MonoisotopicMass.RoundedDouble(6).ToString() + ']');
                }
            }
            return subsequence.ToString();
        }

        public string GetEssentialSequence(IReadOnlyDictionary<string, int> modstoWritePruned)
        {
            string essentialSequence = BaseSequence;
            if (modstoWritePruned != null)
            {
                var sbsequence = new StringBuilder();

                // variable modification on peptide N-terminus
                if (AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
                {
                    if (modstoWritePruned.ContainsKey(pep_n_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_n_term_variable_mod.ModificationType + ":" + pep_n_term_variable_mod.IdWithMotif + ']');
                    }
                }
                for (int r = 0; r < Length; r++)
                {
                    sbsequence.Append(this[r]);
                    // variable modification on this residue
                    if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                    {
                        if (modstoWritePruned.ContainsKey(residue_variable_mod.ModificationType))
                        {
                            sbsequence.Append('[' + residue_variable_mod.ModificationType + ":" + residue_variable_mod.IdWithMotif + ']');
                        }
                    }
                }

                // variable modification on peptide C-terminus
                if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification pep_c_term_variable_mod))
                {
                    if (modstoWritePruned.ContainsKey(pep_c_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_c_term_variable_mod.ModificationType + ":" + pep_c_term_variable_mod.IdWithMotif + ']');
                    }
                }

                essentialSequence = sbsequence.ToString();
            }
            return essentialSequence;
        }


        public static string GetBaseSequenceFromFullSequence(string fullSequence)
        {
            StringBuilder sb = new StringBuilder();
            int bracketCount = 0;
            foreach (char c in fullSequence)
            {
                if (c == '[')
                {
                    bracketCount++;
                }
                else if (c == ']')
                {
                    bracketCount--;
                }
                else if (bracketCount == 0)
                {
                    sb.Append(c);
                }
            }
            return sb.ToString();
        }

    }
}
