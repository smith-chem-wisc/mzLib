using System.Text;
using Chemistry;
using Omics.Modifications;

namespace Omics;

public static class BioPolymerWithSetModsExtensions
{
    /// <summary>
    /// This method returns the full sequence with mass shifts INSTEAD OF PTMs in brackets []
    /// Some external tools cannot parse PTMs, instead requiring a numerical input indicating the mass of a PTM in brackets
    /// after the position of that modification
    /// N-terminal mas shifts are in brackets prior to the first amino acid and apparently missing the + sign
    /// </summary>
    /// <returns></returns>
    public static string FullSequenceWithMassShift(this IBioPolymerWithSetMods withSetMods)
    {
        var subsequence = new StringBuilder();

        // modification on peptide N-terminus
        if (withSetMods.AllModsOneIsNterminus.TryGetValue(1, out Modification? mod))
        {
            subsequence.Append($"[{mod.MonoisotopicMass.RoundedDouble(6)}]");
        }

        for (int r = 0; r < withSetMods.Length; r++)
        {
            subsequence.Append(withSetMods[r]);

            // modification on this residue
            if (withSetMods.AllModsOneIsNterminus.TryGetValue(r + 2, out mod))
            {
                if (mod.MonoisotopicMass > 0)
                {
                    subsequence.Append($"[+{mod.MonoisotopicMass.RoundedDouble(6)}]");
                }
                else
                {
                    subsequence.Append($"[{mod.MonoisotopicMass.RoundedDouble(6)}]");
                }
            }
        }

        // modification on peptide C-terminus
        if (withSetMods.AllModsOneIsNterminus.TryGetValue(withSetMods.Length + 2, out mod))
        {
            if (mod.MonoisotopicMass > 0)
            {
                subsequence.Append($"-[+{mod.MonoisotopicMass.RoundedDouble(6)}]");
            }
            else
            {
                subsequence.Append($"-[{mod.MonoisotopicMass.RoundedDouble(6)}]");
            }
        }
        return subsequence.ToString();
    }

    /// <summary>
    /// This method returns the full sequence only with the specified modifications in the modstoWritePruned dictionary
    /// </summary>
    /// <param name="withSetMods"></param>
    /// <param name="modstoWritePruned"></param>
    /// <returns></returns>
    public static string EssentialSequence(this IBioPolymerWithSetMods withSetMods,
        IReadOnlyDictionary<string, int> modstoWritePruned)
    {
        string essentialSequence = withSetMods.BaseSequence;
        if (modstoWritePruned != null)
        {
            var sbsequence = new StringBuilder(withSetMods.FullSequence.Length);

            // variable modification on peptide N-terminus
            if (withSetMods.AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
            {
                if (modstoWritePruned.ContainsKey(pep_n_term_variable_mod.ModificationType))
                {
                    sbsequence.Append(
                        $"[{pep_n_term_variable_mod.ModificationType}:{pep_n_term_variable_mod.IdWithMotif}]");
                }
            }
            for (int r = 0; r < withSetMods.Length; r++)
            {
                sbsequence.Append(withSetMods[r]);
                // variable modification on this residue
                if (withSetMods.AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                {
                    if (modstoWritePruned.ContainsKey(residue_variable_mod.ModificationType))
                    {
                        sbsequence.Append(
                            $"[{residue_variable_mod.ModificationType}:{residue_variable_mod.IdWithMotif}]");
                    }
                }
            }

            // variable modification on peptide C-terminus
            if (withSetMods.AllModsOneIsNterminus.TryGetValue(withSetMods.Length + 2, out Modification pep_c_term_variable_mod))
            {
                if (modstoWritePruned.ContainsKey(pep_c_term_variable_mod.ModificationType))
                {
                    sbsequence.Append(
                        $"-[{pep_c_term_variable_mod.ModificationType}:{pep_c_term_variable_mod.IdWithMotif}]");
                }
            }

            essentialSequence = sbsequence.ToString();
        }
        return essentialSequence;
    }

    public static string DetermineFullSequence(this IBioPolymerWithSetMods withSetMods) => IBioPolymerWithSetMods
        .DetermineFullSequence(withSetMods.BaseSequence, withSetMods.AllModsOneIsNterminus);
}