using System.Text;
using Chemistry;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Omics;

public static class BioPolymerWithSetModsExtensions
{
    /// <summary>
    /// This method returns the full sequence with mass shifts INSTEAD OF PTMs in brackets []
    /// Some external tools cannot parse PTMs, instead requiring a numerical input indicating the mass of a PTM in brackets
    /// after the position of that modification
    /// N-terminal mas shifts are in brackets prior to the first amino acid and apparently missing the + sign
    /// </summary>
    /// <param name="withSetMods">The biopolymer to serialize.</param>
    /// <param name="lookup">Optional modification lookup to resolve modifications that don't have mass information.</param>
    /// <param name="decimalPlaces">Number of decimal places for mass values (default: 6).</param>
    /// <returns>The sequence with mass shifts in bracket notation.</returns>
    public static string FullSequenceWithMassShift(this IBioPolymerWithSetMods withSetMods, 
        IModificationLookup? lookup = null, 
        int decimalPlaces = 6)
    {
        // Convert to canonical sequence using the extension method
        var canonical = withSetMods.ToCanonicalSequenceBuilder().Build();

        // Create serializer with optional lookup and schema
        var schema = new MassShiftSequenceFormatSchema(decimalPlaces);
        var serializer = new MassShiftSequenceSerializer(schema, lookup);

        // Serialize to mass shift format
        var result = serializer.Serialize(canonical);

        // Should never be null since we're using ThrowException mode by default
        return result ?? throw new InvalidOperationException("Failed to serialize sequence to mass shift format.");
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

    public static IEnumerable<Product> GetMIons(this IBioPolymerWithSetMods withSetMods, IFragmentationParams? fragmentParams)
    {
        // Normal intact molecular ion
        yield return new CustomMProduct("", withSetMods.MonoisotopicMass);

        if (fragmentParams is null)
            yield break;

        // Molecular ion with neutral losses
        foreach (var ionLoss in fragmentParams.MIonLosses)
        {
            yield return new CustomMProduct(ionLoss.Annotation, withSetMods.MonoisotopicMass - ionLoss.MonoisotopicMass);
        }
    }
}