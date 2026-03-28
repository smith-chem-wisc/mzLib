using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications using the mzLib modification database.
/// Searches across all known modifications (MetaMorpheus, UniProt, UNIMOD, RNA mods).
/// Supports mzLib-style identifiers like "Oxidation on M", "Carbamidomethyl on C", etc.
/// </summary>
public class MzLibModificationLookup : ModificationLookupBase
{

    /// <summary>
    /// Singleton instance that searches both protein and RNA modifications.
    /// </summary>
    public static MzLibModificationLookup Instance { get; } = new();

    /// <summary>
    /// Instance that only searches protein modifications.
    /// </summary>
    public static MzLibModificationLookup ProteinOnly { get; } = new(searchProteinMods: true, searchRnaMods: false);

    /// <summary>
    /// Instance that only searches RNA modifications.
    /// </summary>
    public static MzLibModificationLookup RnaOnly { get; } = new(searchProteinMods: false, searchRnaMods: true);

    /// <summary>
    /// Creates a new MzLibModificationLookup.
    /// </summary>
    /// <param name="searchProteinMods">Whether to search protein modifications.</param>
    /// <param name="searchRnaMods">Whether to search RNA modifications.</param>
    /// <param name="massTolerance">Tolerance for mass-based matching in Daltons.</param>
    public MzLibModificationLookup(bool searchProteinMods = true, bool searchRnaMods = true, double massTolerance = 0.001)
        : base(BuildCandidateSet(searchProteinMods, searchRnaMods), massTolerance)
    {
    }

    /// <inheritdoc />
    public override string Name => "mzLib";

    private static IReadOnlyCollection<Modification> BuildCandidateSet(bool searchProteinMods, bool searchRnaMods)
    {
        if (!searchProteinMods && !searchRnaMods)
        {
            throw new ArgumentException("At least one of searchProteinMods or searchRnaMods must be true.");
        }

        if (searchProteinMods && searchRnaMods)
        {
            return Mods.MetaMorpheusModifications;
        }

        return searchProteinMods
            ? Mods.MetaMorpheusProteinModifications
            : Mods.MetaMorpheusRnaModifications;
    }
}
