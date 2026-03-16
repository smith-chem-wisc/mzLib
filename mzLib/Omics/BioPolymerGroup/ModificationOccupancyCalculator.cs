using Omics.BioPolymer;
using Omics.Modifications;

namespace Omics.BioPolymerGroup;

/// <summary>
/// Calculates modification occupancy/stoichiometry from identified peptides.
/// Supports both count-based and intensity-based metrics, at the protein or peptide level.
/// </summary>
public static class ModificationOccupancyCalculator
{
    /// <summary>
    /// Mod types to exclude from occupancy calculations.
    /// </summary>
    private static readonly string[] ExcludedModTypes = ["Common Variable", "Common Fixed"];

    /// <summary>
    /// Location restrictions to exclude (peptide-terminal, not protein-terminal).
    /// </summary>
    private static readonly string[] ExcludedLocations = ["NPep", "PepC"];

    /// <summary>
    /// Calculates per-site modification occupancy mapped to protein coordinates.
    /// </summary>
    /// <param name="bioPolymer">The parent biopolymer whose length defines the coordinate space.</param>
    /// <param name="localizedSequences">Peptides with localized modifications mapped to this biopolymer.</param>
    /// <param name="intensitiesByFullSequence">
    /// Optional map of FullSequence → intensity. When provided, intensity-based stoichiometry is calculated.
    /// When null, only count-based occupancy is populated.
    /// </param>
    /// <returns>
    /// Dictionary keyed by one-based protein position, each value a list of
    /// <see cref="SiteSpecificModificationOccupancy"/> entries for modifications observed at that position.
    /// </returns>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculateProteinLevelOccupancy(
        IBioPolymer bioPolymer,
        IEnumerable<IBioPolymerWithSetMods> localizedSequences,
        Dictionary<string, double>? intensitiesByFullSequence = null)
    {
        var sequences = localizedSequences as IList<IBioPolymerWithSetMods> ?? localizedSequences.ToList();
        // Use an inner dictionary for dedup during construction, then flatten to lists
        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();

        foreach (var sequence in sequences)
        {
            foreach (var mod in sequence.AllModsOneIsNterminus)
            {
                if (!TryGetProteinPosition(mod, sequence, bioPolymer.Length, out int indexInProtein))
                    continue;

                if (!working.TryGetValue(indexInProtein, out var modsAtPosition))
                {
                    modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                    working[indexInProtein] = modsAtPosition;
                }

                if (!modsAtPosition.TryGetValue(mod.Value.IdWithMotif, out var siteOccupancy))
                {
                    siteOccupancy = new SiteSpecificModificationOccupancy(indexInProtein, mod.Value.IdWithMotif);

                    // Count total peptides covering this position
                    foreach (var seq in sequences)
                    {
                        int rangeStart = seq.OneBasedStartResidue - (indexInProtein == 1 ? 1 : 0);
                        if (indexInProtein >= rangeStart && indexInProtein <= seq.OneBasedEndResidue)
                        {
                            siteOccupancy.TotalCount++;
                            if (intensitiesByFullSequence != null &&
                                seq.FullSequence != null &&
                                intensitiesByFullSequence.TryGetValue(seq.FullSequence, out double seqIntensity))
                            {
                                siteOccupancy.TotalIntensity += seqIntensity;
                            }
                        }
                    }

                    modsAtPosition[mod.Value.IdWithMotif] = siteOccupancy;
                }

                siteOccupancy.ModifiedCount++;
                if (intensitiesByFullSequence != null &&
                    sequence.FullSequence != null &&
                    intensitiesByFullSequence.TryGetValue(sequence.FullSequence, out double intensity))
                {
                    siteOccupancy.ModifiedIntensity += intensity;
                }
            }
        }

        return working.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.Values.ToList());
    }

    /// <summary>
    /// Calculates per-site modification occupancy in peptide-local coordinates
    /// for a group of peptides sharing the same base sequence.
    /// Positions use the AllModsOneIsNterminus convention (1 = N-terminus, 2 = first residue, etc.).
    /// </summary>
    /// <param name="peptides">
    /// Peptides sharing the same base sequence. All must have the same BaseSequence.
    /// </param>
    /// <param name="intensitiesByFullSequence">
    /// Optional map of FullSequence → intensity for intensity-based stoichiometry.
    /// </param>
    /// <returns>
    /// Dictionary keyed by peptide-local position (AllModsOneIsNterminus key) containing a list of
    /// <see cref="SiteSpecificModificationOccupancy"/> entries for modifications observed at that position.
    /// </returns>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculatePeptideLevelOccupancy(
        IEnumerable<IBioPolymerWithSetMods> peptides,
        Dictionary<string, double>? intensitiesByFullSequence = null)
    {
        var peptideList = peptides as IList<IBioPolymerWithSetMods> ?? peptides.ToList();
        int totalPeptideCount = peptideList.Count;

        double totalGroupIntensity = 0;
        if (intensitiesByFullSequence != null)
        {
            totalGroupIntensity = peptideList
                .Where(p => p.FullSequence != null && intensitiesByFullSequence.ContainsKey(p.FullSequence))
                .Sum(p => intensitiesByFullSequence[p.FullSequence]);
        }

        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();

        foreach (var peptide in peptideList)
        {
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                if (IsExcludedMod(mod.Value))
                    continue;

                if (!working.TryGetValue(mod.Key, out var modsAtPosition))
                {
                    modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                    working[mod.Key] = modsAtPosition;
                }

                if (!modsAtPosition.TryGetValue(mod.Value.IdWithMotif, out var siteOccupancy))
                {
                    siteOccupancy = new SiteSpecificModificationOccupancy(mod.Key, mod.Value.IdWithMotif)
                    {
                        TotalCount = totalPeptideCount,
                        TotalIntensity = totalGroupIntensity
                    };
                    modsAtPosition[mod.Value.IdWithMotif] = siteOccupancy;
                }

                siteOccupancy.ModifiedCount++;
                if (intensitiesByFullSequence != null &&
                    peptide.FullSequence != null &&
                    intensitiesByFullSequence.TryGetValue(peptide.FullSequence, out double intensity))
                {
                    siteOccupancy.ModifiedIntensity += intensity;
                }
            }
        }

        return working.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.Values.ToList());
    }

    /// <summary>
    /// Maps an AllModsOneIsNterminus entry to a one-based protein position based on the
    /// modification's location restriction. Returns false if the mod should be skipped.
    /// </summary>
    private static bool TryGetProteinPosition(
        KeyValuePair<int, Modification> mod,
        IBioPolymerWithSetMods sequence,
        int bioPolymerLength,
        out int indexInProtein)
    {
        indexInProtein = 0;

        if (IsExcludedMod(mod.Value))
            return false;

        if (mod.Value.LocationRestriction.Equals("N-terminal."))
        {
            indexInProtein = 1;
        }
        else if (mod.Value.LocationRestriction.Equals("Anywhere."))
        {
            indexInProtein = sequence.OneBasedStartResidue + mod.Key - 2;
        }
        else if (mod.Value.LocationRestriction.Equals("C-terminal."))
        {
            indexInProtein = bioPolymerLength;
        }
        else
        {
            return false;
        }

        return true;
    }

    private static bool IsExcludedMod(Modification mod)
    {
        foreach (var excludedType in ExcludedModTypes)
        {
            if (mod.ModificationType.Contains(excludedType))
                return true;
        }

        foreach (var excludedLocation in ExcludedLocations)
        {
            if (mod.LocationRestriction.Equals(excludedLocation))
                return true;
        }

        return false;
    }
}
