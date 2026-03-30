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
    /// <param name="localizedSequences">
    /// All peptide forms from all PSMs mapped to this biopolymer. Used to compute
    /// <see cref="SiteSpecificModificationOccupancy.ModifiedCount"/> (numerator).
    /// </param>
    /// <param name="sequencesForTotalCount">
    /// One representative form per PSM, used to compute
    /// <see cref="SiteSpecificModificationOccupancy.TotalCount"/> (denominator).
    /// Passing a deduplicated list here prevents a single PSM with multiple interpretations
    /// of the same peptide from inflating the denominator.
    /// When null, <paramref name="localizedSequences"/> is used for the denominator as well
    /// (legacy behaviour).
    /// </param>
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
        IEnumerable<IBioPolymerWithSetMods>? sequencesForTotalCount = null,
        Dictionary<string, double>? intensitiesByFullSequence = null)
    {
        var sequences = localizedSequences as IList<IBioPolymerWithSetMods> ?? localizedSequences.ToList();
        // coverageList is used only for TotalCount (denominator): one entry per PSM prevents a
        // single PSM with multiple interpretations of the same peptide from inflating the count.
        var coverageList = sequencesForTotalCount != null
            ? (sequencesForTotalCount as IList<IBioPolymerWithSetMods> ?? sequencesForTotalCount.ToList())
            : sequences;

        // Use an inner dictionary for dedup during construction, then flatten to lists
        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();

        // Cache per-position totals so they are computed once and shared across all mods at the
        // same site. Without this, TotalCount is recalculated for every new mod type at the same
        // position, and occupancies for competing mods at a single site can sum to >1.0.
        var positionTotals = new Dictionary<int, (int count, double intensity)>();

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

                // Compute total coverage for this position once and cache it
                if (!positionTotals.TryGetValue(indexInProtein, out var totals))
                {
                    int totalCount = 0;
                    double totalIntensity = 0;

                    foreach (var seq in coverageList)
                    {
                        int rangeStart = seq.OneBasedStartResidue - (indexInProtein == 1 ? 1 : 0);
                        if (indexInProtein >= rangeStart && indexInProtein <= seq.OneBasedEndResidue)
                        {
                            totalCount++;
                            if (intensitiesByFullSequence != null &&
                                seq.FullSequence != null &&
                                intensitiesByFullSequence.TryGetValue(seq.FullSequence, out double seqIntensity))
                            {
                                totalIntensity += seqIntensity;
                            }
                        }
                    }

                    totals = (totalCount, totalIntensity);
                    positionTotals[indexInProtein] = totals;
                }

                if (!modsAtPosition.TryGetValue(mod.Value.IdWithMotif, out var siteOccupancy))
                {
                    siteOccupancy = new SiteSpecificModificationOccupancy(indexInProtein, mod.Value.IdWithMotif)
                    {
                        TotalCount = totals.count,
                        TotalIntensity = totals.intensity
                    };

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
    /// Provides the forms used for <see cref="SiteSpecificModificationOccupancy.ModifiedCount"/> (numerator).
    /// </param>
    /// <param name="intensitiesByFullSequence">
    /// Optional map of FullSequence → intensity for intensity-based stoichiometry.
    /// </param>
    /// <param name="psmCount">
    /// Optional override for the total PSM count used as the denominator
    /// (<see cref="SiteSpecificModificationOccupancy.TotalCount"/>).
    /// When supplied, this value replaces <c>peptides.Count()</c>, preventing a single PSM
    /// with multiple interpretations of the same base sequence from inflating the denominator.
    /// When null, <c>peptides.Count()</c> is used (legacy behaviour).
    /// </param>
    /// <returns>
    /// Dictionary keyed by peptide-local position (AllModsOneIsNterminus key) containing a list of
    /// <see cref="SiteSpecificModificationOccupancy"/> entries for modifications observed at that position.
    /// </returns>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculatePeptideLevelOccupancy(
        IEnumerable<IBioPolymerWithSetMods> peptides,
        Dictionary<string, double>? intensitiesByFullSequence = null,
        int? psmCount = null)
    {
        var peptideList = peptides as IList<IBioPolymerWithSetMods> ?? peptides.ToList();
        // Use the caller-supplied PSM count when available so that a single PSM with multiple
        // interpretations of the same base sequence does not inflate the denominator.
        int totalPeptideCount = psmCount ?? peptideList.Count;

        double totalGroupIntensity = 0;
        if (intensitiesByFullSequence != null)
        {
            foreach (var p in peptideList)
            {
                if (p.FullSequence != null &&
                    intensitiesByFullSequence.TryGetValue(p.FullSequence, out double val))
                {
                    totalGroupIntensity += val;
                }
            }
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
