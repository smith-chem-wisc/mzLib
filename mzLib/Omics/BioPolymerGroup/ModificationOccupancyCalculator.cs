using Omics.BioPolymer;
using Omics.Modifications;
using Omics.SpectralMatch;

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
    /// Calculates per-site modification occupancy mapped to protein coordinates directly from PSMs.
    /// PSM grouping, form filtering, TotalCount derivation, and intensity lookup are all handled internally.
    /// </summary>
    /// <param name="bioPolymer">The parent biopolymer whose length defines the coordinate space.</param>
    /// <param name="psms">
    /// All PSMs to consider. Forms are filtered to <paramref name="bioPolymer"/> internally.
    /// PSMs whose <see cref="ISpectralMatch.Intensities"/> is a single-element array contribute
    /// to intensity-based stoichiometry; others contribute only to count-based metrics.
    /// </param>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculateProteinLevelOccupancy(
        IBioPolymer bioPolymer,
        IEnumerable<ISpectralMatch> psms)
    {
        var psmList = psms as IList<ISpectralMatch> ?? psms.ToList();

        // Map each PSM to the single form that owns its intensity: matching FullSequence + Accession.
        var psmToBioPolymer = psmList
            .ToDictionary(
                p => p,
                p => p.GetIdentifiedBioPolymersWithSetMods()
                    .FirstOrDefault(s => s.FullSequence != null
                        && s.BaseSequence == p.BaseSequence
                        && s.FullSequence == p.FullSequence
                        && s.Parent.Accession == bioPolymer.Accession));

        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();
        var positionTotals = new Dictionary<int, (int totalCount, double totalIntensity)>();

        foreach (var psm in psmList)
        {
            var sequence = psmToBioPolymer[psm];
            if (sequence is null)  // PSM has no form for this protein, skip 
                continue;

            foreach (var mod in sequence.AllModsOneIsNterminus)
            {
                if (!TryGetProteinPosition(mod, sequence, bioPolymer.Length, out int indexInProtein))
                    continue;

                // Compute TotalCount/TotalIntensity for this position on first encounter only.
                if (!positionTotals.TryGetValue(indexInProtein, out var totals))
                {
                    var totalCount = 0;
                    var totalIntensity = 0.0;
                    foreach (var p in psmList)
                    {
                        var pSeq = psmToBioPolymer[p];
                        if (pSeq is null) 
                            continue;

                        int rangeStart = pSeq.OneBasedStartResidue - (indexInProtein == 1 ? 1 : 0);
                        if (indexInProtein >= rangeStart && indexInProtein <= pSeq.OneBasedEndResidue)
                        {
                            totalCount++;
                            if (p.Intensities is { Length: 1 })
                                totalIntensity += p.Intensities[0];
                        }
                    }
                    totals = (totalCount, totalIntensity);
                    positionTotals[indexInProtein] = totals;
                }

                if (!working.TryGetValue(indexInProtein, out var modsAtPosition))
                {
                    modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                    working[indexInProtein] = modsAtPosition;
                }

                if (!modsAtPosition.ContainsKey(mod.Value.IdWithMotif))
                {
                    modsAtPosition[mod.Value.IdWithMotif] = new SiteSpecificModificationOccupancy(indexInProtein, mod.Value.IdWithMotif)
                    {
                        TotalCount = totals.totalCount,
                        TotalIntensity = totals.totalIntensity
                    };
                }

                var siteOcc = modsAtPosition[mod.Value.IdWithMotif];
                siteOcc.ModifiedCount++;
                if (psm.Intensities is { Length: 1 })
                    siteOcc.ModifiedIntensity += psm.Intensities[0];
            }
        }

        return working.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.Values.ToList());
    }

    /// <summary>
    /// Calculates per-site modification occupancy in peptide-local coordinates directly from PSMs,
    /// returning results for all observed base sequences in a single call.
    /// PSM grouping, intensity derivation, and base-sequence bucketing are all handled internally.
    /// </summary>
    /// <param name="psms">
    /// All PSMs to consider. PSMs are grouped internally by <see cref="ISpectralMatch.BaseSequence"/>.
    /// PSMs whose <see cref="ISpectralMatch.Intensities"/> is a single-element array contribute
    /// to intensity-based stoichiometry; others contribute only to count-based metrics.
    /// </param>
    /// <returns>
    /// Dictionary keyed by base sequence, each value a dictionary keyed by peptide-local position
    /// (AllModsOneIsNterminus convention) containing <see cref="SiteSpecificModificationOccupancy"/> entries.
    /// </returns>
    public static Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>> CalculatePeptideLevelOccupancy(
        IEnumerable<ISpectralMatch> psms)
    {
        var psmList = psms as IList<ISpectralMatch> ?? psms.ToList();
        var result = new Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>>();

        foreach (var baseSeqGroup in psmList
            .Where(p => p.BaseSequence != null)
            .GroupBy(p => p.BaseSequence!))
        {
            // Map each PSM to the single form that owns its intensity: matching FullSequence.
            var psmToForm = baseSeqGroup
                .ToDictionary(
                    p => p,
                    p => p.GetIdentifiedBioPolymersWithSetMods()
                        .FirstOrDefault(s => s.FullSequence != null
                            && s.BaseSequence == baseSeqGroup.Key
                            && s.FullSequence == p.FullSequence));

            // All positions in a peptide share the same denominator: every PSM with this base
            // sequence covers every residue, so TotalCount/TotalIntensity are uniform across sites.
            var totalCount = 0;
            var totalIntensity = 0.0;
            foreach (var p in baseSeqGroup)
            {
                if (psmToForm[p] is null)
                    continue;
                totalCount++;
                if (p.Intensities is { Length: 1 })
                    totalIntensity += p.Intensities[0];
            }

            if (totalCount == 0)
                continue;

            var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();

            foreach (var psm in baseSeqGroup)
            {
                var form = psmToForm[psm];
                if (form is null)
                    continue;

                foreach (var mod in form.AllModsOneIsNterminus)
                {
                    if (IsExcludedMod(mod.Value))
                        continue;

                    if (!working.TryGetValue(mod.Key, out var modsAtPosition))
                    {
                        modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                        working[mod.Key] = modsAtPosition;
                    }

                    if (!modsAtPosition.ContainsKey(mod.Value.IdWithMotif))
                    {
                        modsAtPosition[mod.Value.IdWithMotif] = new SiteSpecificModificationOccupancy(mod.Key, mod.Value.IdWithMotif)
                        {
                            TotalCount = totalCount,
                            TotalIntensity = totalIntensity
                        };
                    }

                    var siteOcc = modsAtPosition[mod.Value.IdWithMotif];
                    siteOcc.ModifiedCount++;
                    if (psm.Intensities is { Length: 1 })
                        siteOcc.ModifiedIntensity += psm.Intensities[0];
                }
            }

            if (working.Count > 0)
                result[baseSeqGroup.Key] = working.ToDictionary(
                    kvp => kvp.Key,
                    kvp => kvp.Value.Values.ToList());
        }

        return result;
    }

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
