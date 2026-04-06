using MzLibUtil;
using Omics.BioPolymer;
using Omics.Modifications;
using Omics.SpectralMatch;
using System.Numerics;

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
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculateParentLevelOccupancy(
        IBioPolymer bioPolymer,
        IEnumerable<ISpectralMatch> psms)
    {
        var psmList = psms as IList<ISpectralMatch> ?? psms.ToList();

        // Map each PSM to the single form that owns its intensity: matching FullSequence + Accession.
        var psmToBioPolymer = psmList
            .ToHashSet() // Ensure distinct PSMs in case of duplicates in input, since we're using PSMs as keys in a dictionary
            .ToDictionary(
                p => p,
                p => p.GetIdentifiedBioPolymersWithSetMods()
                    .FirstOrDefault(s => s.FullSequence != null
                        && s.BaseSequence == p.BaseSequence
                        && s.FullSequence == p.FullSequence
                        && s.Parent.Accession == bioPolymer.Accession));

        var positionTotals = new Dictionary<int, (int totalCount, double totalIntensity)>();
        foreach (var psm in psmList)
        {
            var sequence = psmToBioPolymer[psm];
            if (sequence is null) // PSM for this protein might be ambiguous (e.g. missing full sequence)
            {
                try
                {
                    // Still want to count it toward TotalCount/TotalIntensity for any positions it covers,
                    // so find the best-matching form without the full sequence requirement.
                    sequence = psm.GetIdentifiedBioPolymersWithSetMods()
                        .FirstOrDefault(s => s.BaseSequence == psm.BaseSequence
                            && s.Parent.Accession == bioPolymer.Accession);
                }
                catch (Exception)
                {
                    continue; // If we can't find any form for this PSM, skip it entirely.
                }
            }

            if (sequence is null) // No form found for this PSM, skip it entirely.
                continue;

            int rangeStart = sequence.OneBasedStartResidue - (sequence.OneBasedStartResidue == 1 ? 1 : 0); // Include position 1 if sequence starts at the protein N-terminus
            int rangeEnd = sequence.OneBasedEndResidue + (sequence.OneBasedEndResidue == bioPolymer.Length ? 1 : 0); // Include last position if sequence ends at the protein C-terminus
            for (int i = rangeStart; i <= rangeEnd; i++)
            {
                if (!positionTotals.ContainsKey(i))
                    positionTotals[i] = (0, 0.0);
                var totals = positionTotals[i];
                totals.totalCount++;
                if (psm.Intensities is { Length: 1 })
                    totals.totalIntensity += psm.Intensities[0];
                positionTotals[i] = totals;
            }
        }

        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();
        foreach (var psm in psmList)
        {
            var sequence = psmToBioPolymer[psm];
            if (sequence is null)  // PSM has no form for this protein, skip 
                continue;

            foreach (var mod in sequence.AllModsOneIsNterminus)
            {
                if (IsExcludedMod(mod.Value))
                    continue;

                if (!TryGetProteinPosition(mod, sequence, bioPolymer.Length, out int indexInProtein))
                    continue;

                if (!working.TryGetValue(indexInProtein, out var modsAtPosition))
                {
                    modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                    working[indexInProtein] = modsAtPosition;
                }

                if (!modsAtPosition.ContainsKey(mod.Value.IdWithMotif))
                {
                    modsAtPosition[mod.Value.IdWithMotif] = new SiteSpecificModificationOccupancy(indexInProtein, mod.Value.IdWithMotif)
                    {
                        TotalCount = positionTotals[indexInProtein].totalCount,
                        TotalIntensity = positionTotals[indexInProtein].totalIntensity
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
    /// Dictionary keyed by peptide-local position (AllModsOneIsNterminus convention) containing 
    /// <see cref="SiteSpecificModificationOccupancy"/> entries.
    /// </returns>
    public static Dictionary<int, List<SiteSpecificModificationOccupancy>> CalculateDigestionProductLevelOccupancy(
        IEnumerable<ISpectralMatch> psms)
    {
        var psmList = psms as IList<ISpectralMatch> ?? psms.ToList();
        var result = new Dictionary<string, Dictionary<int, List<SiteSpecificModificationOccupancy>>>();

        var psmsWithBaseSeq = psmList.Where(p => p.BaseSequence != null).ToList();

        if (!psmsWithBaseSeq.Select(p => p.BaseSequence).AllSame())
        {
            throw new ArgumentException("All PSMs must have the same BaseSequence for peptide-level occupancy calculation.");
        }

        // Map each PSM to the single form that owns its intensity: matching FullSequence + Accession.
        // Ambiguous forms (psms without a full sequence match) are filtered out and do not contribute to occupancy.
        var psmToForm = psmsWithBaseSeq
            .ToDictionary(
                p => p,
                p => p.GetIdentifiedBioPolymersWithSetMods()
                    .FirstOrDefault(s => s.FullSequence == p.FullSequence));

        var totalCount = psmsWithBaseSeq.Count;
        var totalIntensity = psmsWithBaseSeq
            .Where(p => p.Intensities is { Length: 1 })
            .Sum(p => p.Intensities[0]);

        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();
        foreach (var psm in psmsWithBaseSeq)
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

        if (working.Count == 0)
            return new Dictionary<int, List<SiteSpecificModificationOccupancy>>(); // Return empty if no mods passed filtering

        return working.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.Values.ToList());
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
        if (ExcludedLocations.Contains(mod.LocationRestriction))
            return true;

        if (ExcludedModTypes.Contains(mod.ModificationType))
            return true;

        return false;
    }
}
