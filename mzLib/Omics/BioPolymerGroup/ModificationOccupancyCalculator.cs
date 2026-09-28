using CsvHelper.Configuration.Attributes;
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

        // Pre-compute the matching form for each PSM by position.
        var psmForms = psmList
            .Select(p => p.GetIdentifiedBioPolymersWithSetMods()
                .FirstOrDefault(s => s.FullSequence != null
                    && s.BaseSequence == p.BaseSequence
                    && s.FullSequence == p.FullSequence
                    && s.Parent.Accession == bioPolymer.Accession))
            .ToArray();

        var positionTotals = new Dictionary<int, (int totalCount, double totalIntensity)>();
        for (int j = 0; j < psmList.Count; j++)
        {
            var psm = psmList[j];
            var sequence = psmForms[j];
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

            // A form that starts after the removed initiator Met covers the protein N-terminus (position 1)
            // but not residue 1 (position 2), so the N-terminus is counted on its own.
            if (StartsAfterInitiatorMethionine(sequence, bioPolymer))
                AddToTotals(positionTotals, 1, psm);

            int rangeStart = sequence.OneBasedStartResidue + (sequence.OneBasedStartResidue == 1 ? 0 : 1); // Include position 1 if sequence starts at the protein N-terminus
            int rangeEnd = sequence.OneBasedEndResidue + (sequence.OneBasedEndResidue == bioPolymer.Length ? 2 : 1); // Include last position if sequence ends at the protein C-terminus
            for (int i = rangeStart; i <= rangeEnd; i++)
                AddToTotals(positionTotals, i, psm);
        }

        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();
        for (int j = 0; j < psmList.Count; j++)
        {
            var psm = psmList[j];
            var sequence = psmForms[j];
            if (sequence is null)  // PSM has no form for this protein, skip
                continue;

            foreach (var mod in sequence.AllModsOneIsNterminus)
            {
                if (IsExcludedMod(mod.Value))
                    continue;

                if (!TryGetProteinPosition(mod, sequence, bioPolymer, out int indexInProtein))
                    continue;

                if (!working.TryGetValue(indexInProtein, out var modsAtPosition))
                {
                    modsAtPosition = new Dictionary<string, SiteSpecificModificationOccupancy>();
                    working[indexInProtein] = modsAtPosition;
                }

                if (!modsAtPosition.ContainsKey(mod.Value.IdWithMotif))
                {
                    if (!positionTotals.TryGetValue(indexInProtein, out var posTotals))
                        continue;

                    modsAtPosition[mod.Value.IdWithMotif] = new SiteSpecificModificationOccupancy(indexInProtein, mod.Value.IdWithMotif)
                    {
                        TotalCount = posTotals.totalCount,
                        TotalIntensity = posTotals.totalIntensity
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

        // Unresolved has two spellings: an implementation may leave BaseSequence null - MetaMorpheus
        // does, for a PSM ambiguous across base sequences - while BaseSpectralMatch coerces that null
        // to empty. Both mean the same thing here, and an empty one left in trips the check below.
        var psmsWithBaseSeq = psmList.Where(p => !string.IsNullOrEmpty(p.BaseSequence)).ToList();

        // Nothing resolved means nothing to attribute a modification to, which is not an error.
        // Returning here also keeps AllSame() off an empty sequence, where its First() would throw.
        if (psmsWithBaseSeq.Count == 0)
            return new Dictionary<int, List<SiteSpecificModificationOccupancy>>();

        if (!psmsWithBaseSeq.Select(p => p.BaseSequence).AllSame())
        {
            throw new ArgumentException("All PSMs must have the same BaseSequence for peptide-level occupancy calculation.");
        }

        var totalCount = psmsWithBaseSeq.Count;
        var totalIntensity = psmsWithBaseSeq
            .Where(p => p.Intensities is { Length: 1 })
            .Sum(p => p.Intensities[0]);

        var working = new Dictionary<int, Dictionary<string, SiteSpecificModificationOccupancy>>();
        foreach (var psm in psmsWithBaseSeq)
        {
            // The form carrying this PSM's modifications. A PSM whose full sequence matches no
            // identified form is ambiguous: it counts toward the denominator but marks no site.
            var form = psm.GetIdentifiedBioPolymersWithSetMods()
                .FirstOrDefault(s => s.FullSequence == psm.FullSequence);

            if (form is null)
                continue;

            foreach (var mod in form.AllModsOneIsNterminus)
            {
                if (IsExcludedMod(mod.Value, ignoreLocation: true))
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

    private static void AddToTotals(Dictionary<int, (int totalCount, double totalIntensity)> positionTotals,
        int position, ISpectralMatch psm)
    {
        positionTotals.TryGetValue(position, out var totals);
        totals.totalCount++;
        if (psm.Intensities is { Length: 1 })
            totals.totalIntensity += psm.Intensities[0];
        positionTotals[position] = totals;
    }

    /// <summary>
    /// True when <paramref name="sequence"/> begins at residue 2 because the initiator Met was removed, which is
    /// the same condition under which digestion produces such a form (Protease: residue 1 must be 'M'). Its
    /// N-terminus is then the protein N-terminus, and ModificationLocalization places "N-terminal." mods there.
    /// </summary>
    /// <remarks>
    /// Not handled: a protease that cleaves C-terminal to Met (e.g. CNBr, "M|") also yields a form starting at
    /// residue 2 from a Met-retained molecule. That form has the same base sequence and span as the Met-removed
    /// N-terminus, so digestion produces one form for both and nothing here can tell them apart; it is counted
    /// as the protein N-terminus. For such proteases, protein N-terminal occupancy may be understated.
    /// </remarks>
    private static bool StartsAfterInitiatorMethionine(IBioPolymerWithSetMods sequence, IBioPolymer bioPolymer)
        => sequence.OneBasedStartResidue == 2
           && bioPolymer.BaseSequence.Length > 0
           && bioPolymer.BaseSequence[0] == 'M';

    private static bool TryGetProteinPosition(
        KeyValuePair<int, Modification> mod,
        IBioPolymerWithSetMods sequence,
        IBioPolymer bioPolymer,
        out int indexInProtein)
    {
        indexInProtein = 0;
        int bioPolymerLength = bioPolymer.Length;

        if (IsExcludedMod(mod.Value))
            return false;

        if (mod.Value.LocationRestriction.Equals("N-terminal."))
        {
            if (sequence.OneBasedStartResidue != 1 && !StartsAfterInitiatorMethionine(sequence, bioPolymer))
                return false;

            indexInProtein = 1;
        }
        else if (mod.Value.LocationRestriction.Equals("Anywhere."))
        {
            indexInProtein = sequence.OneBasedStartResidue + mod.Key - 1;
        }
        else if (mod.Value.LocationRestriction.Equals("C-terminal."))
        {
            if (sequence.OneBasedEndResidue != bioPolymerLength)
                return false;

            indexInProtein = bioPolymerLength + 2;
        }
        else
        {
            return false;
        }

        return true;
    }

    private static bool IsExcludedMod(Modification mod, bool ignoreLocation = false)
    {
        if (ExcludedLocations.Contains(mod.LocationRestriction) && !ignoreLocation)
            return true;

        if (ExcludedModTypes.Contains(mod.ModificationType))
            return true;

        return false;
    }
}
