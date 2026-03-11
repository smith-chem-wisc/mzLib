using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications to their UniProt representations by searching the UniProt PTM database.
/// </summary>
public class UniprotModificationLookup : ModificationLookupBase
{
    public static UniprotModificationLookup Instance { get; } = new();

    public UniprotModificationLookup(IEnumerable<Modification>? candidateSet = null, double massTolerance = 0.001)
        : base(
            conventionForLookup: ModificationNamingConvention.UniProt,
            searchProteinMods: true,
            searchRnaMods: false,
            massTolerance: massTolerance,
            candidateSet: candidateSet ?? Mods.UniprotModifications)
    {
    }

    public override string Name => "UniProt";

    protected override Modification? TryResolvePrimary(CanonicalModification mod)
    {
        if (mod.MzLibModification != null)
        {
            if (IsUniProt(mod.MzLibModification))
            {
                return mod.MzLibModification;
            }

            var residue = mod.TargetResidue ?? mod.MzLibModification.Target?.ToString().FirstOrDefault();
            var equivalent = FindBySourceModification(mod.MzLibModification, residue);
            if (equivalent != null)
            {
                return equivalent;
            }
        }

        return null;
    }

    protected override string NormalizeRepresentation(string representation)
    {
        var normalized = base.NormalizeRepresentation(representation);
        if (string.IsNullOrEmpty(normalized))
        {
            return normalized;
        }

        var colonIndex = normalized.IndexOf(':');
        if (colonIndex >= 0)
        {
            var prefix = normalized[..colonIndex].Trim();
            var remainder = normalized[(colonIndex + 1)..].Trim();

            if (prefix.Equals("UniProt", StringComparison.OrdinalIgnoreCase) || prefix.Contains(' '))
            {
                return string.IsNullOrEmpty(remainder) ? normalized : remainder;
            }
        }

        return normalized;
    }

    private static bool IsUniProt(Modification modification) =>
        modification.ModificationType != null &&
        modification.ModificationType.Equals("UniProt", StringComparison.OrdinalIgnoreCase);

    private Modification? FindBySourceModification(Modification source, char? residue)
    {
        var preferredResidue = residue ?? source.Target?.ToString().FirstOrDefault();

        if (source.ChemicalFormula != null)
        {
            var formulaMatches = FilterCandidates(m => m.ChemicalFormula != null && m.ChemicalFormula.Equals(source.ChemicalFormula));
            var match = SelectWithResiduePreference(formulaMatches, preferredResidue);
            if (match != null)
            {
                return match;
            }
        }

        if (source.MonoisotopicMass.HasValue && MassTolerance.HasValue)
        {
            var sourceMass = source.MonoisotopicMass.Value;
            var tolerance = MassTolerance.Value;
            var massMatches = FilterCandidates(m => m.MonoisotopicMass.HasValue &&
                                                   Math.Abs(m.MonoisotopicMass.Value - sourceMass) <= tolerance);
            var match = SelectWithResiduePreference(massMatches, preferredResidue);
            if (match != null)
            {
                return match;
            }
        }

        var identifiers = new List<string?>
        {
            source.IdWithMotif,
            source.OriginalId,
            source.ModificationType != null ? $"{source.ModificationType}:{source.IdWithMotif}" : null
        };

        foreach (var identifier in identifiers)
        {
            if (string.IsNullOrWhiteSpace(identifier))
            {
                continue;
            }

            var normalized = NormalizeRepresentation(identifier!);
            var resolved = ResolveByIdentifier(normalized);
            if (resolved != null)
            {
                return resolved;
            }
        }

        return null;
    }
}
