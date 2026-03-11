using System;
using System.Collections.Generic;
using Easy.Common.Extensions;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications using UNIMOD identifiers.
/// Supports formats like "UNIMOD:35", "35", or modification names that exist in the UNIMOD database.
/// </summary>
public class UnimodModificationLookup : ModificationLookupBase
{
    /// <summary>
    /// Singleton instance for convenience. Thread-safe due to static initialization.
    /// </summary>
    public static UnimodModificationLookup Instance { get; } = new();

    public UnimodModificationLookup(IEnumerable<Modification>? candidateSet = null)
        : base(
            conventionForLookup: ModificationNamingConvention.Unimod,
            searchProteinMods: true,
            searchRnaMods: false,
            massTolerance: null,
            candidateSet: candidateSet ?? Mods.UnimodModifications)
    {
    }

    /// <inheritdoc />
    public override string Name => "UNIMOD";

    /// <inheritdoc />
    protected override Modification? TryResolvePrimary(CanonicalModification mod)
    {
        // Try to resolve by UNIMOD ID if available
        if (mod.UnimodId.HasValue)
        {
            var withUnimodId = FilterCandidates(p => p.DatabaseReference.Any(p => p.Key.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase) && p.Value.Contains(mod.UnimodId.Value.ToString())));

            if (withUnimodId.IsNotNullOrEmpty())
                return withUnimodId.MaxBy(p => GetOverlapScore(mod.OriginalRepresentation, p.IdWithMotif));

            return ResolveByIdentifier($"UNIMOD:{mod.UnimodId.Value}");
        }

        return null;
    }

    protected override string NormalizeRepresentation(string representation)
    {
        var normalized = base.NormalizeRepresentation(representation);
        if (string.IsNullOrEmpty(normalized))
            return normalized;

        var index = normalized.IndexOf("UNIMOD:", StringComparison.OrdinalIgnoreCase);
        if (index >= 0)
        {
            var id = normalized.Substring(index + 7).Trim();
            return string.IsNullOrEmpty(id) ? normalized : $"UNIMOD:{id}";
        }

        if (int.TryParse(normalized, out _))
        {
            return $"UNIMOD:{normalized}";
        }

        return normalized;
    }
}
