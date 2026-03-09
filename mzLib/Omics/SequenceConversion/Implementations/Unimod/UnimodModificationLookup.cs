using System;
using System.Collections.Generic;
using Omics.Modifications;

namespace Omics.SequenceConversion.Implementations.Unimod;

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

    public UnimodModificationLookup()
        : base(
            conventionForLookup: ModificationNamingConvention.Unimod,
            searchProteinMods: true,
            searchRnaMods: false,
            massTolerance: null,
            candidateSet: Mods.UnimodModifications)
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

    protected override IEnumerable<string> ExpandNameCandidates(string normalizedRepresentation, char? targetResidue)
    {
        if (string.IsNullOrEmpty(normalizedRepresentation))
            yield break;

        yield return normalizedRepresentation;
    }
}
