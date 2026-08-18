using System;
using System.Collections.Generic;
using System.Linq;
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
        : base(candidateSet ?? Mods.UnimodModifications, null)
    {
    }

    /// <inheritdoc />
    public override string Name => "UNIMOD";

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
