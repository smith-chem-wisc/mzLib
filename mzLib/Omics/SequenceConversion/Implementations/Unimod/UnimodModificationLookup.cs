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

    /// <summary>
    /// Adds one step to the base precedence: when the mzLib id names nothing among the Unimod entries
    /// (a UniProt or MetaMorpheus modification, say), the Unimod id the source modification names for
    /// itself, from its "Unimod" database reference or "UNIMOD:n" accession, picks the candidates.
    /// Without it the formula fallback picks any entry sharing the formula: Ethyl (UNIMOD:280) for
    /// N6,N6-dimethyllysine, whose reference is Dimethyl (UNIMOD:36).
    /// The reference only chooses among entries with the source's formula, so it never changes the
    /// mass: N,N-dimethylproline (C2H4) references Delta:H(5)C(2) (C2H5) and keeps the formula path.
    /// </summary>
    protected override IEnumerable<Modification> GetPrimaryCandidates(CanonicalModification mod)
    {
        var primary = base.GetPrimaryCandidates(mod).ToList();
        if (primary.Count > 0)
        {
            return primary;
        }

        var unimodId = mod.UnimodId ?? CanonicalModification.GetUnimodId(mod.MzLibModification);
        if (!unimodId.HasValue)
        {
            return [];
        }

        var byUnimodId = FilterByUnimodId(CandidateSet, unimodId.Value);
        var formula = mod.ChemicalFormula ?? mod.MzLibModification?.ChemicalFormula;
        return formula == null ? byUnimodId : FilterByFormula(byUnimodId, formula);
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
