using Chemistry;
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

    /// <inheritdoc />
    public override string Name => "UNIMOD";

    /// <inheritdoc />
    protected override Modification? TryResolvePrimary(CanonicalModification mod)
    {
        // Try to resolve by UNIMOD ID if available
        if (mod.UnimodId.HasValue)
        {
            var lookupId = $"UNIMOD:{mod.UnimodId.Value}";
            return Mods.GetModification(lookupId, ModificationNamingConvention.Unimod);
        }

        return null;
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByName(string name, char? targetResidue)
    {
        if (string.IsNullOrWhiteSpace(name))
            return null;

        string lookupId;

        if (name.Contains("UNIMOD:", StringComparison.OrdinalIgnoreCase))
        {
            // Already in UNIMOD format
            lookupId = name;
        }
        else if (int.TryParse(name, out _))
        {
            // Numeric ID only - convert to UNIMOD format
            lookupId = $"UNIMOD:{name}";
        }
        else
        {
            // Try as a modification name
            lookupId = name;
        }

        return Mods.GetModification(lookupId, ModificationNamingConvention.Unimod);
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByFormula(ChemicalFormula formula, char? targetResidue)
    {
        // Search UNIMOD modifications for a formula match
        var candidates = Mods.UnimodModifications
            .Where(m => m.ChemicalFormula != null && m.ChemicalFormula.Equals(formula));

        return SelectWithResiduePreference(candidates, targetResidue);
    }
}
