using System.Linq;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications using ALL known modifications from the mzLib modification database.
/// Searches across all modification sources (MetaMorpheus, UniProt, UNIMOD, RNA mods).
/// This is the most comprehensive lookup that searches the entire modification database.
/// </summary>
public class GlobalModificationLookup : ModificationLookupBase
{

    /// <summary>
    /// Singleton instance that searches all known modifications.
    /// </summary>
    public static GlobalModificationLookup Instance { get; } = new();

    /// <summary>
    /// Creates a new GlobalModificationLookup.
    /// </summary>
    /// <param name="massTolerance">Tolerance for mass-based matching in Daltons.</param>
    public GlobalModificationLookup(double massTolerance = 0.001)
        : base(
            conventionForLookup: ModificationNamingConvention.Mixed,
            searchProteinMods: true,
            searchRnaMods: true,
            massTolerance: massTolerance,
            candidateSet: Mods.AllKnownMods)
    {
    }

    /// <inheritdoc />
    public override string Name => "Global (All Mods)";

    /// <inheritdoc />
    protected override Modification? TryResolvePrimary(CanonicalModification mod)
    {
        if (!string.IsNullOrEmpty(mod.MzLibId))
        {
            var resolved = ResolveByIdentifier(mod.MzLibId);
            if (resolved != null)
                return resolved;
        }

        // Try to resolve by UNIMOD ID if available
        // UNIMOD IDs are stored in DatabaseReference["Unimod"] as just the numeric ID
        if (mod.UnimodId.HasValue)
        {
            var unimodIdString = mod.UnimodId.Value.ToString();
            var candidates = Mods.AllKnownMods.Where(m =>
                m.DatabaseReference != null &&
                m.DatabaseReference.TryGetValue("Unimod", out var ids) &&
                ids.Contains(unimodIdString));
            
            // Prefer modification matching the target residue if specified
            return SelectWithResiduePreference(candidates, mod.TargetResidue);
        }

        return null;
    }
}
