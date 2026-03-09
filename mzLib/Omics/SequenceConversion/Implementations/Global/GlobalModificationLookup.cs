using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications using ALL known modifications from the mzLib modification database.
/// Searches across all modification sources (MetaMorpheus, UniProt, UNIMOD, RNA mods).
/// This is the most comprehensive lookup that searches the entire modification database.
/// </summary>
public class GlobalModificationLookup : ModificationLookupBase
{
    private readonly double _massTolerance;

    /// <summary>
    /// Singleton instance that searches all known modifications.
    /// </summary>
    public static GlobalModificationLookup Instance { get; } = new();

    /// <summary>
    /// Creates a new GlobalModificationLookup.
    /// </summary>
    /// <param name="massTolerance">Tolerance for mass-based matching in Daltons.</param>
    public GlobalModificationLookup(double massTolerance = 0.001)
    {
        _massTolerance = massTolerance;
    }

    /// <inheritdoc />
    public override string Name => "Global (All Mods)";

    /// <inheritdoc />
    protected override Modification? TryResolvePrimary(CanonicalModification mod)
    {
        Modification? toReturn = null;

        // Try to resolve by mzLib ID if available
        if (!string.IsNullOrEmpty(mod.MzLibId))
            toReturn = Mods.GetModification(mod.MzLibId, ModificationNamingConvention.Mixed);

        // Try to resolve by UNIMOD ID if available
        // UNIMOD IDs are stored in DatabaseReference["Unimod"] as just the numeric ID
        if (toReturn == null && mod.UnimodId.HasValue)
        {
            var unimodIdString = mod.UnimodId.Value.ToString();
            var candidates = Mods.AllKnownMods.Where(m =>
                m.DatabaseReference != null &&
                m.DatabaseReference.TryGetValue("Unimod", out var ids) &&
                ids.Contains(unimodIdString));
            
            // Prefer modification matching the target residue if specified
            toReturn = SelectWithResiduePreference(candidates, mod.TargetResidue);
        }

        return toReturn;
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByName(string name, char? targetResidue)
    {
        if (string.IsNullOrWhiteSpace(name))
            return null;

        var mod = Mods.GetModification(name, ModificationNamingConvention.Mixed);
        if (mod != null)
            return mod;

        // Try exact match first using AllModsKnownDictionary
        if (Mods.AllModsKnownDictionary.TryGetValue(name, out mod))
            return mod;

        // Try adding motif suffix if target residue is known
        if (targetResidue.HasValue)
        {
            var withMotif = $"{name} on {targetResidue.Value}";
            mod = Mods.GetModification(withMotif, ModificationNamingConvention.Mixed);
            if (mod != null)
                return mod;
        }
        // Try with ModificationType prefix
        mod = Mods.AllKnownMods.FirstOrDefault(m => $"{m.ModificationType}:{m.IdWithMotif}" == name);
        if (mod != null)
            return mod;

        // Try adding motif suffix if target residue is known
        if (targetResidue.HasValue)
        {
            var withMotif = $"{name} on {targetResidue.Value}";
            if (Mods.AllModsKnownDictionary.TryGetValue(withMotif, out mod))
                return mod;
        }

        // Try searching by OriginalId across all modifications
        var candidates = Mods.AllKnownMods
            .Where(m => m.OriginalId == name || m.IdWithMotif == name);

        return SelectWithResiduePreference(candidates, targetResidue);
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByFormula(ChemicalFormula formula, char? targetResidue)
    {
        var candidates = Mods.AllKnownMods
            .Where(m => m.ChemicalFormula != null && m.ChemicalFormula.Equals(formula));

        return SelectWithResiduePreference(candidates, targetResidue);
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByMass(double mass, char? targetResidue)
    {
        var candidates = Mods.AllKnownMods
            .Where(m => m.MonoisotopicMass.HasValue &&
                        Math.Abs(m.MonoisotopicMass.Value - mass) <= _massTolerance);

        return SelectWithResiduePreference(candidates, targetResidue);
    }
}
