using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Resolves modifications using the mzLib modification database.
/// Searches across all known modifications (MetaMorpheus, UniProt, UNIMOD, RNA mods).
/// Supports mzLib-style identifiers like "Oxidation on M", "Carbamidomethyl on C", etc.
/// </summary>
public class MzLibModificationLookup : ModificationLookupBase
{
    private readonly bool _searchProteinMods;
    private readonly bool _searchRnaMods;
    private readonly double _massTolerance;

    /// <summary>
    /// Singleton instance that searches both protein and RNA modifications.
    /// </summary>
    public static MzLibModificationLookup Instance { get; } = new();

    /// <summary>
    /// Instance that only searches protein modifications.
    /// </summary>
    public static MzLibModificationLookup ProteinOnly { get; } = new(searchProteinMods: true, searchRnaMods: false);

    /// <summary>
    /// Instance that only searches RNA modifications.
    /// </summary>
    public static MzLibModificationLookup RnaOnly { get; } = new(searchProteinMods: false, searchRnaMods: true);

    /// <summary>
    /// Creates a new MzLibModificationLookup.
    /// </summary>
    /// <param name="searchProteinMods">Whether to search protein modifications.</param>
    /// <param name="searchRnaMods">Whether to search RNA modifications.</param>
    /// <param name="massTolerance">Tolerance for mass-based matching in Daltons.</param>
    public MzLibModificationLookup(bool searchProteinMods = true, bool searchRnaMods = true, double massTolerance = 0.001)
    {
        _searchProteinMods = searchProteinMods;
        _searchRnaMods = searchRnaMods;
        _massTolerance = massTolerance;
    }

    /// <inheritdoc />
    public override string Name => "mzLib";

    /// <inheritdoc />
    protected override Modification? TryResolvePrimary(CanonicalModification mod)
    {
        // Try to resolve by mzLib ID if available
        if (!string.IsNullOrEmpty(mod.MzLibId))
        {
            return Mods.GetModification(mod.MzLibId, _searchProteinMods, _searchRnaMods);
        }

        return null;
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByName(string name, char? targetResidue)
    {
        if (string.IsNullOrWhiteSpace(name))
            return null;

        // Try exact match first (IdWithMotif)
        var mod = Mods.GetModification(name, _searchProteinMods, _searchRnaMods);
        if (mod != null)
            return mod;

        // Try adding motif suffix if target residue is known
        if (targetResidue.HasValue)
        {
            var withMotif = $"{name} on {targetResidue.Value}";
            mod = Mods.GetModification(withMotif, _searchProteinMods, _searchRnaMods);
            if (mod != null)
                return mod;
        }

        // Try searching by OriginalId
        var candidates = Mods.GetModifications(
            m => m.OriginalId == name || m.IdWithMotif == name,
            proteinOnly: !_searchRnaMods && _searchProteinMods,
            rnaOnly: !_searchProteinMods && _searchRnaMods);

        return SelectWithResiduePreference(candidates, targetResidue);
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByFormula(ChemicalFormula formula, char? targetResidue)
    {
        var candidates = Mods.GetModifications(
            m => m.ChemicalFormula != null && m.ChemicalFormula.Equals(formula),
            proteinOnly: !_searchRnaMods && _searchProteinMods,
            rnaOnly: !_searchProteinMods && _searchRnaMods);

        return SelectWithResiduePreference(candidates, targetResidue);
    }

    /// <inheritdoc />
    protected override Modification? TryResolveByMass(double mass, char? targetResidue)
    {
        var candidates = Mods.GetModifications(
            m => m.MonoisotopicMass.HasValue &&
                 Math.Abs(m.MonoisotopicMass.Value - mass) <= _massTolerance,
            proteinOnly: !_searchRnaMods && _searchProteinMods,
            rnaOnly: !_searchProteinMods && _searchRnaMods);

        return SelectWithResiduePreference(candidates, targetResidue);
    }
}
