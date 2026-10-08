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

    /// <summary>Instance that searches the protein modification catalogs only.</summary>
    public static GlobalModificationLookup ProteinOnly { get; } = new(Mods.AllProteinModsList, 0.001, "Global (Protein Mods)");

    private readonly string _name;

    /// <summary>
    /// Creates a new GlobalModificationLookup.
    /// </summary>
    /// <param name="massTolerance">Tolerance for mass-based matching in Daltons.</param>
    public GlobalModificationLookup(double massTolerance = 0.001)
        : this(Mods.AllKnownMods, massTolerance, "Global (All Mods)")
    {
    }

    private GlobalModificationLookup(IEnumerable<Modification> candidates, double massTolerance, string name)
        : base(candidates, massTolerance)
    {
        _name = name;
    }

    /// <inheritdoc />
    public override string Name => _name;

}
