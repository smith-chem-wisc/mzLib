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
        : base(Mods.AllKnownMods, massTolerance)
    {
    }

    /// <inheritdoc />
    public override string Name => "Global (All Mods)";

}
