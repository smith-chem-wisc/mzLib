namespace Omics.Modifications.IO.Modomics;

/// <summary>
/// Structured report for a MODOMICS load: the modifications that were loaded, and the entries that
/// were tracked instead of loaded, with the reasons why.
/// </summary>
public class ModomicsLoadResult
{
    /// <summary>All loaded MODOMICS modifications (includes terminal caps).</summary>
    public List<Modification> LoadedModifications { get; init; } = [];

    /// <summary>Loaded modifications restricted to the 5' terminus.</summary>
    public List<Modification> TerminalModifications { get; init; } = [];

    /// <summary>Chemistry-equivalent overlaps with curated mods. Both entries are kept; this is reporting only.</summary>
    public List<ModomicsDuplicateEntry> DuplicateModifications { get; init; } = [];

    /// <summary>Entries that cannot yet be represented, with reasons.</summary>
    public List<ModomicsNotYetRepresentableEntry> NotYetRepresentableEntries { get; init; } = [];

    /// <summary>Number of loaded modifications.</summary>
    public int LoadedCount => LoadedModifications.Count;

    /// <summary>Number of loaded terminal modifications.</summary>
    public int TerminalCount => TerminalModifications.Count;

    /// <summary>Number of chemistry-equivalent overlaps with curated mods.</summary>
    public int DuplicateCount => DuplicateModifications.Count;

    /// <summary>Number of entries that cannot yet be represented.</summary>
    public int NotYetRepresentableCount => NotYetRepresentableEntries.Count;
}
