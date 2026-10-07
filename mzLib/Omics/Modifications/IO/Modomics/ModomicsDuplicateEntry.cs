namespace Omics.Modifications.IO.Modomics;

/// <summary>
/// A chemistry-equivalent overlap between a MODOMICS modification and existing curated mods.
/// Both entries are kept; this record exists for reporting only.
/// </summary>
public class ModomicsDuplicateEntry
{
    /// <summary>The curated modifications equivalent to the MODOMICS entry.</summary>
    public required IReadOnlyList<Modification> ExistingModifications { get; init; }

    /// <summary>The MODOMICS-derived modification.</summary>
    public required Modification ModomicsModification { get; init; }
}
