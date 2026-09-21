namespace Omics.Modifications.IO.Modomics;

/// <summary>
/// A MODOMICS entry that cannot yet be represented as an mzLib modification, with the reason it was excluded.
/// </summary>
public class ModomicsNotYetRepresentableEntry
{
    /// <summary>MODOMICS numeric database id.</summary>
    public required int Id { get; init; }

    /// <summary>MODOMICS short name (e.g. "pm1G").</summary>
    public required string ShortName { get; init; }

    /// <summary>MODOMICS full name.</summary>
    public required string Name { get; init; }

    /// <summary>Raw source formula string.</summary>
    public required string Formula { get; init; }

    /// <summary>Moiety type from the MODOMICS CSV ("nucleoside", "nucleotide", "base").</summary>
    public required string MoietyType { get; init; }

    /// <summary>Reference moieties from the MODOMICS JSON.</summary>
    public required IReadOnlyList<string> ReferenceMoieties { get; init; }

    /// <summary>Why the entry could not be represented.</summary>
    public required ModomicsRepresentationFailureReason Reason { get; init; }

    /// <summary>Optional detail (e.g. the unsupported moiety or an exception message).</summary>
    public string? Details { get; init; }
}
