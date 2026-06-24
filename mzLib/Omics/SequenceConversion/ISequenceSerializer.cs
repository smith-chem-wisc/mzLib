using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Serializes a <see cref="CanonicalSequence"/> into a string representation.
/// Each implementation produces output in a specific format (e.g., mzLib, UNIMOD, mass-shift, Chronologer).
/// </summary>
public interface ISequenceSerializer
{
    /// <summary>
    /// Gets the unique name of the format this serializer produces (e.g., "mzLib", "Chronologer").
    /// </summary>
    string FormatName { get; }

    /// <summary>
    /// Gets the schema describing the syntax of this format.
    /// </summary>
    SequenceFormatSchema Schema { get; }

    /// <summary>
    /// Gets the modification lookup used for resolving identifiers during serialization.
    /// </summary>
    IModificationLookup? ModificationLookup { get; }

    /// <summary>
    /// Gets the preferred handling mode for incompatible modifications.
    /// </summary>
    SequenceConversionHandlingMode HandlingMode { get; }

    /// <summary>
    /// Determines whether this serializer can handle the given sequence.
    /// Some serializers may not support certain modification types or terminal modifications.
    /// </summary>
    /// <param name="sequence">The sequence to check.</param>
    /// <returns>True if this serializer can handle the sequence; otherwise, false.</returns>
    bool CanSerialize(CanonicalSequence sequence);

    /// <summary>
    /// Determines whether a modification should be resolved/enriched before serialization.
    /// This allows each serializer to resolve only the missing fields required by its format.
    /// </summary>
    /// <param name="mod">The modification to evaluate.</param>
    /// <returns>True if this modification should be resolved; otherwise, false.</returns>
    bool ShouldResolveMod(CanonicalModification mod);

    /// <summary>
    /// Serializes the sequence into a string representation.
    /// </summary>
    /// <param name="sequence">The sequence to serialize.</param>
    /// <param name="warnings">Optional. Accumulates warnings and errors during serialization.
    /// If null, a new instance will be created internally if issues occur.</param>
    /// <param name="mode">Specifies how to handle serialization issues (throw, return null, etc.).</param>
    /// <returns>The serialized string, or null if serialization failed and mode allows null returns.</returns>
    /// <exception cref="SequenceConversionException">Thrown when serialization fails and mode is ThrowException.</exception>
    string? Serialize(
        CanonicalSequence sequence,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Projects a canonical sequence into a OneIsNterminus modification dictionary.
    /// Indexing convention: 1 = N-terminus, 2 = first residue, Length + 2 = C-terminus.
    /// </summary>
    /// <param name="sequence">The canonical sequence to project.</param>
    /// <param name="knownMods">Optional override dictionary used for final resolution.
    /// If null, serializer lookup and canonical metadata are used.</param>
    /// <param name="warnings">Optional warnings collection.</param>
    /// <param name="mode">Specifies how to handle incompatible/ambiguous modifications.</param>
    /// <returns>A modification dictionary keyed by OneIsNterminus positions.</returns>
    Dictionary<int, Modification> ToOneIsNterminusModificationDictionary(
        CanonicalSequence sequence,
        Dictionary<string, Modification>? knownMods = null,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException);
}
