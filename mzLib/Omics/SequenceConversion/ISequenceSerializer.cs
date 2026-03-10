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
}
