namespace Omics.SequenceConversion;

/// <summary>
/// Defines how converters should respond when encountering incompatible details.
/// </summary>
public enum SequenceConversionHandlingMode
{
    /// <summary>
    /// Converter throws an exception immediately.
    /// </summary>
    ThrowException,

    /// <summary>
    /// Converter returns null without throwing.
    /// </summary>
    ReturnNull,

    /// <summary>
    /// Converter removes incompatible annotations before returning such as modifications it cannot handle, with warnings about what was removed.
    /// </summary>
    RemoveIncompatibleElements,

    /// <summary>
    /// Converter falls back to the underlying primary sequence when available.
    /// </summary>
    UsePrimarySequence
}
