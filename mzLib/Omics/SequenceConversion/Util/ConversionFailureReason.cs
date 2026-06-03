namespace Omics.SequenceConversion;

/// <summary>
/// Describes why a sequence conversion failed.
/// </summary>
public enum ConversionFailureReason
{
    /// <summary>
    /// The provided sequence was null, empty, or otherwise invalid.
    /// </summary>
    InvalidSequence,

    /// <summary>
    /// The converter could not handle one or more modifications.
    /// </summary>
    IncompatibleModifications,

    /// <summary>
    /// The converter does not support the requested direction.
    /// </summary>
    UnsupportedDirection,

    /// <summary>
    /// The converter encountered an unexpected format or token.
    /// </summary>
    UnknownFormat
}
