namespace Omics.SequenceConversion;

/// <summary>
/// Exception thrown when a converter cannot satisfy the requested handling mode.
/// </summary>
public class SequenceConversionException : Exception
{
    public SequenceConversionException(string message, ConversionFailureReason failureReason,
        IReadOnlyList<string>? incompatibleItems = null, IReadOnlyList<string>? warnings = null, Exception? innerException = null)
        : base(message, innerException)
    {
        FailureReason = failureReason;
        IncompatibleItems = incompatibleItems ?? Array.Empty<string>();
        Warnings = warnings ?? Array.Empty<string>();
    }

    public ConversionFailureReason FailureReason { get; }

    public IReadOnlyList<string> IncompatibleItems { get; }

    public IReadOnlyList<string> Warnings { get; }
}
