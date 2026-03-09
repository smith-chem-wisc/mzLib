namespace Omics.SequenceConversion;

/// <summary>
/// Shared helper methods for sequence parsers and serializers.
/// </summary>
public static class SequenceConversionHelpers
{
    /// <summary>
    /// Handles an error for parsers (returns nullable CanonicalSequence).
    /// </summary>
    /// <param name="warnings">Warnings collection to add the failure to.</param>
    /// <param name="mode">How to handle the error.</param>
    /// <param name="reason">The reason for the failure.</param>
    /// <param name="message">The error message.</param>
    /// <returns>Null in most cases, or throws an exception if mode is ThrowException.</returns>
    public static CanonicalSequence? HandleParserError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        warnings.SetFailure(reason, message);

        return mode switch
        {
            SequenceConversionHandlingMode.ThrowException =>
                throw warnings.ToException(message),
            SequenceConversionHandlingMode.ReturnNull => null,
            _ => null
        };
    }

    /// <summary>
    /// Handles an error for serializers (returns nullable string).
    /// </summary>
    /// <param name="warnings">Warnings collection to add the failure to.</param>
    /// <param name="mode">How to handle the error.</param>
    /// <param name="reason">The reason for the failure.</param>
    /// <param name="message">The error message.</param>
    /// <returns>Null in most cases, or throws an exception if mode is ThrowException.</returns>
    public static string? HandleSerializerError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        warnings.SetFailure(reason, message);

        return mode switch
        {
            SequenceConversionHandlingMode.ThrowException =>
                throw warnings.ToException(message),
            SequenceConversionHandlingMode.ReturnNull => null,
            SequenceConversionHandlingMode.UsePrimarySequence => null, // Caller should handle
            _ => null
        };
    }
}
