namespace Omics.SequenceConversion;

/// <summary>
/// Represents the outcome of a sequence conversion.
/// </summary>
/// <param name="Success">True when conversion completed without errors.</param>
/// <param name="FailureReason">Optional failure reason when unsuccessful.</param>
/// <param name="Warnings">Optional warnings produced during conversion.</param>
/// <param name="IncompatibleItems">Optional details about incompatible annotations.</param>
public sealed record ConversionResult(
    bool Success,
    ConversionFailureReason? FailureReason = null,
    IReadOnlyList<string>? Warnings = null,
    IReadOnlyList<string>? IncompatibleItems = null)
{
    public static ConversionResult Successful(IReadOnlyList<string>? warnings = null) =>
        new(true, null, warnings, null);

    public static ConversionResult Failure(ConversionFailureReason reason,
        IReadOnlyList<string>? warnings = null,
        IReadOnlyList<string>? incompatibleItems = null) =>
        new(false, reason, warnings, incompatibleItems);
}
