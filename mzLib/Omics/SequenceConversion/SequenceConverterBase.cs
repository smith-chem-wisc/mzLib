namespace Omics.SequenceConversion;

/// <summary>
/// Provides common plumbing for converters, including handling modes and result creation.
/// </summary>
public abstract class SequenceConverterBase : ISequenceConverter
{
    protected SequenceConverterBase()
    {
        HandlingMode = DefaultHandlingMode;
    }

    public abstract string ConverterName { get; }

    public virtual bool IsBidirectional => false;

    public virtual SequenceConversionHandlingMode DefaultHandlingMode => SequenceConversionHandlingMode.RemoveIncompatibleElements;

    public SequenceConversionHandlingMode HandlingMode { get; set; }

    public virtual string? ConvertTo(string fullSequence, out ConversionResult result)
    {
        if (string.IsNullOrWhiteSpace(fullSequence))
        {
            result = ConversionResult.Failure(ConversionFailureReason.InvalidSequence);
            return null;
        }

        try
        {
            var warnings = new List<string>();
            string converted = ConvertSequenceCore(fullSequence, HandlingMode, warnings);
            result = ConversionResult.Successful(warnings.Count == 0 ? Array.Empty<string>() : warnings);
            return converted;
        }
        catch (SequenceConversionException ex)
        {
            if (HandlingMode == SequenceConversionHandlingMode.ThrowException)
            {
                throw;
            }

            result = ConversionResult.Failure(ex.FailureReason, ex.Warnings, ex.IncompatibleItems);
            return null;
        }
    }

    public virtual string? ConvertFrom(string targetSequence, out ConversionResult result)
    {
        if (!IsBidirectional)
        {
            result = ConversionResult.Failure(ConversionFailureReason.UnsupportedDirection);
            return null;
        }

        if (string.IsNullOrWhiteSpace(targetSequence))
        {
            result = ConversionResult.Failure(ConversionFailureReason.InvalidSequence);
            return null;
        }

        try
        {
            var warnings = new List<string>();
            var converted = ConvertFromCore(targetSequence, HandlingMode, warnings);
            result = converted is null
                ? ConversionResult.Failure(ConversionFailureReason.UnknownFormat, warnings)
                : ConversionResult.Successful(warnings.Count == 0 ? Array.Empty<string>() : warnings);
            return converted;
        }
        catch (SequenceConversionException ex)
        {
            if (HandlingMode == SequenceConversionHandlingMode.ThrowException)
            {
                throw;
            }

            result = ConversionResult.Failure(ex.FailureReason, ex.Warnings, ex.IncompatibleItems);
            return null;
        }
    }

    /// <summary>
    /// Performs the forward conversion to the target representation.
    /// </summary>
    /// <param name="fullSequence">Full mzLib sequence (with modification annotations).
    /// </param>
    /// <param name="handlingMode">Current handling mode.</param>
    /// <param name="warnings">A collection to append warning messages.</param>
    protected abstract string ConvertSequenceCore(string fullSequence, SequenceConversionHandlingMode handlingMode,
        IList<string> warnings);

    /// <summary>
    /// Optional hook for reverse conversion when <see cref="IsBidirectional"/> is true.
    /// </summary>
    protected virtual string? ConvertFromCore(string targetSequence,
        SequenceConversionHandlingMode handlingMode, IList<string> warnings) => null;
}
