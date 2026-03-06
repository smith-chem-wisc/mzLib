namespace Omics.SequenceConversion;

/// <summary>
/// Converts sequences represented as mzLib full-sequence strings to and from alternative representations.
/// </summary>
public interface ISequenceConverter
{
    string ConverterName { get; }

    bool IsBidirectional { get; }

    SequenceConversionHandlingMode HandlingMode { get; set; }

    SequenceConversionHandlingMode DefaultHandlingMode { get; }

    string? ConvertTo(string fullSequence, out ConversionResult result);

    string? ConvertFrom(string targetSequence, out ConversionResult result);
}
