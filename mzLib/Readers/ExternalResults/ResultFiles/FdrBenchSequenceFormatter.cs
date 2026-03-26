using Omics.SequenceConversion;

namespace Readers;

internal static class FdrBenchSequenceFormatter
{
    private static readonly ISequenceSerializer SequenceSerializer = new UnimodSequenceSerializer(
        new UnimodSequenceFormatSchema(UnimodLabelStyle.CamelCase, '(', ')'));

    public static string FormatSequence(string? fullSequence, string fallbackBaseSequence, FdrBenchWriterSettings settings)
    {
        if (string.IsNullOrWhiteSpace(fullSequence))
        {
            return fallbackBaseSequence;
        }

        settings ??= FdrBenchWriterSettings.Default;
        var warnings = settings.ConversionWarnings ?? new ConversionWarnings();

        var canonical = SequenceConversionService.Default.Parse(
            fullSequence,
            settings.SourceFormat,
            warnings,
            settings.ConversionHandlingMode);

        if (canonical == null)
        {
            return fallbackBaseSequence;
        }

        var serialized = SequenceSerializer.Serialize(
            canonical.Value,
            warnings,
            settings.ConversionHandlingMode);

        return string.IsNullOrEmpty(serialized)
            ? fallbackBaseSequence
            : serialized;
    }
}
