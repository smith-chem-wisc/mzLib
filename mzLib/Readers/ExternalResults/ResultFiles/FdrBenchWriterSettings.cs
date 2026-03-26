using Omics.SequenceConversion;

namespace Readers;

public class FdrBenchWriterSettings
{
    public static FdrBenchWriterSettings Default { get; } = new();

    /// <summary>
    /// Determines if decoy identifications should be included in the export.
    /// Defaults to false because FDRBench expects only target hits for FDP evaluation.
    /// </summary>
    public bool IncludeDecoys { get; init; }

    /// <summary>
    /// Specifies how sequence conversion issues should be handled when serializing
    /// mzLib formatted sequences into the UniMod syntax expected by FDRBench.
    /// </summary>
    public SequenceConversionHandlingMode ConversionHandlingMode { get; init; } = SequenceConversionHandlingMode.RemoveIncompatibleElements;

    /// <summary>
    /// Optional container for warnings emitted during sequence conversion.
    /// When not provided, conversions collect warnings in a throwaway instance.
    /// </summary>
    public ConversionWarnings? ConversionWarnings { get; init; }

    /// <summary>
    /// Source sequence format used when parsing the mzLib full sequence string.
    /// Defaults to the native mzLib schema.
    /// </summary>
    public string SourceFormat { get; init; } = MzLibSequenceFormatSchema.Instance.FormatName;
}
