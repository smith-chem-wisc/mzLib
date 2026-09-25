namespace Omics.SequenceConversion;

/// <summary>
/// Syntax metadata for MODOMICS one-letter RNA sequences.
/// </summary>
public sealed class ModomicsSequenceFormatSchema : SequenceFormatSchema
{
    public static ModomicsSequenceFormatSchema Instance { get; } = new();

    private ModomicsSequenceFormatSchema()
        : base(modOpen: '\0', modClosed: '\0', nTermSeparator: null, cTermSeparator: null)
    {
    }

    public override string FormatName => "Modomics";
}
