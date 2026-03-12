namespace Omics.SequenceConversion;

public class UniProtSequenceSchema : MzLibSequenceFormatSchema
{
    /// <summary>
    /// Singleton instance of the UniProt schema.
    /// Uses OriginalId (e.g., "Phosphoserine") without modification type prefix.
    /// </summary>
    public static readonly UniProtSequenceSchema Instance = new(false);

    /// <summary>
    /// Private constructor to enforce singleton pattern.
    /// </summary>
    private UniProtSequenceSchema(bool writeModTypes, char modOpen = '[', char modClosed = ']', string? nTermSeparator = "", string? cTermSeparator = "-")
        : base(modOpen, modClosed, nTermSeparator, cTermSeparator) 
    {
        WriteModType = writeModTypes;
    }

    /// <inheritdoc />
    public override string FormatName => "UniProt";

    public bool WriteModType { get; set; }
}