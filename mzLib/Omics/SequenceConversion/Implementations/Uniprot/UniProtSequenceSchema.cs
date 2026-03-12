namespace Omics.SequenceConversion;

public class UniProtSequenceSchema : MzLibSequenceFormatSchema
{
    /// <summary>
    /// Singleton instance of the UniProt schema.
    /// Uses OriginalId (e.g., "Phosphoserine") without modification type prefix.
    /// </summary>
    public new static readonly UniProtSequenceSchema Instance = new(false, false);

    /// <summary>
    /// Public constructor for UniProtSequenceSchema.
    /// By default, does not include modification type prefixes (e.g., "Phosphoserine" instead of "UniProt:Phosphoserine").
    /// </summary>
    public UniProtSequenceSchema(bool writeModTypes, bool writeMotifs, char modOpen = '[', char modClosed = ']', string? nTermSeparator = "", string? cTermSeparator = "-")
        : base(modOpen, modClosed, nTermSeparator, cTermSeparator) 
    {
        WriteModType = writeModTypes;
        WriteMotifs = writeMotifs;
    }

    /// <inheritdoc />
    public override string FormatName => "UniProt";

    public bool WriteModType { get; set; }
    public bool WriteMotifs { get; set; }
}