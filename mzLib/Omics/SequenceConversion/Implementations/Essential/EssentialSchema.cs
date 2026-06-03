namespace Omics.SequenceConversion;

public sealed class EssentialSequenceFormatSchema : SequenceFormatSchema
{
    public static readonly EssentialSequenceFormatSchema Instance = new();

    private EssentialSequenceFormatSchema(char modOpen = '[', char modClosed = ']', string? nTermSeparator = "", string? cTermSeparator = "-")
        : base(modOpen, modClosed, nTermSeparator, cTermSeparator)
    {
    }

    public override string FormatName => "Essential";
}
