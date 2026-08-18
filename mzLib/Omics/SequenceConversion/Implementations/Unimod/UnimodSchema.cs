namespace Omics.SequenceConversion;

public enum UnimodLabelStyle
{
    UpperCase,
    CamelCase,
    LowerCase,
    NoLabel
}

public sealed class UnimodSequenceFormatSchema : SequenceFormatSchema
{
    public static readonly UnimodSequenceFormatSchema Instance = new();

    public UnimodSequenceFormatSchema(UnimodLabelStyle labelStyle = UnimodLabelStyle.UpperCase, char modOpen = '[', char modClosed = ']', string nTermSeparator = "", string cTermSeparator = "-")
        : base(modOpen, modClosed, nTermSeparator, cTermSeparator)
    {
        LabelStyle = labelStyle;
    }

    public UnimodLabelStyle LabelStyle { get; init; }
    public override string FormatName => "Unimod";

    public string GetLabel()
    {
        return LabelStyle switch
        {
            UnimodLabelStyle.UpperCase => "UNIMOD",
            UnimodLabelStyle.CamelCase => "Unimod",
            UnimodLabelStyle.LowerCase => "unimod",
            UnimodLabelStyle.NoLabel => string.Empty,
            _ => "UNIMOD"
        };
    }
}
