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

    /// <summary>
    /// Makes <see cref="PlaceModification"/> put an N-terminal pyro-Glu (UNIMOD:27 before E, UNIMOD:28 before Q) on the first
    /// residue, when that residue carries no modification of its own. The serializer doesn't apply it; a caller that wants the
    /// placement calls PlaceModification before serializing. Off by default: mzLib keeps pyro-Glu N-terminal, where it formed.
    /// </summary>
    public bool PyroGluOnFirstResidue { get; init; }
    public override string FormatName => "Unimod";

    /// <summary>
    /// Where this format places the modification in <paramref name="sequence"/>: moved by <see cref="PyroGluOnFirstResidue"/>, or unchanged.
    /// </summary>
    public override CanonicalModification PlaceModification(CanonicalModification modification, CanonicalSequence sequence)
    {
        if (!PyroGluOnFirstResidue || modification.PositionType != ModificationPositionType.NTerminus || sequence.HasModificationAt(0))
            return modification;

        char first = sequence.BaseSequence[0];
        return (modification.UnimodId, first) is (27, 'E') or (28, 'Q')
            ? modification with { PositionType = ModificationPositionType.Residue, ResidueIndex = 0, TargetResidue = first }
            : modification;
    }

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
