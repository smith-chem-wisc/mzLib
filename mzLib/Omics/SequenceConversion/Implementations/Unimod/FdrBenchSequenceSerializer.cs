namespace Omics.SequenceConversion;

/// <summary>
/// Sequence serializer for FDRBench format.
/// FDRBench uses UniMod modification notation with parentheses instead of brackets:
/// e.g., "PEPC(UniMod:4)IDE" instead of "PEPC[UNIMOD:4]IDE"
/// </summary>
public class FdrBenchSequenceSerializer : UnimodSequenceSerializer
{
    /// <summary>
    /// Default instance using FDRBench format: CamelCase UniMod labels with parentheses.
    /// Output format: "PEPC(UniMod:4)IDE"
    /// </summary>
    public static FdrBenchSequenceSerializer Instance { get; } = new();

    /// <summary>
    /// Creates a new FDRBench serializer with the specified label style.
    /// </summary>
    /// <param name="labelStyle">UniMod label style (default: CamelCase)</param>
    public FdrBenchSequenceSerializer(UnimodLabelStyle labelStyle = UnimodLabelStyle.CamelCase)
        : base(new UnimodSequenceFormatSchema(labelStyle, '(', ')'))
    {
    }

    public override string FormatName => "FdrBench";
}
