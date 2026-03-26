namespace Omics.SequenceConversion;

public class UnimodSequenceSerializer : SequenceSerializerBase
{
    public static UnimodSequenceSerializer Instance { get; } = new();

    public UnimodSequenceSerializer(UnimodSequenceFormatSchema? schema = null, IModificationLookup? lookup = null)
        : base(lookup ?? GlobalModificationLookup.Instance)
    {
        Schema = schema ?? UnimodSequenceFormatSchema.Instance;
    }

    public override string FormatName => Schema.FormatName;

    public override SequenceFormatSchema Schema { get; }

    public override bool CanSerialize(CanonicalSequence sequence)
    {
        return !string.IsNullOrEmpty(sequence.BaseSequence);
    }

    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        return !mod.UnimodId.HasValue;
    }

    protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        if (mod.UnimodId.HasValue)
        {
            var schema = (UnimodSequenceFormatSchema)Schema;
            var label = schema.GetLabel();
            return string.IsNullOrEmpty(label)
                ? mod.UnimodId.Value.ToString()
                : $"{label}:{mod.UnimodId.Value}";
        }

        warnings.AddIncompatibleItem(mod.ToString());

        if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
        {
            warnings.AddWarning($"Removing incompatible modification without UNIMOD ID: {mod}");
            return null;
        }

        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new SequenceConversionException(
                $"Cannot serialize modification in Unimod format - no UNIMOD ID available: {mod}",
                ConversionFailureReason.IncompatibleModifications,
                new[] { mod.ToString() });
        }

        return null;
    }
}
