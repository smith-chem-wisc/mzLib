using System;

namespace Omics.SequenceConversion;

/// <summary>
/// Serializes canonical sequences to the UniProt modification notation (e.g., [UniProt:Phosphoserine on S]).
/// </summary>
public class UniprotSequenceSerializer : SequenceSerializerBase
{
    public static UniprotSequenceSerializer Instance { get; } = new();

    public UniprotSequenceSerializer(SequenceFormatSchema? schema = null, IModificationLookup? lookup = null)
        : base(lookup ?? UniprotModificationLookup.Instance)
    {
        Schema = schema ?? MzLibSequenceFormatSchema.Instance;
    }

    public override string FormatName => Schema.FormatName;

    public override SequenceFormatSchema Schema { get; }

    public override bool CanSerialize(CanonicalSequence sequence) => !string.IsNullOrEmpty(sequence.BaseSequence);

    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        var resolved = mod.MzLibModification;
        return resolved == null || !string.Equals(resolved.ModificationType, "UniProt", StringComparison.OrdinalIgnoreCase);
    }

    protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        var resolved = mod.MzLibModification;
        if (resolved != null &&
            !string.IsNullOrWhiteSpace(resolved.IdWithMotif) &&
            string.Equals(resolved.ModificationType, "UniProt", StringComparison.OrdinalIgnoreCase))
        {
            return $"UniProt:{resolved.IdWithMotif}";
        }

        warnings.AddIncompatibleItem(mod.ToString());

        if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
        {
            warnings.AddWarning($"Removing incompatible modification without UniProt mapping: {mod}");
            return null;
        }

        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new SequenceConversionException(
                $"Cannot serialize modification in UniProt format - mapping unavailable: {mod}",
                ConversionFailureReason.IncompatibleModifications,
                new[] { mod.ToString() });
        }

        return null;
    }
}
