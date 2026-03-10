namespace Omics.SequenceConversion;

/// <summary>
/// Serializes a <see cref="CanonicalSequence"/> into mzLib format strings.
/// 
/// Output format examples:
/// - Simple: "PEPTIDE"
/// - With residue modification: "PEP[Oxidation on M]TIDE"
/// - With N-terminal modification: "[Acetyl]PEPTIDE" (no separator)
/// - With C-terminal modification: "PEPTIDE-[Amidated]"
/// 
/// Note: For mass shift notation output, use a separate MassShiftSequenceSerializer instead.
/// </summary>
public class MzLibSequenceSerializer : SequenceSerializerBase
{
    /// <summary>
    /// Singleton instance.
    /// </summary>
    public static MzLibSequenceSerializer Instance { get; } = new();

    /// <summary>
    /// Creates a new MzLibSequenceSerializer.
    /// </summary>
    /// <param name="lookup">Optional modification lookup to resolve modifications.</param>
    public MzLibSequenceSerializer(IModificationLookup? lookup = null)
        : base(lookup ?? GlobalModificationLookup.Instance)
    {
    }

    /// <inheritdoc />
    public override string FormatName => MzLibSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public override SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public override bool CanSerialize(CanonicalSequence sequence)
    {
        return !string.IsNullOrEmpty(sequence.BaseSequence);
    }

    /// <inheritdoc />
    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        return !TryGetStrictMzLibToken(mod, out _);
    }

    /// <summary>
    /// Gets the string representation of a modification for serialization.
    /// </summary>
    protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        if (TryGetStrictMzLibToken(mod, out var token))
        {
            return token;
        }

        // Cannot serialize this modification in mzLib format
        warnings.AddIncompatibleItem(mod.ToString());
        
        if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
        {
            warnings.AddWarning($"Removing incompatible modification: {mod}");
            return null; // Signal to skip this modification
        }

        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new SequenceConversionException(
                $"Cannot serialize modification in mzLib format: {mod}. Consider using MassShiftSequenceSerializer instead.",
                ConversionFailureReason.IncompatibleModifications,
                new[] { mod.ToString() });
        }

        return null;
    }

    private static bool TryGetStrictMzLibToken(CanonicalModification mod, out string token)
    {
        if (mod.MzLibModification != null)
        {
            var modificationType = mod.MzLibModification.ModificationType;
            var idWithMotif = mod.MzLibModification.IdWithMotif;

            if (!string.IsNullOrWhiteSpace(modificationType) && !string.IsNullOrWhiteSpace(idWithMotif))
            {
                token = $"{modificationType}:{idWithMotif}";
                return true;
            }
        }

        if (TryNormalizeStrictToken(mod.MzLibId, out token))
        {
            return true;
        }

        if (TryNormalizeStrictToken(mod.OriginalRepresentation, out token))
        {
            return true;
        }

        token = string.Empty;
        return false;
    }

    private static bool TryNormalizeStrictToken(string? value, out string token)
    {
        token = string.Empty;

        if (string.IsNullOrWhiteSpace(value))
        {
            return false;
        }

        var trimmed = value.Trim();
        var separatorIndex = trimmed.IndexOf(':');
        if (separatorIndex <= 0 || separatorIndex >= trimmed.Length - 1)
        {
            return false;
        }

        var modificationType = trimmed.Substring(0, separatorIndex).Trim();
        var idWithMotif = trimmed.Substring(separatorIndex + 1).Trim();
        if (string.IsNullOrEmpty(modificationType) || string.IsNullOrEmpty(idWithMotif))
        {
            return false;
        }

        token = $"{modificationType}:{idWithMotif}";
        return true;
    }
}
