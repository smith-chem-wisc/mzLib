using MathNet.Numerics;

namespace Omics.SequenceConversion;

/// <summary>
/// Serializes a <see cref="CanonicalSequence"/> into mass shift format strings.
/// 
/// Output format examples:
/// - Simple: "PEPTIDE"
/// - With residue modification: "PEPTM[+15.9949]IDE" (oxidation)
/// - With N-terminal modification: "[+42.0106]PEPTIDE" (acetylation, no separator)
/// - With C-terminal modification: "PEPTIDE-[+0.9840]" (amidation)
/// 
/// Mass shifts are formatted with:
/// - Always include sign (+ or -)
/// - 4 decimal places precision
/// - Format: [±X.XXXX]
/// </summary>
public class MassShiftSequenceSerializer : SequenceSerializerBase
{
    /// <summary>
    /// Singleton instance.
    /// </summary>
    public static MassShiftSequenceSerializer Instance { get; } = new();

    /// <summary>
    /// Creates a new MassShiftSequenceSerializer.
    /// </summary>
    /// <param name="schema">Optional schema to use for formatting.</param>
    /// <param name="lookup">Optional modification lookup to resolve modifications before serialization.</param>
    public MassShiftSequenceSerializer(MassShiftSequenceFormatSchema? schema = null, IModificationLookup? lookup = null)
        : base(lookup ?? GlobalModificationLookup.Instance)
    {
        Schema = schema ?? MassShiftSequenceFormatSchema.Instance;
    }

    /// <inheritdoc />
    public override string FormatName => MassShiftSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public override SequenceFormatSchema Schema { get; }

    /// <inheritdoc />
    public override bool CanSerialize(CanonicalSequence sequence)
    {
        return !string.IsNullOrEmpty(sequence.BaseSequence);
    }

    /// <inheritdoc />
    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        return !mod.EffectiveMass.HasValue;
    }

    /// <summary>
    /// Gets the string representation of a modification as a mass shift.
    /// </summary>
    protected override string? GetModificationString(
        CanonicalModification mod, 
        ConversionWarnings warnings, 
        SequenceConversionHandlingMode mode)
    {
        // Get the effective mass (from MonoisotopicMass, ChemicalFormula, or MzLibModification)
        var mass = mod.EffectiveMass;

        if (!mass.HasValue)
        {
            // Cannot serialize this modification without mass information
            warnings.AddIncompatibleItem(mod.ToString());
            
            if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
            {
                warnings.AddWarning($"Removing modification without mass information: {mod}");
                return null; // Signal to skip this modification
            }

            if (mode == SequenceConversionHandlingMode.ThrowException)
            {
                throw new SequenceConversionException(
                    $"Cannot serialize modification in mass shift format - no mass information available: {mod}",
                    ConversionFailureReason.IncompatibleModifications,
                    new[] { mod.ToString() });
            }

            return null;
        }

        // Format the mass with sign and configurable decimal places
        // Use explicit + sign for positive values
        var sign = mass.Value >= 0 ? "+" : "";
        return $"{sign}{mass.Value.Round((Schema as MassShiftSequenceFormatSchema)?.DecimalPlaces ?? 6)}";
    }
}
