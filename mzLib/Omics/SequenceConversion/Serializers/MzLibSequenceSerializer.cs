using System.Text;

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
public class MzLibSequenceSerializer : ISequenceSerializer
{
    private readonly IModificationLookup? _lookup;

    /// <summary>
    /// Singleton instance.
    /// </summary>
    public static MzLibSequenceSerializer Instance { get; } = new();

    /// <summary>
    /// Creates a new MzLibSequenceSerializer.
    /// </summary>
    /// <param name="lookup">Optional modification lookup to resolve modifications.</param>
    public MzLibSequenceSerializer(IModificationLookup? lookup = null)
    {
        _lookup = lookup;
    }

    /// <inheritdoc />
    public string FormatName => MzLibSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public bool CanSerialize(CanonicalSequence sequence)
    {
        // All modifications must have mzLib IDs or be resolvable
        return sequence.Modifications.All(m => !string.IsNullOrEmpty(m.MzLibId) || m.IsResolved);
    }

    /// <inheritdoc />
    public string? Serialize(
        CanonicalSequence sequence,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        if (string.IsNullOrEmpty(sequence.BaseSequence))
        {
            return HandleError(warnings, mode, ConversionFailureReason.InvalidSequence,
                "Sequence has no base sequence.");
        }

        try
        {
            return SerializeInternal(sequence, warnings, mode);
        }
        catch (SequenceConversionException)
        {
            throw;
        }
        catch (Exception ex)
        {
            return HandleError(warnings, mode, ConversionFailureReason.UnknownFormat,
                $"Unexpected error serializing sequence: {ex.Message}");
        }
    }

    private string? SerializeInternal(CanonicalSequence sequence, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        var sb = new StringBuilder();

        // Handle N-terminal modification (no separator in mzLib format)
        var nTermMod = sequence.NTerminalModification;
        if (nTermMod.HasValue)
        {
            var modString = GetModificationString(nTermMod.Value, warnings, mode);
            if (modString == null && mode == SequenceConversionHandlingMode.ReturnNull)
                return null;
            
            if (modString != null)
            {
                sb.Append('[');
                sb.Append(modString);
                sb.Append(']');
            }
        }

        // Handle residue modifications - build a lookup for quick access
        var residueMods = sequence.ResidueModifications
            .Where(m => m.ResidueIndex.HasValue)
            .ToDictionary(m => m.ResidueIndex!.Value, m => m);

        // Write sequence with modifications
        for (int i = 0; i < sequence.BaseSequence.Length; i++)
        {
            sb.Append(sequence.BaseSequence[i]);

            if (residueMods.TryGetValue(i, out var mod))
            {
                var modString = GetModificationString(mod, warnings, mode);
                if (modString == null && mode == SequenceConversionHandlingMode.ReturnNull)
                    return null;

                if (modString != null)
                {
                    sb.Append('[');
                    sb.Append(modString);
                    sb.Append(']');
                }
            }
        }

        // Handle C-terminal modification
        var cTermMod = sequence.CTerminalModification;
        if (cTermMod.HasValue)
        {
            var modString = GetModificationString(cTermMod.Value, warnings, mode);
            if (modString == null && mode == SequenceConversionHandlingMode.ReturnNull)
                return null;

            if (modString != null)
            {
                sb.Append("-[");
                sb.Append(modString);
                sb.Append(']');
            }
        }

        return sb.ToString();
    }

    /// <summary>
    /// Gets the string representation of a modification for serialization.
    /// </summary>
    private string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        // Try to use mzLib ID first
        if (!string.IsNullOrEmpty(mod.MzLibId))
        {
            return mod.MzLibId;
        }

        // Try to get IdWithMotif from resolved modification
        if (mod.IsResolved && mod.MzLibModification != null)
        {
            return mod.MzLibModification.IdWithMotif;
        }

        // Try to resolve via lookup
        if (_lookup != null)
        {
            var resolved = _lookup.TryResolve(mod);
            if (resolved.HasValue && resolved.Value.IsResolved)
            {
                return resolved.Value.MzLibModification?.IdWithMotif ?? resolved.Value.MzLibId;
            }
        }

        // Fallback to original representation if it's mzLib-style (not a mass shift)
        if (!string.IsNullOrEmpty(mod.OriginalRepresentation) && !IsNumericMassShift(mod.OriginalRepresentation))
        {
            warnings.AddWarning($"Using original representation for modification: {mod.OriginalRepresentation}");
            return mod.OriginalRepresentation;
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

    /// <summary>
    /// Checks if a string looks like a numeric mass shift (e.g., "+15.995" or "-18.011").
    /// </summary>
    private static bool IsNumericMassShift(string value)
    {
        if (string.IsNullOrEmpty(value))
            return false;

        // Check if it matches the pattern for mass shifts: optional sign, digits, optional decimal point and more digits
        return System.Text.RegularExpressions.Regex.IsMatch(value, @"^[+-]?\d+(\.\d+)?$");
    }

    private static string? HandleError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        warnings.SetFailure(reason, message);

        return mode switch
        {
            SequenceConversionHandlingMode.ThrowException =>
                throw warnings.ToException(message),
            SequenceConversionHandlingMode.ReturnNull => null,
            SequenceConversionHandlingMode.UsePrimarySequence => null, // Caller should handle
            _ => null
        };
    }
}
