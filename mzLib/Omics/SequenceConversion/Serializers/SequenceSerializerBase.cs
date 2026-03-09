using System.Text;

namespace Omics.SequenceConversion;

/// <summary>
/// Base class for sequence serializers that use bracket notation for modifications.
/// Provides common serialization logic for formats like mzLib and mass-shift.
/// </summary>
public abstract class SequenceSerializerBase : ISequenceSerializer
{
    /// <inheritdoc />
    public abstract string FormatName { get; }

    /// <inheritdoc />
    public abstract SequenceFormatSchema Schema { get; }

    /// <inheritdoc />
    public abstract bool CanSerialize(CanonicalSequence sequence);

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

    /// <summary>
    /// Internal serialization implementation shared across bracket-based serializers.
    /// </summary>
    protected virtual string? SerializeInternal(
        CanonicalSequence sequence, 
        ConversionWarnings warnings, 
        SequenceConversionHandlingMode mode)
    {
        var sb = new StringBuilder();

        // Handle N-terminal modification
        var nTermMod = sequence.NTerminalModification;
        if (nTermMod.HasValue)
        {
            var modString = GetModificationString(nTermMod.Value, warnings, mode);
            if (modString == null && mode == SequenceConversionHandlingMode.ReturnNull)
                return null;
            
            if (modString != null)
            {
                sb.Append(Schema.ModOpenBracket);
                sb.Append(modString);
                sb.Append(Schema.ModCloseBracket);
                
                // Add N-terminal separator if defined and not empty
                if (Schema.NTermSeparator != null && !string.IsNullOrEmpty(Schema.NTermSeparator))
                {
                    sb.Append(Schema.NTermSeparator);
                }
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
                    sb.Append(Schema.ModOpenBracket);
                    sb.Append(modString);
                    sb.Append(Schema.ModCloseBracket);
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
                // Add C-terminal separator if defined
                if (Schema.CTermSeparator != null && !string.IsNullOrEmpty(Schema.CTermSeparator))
                {
                    sb.Append(Schema.CTermSeparator);
                }
                
                sb.Append(Schema.ModOpenBracket);
                sb.Append(modString);
                sb.Append(Schema.ModCloseBracket);
            }
        }

        return sb.ToString();
    }

    /// <summary>
    /// Gets the string representation of a modification for serialization.
    /// Must be implemented by derived classes to handle format-specific modification syntax.
    /// </summary>
    /// <param name="mod">The modification to serialize.</param>
    /// <param name="warnings">Warnings collection to add any warnings.</param>
    /// <param name="mode">Handling mode for errors.</param>
    /// <returns>The modification string (without brackets), or null if the modification should be skipped.</returns>
    protected abstract string? GetModificationString(
        CanonicalModification mod, 
        ConversionWarnings warnings, 
        SequenceConversionHandlingMode mode);

    /// <summary>
    /// Handles an error based on the handling mode.
    /// </summary>
    protected static string? HandleError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        return SequenceConversionHelpers.HandleSerializerError(warnings, mode, reason, message);
    }
}
