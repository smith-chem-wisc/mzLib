using System.Text;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Base class for sequence serializers that use bracket notation for modifications.
/// Provides common serialization logic for formats like mzLib and mass-shift.
/// </summary>
public abstract class SequenceSerializerBase : ISequenceSerializer
{
    private readonly IModificationLookup? _lookup;

    protected SequenceSerializerBase(IModificationLookup? lookup = null)
    {
        _lookup = lookup;
    }

    /// <inheritdoc />
    public abstract string FormatName { get; }

    /// <inheritdoc />
    public abstract SequenceFormatSchema Schema { get; }

    /// <inheritdoc />
    public IModificationLookup? ModificationLookup => _lookup;

    /// <inheritdoc />
    public virtual SequenceConversionHandlingMode HandlingMode => SequenceConversionHandlingMode.ThrowException;

    /// <inheritdoc />
    public abstract bool CanSerialize(CanonicalSequence sequence);

    /// <inheritdoc />
    public abstract bool ShouldResolveMod(CanonicalModification mod);

    /// <inheritdoc />
    public virtual string? Serialize(
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
            sequence = EnrichModificationsIfNeeded(sequence);
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

    /// <inheritdoc />
    public virtual Dictionary<int, Modification> ToOneIsNterminusModificationDictionary(
        CanonicalSequence sequence,
        Dictionary<string, Modification>? knownMods = null,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        if (string.IsNullOrEmpty(sequence.BaseSequence))
        {
            HandleDictionaryError(
                warnings,
                mode,
                ConversionFailureReason.InvalidSequence,
                "Sequence has no base sequence.");
            return new Dictionary<int, Modification>();
        }

        try
        {
            sequence = ResolveModificationsForProjection(sequence, knownMods);
            var allModsOneIsNterminus = new Dictionary<int, Modification>();

            foreach (var canonicalModification in sequence.Modifications)
            {
                var index = ResolveOneIsNterminusIndex(canonicalModification, sequence.BaseSequence.Length);
                var resolved = canonicalModification.MzLibModification;
                if (resolved == null)
                {
                    if (mode == SequenceConversionHandlingMode.ThrowException)
                    {
                        throw new SequenceConversionException(
                            $"Unable to resolve projected modification {canonicalModification}",
                            ConversionFailureReason.IncompatibleModifications,
                            new[] { canonicalModification.ToString() });
                    }

                    warnings.AddWarning($"Unable to resolve projected modification {canonicalModification}; skipping.");
                    warnings.AddIncompatibleItem(canonicalModification.ToString());
                    continue;
                }

                if (allModsOneIsNterminus.ContainsKey(index))
                {
                    if (mode == SequenceConversionHandlingMode.ThrowException)
                    {
                        throw new SequenceConversionException(
                            $"Multiple modifications map to OneIsNterminus index {index}.",
                            ConversionFailureReason.InvalidSequence,
                            new[] { canonicalModification.ToString() });
                    }

                    warnings.AddWarning($"Multiple modifications map to OneIsNterminus index {index}; keeping the first entry.");
                    continue;
                }

                allModsOneIsNterminus.Add(index, resolved);
            }

            return allModsOneIsNterminus;
        }
        catch (SequenceConversionException)
        {
            throw;
        }
        catch (Exception ex)
        {
            HandleDictionaryError(
                warnings,
                mode,
                ConversionFailureReason.UnknownFormat,
                $"Unexpected error projecting sequence modifications: {ex.Message}");
            return new Dictionary<int, Modification>();
        }
    }

    private CanonicalSequence ResolveModificationsForProjection(
        CanonicalSequence sequence,
        Dictionary<string, Modification>? knownMods)
    {
        if (!sequence.HasModifications || (_lookup == null && knownMods == null))
        {
            return sequence;
        }

        var modifications = sequence.Modifications;
        var updated = new CanonicalModification[modifications.Length];
        var changed = false;
        var fallbackLookup = knownMods != null ? new KnownModsLookup(knownMods) : null;

        for (int i = 0; i < modifications.Length; i++)
        {
            var mod = modifications[i];
            var enriched = mod;

            if (mod.MzLibModification == null)
            {
                if (fallbackLookup != null)
                {
                    var resolved = fallbackLookup.TryResolve(mod);
                    if (resolved.HasValue)
                    {
                        enriched = resolved.Value with
                        {
                            PositionType = mod.PositionType,
                            ResidueIndex = mod.ResidueIndex,
                            TargetResidue = mod.TargetResidue ?? resolved.Value.TargetResidue,
                            OriginalRepresentation = mod.OriginalRepresentation
                        };
                    }
                }

                if (enriched.MzLibModification == null && _lookup != null)
                {
                    var resolved = _lookup.TryResolve(mod);
                    if (resolved.HasValue)
                    {
                        enriched = resolved.Value with
                        {
                            PositionType = mod.PositionType,
                            ResidueIndex = mod.ResidueIndex,
                            TargetResidue = mod.TargetResidue ?? resolved.Value.TargetResidue,
                            OriginalRepresentation = mod.OriginalRepresentation
                        };
                    }
                }
            }

            if (!enriched.Equals(mod))
            {
                changed = true;
            }

            updated[i] = enriched;
        }

        return changed ? sequence.WithModifications(updated) : sequence;
    }

    /// <summary>
    /// Resolves and enriches modifications only when required by <see cref="ShouldResolveMod"/>.
    /// </summary>
    protected virtual CanonicalSequence EnrichModificationsIfNeeded(CanonicalSequence sequence)
    {
        if (!sequence.HasModifications || _lookup == null)
        {
            return sequence;
        }

        var modifications = sequence.Modifications;
        var updated = new CanonicalModification[modifications.Length];
        var changed = false;

        for (int i = 0; i < modifications.Length; i++)
        {
            var mod = modifications[i];
            var enriched = mod;

            if (ShouldResolveMod(mod))
            {
                var resolved = _lookup.TryResolve(mod);
                if (resolved.HasValue)
                {
                    enriched = resolved.Value with
                    {
                        PositionType = mod.PositionType,
                        ResidueIndex = mod.ResidueIndex,
                        TargetResidue = mod.TargetResidue ?? resolved.Value.TargetResidue,
                        OriginalRepresentation = mod.OriginalRepresentation
                    };

                    if (!enriched.Equals(mod))
                    {
                        changed = true;
                    }
                }
            }

            updated[i] = enriched;
        }

        return changed ? sequence.WithModifications(updated) : sequence;
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
                if (!string.IsNullOrEmpty(Schema.NTermSeparator))
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
                if (!string.IsNullOrEmpty(Schema.CTermSeparator))
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

    private static int ResolveOneIsNterminusIndex(CanonicalModification canonicalModification, int baseSequenceLength)
    {
        switch (canonicalModification.PositionType)
        {
            case ModificationPositionType.NTerminus:
                return 1;

            case ModificationPositionType.CTerminus:
                return baseSequenceLength + 2;

            case ModificationPositionType.Residue:
                if (!canonicalModification.ResidueIndex.HasValue ||
                    canonicalModification.ResidueIndex.Value < 0 ||
                    canonicalModification.ResidueIndex.Value >= baseSequenceLength)
                {
                    throw new SequenceConversionException(
                        $"Residue index {canonicalModification.ResidueIndex?.ToString() ?? "null"} is invalid for base sequence length {baseSequenceLength}.",
                        ConversionFailureReason.InvalidSequence,
                        new[] { canonicalModification.ToString() });
                }

                return canonicalModification.ResidueIndex.Value + 2;

            default:
                throw new SequenceConversionException(
                    $"Unsupported modification position type: {canonicalModification.PositionType}",
                    ConversionFailureReason.InvalidSequence,
                    new[] { canonicalModification.ToString() });
        }
    }

    private static void HandleDictionaryError(
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        ConversionFailureReason reason,
        string message)
    {
        warnings.SetFailure(reason, message);

        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw warnings.ToException(message);
        }
    }
}
