using System.Globalization;
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
            sequence = EnrichModificationsIfNeeded(sequence);
            var allModsOneIsNterminus = new Dictionary<int, Modification>();

            foreach (var canonicalModification in sequence.Modifications)
            {
                var index = ResolveOneIsNterminusIndex(canonicalModification, sequence.BaseSequence.Length);
                var targetResidue = ResolveTargetResidue(canonicalModification, sequence.BaseSequence);

                if (!TryResolveKnownModification(canonicalModification, targetResidue, knownMods, out var resolved))
                {
                    if (mode == SequenceConversionHandlingMode.ThrowException)
                    {
                        throw new SequenceConversionException(
                            $"Unable to resolve modification {canonicalModification}",
                            ConversionFailureReason.IncompatibleModifications,
                            new[] { canonicalModification.ToString() });
                    }

                    warnings.AddWarning($"Unable to resolve modification {canonicalModification}; skipping.");
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

                allModsOneIsNterminus.Add(index, resolved!);
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

    private bool TryResolveKnownModification(
        CanonicalModification canonicalModification,
        char? targetResidue,
        Dictionary<string, Modification>? knownMods,
        out Modification? resolved)
    {
        resolved = null;

        if (knownMods != null)
        {
            if (canonicalModification.MzLibModification != null &&
                !string.IsNullOrWhiteSpace(canonicalModification.MzLibModification.IdWithMotif) &&
                knownMods.TryGetValue(canonicalModification.MzLibModification.IdWithMotif, out var byIdWithMotif))
            {
                resolved = byIdWithMotif;
                return true;
            }

            if (TryResolveByToken(canonicalModification.MzLibId, targetResidue, knownMods, out resolved))
            {
                return true;
            }

            if (TryResolveByToken(canonicalModification.OriginalRepresentation, targetResidue, knownMods, out resolved))
            {
                return true;
            }

            if (canonicalModification.UnimodId.HasValue &&
                TryResolveByUnimodId(canonicalModification.UnimodId.Value, targetResidue, knownMods, out resolved))
            {
                return true;
            }
        }

        if (canonicalModification.MzLibModification != null)
        {
            resolved = canonicalModification.MzLibModification;
            return true;
        }

        if (_lookup != null)
        {
            var resolvedCanonical = _lookup.TryResolve(canonicalModification);
            if (resolvedCanonical.HasValue && resolvedCanonical.Value.MzLibModification != null)
            {
                var lookupResolved = resolvedCanonical.Value.MzLibModification;

                if (knownMods != null &&
                    !string.IsNullOrWhiteSpace(lookupResolved.IdWithMotif) &&
                    knownMods.TryGetValue(lookupResolved.IdWithMotif, out var remapped))
                {
                    resolved = remapped;
                    return true;
                }

                resolved = lookupResolved;
                return true;
            }
        }

        return false;
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

    private static char? ResolveTargetResidue(CanonicalModification canonicalModification, string baseSequence)
    {
        if (canonicalModification.TargetResidue.HasValue)
        {
            return canonicalModification.TargetResidue.Value;
        }

        if (canonicalModification.PositionType == ModificationPositionType.NTerminus)
        {
            return baseSequence.Length > 0 ? baseSequence[0] : (char?)null;
        }

        if (canonicalModification.PositionType == ModificationPositionType.CTerminus)
        {
            return baseSequence.Length > 0 ? baseSequence[^1] : (char?)null;
        }

        if (canonicalModification.ResidueIndex.HasValue &&
            canonicalModification.ResidueIndex.Value >= 0 &&
            canonicalModification.ResidueIndex.Value < baseSequence.Length)
        {
            return baseSequence[canonicalModification.ResidueIndex.Value];
        }

        return null;
    }

    private static bool TryResolveByToken(
        string? modToken,
        char? targetResidue,
        Dictionary<string, Modification> knownMods,
        out Modification? resolved)
    {
        resolved = null;
        if (string.IsNullOrWhiteSpace(modToken))
        {
            return false;
        }

        var trimmed = modToken.Trim();
        if (knownMods.TryGetValue(trimmed, out var exactMatch))
        {
            resolved = exactMatch;
            return true;
        }

        var splitIndex = trimmed.IndexOf(':');
        if (splitIndex > 0 && splitIndex < trimmed.Length - 1)
        {
            var tokenWithoutPrefix = trimmed.Substring(splitIndex + 1).Trim();
            if (knownMods.TryGetValue(tokenWithoutPrefix, out var withoutPrefixMatch))
            {
                resolved = withoutPrefixMatch;
                return true;
            }
        }

        if (TryParseUnimodToken(trimmed, out var unimodId) &&
            TryResolveByUnimodId(unimodId, targetResidue, knownMods, out resolved))
        {
            return true;
        }

        return false;
    }

    private static bool TryResolveByUnimodId(
        int unimodId,
        char? targetResidue,
        Dictionary<string, Modification> knownMods,
        out Modification? resolved)
    {
        resolved = null;

        var candidates = knownMods.Values
            .Where(modification => modification != null && MatchesUnimodId(modification, unimodId))
            .Distinct()
            .ToList();

        if (candidates.Count == 0)
        {
            return false;
        }

        if (candidates.Count == 1)
        {
            resolved = candidates[0];
            return true;
        }

        if (!targetResidue.HasValue)
        {
            return false;
        }

        var residueMatches = candidates
            .Where(candidate =>
            {
                var target = candidate.Target?.ToString();
                return !string.IsNullOrWhiteSpace(target) &&
                       target.Contains(targetResidue.Value.ToString(), StringComparison.OrdinalIgnoreCase);
            })
            .Distinct()
            .ToList();

        if (residueMatches.Count != 1)
        {
            return false;
        }

        resolved = residueMatches[0];
        return true;
    }

    private static bool MatchesUnimodId(Modification modification, int unimodId)
    {
        if (modification.DatabaseReference != null)
        {
            foreach (var databaseReference in modification.DatabaseReference)
            {
                if (!databaseReference.Key.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase))
                {
                    continue;
                }

                foreach (var token in databaseReference.Value)
                {
                    if (TryParseUnimodToken(token, out var parsedId) && parsedId == unimodId)
                    {
                        return true;
                    }
                }
            }
        }

        return TryParseUnimodToken(modification.Accession, out var accessionId) && accessionId == unimodId;
    }

    private static bool TryParseUnimodToken(string? token, out int unimodId)
    {
        unimodId = 0;
        if (string.IsNullOrWhiteSpace(token))
        {
            return false;
        }

        var trimmed = token.Trim();
        var lastColon = trimmed.LastIndexOf(':');
        if (lastColon >= 0 && lastColon < trimmed.Length - 1)
        {
            var prefix = trimmed.Substring(0, lastColon).Trim();
            if (!string.IsNullOrWhiteSpace(prefix) &&
                !prefix.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase) &&
                !prefix.EndsWith("UNIMOD", StringComparison.OrdinalIgnoreCase))
            {
                return false;
            }

            trimmed = trimmed.Substring(lastColon + 1).Trim();
        }

        return int.TryParse(trimmed, NumberStyles.Integer, CultureInfo.InvariantCulture, out unimodId);
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
