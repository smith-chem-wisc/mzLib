using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Chemistry;

namespace Omics.Modifications.Conversion;

/// <summary>
/// Defines how to handle modifications that cannot be converted to the target convention.
/// </summary>
public enum SequenceConversionHandlingMode
{
    ThrowException,
    RemoveIncompatibleMods,
    UsePrimarySequence,
    UseMassShifts,
    KeepOriginalAnnotation,
    ReturnNull
}

/// <summary>
/// Describes why a conversion succeeded with fallback or failed.
/// </summary>
public enum SequenceConversionFailureReason
{
    None,
    NoEquivalent,
    AmbiguousEquivalent,
    IncompatibleModifications,
    UsedPrimarySequence,
    UsedMassShifts,
    ModificationsRemoved,
    ReturnedNull,
    InvalidTargetConvention
}

/// <summary>
/// Exception thrown when a sequence conversion cannot be completed.
/// </summary>
public class SequenceConversionException : Exception
{
    /// <summary>Reason for the conversion failure.</summary>
    public SequenceConversionFailureReason Reason { get; }
    /// <summary>The modification that could not be converted, if applicable.</summary>
    public Modification? SourceModification { get; }
    /// <summary>Target naming convention requested by the caller.</summary>
    public ModificationNamingConvention? TargetConvention { get; }

    /// <summary>
    /// Initializes a new instance of the <see cref="SequenceConversionException"/> class.
    /// </summary>
    public SequenceConversionException(string message, SequenceConversionFailureReason reason, Modification? sourceModification = null, ModificationNamingConvention? targetConvention = null)
        : base(message)
    {
        Reason = reason;
        SourceModification = sourceModification;
        TargetConvention = targetConvention;
    }
}

/// <summary>
/// Converts sequences and modification dictionaries between naming conventions.
/// </summary>
public sealed class SequenceConverter : ISequenceConverter
{
    /// <summary>Default number of decimal places when formatting mass shifts.</summary>
    public const int DefaultMassShiftDecimalPlaces = 4;
    /// <summary>Default sign behavior for formatted mass shifts.</summary>
    public const bool DefaultMassShiftSignedNotation = true;

    private static readonly Lazy<SequenceConverter> _default = new(() => new SequenceConverter());
    private static readonly Lazy<IReadOnlyList<Modification>> _defaultMassShiftCandidates = new(() =>
    {
        Mods.ModsByConvention.TryGetValue(ModificationNamingConvention.MetaMorpheus, out var proteinMods);
        Mods.ModsByConvention.TryGetValue(ModificationNamingConvention.MetaMorpheus_Rna, out var rnaMods);
        proteinMods ??= new List<Modification>();
        rnaMods ??= new List<Modification>();
        return proteinMods.Concat(rnaMods).ToList();
    });

    private ModificationCrossRefIndex _crossRefIndex;
    private readonly object _cacheLock = new();
    private Dictionary<ModificationNamingConvention, HashSet<Modification>> _conventionSets;
    private Dictionary<ModificationNamingConvention, Dictionary<string, Modification>> _conventionById;
    private int _cachedVersion;
    private int _indexVersion;

    /// <summary>Gets the shared default instance.</summary>
    public static SequenceConverter Default => _default.Value;

    /// <summary>
    /// Initializes a new instance of <see cref="SequenceConverter"/>.
    /// </summary>
    public SequenceConverter(ModificationCrossRefIndex? crossRefIndex = null)
    {
        _crossRefIndex = crossRefIndex ?? ModificationCrossRefIndex.Global;
        _conventionSets = new Dictionary<ModificationNamingConvention, HashSet<Modification>>();
        _conventionById = new Dictionary<ModificationNamingConvention, Dictionary<string, Modification>>();
        _cachedVersion = -1;
        _indexVersion = Mods.RegistryVersion;
        EnsureConventionCaches();
    }

    /// <summary>
    /// Formats a modification mass shift in bracket notation.
    /// </summary>
    public static string GetMassShiftNotation(
        Modification modification,
        int decimalPlaces = DefaultMassShiftDecimalPlaces,
        bool signed = DefaultMassShiftSignedNotation)
    {
        ArgumentNullException.ThrowIfNull(modification);
        if (decimalPlaces < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(decimalPlaces));
        }

        if (!modification.MonoisotopicMass.HasValue)
        {
            return "[]";
        }

        var rounded = modification.MonoisotopicMass.Value.RoundedDouble(decimalPlaces) ?? modification.MonoisotopicMass.Value;
        var formatted = rounded.ToString($"F{decimalPlaces}", CultureInfo.InvariantCulture);

        if (!signed || rounded < 0)
        {
            return $"[{formatted}]";
        }

        return $"[+{formatted}]";
    }

    /// <summary>
    /// Builds a mass-shift sequence from a base sequence and modification dictionary.
    /// </summary>
    public static string BuildMassShiftSequence(
        string baseSequence,
        IReadOnlyDictionary<int, Modification> modifications,
        int decimalPlaces = DefaultMassShiftDecimalPlaces,
        bool signed = DefaultMassShiftSignedNotation,
        Func<int, bool>? includeModificationPredicate = null)
    {
        ArgumentNullException.ThrowIfNull(baseSequence);
        ArgumentNullException.ThrowIfNull(modifications);

        if (modifications.Count == 0)
        {
            return baseSequence;
        }

        var builder = new StringBuilder(baseSequence.Length + modifications.Count * (decimalPlaces + 6));

        if (modifications.TryGetValue(1, out var nTerm) && ShouldInclude(1))
        {
            builder.Append(GetMassShiftNotation(nTerm, decimalPlaces, signed));
        }

        for (var i = 0; i < baseSequence.Length; i++)
        {
            builder.Append(baseSequence[i]);
            var key = i + 2;
            if (modifications.TryGetValue(key, out var residueMod) && ShouldInclude(key))
            {
                builder.Append(GetMassShiftNotation(residueMod, decimalPlaces, signed));
            }
        }

        var cTermKey = baseSequence.Length + 2;
        if (modifications.TryGetValue(cTermKey, out var cTerm) && ShouldInclude(cTermKey))
        {
            builder.Append('-');
            builder.Append(GetMassShiftNotation(cTerm, decimalPlaces, signed));
        }

        return builder.ToString();

        bool ShouldInclude(int index) => includeModificationPredicate?.Invoke(index) ?? true;
    }

    /// <summary>
    /// Converts a full annotated sequence to mass-shift notation.
    /// </summary>
    public static string ToMassShiftNotation(
        string fullSequence,
        Dictionary<string, Modification>? modDictionary = null,
        int decimalPlaces = DefaultMassShiftDecimalPlaces,
        bool signed = DefaultMassShiftSignedNotation,
        Func<int, bool>? includeModificationPredicate = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(fullSequence);
        var lookup = modDictionary ?? Mods.AllModsKnownDictionary;
        var baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);
        var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(fullSequence, lookup);
        return BuildMassShiftSequence(baseSequence, mods, decimalPlaces, signed, includeModificationPredicate);
    }

    /// <summary>
    /// Maps a mass-shift sequence back to modifications using candidate mods.
    /// </summary>
    public static Dictionary<int, Modification> FromMassShiftNotation(
        string massShiftSequence,
        string? baseSequence = null,
        IReadOnlyList<Modification>? candidateMods = null,
        double tolerance = 0.001)
    {
        ArgumentException.ThrowIfNullOrEmpty(massShiftSequence);
        if (tolerance < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(tolerance));
        }

        var parsed = ParseMassShiftSequence(massShiftSequence, out var parsedBaseSequence);
        if (baseSequence != null && !string.Equals(parsedBaseSequence, baseSequence, StringComparison.Ordinal))
        {
            throw new ArgumentException("Provided base sequence does not match mass shift notation.", nameof(baseSequence));
        }

        var resolvedBaseSequence = baseSequence ?? parsedBaseSequence;
        var candidates = candidateMods ?? _defaultMassShiftCandidates.Value;
        return MapMassShiftsToModifications(parsed, resolvedBaseSequence, candidates, tolerance);
    }

    /// <summary>
    /// Attempts to convert a bio-polymer full sequence to the target convention.
    /// </summary>
    public bool TryConvertFullSequence(
        IBioPolymerWithSetMods bioPolymer,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out string? converted,
        out SequenceConversionFailureReason? failureReason)
    {
        ArgumentNullException.ThrowIfNull(bioPolymer);
        failureReason = null;
        converted = null;

        if (targetConvention == ModificationNamingConvention.Mixed)
        {
            failureReason = SequenceConversionFailureReason.InvalidTargetConvention;
            return false;
        }

        EnsureConventionCaches();
        EnsureCrossRefIndex();

        var conversionSucceeded = TryConvertModifications(
            bioPolymer.AllModsOneIsNterminus,
            bioPolymer.BaseSequence,
            targetConvention,
            handlingMode,
            out var convertedModDictionary,
            out var reason);

        if (!conversionSucceeded || convertedModDictionary == null)
        {
            return HandleBioConversionFailure(bioPolymer, handlingMode, reason, out converted, out failureReason);
        }

        converted = IBioPolymerWithSetMods.DetermineFullSequence(bioPolymer.BaseSequence, convertedModDictionary);
        failureReason = reason is SequenceConversionFailureReason.None ? null : reason;
        return true;
    }

    /// <summary>
    /// Converts a bio-polymer full sequence to the target convention.
    /// </summary>
    public string ConvertFullSequence(
        IBioPolymerWithSetMods bioPolymer,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        if (!TryConvertFullSequence(bioPolymer, targetConvention, handlingMode, out var converted, out var reason) || converted == null)
        {
            throw new SequenceConversionException(
                $"Unable to convert sequence to {targetConvention} (reason: {reason})",
                reason ?? SequenceConversionFailureReason.NoEquivalent);
        }

        return converted;
    }

    /// <summary>
    /// Attempts to convert a full sequence string between conventions.
    /// </summary>
    public bool TryConvertFullSequence(
        string fullSequence,
        ModificationNamingConvention sourceConvention,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out string? converted,
        out SequenceConversionFailureReason? failureReason)
    {
        ArgumentNullException.ThrowIfNull(fullSequence);
        EnsureConventionCaches();

        failureReason = null;
        var baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);
        var sourceLookup = GetSourceLookup(sourceConvention);
        Dictionary<int, Modification> sourceMods;
        try
        {
            sourceMods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(fullSequence, sourceLookup);
        }
        catch
        {
            converted = null;
            failureReason = SequenceConversionFailureReason.NoEquivalent;
            return HandleStringConversionFailure(fullSequence, baseSequence, null, handlingMode, ref failureReason, out converted);
        }

        if (!TryConvertModifications(sourceMods, baseSequence, targetConvention, handlingMode, out var convertedMods, out var conversionReason) || convertedMods == null)
        {
            failureReason = conversionReason;
            return HandleStringConversionFailure(fullSequence, baseSequence, sourceMods, handlingMode, ref failureReason, out converted);
        }

        converted = IBioPolymerWithSetMods.DetermineFullSequence(baseSequence, convertedMods);
        failureReason = conversionReason is SequenceConversionFailureReason.None ? null : conversionReason;
        return true;
    }

    /// <summary>
    /// Converts a full sequence string between conventions.
    /// </summary>
    public string ConvertFullSequence(
        string fullSequence,
        ModificationNamingConvention sourceConvention,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        if (!TryConvertFullSequence(fullSequence, sourceConvention, targetConvention, handlingMode, out var converted, out var reason) || converted == null)
        {
            throw new SequenceConversionException(
                $"Unable to convert sequence to {targetConvention} (reason: {reason})",
                reason ?? SequenceConversionFailureReason.NoEquivalent);
        }

        return converted;
    }

    /// <summary>
    /// Attempts to build a Chronologer-compatible sequence.
    /// </summary>
    public bool TryGetChronologerSequence(
        IBioPolymerWithSetMods bioPolymer,
        SequenceConversionHandlingMode handlingMode,
        out string? massShiftSequence,
        out SequenceConversionFailureReason? reason)
    {
        ArgumentNullException.ThrowIfNull(bioPolymer);
        reason = null;
        massShiftSequence = null;

        var incompatiblePositions = bioPolymer.AllModsOneIsNterminus
            .Where(kvp => !kvp.Value.MonoisotopicMass.HasValue)
            .Select(kvp => kvp.Key)
            .ToHashSet();

        if (incompatiblePositions.Count == 0)
        {
            massShiftSequence = bioPolymer.FullSequenceWithMassShift();
            return true;
        }

        switch (handlingMode)
        {
            case SequenceConversionHandlingMode.RemoveIncompatibleMods:
                massShiftSequence = bioPolymer.FullSequenceWithMassShift(index => !incompatiblePositions.Contains(index));
                reason = SequenceConversionFailureReason.ModificationsRemoved;
                return true;

            case SequenceConversionHandlingMode.UsePrimarySequence:
                massShiftSequence = bioPolymer.BaseSequence;
                reason = SequenceConversionFailureReason.UsedPrimarySequence;
                return true;

            case SequenceConversionHandlingMode.UseMassShifts:
                reason = SequenceConversionFailureReason.IncompatibleModifications;
                return false;

            case SequenceConversionHandlingMode.KeepOriginalAnnotation:
                massShiftSequence = bioPolymer.FullSequenceWithMassShift();
                reason = SequenceConversionFailureReason.IncompatibleModifications;
                return true;

            case SequenceConversionHandlingMode.ReturnNull:
                reason = SequenceConversionFailureReason.ReturnedNull;
                return false;

            case SequenceConversionHandlingMode.ThrowException:
            default:
                throw new SequenceConversionException(
                    "Mass shift conversion failed due to incompatible modifications.",
                    SequenceConversionFailureReason.IncompatibleModifications);
        }
    }

    /// <summary>
    /// Converts a modification dictionary to the target convention.
    /// </summary>
    public Dictionary<int, Modification> ConvertModifications(
        IReadOnlyDictionary<int, Modification> source,
        string baseSequence,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        var converted = ConvertModificationsInternal(source, baseSequence, targetConvention, handlingMode, out var reason);

        if (converted == null)
        {
            throw new SequenceConversionException(
                $"Unable to convert modifications to {targetConvention} (reason: {reason})",
                reason ?? SequenceConversionFailureReason.NoEquivalent);
        }

        return converted;
    }

    /// <summary>
    /// Attempts to convert a modification dictionary to the target convention.
    /// </summary>
    public bool TryConvertModifications(
        IReadOnlyDictionary<int, Modification> source,
        string baseSequence,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out Dictionary<int, Modification>? converted,
        out SequenceConversionFailureReason? reason)
    {
        converted = ConvertModificationsInternal(source, baseSequence, targetConvention, handlingMode, out reason);
        return converted != null;
    }

    /// <summary>
    /// Attempts to convert a single modification at a specific position.
    /// </summary>
    public bool TryConvertModification(
        Modification source,
        string baseSequence,
        int modificationKey,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out Modification? converted,
        out SequenceConversionFailureReason? reason)
    {
        ArgumentNullException.ThrowIfNull(source);
        ArgumentNullException.ThrowIfNull(baseSequence);

        var single = new Dictionary<int, Modification>(1)
        {
            { modificationKey, source }
        };

        var success = TryConvertModifications(single, baseSequence, targetConvention, handlingMode, out var convertedDict, out reason);

        if (!success || convertedDict == null)
        {
            converted = null;
            return success;
        }

        converted = convertedDict.TryGetValue(modificationKey, out var mapped)
            ? mapped
            : null;

        return true;
    }

    /// <summary>
    /// Converts a single modification at a specific position.
    /// </summary>
    public Modification ConvertModification(
        Modification source,
        string baseSequence,
        int modificationKey,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        if (!TryConvertModification(source, baseSequence, modificationKey, targetConvention, handlingMode, out var converted, out var reason) || converted == null)
        {
            throw new SequenceConversionException(
                $"Unable to convert modification {source.IdWithMotif} (reason: {reason}).",
                reason ?? SequenceConversionFailureReason.NoEquivalent,
                source,
                targetConvention);
        }

        return converted;
    }

    /// <summary>
    /// Attempts to convert a modification definition to the target convention.
    /// </summary>
    public bool TryConvertModificationDefinition(
        Modification source,
        ModificationNamingConvention targetConvention,
        out Modification? converted)
    {
        ArgumentNullException.ThrowIfNull(source);
        EnsureConventionCaches();
        EnsureCrossRefIndex();

        if (_conventionSets.TryGetValue(targetConvention, out var targetSet) && targetSet.Contains(source))
        {
            converted = source;
            return true;
        }

        if (_conventionById.TryGetValue(targetConvention, out var targetLookup) && targetLookup.TryGetValue(source.IdWithMotif, out var existing))
        {
            converted = existing;
            return true;
        }

        if (!_conventionSets.TryGetValue(targetConvention, out targetSet))
        {
            converted = null;
            return false;
        }

        var candidates = _crossRefIndex.GetCandidates(source, targetSet.Contains);

        if (candidates.Count == 0)
        {
            converted = null;
            return false;
        }

        var filtered = FilterByChemicalFormula(source, candidates);

        if (filtered.Count == 0)
        {
            converted = null;
            return false;
        }

        if (filtered.Count > 1)
        {
            throw new SequenceConversionException(
                $"Multiple equivalent modifications found for {source.IdWithMotif}.",
                SequenceConversionFailureReason.AmbiguousEquivalent,
                source,
                targetConvention);
        }

        converted = filtered[0];
        return true;
    }

    /// <summary>
    /// Converts a modification definition to the target convention.
    /// </summary>
    public Modification ConvertModificationDefinition(
        Modification source,
        ModificationNamingConvention targetConvention)
    {
        if (!TryConvertModificationDefinition(source, targetConvention, out var converted) || converted == null)
        {
            throw new SequenceConversionException(
                $"Unable to convert modification {source.IdWithMotif} to {targetConvention}.",
                SequenceConversionFailureReason.NoEquivalent,
                source,
                targetConvention);
        }

        return converted;
    }

    private Dictionary<int, Modification>? ConvertModificationsInternal(
        IReadOnlyDictionary<int, Modification> source,
        string baseSequence,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out SequenceConversionFailureReason? reason)
    {
        reason = SequenceConversionFailureReason.None;
        if (source.Count == 0)
        {
            return new Dictionary<int, Modification>();
        }

        var targetDictionary = new Dictionary<int, Modification>(source.Count);
        var targetSet = _conventionSets[targetConvention];
        var removedAny = false;
        var keptOriginal = false;

        foreach (var kvp in source)
        {
            var key = kvp.Key;
            var originalMod = kvp.Value;

            if (targetSet.Contains(originalMod))
            {
                targetDictionary[key] = originalMod;
                continue;
            }

            if (_conventionById[targetConvention].TryGetValue(originalMod.IdWithMotif, out var existing))
            {
                targetDictionary[key] = existing;
                continue;
            }

            var candidates = _crossRefIndex.GetCandidates(originalMod, targetSet.Contains);

            if (candidates.Count == 0)
            {
                if (!HandleConversionFallback(originalMod, key, handlingMode, targetDictionary, ref removedAny, ref keptOriginal, out reason))
                {
                    return null;
                }

                continue;
            }

            var filtered = FilterByChemicalFormula(originalMod, candidates);
            filtered = FilterByResidueCompatibility(filtered, baseSequence, key);

            if (filtered.Count == 0)
            {
                if (!HandleConversionFallback(originalMod, key, handlingMode, targetDictionary, ref removedAny, ref keptOriginal, out reason))
                {
                    return null;
                }

                continue;
            }

            if (filtered.Count > 1)
            {
                throw new SequenceConversionException(
                    $"Multiple equivalent modifications found for {originalMod.IdWithMotif}.",
                    SequenceConversionFailureReason.AmbiguousEquivalent,
                    originalMod,
                    targetConvention);
            }

            targetDictionary[key] = filtered[0];
        }

        if (reason == SequenceConversionFailureReason.None)
        {
            if (removedAny)
            {
                reason = SequenceConversionFailureReason.ModificationsRemoved;
            }
            else if (keptOriginal)
            {
                reason = SequenceConversionFailureReason.IncompatibleModifications;
            }
        }

        return targetDictionary;
    }

    private static bool HandleConversionFallback(
        Modification source,
        int key,
        SequenceConversionHandlingMode handlingMode,
        Dictionary<int, Modification> targetDictionary,
        ref bool removedAny,
        ref bool keptOriginal,
        out SequenceConversionFailureReason? reason)
    {
        reason = SequenceConversionFailureReason.None;

        switch (handlingMode)
        {
            case SequenceConversionHandlingMode.RemoveIncompatibleMods:
                removedAny = true;
                return true;

            case SequenceConversionHandlingMode.KeepOriginalAnnotation:
                keptOriginal = true;
                targetDictionary[key] = source;
                return true;

            case SequenceConversionHandlingMode.UsePrimarySequence:
                reason = SequenceConversionFailureReason.UsedPrimarySequence;
                return false;

            case SequenceConversionHandlingMode.ReturnNull:
                reason = SequenceConversionFailureReason.ReturnedNull;
                return false;

            case SequenceConversionHandlingMode.UseMassShifts:
                reason = SequenceConversionFailureReason.UsedMassShifts;
                return false;

            case SequenceConversionHandlingMode.ThrowException:
            default:
                throw new SequenceConversionException(
                    $"Unable to convert modification {source.IdWithMotif} at position {key}.",
                    SequenceConversionFailureReason.NoEquivalent,
                    source);
        }
    }

    private static List<Modification> FilterByChemicalFormula(Modification source, IReadOnlyList<Modification> candidates)
    {
        if (source.ChemicalFormula == null || candidates.Count <= 1)
        {
            return candidates.ToList();
        }

        var filtered = candidates
            .Where(c => c.ChemicalFormula != null && c.ChemicalFormula.Equals(source.ChemicalFormula))
            .ToList();

        return filtered.Count > 0 ? filtered : candidates.ToList();
    }

    private static List<Modification> FilterByResidueCompatibility(IReadOnlyList<Modification> candidates, string baseSequence, int modKey)
    {
        if (candidates.Count <= 1)
        {
            return candidates.ToList();
        }

        var isNTerm = modKey == 1;
        var isCTerm = modKey == baseSequence.Length + 2;
        char? residue = null;

        if (!isNTerm && !isCTerm)
        {
            var residueIndex = modKey - 2;
            if (residueIndex >= 0 && residueIndex < baseSequence.Length)
            {
                residue = baseSequence[residueIndex];
            }
        }

        var filtered = candidates.Where(mod => ResidueMatches(mod, residue, isNTerm, isCTerm)).ToList();

        return filtered.Count > 0 ? filtered : candidates.ToList();
    }

    private static bool ResidueMatches(Modification mod, char? residue, bool isNTerm, bool isCTerm)
    {
        var restriction = mod.LocationRestriction ?? string.Empty;

        if (isNTerm)
        {
            if (restriction.Contains("C-terminal", StringComparison.OrdinalIgnoreCase) ||
                restriction.Contains("3'-terminal", StringComparison.OrdinalIgnoreCase))
            {
                return false;
            }
        }
        else if (isCTerm)
        {
            if (restriction.Contains("N-terminal", StringComparison.OrdinalIgnoreCase) ||
                restriction.Contains("5'-terminal", StringComparison.OrdinalIgnoreCase))
            {
                return false;
            }
        }
        else if (restriction.Contains("terminal", StringComparison.OrdinalIgnoreCase))
        {
            return false;
        }

        if (residue == null)
        {
            return true;
        }

        if (mod.Target == null)
        {
            return true;
        }

        var motif = mod.Target.ToString();
        foreach (var ch in motif)
        {
            if (char.IsUpper(ch))
            {
                if (MotifCharMatches(ch, residue.Value))
                {
                    return true;
                }

                return false;
            }
        }

        return true;
    }

    private static bool MotifCharMatches(char motifChar, char residue)
    {
        var upperMotif = char.ToUpperInvariant(motifChar);
        var upperResidue = char.ToUpperInvariant(residue);
        return upperMotif switch
        {
            'X' => true,
            'B' => upperResidue is 'D' or 'N',
            'J' => upperResidue is 'I' or 'L',
            'Z' => upperResidue is 'E' or 'Q',
            _ => upperMotif == upperResidue
        };
    }

    private void EnsureConventionCaches()
    {
        if (_cachedVersion == Mods.RegistryVersion)
        {
            return;
        }

        lock (_cacheLock)
        {
            if (_cachedVersion == Mods.RegistryVersion)
            {
                return;
            }

            _conventionSets = Mods.ModsByConvention
                .ToDictionary(kvp => kvp.Key, kvp => kvp.Value.ToHashSet());

            _conventionById = Mods.ModsByConvention
                .ToDictionary(
                    kvp => kvp.Key,
                    kvp => kvp.Value
                        .DistinctBy(mod => mod.IdWithMotif, StringComparer.Ordinal)
                        .ToDictionary(mod => mod.IdWithMotif, StringComparer.Ordinal));

            _cachedVersion = Mods.RegistryVersion;
        }
    }

    private void EnsureCrossRefIndex()
    {
        if (_indexVersion == Mods.RegistryVersion)
        {
            return;
        }

        lock (_cacheLock)
        {
            if (_indexVersion == Mods.RegistryVersion)
            {
                return;
            }

            _crossRefIndex = ModificationCrossRefIndex.Global;
            _indexVersion = Mods.RegistryVersion;
        }
    }

    private Dictionary<string, Modification> GetSourceLookup(ModificationNamingConvention convention)
    {
        if (_conventionById.TryGetValue(convention, out var lookup))
        {
            return lookup;
        }

        return Mods.AllModsKnownDictionary;
    }

    private static bool HandleBioConversionFailure(
        IBioPolymerWithSetMods bioPolymer,
        SequenceConversionHandlingMode handlingMode,
        SequenceConversionFailureReason? conversionReason,
        out string? converted,
        out SequenceConversionFailureReason? failureReason)
    {
        failureReason = conversionReason;

        switch (handlingMode)
        {
            case SequenceConversionHandlingMode.UsePrimarySequence:
                converted = bioPolymer.BaseSequence;
                failureReason = SequenceConversionFailureReason.UsedPrimarySequence;
                return true;
            case SequenceConversionHandlingMode.UseMassShifts:
                if (TryBuildMassShiftSequence(bioPolymer, out converted))
                {
                    failureReason = SequenceConversionFailureReason.UsedMassShifts;
                    return true;
                }

                failureReason = SequenceConversionFailureReason.IncompatibleModifications;
                converted = null;
                return false;
            case SequenceConversionHandlingMode.KeepOriginalAnnotation:
                converted = bioPolymer.FullSequence;
                failureReason = SequenceConversionFailureReason.IncompatibleModifications;
                return true;
            case SequenceConversionHandlingMode.RemoveIncompatibleMods:
                converted = bioPolymer.BaseSequence;
                failureReason = SequenceConversionFailureReason.ModificationsRemoved;
                return true;
            case SequenceConversionHandlingMode.ReturnNull:
                failureReason = SequenceConversionFailureReason.ReturnedNull;
                converted = null;
                return false;
            case SequenceConversionHandlingMode.ThrowException:
            default:
                throw new SequenceConversionException(
                    "Sequence conversion failed and no fallback was allowed by the selected handling mode.",
                    conversionReason ?? SequenceConversionFailureReason.NoEquivalent);
        }
    }

    private static bool HandleStringConversionFailure(
        string originalSequence,
        string baseSequence,
        IReadOnlyDictionary<int, Modification>? sourceMods,
        SequenceConversionHandlingMode handlingMode,
        ref SequenceConversionFailureReason? failureReason,
        out string? converted)
    {
        switch (handlingMode)
        {
            case SequenceConversionHandlingMode.UsePrimarySequence:
                converted = baseSequence;
                failureReason = SequenceConversionFailureReason.UsedPrimarySequence;
                return true;
            case SequenceConversionHandlingMode.UseMassShifts:
                if (sourceMods != null && TryBuildMassShiftSequence(baseSequence, sourceMods, out converted))
                {
                    failureReason = SequenceConversionFailureReason.UsedMassShifts;
                    return true;
                }

                failureReason = SequenceConversionFailureReason.IncompatibleModifications;
                converted = null;
                return false;
            case SequenceConversionHandlingMode.KeepOriginalAnnotation:
                converted = originalSequence;
                failureReason = SequenceConversionFailureReason.IncompatibleModifications;
                return true;
            case SequenceConversionHandlingMode.RemoveIncompatibleMods:
                converted = baseSequence;
                failureReason = SequenceConversionFailureReason.ModificationsRemoved;
                return true;
            case SequenceConversionHandlingMode.ReturnNull:
                failureReason = SequenceConversionFailureReason.ReturnedNull;
                converted = null;
                return false;
            default:
                converted = null;
                return false;
        }
    }

    private static Dictionary<int, double> ParseMassShiftSequence(string massShiftSequence, out string baseSequence)
    {
        var baseBuilder = new StringBuilder(massShiftSequence.Length);
        var shifts = new Dictionary<int, double>();
        var lastResidueIndex = -1;
        var parsingCTerm = false;

        for (var i = 0; i < massShiftSequence.Length;)
        {
            var current = massShiftSequence[i];

            if (current == '-' && i + 1 < massShiftSequence.Length && massShiftSequence[i + 1] == '[')
            {
                parsingCTerm = true;
                i++;
                continue;
            }

            if (current == '[')
            {
                var closeIndex = massShiftSequence.IndexOf(']', i);
                if (closeIndex < 0)
                {
                    throw new FormatException("Unterminated mass shift token detected.");
                }

                var token = massShiftSequence.Substring(i + 1, closeIndex - i - 1);
                if (!double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out var mass))
                {
                    throw new FormatException($"Unable to parse mass shift token '{token}'.");
                }

                var key = parsingCTerm
                    ? baseBuilder.Length + 2
                    : baseBuilder.Length == 0 && lastResidueIndex == -1
                        ? 1
                        : lastResidueIndex + 2;

                shifts[key] = mass;
                parsingCTerm = false;
                i = closeIndex + 1;
                continue;
            }

            baseBuilder.Append(current);
            lastResidueIndex = baseBuilder.Length - 1;
            i++;
        }

        baseSequence = baseBuilder.ToString();
        return shifts;
    }

    private static Dictionary<int, Modification> MapMassShiftsToModifications(
        IReadOnlyDictionary<int, double> massShifts,
        string baseSequence,
        IReadOnlyList<Modification> candidateMods,
        double tolerance)
    {
        var mapped = new Dictionary<int, Modification>(massShifts.Count);

        foreach (var (key, mass) in massShifts)
        {
            var matches = candidateMods
                .Where(mod => mod.MonoisotopicMass.HasValue && Math.Abs(mod.MonoisotopicMass.Value - mass) <= tolerance)
                .ToList();

            if (matches.Count == 0)
            {
                throw new SequenceConversionException(
                    $"No modification matches mass shift {mass} at position {key}.",
                    SequenceConversionFailureReason.NoEquivalent);
            }

            var filtered = FilterByResidueCompatibility(matches, baseSequence, key);

            if (filtered.Count == 0)
            {
                throw new SequenceConversionException(
                    $"No residue-compatible modification matches mass shift {mass} at position {key}.",
                    SequenceConversionFailureReason.NoEquivalent);
            }

            if (filtered.Count > 1)
            {
                throw new SequenceConversionException(
                    $"Multiple modifications match mass shift {mass} at position {key}.",
                    SequenceConversionFailureReason.AmbiguousEquivalent);
            }

            mapped[key] = filtered[0];
        }

        return mapped;
    }

    private static bool TryBuildMassShiftSequence(IBioPolymerWithSetMods bioPolymer, out string? massShiftSequence)
    {
        if (bioPolymer.AllModsOneIsNterminus.Values.Any(mod => !mod.MonoisotopicMass.HasValue))
        {
            massShiftSequence = null;
            return false;
        }

        massShiftSequence = bioPolymer.FullSequenceWithMassShift(decimalPlaces: DefaultMassShiftDecimalPlaces);
        return true;
    }

    private static bool TryBuildMassShiftSequence(
        string baseSequence,
        IReadOnlyDictionary<int, Modification> sourceMods,
        out string? massShiftSequence)
    {
        if (sourceMods.Count == 0)
        {
            massShiftSequence = baseSequence;
            return true;
        }

        if (sourceMods.Values.Any(mod => !mod.MonoisotopicMass.HasValue))
        {
            massShiftSequence = null;
            return false;
        }

        massShiftSequence = BuildMassShiftSequence(
            baseSequence,
            sourceMods,
            DefaultMassShiftDecimalPlaces,
            DefaultMassShiftSignedNotation);
        return true;
    }
}
