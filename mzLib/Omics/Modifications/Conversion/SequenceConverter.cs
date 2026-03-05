using System;
using System.Collections.Generic;
using System.Linq;

namespace Omics.Modifications.Conversion;

public enum SequenceConversionHandlingMode
{
    ThrowException,
    RemoveIncompatibleMods,
    UsePrimarySequence,
    KeepOriginalAnnotation,
    ReturnNull
}

public enum SequenceConversionFailureReason
{
    None,
    NoEquivalent,
    AmbiguousEquivalent,
    IncompatibleModifications,
    UsedPrimarySequence,
    ModificationsRemoved,
    ReturnedNull,
    InvalidTargetConvention
}

public class SequenceConversionException : Exception
{
    public SequenceConversionFailureReason Reason { get; }
    public Modification? SourceModification { get; }
    public ModificationNamingConvention? TargetConvention { get; }

    public SequenceConversionException(string message, SequenceConversionFailureReason reason, Modification? sourceModification = null, ModificationNamingConvention? targetConvention = null)
        : base(message)
    {
        Reason = reason;
        SourceModification = sourceModification;
        TargetConvention = targetConvention;
    }
}

public sealed class SequenceConverter
{
    private static readonly Lazy<SequenceConverter> _default = new(() => new SequenceConverter());

    private ModificationCrossRefIndex _crossRefIndex;
    private readonly object _cacheLock = new();
    private Dictionary<ModificationNamingConvention, HashSet<Modification>> _conventionSets;
    private Dictionary<ModificationNamingConvention, Dictionary<string, Modification>> _conventionById;
    private int _cachedVersion;
    private int _indexVersion;

    public static SequenceConverter Default => _default.Value;

    public SequenceConverter(ModificationCrossRefIndex? crossRefIndex = null)
    {
        _crossRefIndex = crossRefIndex ?? ModificationCrossRefIndex.Global;
        _conventionSets = new Dictionary<ModificationNamingConvention, HashSet<Modification>>();
        _conventionById = new Dictionary<ModificationNamingConvention, Dictionary<string, Modification>>();
        _cachedVersion = -1;
        _indexVersion = Mods.RegistryVersion;
        EnsureConventionCaches();
    }

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

        var convertedModDictionary = ConvertModifications(
            bioPolymer.AllModsOneIsNterminus,
            bioPolymer.BaseSequence,
            targetConvention,
            handlingMode,
            out var reason);

        if (reason.HasValue && reason.Value != SequenceConversionFailureReason.None)
        {
            switch (handlingMode)
            {
                case SequenceConversionHandlingMode.ReturnNull:
                    failureReason = SequenceConversionFailureReason.ReturnedNull;
                    return false;
                case SequenceConversionHandlingMode.UsePrimarySequence when reason.Value == SequenceConversionFailureReason.UsedPrimarySequence:
                    converted = bioPolymer.BaseSequence;
                    failureReason = reason;
                    return true;
            }
        }

        if (convertedModDictionary == null)
        {
            converted = bioPolymer.BaseSequence;
        }
        else
        {
            converted = IBioPolymerWithSetMods.DetermineFullSequence(bioPolymer.BaseSequence, convertedModDictionary);
        }

        return true;
    }

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

    private Dictionary<int, Modification>? ConvertModifications(
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

        foreach (var kvp in source)
        {
            var key = kvp.Key;
            var originalMod = kvp.Value;

            try
            {
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
                    if (!HandleMissingEquivalent(originalMod, key, handlingMode, ref targetDictionary, baseSequence))
                    {
                        reason = handlingMode switch
                        {
                            SequenceConversionHandlingMode.ReturnNull => SequenceConversionFailureReason.ReturnedNull,
                            SequenceConversionHandlingMode.UsePrimarySequence => SequenceConversionFailureReason.UsedPrimarySequence,
                            _ => SequenceConversionFailureReason.NoEquivalent
                        };

                        if (handlingMode == SequenceConversionHandlingMode.UsePrimarySequence)
                        {
                            return null;
                        }

                        if (handlingMode == SequenceConversionHandlingMode.ReturnNull)
                        {
                            return null;
                        }
                    }

                    continue;
                }

                var filtered = FilterByChemicalFormula(originalMod, candidates);
                filtered = FilterByResidueCompatibility(filtered, baseSequence, key);

                if (filtered.Count == 0)
                {
                    if (!HandleMissingEquivalent(originalMod, key, handlingMode, ref targetDictionary, baseSequence))
                    {
                        reason = SequenceConversionFailureReason.NoEquivalent;
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
            catch (SequenceConversionException)
            {
                throw;
            }
        }

        return targetDictionary;
    }

    private bool HandleMissingEquivalent(
        Modification source,
        int key,
        SequenceConversionHandlingMode handlingMode,
        ref Dictionary<int, Modification> targetDictionary,
        string baseSequence)
    {
        switch (handlingMode)
        {
            case SequenceConversionHandlingMode.RemoveIncompatibleMods:
                return true;

            case SequenceConversionHandlingMode.UsePrimarySequence:
            case SequenceConversionHandlingMode.ReturnNull:
                return false;

            case SequenceConversionHandlingMode.KeepOriginalAnnotation:
                targetDictionary[key] = source;
                return true;

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

    private List<Modification> FilterByResidueCompatibility(IReadOnlyList<Modification> candidates, string baseSequence, int modKey)
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
                    kvp => kvp.Value.ToDictionary(mod => mod.IdWithMotif, StringComparer.Ordinal));

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
}
