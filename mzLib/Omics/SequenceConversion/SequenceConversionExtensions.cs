using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Omics.SequenceConversion;

public static class SequenceConversionExtensions
{
    public static CanonicalSequence ToCanonicalSequence(
        this IBioPolymerWithSetMods bioPolymer,
        string? sourceFormat = null)
    {
        var format = string.IsNullOrWhiteSpace(sourceFormat)
            ? MzLibSequenceFormatSchema.Instance.FormatName
            : sourceFormat;

        var builder = new CanonicalSequenceBuilder(bioPolymer.BaseSequence)
            .WithSourceFormat(format);

        var baseSequence = bioPolymer.BaseSequence;
        foreach (var kvp in bioPolymer.AllModsOneIsNterminus)
        {
            var index = kvp.Key;
            var mod = kvp.Value;
            var representation = BuildOriginalRepresentation(mod);

            if (index == 1)
            {
                builder.AddNTerminalModification(
                    representation,
                    mod.MonoisotopicMass,
                    mod.ChemicalFormula,
                    mzLibId: mod.IdWithMotif,
                    mzLibModification: mod);
            }
            else if (index == baseSequence.Length + 2)
            {
                builder.AddCTerminalModification(
                    representation,
                    mod.MonoisotopicMass,
                    mod.ChemicalFormula,
                    mzLibId: mod.IdWithMotif,
                    mzLibModification: mod);
            }
            else
            {
                var residueIndex = index - 2;
                if (residueIndex < 0 || residueIndex >= baseSequence.Length)
                {
                    continue;
                }

                builder.AddResidueModification(
                    residueIndex,
                    representation,
                    mod.MonoisotopicMass,
                    mod.ChemicalFormula,
                    mzLibId: mod.IdWithMotif,
                    mzLibModification: mod);
            }
        }

        return builder.Build();
    }

    public static CanonicalSequence ToCanonicalSequence(
        this IBioPolymer bioPolymer,
        ConversionWarnings? warnings = null,
        string? sourceFormat = null)
    {
        warnings ??= new ConversionWarnings();

        var format = string.IsNullOrWhiteSpace(sourceFormat)
            ? MzLibSequenceFormatSchema.Instance.FormatName
            : sourceFormat;

        var builder = new CanonicalSequenceBuilder(bioPolymer.BaseSequence)
            .WithSourceFormat(format);

        var baseSequence = bioPolymer.BaseSequence;
        var mods = SelectModDictionary(bioPolymer, warnings);

        foreach (var kvp in mods)
        {
            var index = kvp.Key;
            var mod = SelectFirstModification(kvp.Value, warnings, index);
            if (mod == null)
            {
                continue;
            }

            var representation = BuildOriginalRepresentation(mod);
            var isNTerm = IsNTerminal(mod.LocationRestriction, index);
            var isCTerm = IsCTerminal(mod.LocationRestriction, index, baseSequence.Length);

            if (isNTerm)
            {
                builder.AddNTerminalModification(
                    representation,
                    mod.MonoisotopicMass,
                    mod.ChemicalFormula,
                    mzLibId: mod.IdWithMotif,
                    mzLibModification: mod);
                continue;
            }

            if (isCTerm)
            {
                builder.AddCTerminalModification(
                    representation,
                    mod.MonoisotopicMass,
                    mod.ChemicalFormula,
                    mzLibId: mod.IdWithMotif,
                    mzLibModification: mod);
                continue;
            }

            var residueIndex = index - 1;
            if (residueIndex < 0 || residueIndex >= baseSequence.Length)
            {
                continue;
            }

            builder.AddResidueModification(
                residueIndex,
                representation,
                mod.MonoisotopicMass,
                mod.ChemicalFormula,
                mzLibId: mod.IdWithMotif,
                mzLibModification: mod);
        }

        return builder.Build();
    }

    public static string? Serialize(
        this IBioPolymerWithSetMods bioPolymer,
        string targetFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException,
        SequenceConversionService? service = null)
    {
        var converter = bioPolymer.ToCanonicalSequence();
        return (service ?? SequenceConversionService.Default).Serialize(converter, targetFormat, warnings, mode);
    }

    public static string? Serialize(
        this IBioPolymer bioPolymer,
        string targetFormat,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException,
        SequenceConversionService? service = null)
    {
        var canonical = bioPolymer.ToCanonicalSequence(warnings);
        return (service ?? SequenceConversionService.Default).Serialize(canonical, targetFormat, warnings, mode);
    }

    public static string? Serialize(
        this IBioPolymerWithSetMods bioPolymer,
        ISequenceSerializer serializer,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        var canonical = bioPolymer.ToCanonicalSequence();
        return serializer.Serialize(canonical, warnings, mode);
    }

    public static string? Serialize(
        this IBioPolymer bioPolymer,
        ISequenceSerializer serializer,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        var canonical = bioPolymer.ToCanonicalSequence(warnings);
        return serializer.Serialize(canonical, warnings, mode);
    }

    public static Dictionary<int, Modification> ParseToOneIsNterminusModificationDictionary(
        this SequenceConversionService service,
        string sequence,
        Dictionary<string, Modification> knownMods,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        ArgumentNullException.ThrowIfNull(service);
        ArgumentNullException.ThrowIfNull(knownMods);

        warnings ??= new ConversionWarnings();
        var canonical = service.ParseAutoDetect(sequence, warnings, mode);
        if (!canonical.HasValue)
        {
            return new Dictionary<int, Modification>();
        }

        var serializer = ResolveSerializerForCanonical(service, canonical.Value);
        return serializer.ToOneIsNterminusModificationDictionary(canonical.Value, knownMods, warnings, mode);
    }

    public static Dictionary<int, Modification> ToOneIsNterminusModificationDictionary(
        this CanonicalSequence canonical,
        Dictionary<string, Modification> knownMods,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException,
        ConversionWarnings? warnings = null)
    {
        ArgumentNullException.ThrowIfNull(knownMods);
        var serializer = ResolveSerializerForCanonical(SequenceConversionService.Default, canonical);
        return serializer.ToOneIsNterminusModificationDictionary(canonical, knownMods, warnings, mode);
    }

    private static ISequenceSerializer ResolveSerializerForCanonical(SequenceConversionService service, CanonicalSequence canonical)
    {
        var sourceFormat = canonical.SourceFormat;
        if (!string.IsNullOrWhiteSpace(sourceFormat))
        {
            var serializer = service.GetSerializer(sourceFormat);
            if (serializer != null)
            {
                return serializer;
            }
        }

        return MzLibSequenceSerializer.Instance;
    }

    private static IDictionary<int, List<Modification>> SelectModDictionary(
        IBioPolymer bioPolymer,
        ConversionWarnings warnings)
    {
        var originals = bioPolymer.OriginalNonVariantModifications;
        if (originals != null && originals.Count > 0)
        {
            return originals;
        }

        var localized = bioPolymer.OneBasedPossibleLocalizedModifications;
        if (localized != null && localized.Count > 0)
        {
            warnings.AddWarning("Falling back to possible localized modifications for canonical conversion.");
            return localized;
        }

        return new Dictionary<int, List<Modification>>();
    }

    private static Modification? SelectFirstModification(
        IList<Modification> modifications,
        ConversionWarnings warnings,
        int index)
    {
        if (modifications == null || modifications.Count == 0)
        {
            return null;
        }

        if (modifications.Count > 1)
        {
            warnings.AddWarning($"Multiple modifications at position {index}; using the first entry for conversion.");
        }

        return modifications[0];
    }

    private static bool IsNTerminal(string? locationRestriction, int index)
    {
        if (string.IsNullOrWhiteSpace(locationRestriction))
        {
            return false;
        }

        return index == 1 && locationRestriction.Contains("N-terminal", StringComparison.OrdinalIgnoreCase) ||
               index == 1 && locationRestriction.Contains("5'-terminal", StringComparison.OrdinalIgnoreCase);
    }

    private static bool IsCTerminal(string? locationRestriction, int index, int sequenceLength)
    {
        if (string.IsNullOrWhiteSpace(locationRestriction))
        {
            return false;
        }

        return index == sequenceLength && locationRestriction.Contains("C-terminal", StringComparison.OrdinalIgnoreCase) ||
               index == sequenceLength && locationRestriction.Contains("3'-terminal", StringComparison.OrdinalIgnoreCase);
    }

    #region Peptides/Oligos
    public static void ConvertModifications(this IBioPolymerWithSetMods withSetMods, IModificationLookup modificationLookup)
    {
        ConvertModifications(withSetMods, modificationLookup, SequenceConversionHandlingMode.RemoveIncompatibleElements);
    }

    public static void ConvertModifications(this IBioPolymerWithSetMods withSetMods, ISequenceConverter sequenceConverter)
    {
        ArgumentNullException.ThrowIfNull(sequenceConverter);
        ConvertModifications(withSetMods, sequenceConverter.Serializer);
    }

    public static void ConvertModifications(this IBioPolymerWithSetMods withSetMods, ISequenceSerializer sequenceSerializer)
    {
        ArgumentNullException.ThrowIfNull(withSetMods);
        ArgumentNullException.ThrowIfNull(sequenceSerializer);

        if (sequenceSerializer.HandlingMode == SequenceConversionHandlingMode.UsePrimarySequence)
        {
            withSetMods.AllModsOneIsNterminus.Clear();
            return;
        }

        var canonical = withSetMods.ToCanonicalSequence();
        var projected = sequenceSerializer.ToOneIsNterminusModificationDictionary(
            canonical,
            knownMods: null,
            warnings: null,
            mode: sequenceSerializer.HandlingMode);

        withSetMods.AllModsOneIsNterminus.Clear();
        foreach (var kvp in projected)
        {
            withSetMods.AllModsOneIsNterminus[kvp.Key] = kvp.Value;
        }
    }

    #endregion

    #region Proteins/Transcripts

    public static void ConvertModifications(this IBioPolymer withSetMods, ISequenceConverter sequenceConverter)
    {
        ArgumentNullException.ThrowIfNull(sequenceConverter);
        ConvertModifications(withSetMods, sequenceConverter.Serializer);
    }

    public static IBioPolymer ConvertModifications(this IBioPolymer withSetMods, IModificationLookup modificationLookup)
    {
        ConvertModifications(withSetMods, modificationLookup, SequenceConversionHandlingMode.RemoveIncompatibleElements);
        return withSetMods;
    }

    public static void ConvertModifications(this IBioPolymer withSetMods, ISequenceSerializer sequenceSerializer)
    {
        ArgumentNullException.ThrowIfNull(sequenceSerializer);
        var lookup = sequenceSerializer.ModificationLookup ?? GlobalModificationLookup.Instance;
        ConvertModifications(withSetMods, lookup, sequenceSerializer.HandlingMode);
    }

    #endregion

    private static void ConvertModifications(
        IBioPolymerWithSetMods withSetMods,
        IModificationLookup lookup,
        SequenceConversionHandlingMode mode)
    {
        ArgumentNullException.ThrowIfNull(withSetMods);
        ArgumentNullException.ThrowIfNull(lookup);

        if (mode == SequenceConversionHandlingMode.UsePrimarySequence)
        {
            withSetMods.AllModsOneIsNterminus.Clear();
            return;
        }

        var baseSequence = withSetMods.BaseSequence;
        var keys = withSetMods.AllModsOneIsNterminus.Keys.ToList();
        foreach (var key in keys)
        {
            if (!withSetMods.AllModsOneIsNterminus.TryGetValue(key, out var mod))
            {
                continue;
            }

            var position = ResolveWithSetModsPosition(baseSequence, key, out var residueIndex, out var targetResidue);
            var canonical = BuildCanonicalModification(mod, position, residueIndex, targetResidue);
            var resolved = lookup.TryResolve(canonical);
            if (!resolved.HasValue || resolved.Value.MzLibModification == null)
            {
                HandleUnresolvedModification(mode, withSetMods.AllModsOneIsNterminus, key, canonical);
                continue;
            }

            withSetMods.AllModsOneIsNterminus[key] = resolved.Value.MzLibModification;
        }
    }

    private static void ConvertModifications(
        IBioPolymer bioPolymer,
        IModificationLookup lookup,
        SequenceConversionHandlingMode mode)
    {
        ArgumentNullException.ThrowIfNull(bioPolymer);
        ArgumentNullException.ThrowIfNull(lookup);

        if (mode == SequenceConversionHandlingMode.UsePrimarySequence)
        {
            bioPolymer.OneBasedPossibleLocalizedModifications.Clear();
            bioPolymer.OriginalNonVariantModifications.Clear();
            foreach (var variant in bioPolymer.SequenceVariations)
            {
                variant.OneBasedModifications.Clear();
            }

            foreach (var variant in bioPolymer.AppliedSequenceVariations)
            {
                variant.OneBasedModifications.Clear();
            }

            return;
        }

        ConvertModificationLists(
            bioPolymer.OneBasedPossibleLocalizedModifications,
            index => ResolveBioPolymerPosition(bioPolymer.BaseSequence, index),
            lookup,
            mode);

        ConvertModificationLists(
            bioPolymer.OriginalNonVariantModifications,
            index => ResolveBioPolymerPosition(bioPolymer.BaseSequence, index),
            lookup,
            mode);

        foreach (var variant in bioPolymer.SequenceVariations)
        {
            ConvertModificationLists(
                variant.OneBasedModifications,
                index => ResolveVariantPosition(variant.VariantSequence, index),
                lookup,
                mode);
        }

        foreach (var variant in bioPolymer.AppliedSequenceVariations)
        {
            ConvertModificationLists(
                variant.OneBasedModifications,
                index => ResolveVariantPosition(variant.VariantSequence, index),
                lookup,
                mode);
        }

        if (bioPolymer.ConsensusVariant is not null && !ReferenceEquals(bioPolymer, bioPolymer.ConsensusVariant))
        {
            ConvertModifications(bioPolymer.ConsensusVariant, lookup, mode);
        }
    }

    private static void ConvertModificationLists(
        IDictionary<int, List<Modification>> modifications,
        Func<int, (ModificationPositionType positionType, int? residueIndex, char? targetResidue)> positionResolver,
        IModificationLookup lookup,
        SequenceConversionHandlingMode mode)
    {
        var keys = modifications.Keys.ToList();
        foreach (var key in keys)
        {
            if (!modifications.TryGetValue(key, out var list) || list == null || list.Count == 0)
            {
                modifications.Remove(key);
                continue;
            }

            var (positionType, residueIndex, targetResidue) = positionResolver(key);
            for (var index = 0; index < list.Count; index++)
            {
                var mod = list[index];
                var canonical = BuildCanonicalModification(mod, positionType, residueIndex, targetResidue);
                var resolved = lookup.TryResolve(canonical);
                if (!resolved.HasValue || resolved.Value.MzLibModification == null)
                {
                    HandleUnresolvedModification(mode, list, index, canonical);
                    if (mode != SequenceConversionHandlingMode.ThrowException)
                    {
                        index--;
                    }
                    continue;
                }

                list[index] = resolved.Value.MzLibModification;
            }

            if (list.Count == 0)
            {
                modifications.Remove(key);
            }
        }
    }

    private static ModificationPositionType ResolveWithSetModsPosition(
        string baseSequence,
        int index,
        out int? residueIndex,
        out char? targetResidue)
    {
        residueIndex = null;
        targetResidue = null;

        if (index == 1)
        {
            targetResidue = baseSequence.Length > 0 ? baseSequence[0] : (char?)null;
            return ModificationPositionType.NTerminus;
        }

        if (index == baseSequence.Length + 2)
        {
            targetResidue = baseSequence.Length > 0 ? baseSequence[^1] : (char?)null;
            return ModificationPositionType.CTerminus;
        }

        residueIndex = index - 2;
        if (residueIndex >= 0 && residueIndex < baseSequence.Length)
        {
            targetResidue = baseSequence[residueIndex.Value];
        }

        return ModificationPositionType.Residue;
    }

    private static (ModificationPositionType positionType, int? residueIndex, char? targetResidue) ResolveBioPolymerPosition(
        string baseSequence,
        int index)
    {
        if (index == 1)
        {
            return (ModificationPositionType.NTerminus, null, baseSequence.Length > 0 ? baseSequence[0] : (char?)null);
        }

        if (index == baseSequence.Length + 2)
        {
            return (ModificationPositionType.CTerminus, null, baseSequence.Length > 0 ? baseSequence[^1] : (char?)null);
        }

        var residueIndex = index - 1;
        var residue = residueIndex >= 0 && residueIndex < baseSequence.Length
            ? baseSequence[residueIndex]
            : (char?)null;
        return (ModificationPositionType.Residue, residueIndex, residue);
    }

    private static (ModificationPositionType positionType, int? residueIndex, char? targetResidue) ResolveVariantPosition(
        string? variantSequence,
        int index)
    {
        if (string.IsNullOrEmpty(variantSequence))
        {
            return (ModificationPositionType.Residue, null, 'X');
        }

        if (index == 1)
        {
            return (ModificationPositionType.Residue, 0, variantSequence[0]);
        }

        if (index <= variantSequence.Length)
        {
            var residueIndex = index - 1;
            return (ModificationPositionType.Residue, residueIndex, variantSequence[residueIndex]);
        }

        return (ModificationPositionType.Residue, null, 'X');
    }

    private static CanonicalModification BuildCanonicalModification(
        Modification mod,
        ModificationPositionType positionType,
        int? residueIndex,
        char? targetResidue)
    {
        var originalRepresentation = BuildOriginalRepresentation(mod);
        var modTarget = mod.Target?.ToString();
        var effectiveTargetResidue = !string.IsNullOrWhiteSpace(modTarget) && !string.Equals(modTarget, "X", StringComparison.OrdinalIgnoreCase)
            ? modTarget[0]
            : targetResidue;
        return positionType switch
        {
            ModificationPositionType.NTerminus => CanonicalModification.AtNTerminus(
                originalRepresentation,
                effectiveTargetResidue,
                mod.MonoisotopicMass,
                mod.ChemicalFormula,
                mzLibId: mod.IdWithMotif,
                mzLibModification: mod),
            ModificationPositionType.CTerminus => CanonicalModification.AtCTerminus(
                originalRepresentation,
                effectiveTargetResidue,
                mod.MonoisotopicMass,
                mod.ChemicalFormula,
                mzLibId: mod.IdWithMotif,
                mzLibModification: mod),
            _ => CanonicalModification.AtResidue(
                residueIndex ?? 0,
                effectiveTargetResidue ?? 'X',
                originalRepresentation,
                mod.MonoisotopicMass,
                mod.ChemicalFormula,
                mzLibId: mod.IdWithMotif,
                mzLibModification: mod)
        };
    }

    private static void HandleUnresolvedModification(
        SequenceConversionHandlingMode mode,
        IDictionary<int, Modification> modifications,
        int key,
        CanonicalModification canonical)
    {
        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new SequenceConversionException(
                $"Unable to resolve modification {canonical}",
                ConversionFailureReason.IncompatibleModifications,
                new[] { canonical.ToString() });
        }

        modifications.Remove(key);
    }

    private static void HandleUnresolvedModification(
        SequenceConversionHandlingMode mode,
        IList<Modification> modifications,
        int index,
        CanonicalModification canonical)
    {
        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new SequenceConversionException(
                $"Unable to resolve modification {canonical}",
                ConversionFailureReason.IncompatibleModifications,
                new[] { canonical.ToString() });
        }

        modifications.RemoveAt(index);
    }

    private static string BuildOriginalRepresentation(Modification mod)
    {
        if (!string.IsNullOrWhiteSpace(mod.ModificationType) && !string.IsNullOrWhiteSpace(mod.IdWithMotif))
        {
            return $"{mod.ModificationType}:{mod.IdWithMotif}";
        }

        if (!string.IsNullOrWhiteSpace(mod.IdWithMotif))
        {
            return mod.IdWithMotif;
        }

        if (!string.IsNullOrWhiteSpace(mod.OriginalId))
        {
            return mod.OriginalId;
        }

        return "Unknown";
    }
}
