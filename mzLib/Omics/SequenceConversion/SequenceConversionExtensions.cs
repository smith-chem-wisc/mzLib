using Omics.Modifications;

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
