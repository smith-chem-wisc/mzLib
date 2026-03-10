namespace Omics.SequenceConversion;

public class EssentialSequenceSerializer : SequenceSerializerBase
{    
    /// <summary>
    /// Singleton instance.
    /// </summary>
    public static EssentialSequenceSerializer Instance { get; } = new();

    private readonly IReadOnlyDictionary<string, int> _modsToWritePruned;

    public EssentialSequenceSerializer(
        IReadOnlyDictionary<string, int>? modsToWritePruned = null,
        IModificationLookup? lookup = null)
        : base(lookup ?? GlobalModificationLookup.Instance)
    {
        // Copied straight from #MeatMorpheus.TaskLayer.SearchTask.SearchParameters
        _modsToWritePruned = modsToWritePruned ?? new Dictionary<string, int>
        {
            //Key is modification type.

            //Value is integer 0, 1, 2 and 3 interpreted as:
            //   0:   Do not Write
            //   1:   Write if in DB and Observed
            //   2:   Write if in DB
            //   3:   Write if Observed

            {"N-linked glycosylation", 3},
            {"O-linked glycosylation", 3},
            {"Other glycosylation", 3},
            {"Common Biological", 3},
            {"Less Common", 3},
            {"Metal", 3},
            {"2+ nucleotide substitution", 3},
            {"1 nucleotide substitution", 3},
            {"UniProt", 2},
        };
    }

    public override string FormatName => EssentialSequenceFormatSchema.Instance.FormatName;

    public override SequenceFormatSchema Schema => EssentialSequenceFormatSchema.Instance;

    public override bool CanSerialize(CanonicalSequence sequence)
    {
        return !string.IsNullOrEmpty(sequence.BaseSequence);
    }

    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        if (mod.MzLibModification != null)
        {
            return false;
        }

        return !TryParseTypeAndId(mod.OriginalRepresentation, out _, out _)
               && !TryParseTypeAndId(mod.MzLibId, out _, out _);
    }

    protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        if (!TryGetTypeAndId(mod, out var modificationType, out var idWithMotif))
        {
            warnings.AddIncompatibleItem(mod.ToString());

            if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
            {
                warnings.AddWarning($"Removing incompatible modification for essential sequence: {mod}");
                return null;
            }

            if (mode == SequenceConversionHandlingMode.ThrowException)
            {
                throw new SequenceConversionException(
                    $"Cannot serialize modification in essential format - missing modification type or ID: {mod}",
                    ConversionFailureReason.IncompatibleModifications,
                    new[] { mod.ToString() });
            }

            return null;
        }

        if (!_modsToWritePruned.ContainsKey(modificationType))
        {
            return null;
        }

        return $"{modificationType}:{idWithMotif}";
    }

    private static bool TryGetTypeAndId(CanonicalModification mod, out string modificationType, out string idWithMotif)
    {
        if (mod.MzLibModification != null
            && !string.IsNullOrWhiteSpace(mod.MzLibModification.ModificationType)
            && !string.IsNullOrWhiteSpace(mod.MzLibModification.IdWithMotif))
        {
            modificationType = mod.MzLibModification.ModificationType;
            idWithMotif = mod.MzLibModification.IdWithMotif;
            return true;
        }

        if (TryParseTypeAndId(mod.OriginalRepresentation, out modificationType, out idWithMotif))
        {
            return true;
        }

        if (TryParseTypeAndId(mod.MzLibId, out modificationType, out idWithMotif))
        {
            return true;
        }

        modificationType = string.Empty;
        idWithMotif = string.Empty;
        return false;
    }

    private static bool TryParseTypeAndId(string? value, out string modificationType, out string idWithMotif)
    {
        modificationType = string.Empty;
        idWithMotif = string.Empty;

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

        modificationType = trimmed.Substring(0, separatorIndex).Trim();
        idWithMotif = trimmed.Substring(separatorIndex + 1).Trim();

        return !string.IsNullOrEmpty(modificationType) && !string.IsNullOrEmpty(idWithMotif);
    }
}
