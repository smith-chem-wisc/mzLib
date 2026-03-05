namespace Omics.Modifications.Conversion;

public static class SequenceTargetConverter
{
    public static bool TryConvert(
        IBioPolymerWithSetMods bioPolymer,
        SequenceConversionTarget target,
        SequenceConversionHandlingMode handlingMode,
        out string? converted,
        out SequenceConversionFailureReason? reason,
        ISequenceConverter? converter = null)
    {
        ArgumentNullException.ThrowIfNull(bioPolymer);
        converter ??= SequenceConverter.Default;

        if (target == SequenceConversionTarget.Chronologer)
        {
            if (!converter.TryGetChronologerSequence(bioPolymer, handlingMode, out var massShiftSequence, out reason) || massShiftSequence == null)
            {
                converted = null;
                return false;
            }

            return ChronologerSequenceFormatter.TryFormatChronologerSequence(
                bioPolymer,
                massShiftSequence,
                handlingMode,
                out converted,
                out reason);
        }

        var convention = target switch
        {
            SequenceConversionTarget.MetaMorpheus => ModificationNamingConvention.MetaMorpheus,
            SequenceConversionTarget.UniProt => ModificationNamingConvention.UniProt,
            SequenceConversionTarget.Unimod => ModificationNamingConvention.Unimod,
            _ => ModificationNamingConvention.MetaMorpheus
        };

        return converter.TryConvertFullSequence(
            bioPolymer,
            convention,
            handlingMode,
            out converted,
            out reason);
    }

    public static string Convert(
        IBioPolymerWithSetMods bioPolymer,
        SequenceConversionTarget target,
        SequenceConversionHandlingMode handlingMode,
        ISequenceConverter? converter = null)
    {
        if (!TryConvert(bioPolymer, target, handlingMode, out var converted, out var reason, converter) || converted == null)
        {
            throw new SequenceConversionException(
                $"Unable to convert sequence to {target} (reason: {reason})",
                reason ?? SequenceConversionFailureReason.NoEquivalent);
        }

        return converted;
    }
}
