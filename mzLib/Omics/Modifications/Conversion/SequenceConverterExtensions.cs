namespace Omics.Modifications.Conversion;

public static class SequenceConverterExtensions
{
    public static string ToMetaMorpheusSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        return SequenceConverter.Default.ConvertFullSequence(bioPolymer, ModificationNamingConvention.MetaMorpheus, handlingMode);
    }

    public static string ToUnimodSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        return SequenceConverter.Default.ConvertFullSequence(bioPolymer, ModificationNamingConvention.Unimod, handlingMode);
    }

    public static string ToUniProtSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        return SequenceConverter.Default.ConvertFullSequence(bioPolymer, ModificationNamingConvention.UniProt, handlingMode);
    }

    public static string ToChronologerSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.RemoveIncompatibleMods)
    {
        return SequenceTargetConverter.Convert(bioPolymer, SequenceConversionTarget.Chronologer, handlingMode);
    }
}
