namespace Omics.Modifications.Conversion;

/// <summary>
/// Convenience extensions for sequence conversion.
/// </summary>
public static class SequenceConverterExtensions
{
    /// <summary>
    /// Converts to MetaMorpheus naming convention.
    /// </summary>
    public static string ToMetaMorpheusSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        return SequenceConverter.Default.ConvertFullSequence(bioPolymer, ModificationNamingConvention.MetaMorpheus, handlingMode);
    }

    /// <summary>
    /// Converts to Unimod naming convention.
    /// </summary>
    public static string ToUnimodSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        return SequenceConverter.Default.ConvertFullSequence(bioPolymer, ModificationNamingConvention.Unimod, handlingMode);
    }

    /// <summary>
    /// Converts to UniProt naming convention.
    /// </summary>
    public static string ToUniProtSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
    {
        return SequenceConverter.Default.ConvertFullSequence(bioPolymer, ModificationNamingConvention.UniProt, handlingMode);
    }

    /// <summary>
    /// Converts to Chronologer-compatible sequence encoding.
    /// </summary>
    public static string ToChronologerSequence(this IBioPolymerWithSetMods bioPolymer, SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.RemoveIncompatibleMods)
    {
        return SequenceTargetConverter.Convert(bioPolymer, SequenceConversionTarget.Chronologer, handlingMode);
    }
}
