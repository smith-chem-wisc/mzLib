namespace Omics.Modifications.Conversion;

public static class IncompatibleModHandlingModeExtensions
{
    public static SequenceConversionHandlingMode ToSequenceConversionHandlingMode(this IncompatibleModHandlingMode handlingMode)
    {
        return handlingMode switch
        {
            IncompatibleModHandlingMode.RemoveIncompatibleMods => SequenceConversionHandlingMode.RemoveIncompatibleMods,
            IncompatibleModHandlingMode.UsePrimarySequence => SequenceConversionHandlingMode.UsePrimarySequence,
            IncompatibleModHandlingMode.ReturnNull => SequenceConversionHandlingMode.ReturnNull,
            _ => SequenceConversionHandlingMode.ThrowException
        };
    }
}
