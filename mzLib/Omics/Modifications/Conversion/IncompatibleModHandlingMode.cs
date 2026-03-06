namespace Omics.Modifications.Conversion;

/// <summary>
/// Mode for handling incompatible modifications when converting sequences or formatting
/// for downstream consumers such as retention time predictors.
/// </summary>
public enum IncompatibleModHandlingMode
{
    RemoveIncompatibleMods,
    UsePrimarySequence,
    ThrowException,
    ReturnNull
}
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