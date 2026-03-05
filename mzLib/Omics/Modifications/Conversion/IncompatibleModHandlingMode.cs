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
