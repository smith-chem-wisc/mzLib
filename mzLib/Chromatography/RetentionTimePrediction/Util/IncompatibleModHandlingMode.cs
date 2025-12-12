namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Mode for handling incompatible modifications during RT prediction.
/// Simple enum - no object allocation needed.
/// </summary>
public enum IncompatibleModHandlingMode
{
    /// <summary>
    /// Strip unsupported modifications and predict with remaining mods
    /// </summary>
    RemoveIncompatibleMods,

    /// <summary>
    /// Ignore all modifications and use only base sequence
    /// </summary>
    UsePrimarySequence,

    /// <summary>
    /// Throw exception if incompatible modifications present
    /// </summary>
    ThrowException,

    /// <summary>
    /// Return null if incompatible modifications present
    /// </summary>
    ReturnNull
}