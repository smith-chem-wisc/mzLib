namespace Chromatography.RetentionTimePrediction.Util;

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

/// <summary>
/// Exception thrown when incompatible modifications are encountered with ThrowException mode
/// </summary>
public class IncompatibleModificationException : Exception
{
    public string PeptideFullSequence { get; }
    public string WorkingSequence { get; }

    public IncompatibleModificationException(
        string peptideFullSequence,
        string workingSequence,
        string predictorName)
        : base($"Peptide '{peptideFullSequence}' contains incompatible modification(s) " +
               $"sequence was converted to {workingSequence} before {predictorName} gave up")
    {
        PeptideFullSequence = peptideFullSequence;
        WorkingSequence = workingSequence;
    }
}