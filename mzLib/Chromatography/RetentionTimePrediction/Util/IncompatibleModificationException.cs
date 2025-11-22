namespace Chromatography.RetentionTimePrediction.Util;

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