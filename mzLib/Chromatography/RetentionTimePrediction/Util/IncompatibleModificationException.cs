using MzLibUtil;

namespace Chromatography.RetentionTimePrediction.Util;

/// <summary>
/// Exception thrown when incompatible modifications are encountered with ThrowException mode
/// </summary>
public class IncompatibleModificationException(string peptideFullSequence, string workingSequence, string predictorName)
    : MzLibException($"Peptide '{peptideFullSequence}' contains incompatible modification(s) " +
                     $"sequence was converted to {workingSequence} before {predictorName} gave up")
{
    public string OriginalSequence { get; } = peptideFullSequence;
    public string WorkingSequence { get; } = workingSequence;
    public string PredictorName { get; } = predictorName;
}