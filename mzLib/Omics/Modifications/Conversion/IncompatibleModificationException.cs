using MzLibUtil;

namespace Omics.Modifications.Conversion;

/// <summary>
/// Exception thrown when incompatible modifications are encountered and the selected handling mode requires an error.
/// </summary>
public class IncompatibleModificationException(string? peptideFullSequence, string? workingSequence, string? predictorName)
    : MzLibException(
        $"Peptide '{peptideFullSequence}' contains incompatible modification(s); sequence was converted to {workingSequence} before {predictorName} aborted")
{
    public string? OriginalSequence { get; } = peptideFullSequence;
    public string? WorkingSequence { get; } = workingSequence;
    public string? PredictorName { get; } = predictorName;
}
