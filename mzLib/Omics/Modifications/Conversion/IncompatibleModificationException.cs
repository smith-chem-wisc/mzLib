using MzLibUtil;

namespace Omics.Modifications.Conversion;

/// <summary>
/// Exception thrown when incompatible modifications are encountered and the selected handling mode requires an error.
/// </summary>
public class IncompatibleModificationException(string? peptideFullSequence, string? workingSequence, string? contextName)
    : MzLibException(
        $"Peptide '{peptideFullSequence}' contains incompatible modification(s); sequence was converted to {workingSequence} before {contextName} aborted")
{
    public string? OriginalSequence { get; } = peptideFullSequence;
    public string? WorkingSequence { get; } = workingSequence;
    public string? ContextName { get; } = contextName;
}
