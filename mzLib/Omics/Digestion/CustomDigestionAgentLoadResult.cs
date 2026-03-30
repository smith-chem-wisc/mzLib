namespace Omics.Digestion
{
    /// <summary>
    /// Represents the outcome of a custom digestion agent (protease or RNase) load operation.
    /// Returned by <c>LoadAndMergeCustomProteases</c> and <c>LoadAndMergeCustomRnases</c> to
    /// give callers full visibility into what was accepted and what was rejected.
    /// </summary>
    /// <param name="Added">
    /// Names of digestion agents that were successfully added to the dictionary from the custom file(s).
    /// Only agents whose names did not already exist in the embedded resource are included here.
    /// </param>
    /// <param name="Skipped">
    /// Names of digestion agents that were rejected because a same-named agent already exists in the
    /// embedded resource or was already loaded from an earlier custom file in the same call.
    /// The embedded definition is always retained; the custom definition is discarded.
    /// </param>
    public record CustomDigestionAgentLoadResult(
        List<string> Added,
        List<string> Skipped);
}