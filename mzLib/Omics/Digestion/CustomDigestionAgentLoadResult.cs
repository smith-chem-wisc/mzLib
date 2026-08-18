using System.Collections.ObjectModel;

namespace Omics.Digestion
{
    /// <summary>
    /// Represents the outcome of a custom digestion agent (protease or RNase) load operation.
    /// Returned by <c>LoadAndMergeCustomProteases</c> and <c>LoadAndMergeCustomRnases</c> to
    /// give callers full visibility into what was accepted and what was rejected.
    /// Both collections are exposed as <see cref="IReadOnlyList{T}"/> backed by
    /// <see cref="ReadOnlyCollection{T}"/> wrapping a defensive copy of the constructor input,
    /// so callers cannot mutate the reported results even if they retain a reference to the
    /// original list passed in.
    /// </summary>
    public record CustomDigestionAgentLoadResult
    {
        /// <summary>
        /// Names of digestion agents that were successfully added to the dictionary from the custom file(s).
        /// Only agents whose names did not already exist in the embedded resource are included here.
        /// This collection is read-only and cannot be modified by the caller.
        /// </summary>
        public IReadOnlyList<string> Added { get; }

        /// <summary>
        /// Names of digestion agents that were rejected because a same-named agent already exists in the
        /// embedded resource, was previously loaded by an earlier call, or was already added by an earlier
        /// file in the same batch. These three cases are intentionally not distinguished; callers that need
        /// to differentiate "collides with built-in" from "duplicate across your custom files" must compare
        /// against the embedded baseline themselves. The embedded definition is always retained; the custom
        /// definition is discarded. This collection is read-only and cannot be modified by the caller.
        /// </summary>
        public IReadOnlyList<string> Skipped { get; }

        public CustomDigestionAgentLoadResult(IEnumerable<string> Added, IEnumerable<string> Skipped)
        {
            this.Added = new ReadOnlyCollection<string>(Added?.ToList() ?? new List<string>());
            this.Skipped = new ReadOnlyCollection<string>(Skipped?.ToList() ?? new List<string>());
        }
    }
}
