namespace Omics
{
    /// <summary>
    /// Interface for objects that can calculate and provide fragment-level sequence coverage.
    /// Typically implemented by spectral match types that have matched fragment ion data
    /// which can be used to determine which residues in the sequence are covered by fragments.
    /// </summary>
    public interface IHasSequenceCoverageFromFragments
    {
        /// <summary>
        /// Positions in the biopolymer sequence (one-based) that are covered by fragment ions.
        /// Populated by <see cref="GetSequenceCoverage"/> when fragment coverage data is available.
        /// Returns null if coverage has not been calculated or no fragment data is available.
        /// </summary>
        HashSet<int>? FragmentCoveragePositionInPeptide { get; }

        /// <summary>
        /// Calculates sequence coverage from fragment ions and populates 
        /// <see cref="FragmentCoveragePositionInPeptide"/> with one-based positions
        /// of residues that are covered by matched fragment ions.
        /// 
        /// Coverage is typically determined by identifying which residues have 
        /// fragments on both sides (N-terminal and C-terminal for peptides, 
        /// or 5' and 3' for nucleic acids).
        /// 
        /// Implementations without fragment ion data should leave the property null.
        /// This method may be called multiple times; implementations should 
        /// recalculate or return cached results as appropriate.
        /// </summary>
        void GetSequenceCoverage();
    }
}