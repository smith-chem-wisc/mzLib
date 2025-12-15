
namespace Omics
{
    public interface IHasSequenceCoverageFromFragments
    {
        /// <summary>
        /// Positions in the biopolymer sequence (one-based) that are covered by fragment ions.
        /// Populated by <see cref="GetSequenceCoverage"/> when fragment coverage data is available.
        /// May be null if coverage has not been calculated or is not available.
        /// </summary>
        HashSet<int>? FragmentCoveragePositionInPeptide { get; }

        /// <summary>
        /// Calculates sequence coverage from fragment ions for this spectral match.
        /// Populates <see cref="FragmentCoveragePositionInPeptide"/> with one-based positions
        /// of residues that are covered by matched fragment ions.
        /// Works for any biopolymer type (proteins, nucleic acids, etc.).
        /// Implementations without fragment ion data may leave the property null or empty.
        /// </summary>
        void GetSequenceCoverage();
    }
}
