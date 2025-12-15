
namespace Omics
{
    public interface IHasSequenceCoverageFromFragments
    {
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
