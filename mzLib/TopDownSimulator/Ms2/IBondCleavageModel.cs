#nullable enable
namespace TopDownSimulator.Ms2;

/// <summary>
/// Supplies the probability that a single backbone peptide bond cleaves during MS2 activation.
/// </summary>
/// <remarks>
/// This is the innermost, most-likely-to-be-replaced seam of the MS2 model: swap this to change
/// <em>where</em> the backbone breaks without touching how intensity is apportioned or how spectra
/// are built. <see cref="ResiduePairCleavagePropensityModel"/> is the shipped implementation, which
/// reads a published amino-acid-pair propensity table.
/// </remarks>
public interface IBondCleavageModel
{
    /// <summary>
    /// Probability in [0, 1] that the peptide bond at <paramref name="zeroBasedBondIndex"/> cleaves.
    /// </summary>
    /// <param name="baseSequence">
    /// Unmodified amino-acid sequence written N-terminus to C-terminus.
    /// </param>
    /// <param name="zeroBasedBondIndex">
    /// Index of the bond joining <c>baseSequence[zeroBasedBondIndex]</c> (N-terminal side) to
    /// <c>baseSequence[zeroBasedBondIndex + 1]</c> (C-terminal side). Valid range is
    /// <c>0 .. baseSequence.Length - 2</c>.
    /// </param>
    double GetCleavageProbability(string baseSequence, int zeroBasedBondIndex);
}
