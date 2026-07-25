#nullable enable
using System.Collections.Generic;
using Omics;
using Omics.Fragmentation;

namespace TopDownSimulator.Ms2;

/// <summary>
/// One simulated MS2 peak: an mzLib <see cref="Product"/> observed at a particular charge and
/// isotopologue, with the ion count assigned to it by the fragmentation model.
/// </summary>
/// <param name="Product">
/// The mzLib fragment. All mass chemistry (terminal caps, product-type shifts, modifications)
/// comes from <c>IBioPolymerWithSetMods.Fragment</c> — this project computes none of it.
/// </param>
/// <param name="OneBasedBondNumber">
/// Which backbone bond produced the ion, counting from the N-terminus. Bond <c>k</c> joins residue
/// <c>k</c> to residue <c>k + 1</c>. For an N-terminal series this equals
/// <c>Product.FragmentNumber</c>; for the complementary C-terminal series it is
/// <c>sequenceLength - Product.FragmentNumber</c>.
/// </param>
/// <param name="Charge">Charge state at which the fragment is observed.</param>
/// <param name="IsotopologueIndex">0 for the monoisotopic peak, 1 for M+1, and so on.</param>
/// <param name="Mz">Observed m/z, from <c>Chemistry.ClassExtensions.ToMz</c>.</param>
/// <param name="Intensity">Ion count assigned to this peak.</param>
public sealed record SimulatedFragmentIon(
    Product Product,
    int OneBasedBondNumber,
    int Charge,
    int IsotopologueIndex,
    double Mz,
    double Intensity);

/// <summary>
/// Output of an <see cref="IMs2FragmentationModel"/> for one precursor.
/// </summary>
/// <param name="BaseSequence">Unmodified sequence that was fragmented.</param>
/// <param name="StartingIonCount">
/// Ion count the model started from, before any bond consumed intensity.
/// </param>
/// <param name="Fragments">Every simulated fragment peak, in no particular m/z order.</param>
/// <param name="SurvivingPrecursorIntensity">
/// Ion count still in the intact precursor after the whole backbone has been walked.
/// <c>Fragments.Sum(Intensity) + SurvivingPrecursorIntensity == StartingIonCount</c>.
/// </param>
public sealed record Ms2FragmentationResult(
    string BaseSequence,
    double StartingIonCount,
    IReadOnlyList<SimulatedFragmentIon> Fragments,
    double SurvivingPrecursorIntensity);

/// <summary>
/// Turns a precursor into a list of MS2 fragment peaks with intensities.
/// </summary>
/// <remarks>
/// This is the swap point for the whole MS2 spectrum model. <see cref="Ms2Simulator"/> and
/// everything downstream of it depend only on this interface, never on
/// <see cref="PropensityCascadeFragmentationModel"/>; replacing the model means passing a different
/// implementation to the <see cref="Ms2Simulator"/> constructor and changing nothing else.
/// </remarks>
public interface IMs2FragmentationModel
{
    /// <summary>Fragments a precursor and assigns an intensity to every product ion.</summary>
    Ms2FragmentationResult Fragment(IBioPolymerWithSetMods precursor);
}
