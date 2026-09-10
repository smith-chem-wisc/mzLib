using Chemistry;
using System;

namespace MassSpectrometry;

/// <summary>
/// A decoy isotope model that permutes envelope <i>intensities</i> while leaving peak
/// <i>positions</i> on the real <sup>13</sup>C lattice.
///
/// <para><b>Why this exists alongside <see cref="DecoyAveragine"/></b></para>
/// <para>
/// <see cref="DecoyAveragine"/> builds its decoy by moving peaks: the n-th isotope peak is
/// displaced by <c>n * (decoySpacing - C13MinusC12)</c>, a fixed offset in <b>neutral mass</b>
/// (58.955 mDa per isotope step at the default 0.9444 Da spacing).
/// </para>
/// <para>
/// Deconvolution, however, matches peaks in <b>m/z</b>, not in neutral mass. A neutral-mass
/// displacement of <c>n * 58.955</c> mDa is seen by the matcher as <c>n * 58.955 / z</c> mDa,
/// so the decoy's distance from the real lattice shrinks as <c>1/z</c>. At m/z ≈ 900 the first
/// decoy tooth sits:
/// </para>
/// <list type="table">
///   <listheader><term>charge</term><description>offset of the first decoy tooth</description></listheader>
///   <item><term>z = 2</term><description>29.5 mDa ≈ 33 ppm — comfortably outside any matching tolerance</description></item>
///   <item><term>z = 5</term><description>11.8 mDa ≈ 13 ppm — inside a 20 ppm tolerance</description></item>
///   <item><term>z = 15</term><description>3.9 mDa ≈ 4.4 ppm</description></item>
///   <item><term>z = 30</term><description>2.0 mDa ≈ 2.2 ppm — inside even a 4 ppm tolerance</description></item>
/// </list>
/// <para>
/// Past roughly z = 4 at a 20 ppm tolerance (z = 17 at 4 ppm), the "impossible" comb lands on the
/// real peaks and the decoy stops being a decoy: it is scored as though it were the target, the
/// decoy score distribution creeps toward the target distribution, and the estimated FDR is
/// biased low. Bottom-up precursors are mostly z = 2–4, which is where the 0.9444 Da constant
/// comes from (OpenMS FLASHDeconv <c>noise_iso_delta_</c>) and where it behaves correctly.
/// Top-down proteoforms routinely carry z = 10–40, which is exactly where it degenerates.
/// </para>
///
/// <para><b>How this model avoids the problem</b></para>
/// <para>
/// Peak <i>positions</i> are left untouched, so there is no offset to be divided by charge —
/// the decoy is charge-invariant by construction. What is falsified instead is the envelope's
/// <i>shape</i>: the theoretical intensities are permuted, so a real isotope envelope and this
/// model disagree about which peaks should be tall. A scorer that rewards fitting the observed
/// intensity pattern (rather than merely finding peaks at the right m/z) separates the two at
/// any charge state.
/// </para>
/// <para>
/// The apex intensity is deliberately held at its own index rather than permuted with the rest.
/// Callers and the deconvolution algorithms treat index 0 of the intensity-descending arrays as
/// the apex; permuting it would change envelope selection behaviour rather than just the decoy's
/// shape, and would confound "this decoy is hard to match" with "this decoy is indexed oddly".
/// </para>
///
/// <para><b>Usage</b></para>
/// <code>
/// var decoyModel = new ShuffledAveragine(new Averagine());
/// var decoyParams = new ClassicDeconvolutionParameters(
///     minCharge, maxCharge, ppm, ratio, averageResidueModel: decoyModel);
/// var decoys = Deconvoluter.Deconvolute(spectrum, decoyParams);
/// </code>
///
/// <para>
/// Which decoy is stronger is an empirical question and is expected to depend on the regime;
/// this class exists so the comparison can be made rather than assumed. Nothing here changes
/// the default: <see cref="DecoyAveragine"/> remains what
/// <see cref="DeconvolutionParameters.ToDecoyParameters"/> produces.
/// </para>
/// </summary>
public sealed class ShuffledAveragine : AverageResidue
{
    /// <summary>
    /// Default permutation seed. Fixed so that a decoy run is reproducible; vary it to check
    /// that a result does not depend on one particular shuffle.
    /// </summary>
    public const int DefaultShuffleSeed = 42;

    private readonly AverageResidue _real;
    private readonly double[][] _shuffledIntensities;

    /// <summary>The seed used to build this model's permutations.</summary>
    public int ShuffleSeed { get; }

    /// <summary>
    /// Builds the permuted intensity table once, mirroring <see cref="DecoyAveragine"/>'s
    /// precompute-at-construction pattern so that deconvolution stays an O(1) lookup.
    /// </summary>
    /// <param name="realModel">The real model whose envelopes are being falsified.</param>
    /// <param name="shuffleSeed">Permutation seed; see <see cref="DefaultShuffleSeed"/>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="realModel"/> is null.</exception>
    public ShuffledAveragine(AverageResidue realModel, int shuffleSeed = DefaultShuffleSeed)
    {
        _real = realModel ?? throw new ArgumentNullException(nameof(realModel));
        ShuffleSeed = shuffleSeed;

        _shuffledIntensities = new double[NumAveraginesToGenerate][];
        var rng = new Random(shuffleSeed);

        for (int i = 0; i < NumAveraginesToGenerate; i++)
        {
            double[] real = realModel.GetAllTheoreticalIntensities(i);
            double[] shuffled = (double[])real.Clone();

            // Fisher-Yates over indices 1..n-1, leaving the apex at index 0 in place.
            // An envelope with fewer than three peaks has nothing to permute below the apex,
            // so it is left identical to the real model rather than silently passed off as a
            // decoy -- see IsDegenerate.
            for (int j = shuffled.Length - 1; j > 1; j--)
            {
                int k = rng.Next(1, j + 1);
                (shuffled[j], shuffled[k]) = (shuffled[k], shuffled[j]);
            }

            _shuffledIntensities[i] = shuffled;
        }
    }

    /// <summary>
    /// True when the envelope at <paramref name="index"/> is too short to permute (fewer than
    /// three peaks), so this model returns the real intensities and is not a decoy there.
    /// Exposed so that callers computing FDR can exclude these rather than count them as decoys.
    /// </summary>
    public bool IsDegenerate(int index) => _real.GetAllTheoreticalIntensities(index).Length < 3;

    /// <inheritdoc />
    public override int GetMostIntenseMassIndex(double testMass)
        => _real.GetMostIntenseMassIndex(testMass);

    /// <summary>
    /// Real masses, unchanged. Positions stay on the <sup>13</sup>C lattice, which is what makes
    /// this decoy charge-invariant.
    /// </summary>
    public override double[] GetAllTheoreticalMasses(int index)
        => _real.GetAllTheoreticalMasses(index);

    /// <summary>Permuted intensities: the falsified envelope shape.</summary>
    public override double[] GetAllTheoreticalIntensities(int index)
        => _shuffledIntensities[index];

    /// <inheritdoc />
    public override double GetDiffToMonoisotopic(int index)
        => _real.GetDiffToMonoisotopic(index);

    #region IEquatable

    /// <inheritdoc />
    protected override bool EqualProperties(AverageResidue other)
        => other is ShuffledAveragine o && ShuffleSeed == o.ShuffleSeed && _real.Equals(o._real);

    /// <inheritdoc />
    protected override void AddHashCodes(HashCode hash)
    {
        hash.Add(ShuffleSeed);
        hash.Add(_real);
    }

    #endregion
}
