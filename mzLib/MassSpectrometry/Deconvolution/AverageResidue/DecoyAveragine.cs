using Chemistry;
using System;
using System.Linq;

namespace MassSpectrometry;

/// <summary>
/// A pre-computed lookup table of isotopic envelopes with a physically impossible isotope
/// spacing, used to generate decoy deconvolution results for target-decoy FDR estimation.
///
/// <para><b>Design rationale</b></para>
/// <para>
/// <see cref="Averagine"/> pre-computes 1500 theoretical isotope envelopes at static
/// construction time and stores them as lookup arrays indexed by mass. Every deconvolution
/// call is then just an array lookup — O(1) per spectrum peak.  This class follows the
/// identical pattern: all shifted envelopes are computed once in the instance constructor
/// and stored in private arrays of the same shape.
/// </para>
///
/// <para><b>How the decoy shift works</b></para>
/// <para>
/// The real C13-C12 isotope spacing is <see cref="Constants.C13MinusC12"/> ≈ 1.003355 Da.
/// Each isotope peak n steps above the monoisotopic peak sits at approximately
/// <c>monoisotopicMass + n * 1.003355 Da</c>.
/// </para>
/// <para>
/// The decoy replaces this spacing with <see cref="DecoyIsotopeSpacing"/> (default 0.9444 Da,
/// taken from OpenMS FLASHDeconv <c>noise_iso_delta_</c>).  This value is:
/// <list type="bullet">
///   <item>Close enough to the real spacing that the classic charge-detection heuristic
///   (<c>charge ≈ 1 / deltaMass</c>) still assigns plausible charge states.</item>
///   <item>Far enough that no real molecule produces this pattern, so surviving decoy
///   envelopes represent false-positive detections.</item>
/// </list>
/// </para>
/// <para>
/// For each envelope in the table, the n-th isotope peak is shifted by
/// <c>n * (decoySpacing − realSpacing)</c> relative to its real position, where n is the
/// peak's 0-based index counting up from the monoisotopic (lightest) peak.  The apex and
/// intensity values are unchanged — the envelope shape looks the same, only the positions
/// of the non-monoisotopic peaks move.
/// </para>
///
/// <para><b>Usage</b></para>
/// <code>
/// // Create once per decoy run — construction is O(N) where N = NumAveraginesToGenerate
/// var decoyModel = new DecoyAveragine(new Averagine());
///
/// // Plug into parameters — all subsequent deconvolution calls are O(1) lookups
/// var decoyParams = new ClassicDeconvolutionParameters(
///     minCharge, maxCharge, ppm, ratio,
///     averageResidueModel: decoyModel);
///
/// var decoys = Deconvoluter.Deconvolute(spectrum, decoyParams);
/// </code>
///
/// <para>Reference: Käll et al. (2008) target-decoy approach; OpenMS FLASHDeconv
/// <c>SpectralDeconvolution.cpp</c> <c>noise_iso_delta_ = 0.9444</c>.</para>
/// </summary>
public sealed class DecoyAveragine : AverageResidue
{
    // ── Constants ─────────────────────────────────────────────────────────────

    /// <summary>
    /// The default decoy isotope spacing (Da), taken from OpenMS FLASHDeconv
    /// <c>noise_iso_delta_</c>. Close to the real C13-C12 spacing (1.003355 Da)
    /// but physically impossible for any real molecule.
    /// </summary>
    public const double DefaultDecoyIsotopeSpacing = 0.9444;

    // ── Pre-computed lookup tables (same shape as Averagine's static arrays) ──

    /// <summary>Shifted mass arrays, one per averagine, in intensity-descending order.</summary>
    private readonly double[][] _allMasses;

    /// <summary>Mass of the apex (most intense) isotope peak for each averagine entry.</summary>
    private readonly double[] _mostIntenseMasses;

    // ── Delegates to the real model ───────────────────────────────────────────

    /// <summary>The underlying real Averagine model (or other AverageResidue).</summary>
    private readonly AverageResidue _real;

    // ── Public properties ─────────────────────────────────────────────────────

    /// <summary>The decoy isotope spacing used to build this table (Da).</summary>
    public double DecoyIsotopeSpacing { get; }

    // ── AverageResidue overrides ──────────────────────────────────────────────

    /// <inheritdoc />
    public override int GetMostIntenseMassIndex(double testMass)
        => _mostIntenseMasses.GetClosestIndex(testMass);

    /// <summary>
    /// Returns the pre-computed shifted mass array for this averagine entry.
    /// Array is in intensity-descending order (apex at index 0), matching the
    /// convention of <see cref="Averagine.AllMasses"/>.
    /// </summary>
    public override double[] GetAllTheoreticalMasses(int index)
        => _allMasses[index];

    /// <summary>
    /// Returns the real Averagine intensity array unchanged.
    /// The decoy has the same envelope shape — only peak positions are shifted.
    /// </summary>
    public override double[] GetAllTheoreticalIntensities(int index)
        => _real.GetAllTheoreticalIntensities(index);

    /// <summary>
    /// Returns the real DiffToMonoisotopic unchanged.
    /// The monoisotopic peak (index 0 by convention from the real model's perspective)
    /// is not shifted, so this offset remains valid for monoisotopic mass recovery.
    /// </summary>
    public override double GetDiffToMonoisotopic(int index)
        => _real.GetDiffToMonoisotopic(index);

    // ── Constructor ───────────────────────────────────────────────────────────

    /// <summary>
    /// Constructs a <see cref="DecoyAveragine"/> by pre-computing shifted mass tables
    /// from the supplied real model.
    /// </summary>
    /// <param name="realModel">
    ///   The real averagine model to derive intensities and DiffToMonoisotopic from.
    ///   Typically <c>new Averagine()</c>.
    /// </param>
    /// <param name="decoyIsotopeSpacing">
    ///   The isotope spacing (Da) to use for decoy peaks. Defaults to
    ///   <see cref="DefaultDecoyIsotopeSpacing"/> (0.9444 Da).
    /// </param>
    /// <exception cref="ArgumentException">
    ///   Thrown if <paramref name="decoyIsotopeSpacing"/> equals
    ///   <see cref="Constants.C13MinusC12"/> — that would produce identical
    ///   envelopes to the real model, which is not a valid decoy.
    /// </exception>
    public DecoyAveragine(AverageResidue realModel,
        double decoyIsotopeSpacing = DefaultDecoyIsotopeSpacing, double targetIsotopeSpacing = Constants.C13MinusC12)
    {
        if (Math.Abs(decoyIsotopeSpacing - Constants.C13MinusC12) < 1e-6)
            throw new ArgumentException(
                $"{nameof(decoyIsotopeSpacing)} must differ from the real C13-C12 " +
                $"spacing ({Constants.C13MinusC12} Da). Decoy envelopes would be " +
                $"identical to target envelopes.", nameof(decoyIsotopeSpacing));

        _real = realModel;
        DecoyIsotopeSpacing = decoyIsotopeSpacing;

        // Per-isotope shift relative to the real spacing.
        // Positive offset → peaks spread further apart.
        // Negative offset (default) → peaks compressed together.
        double perPeakOffset = decoyIsotopeSpacing - targetIsotopeSpacing;

        _allMasses = new double[NumAveraginesToGenerate][];
        _mostIntenseMasses = new double[NumAveraginesToGenerate];

        for (int i = 0; i < NumAveraginesToGenerate; i++)
        {
            double[] realMasses = realModel.GetAllTheoreticalMasses(i);

            // realMasses is intensity-descending (apex at [0]).
            // We need to know each peak's isotope index, defined as how many
            // neutrons above the monoisotopic peak it sits.
            //
            // The monoisotopic peak is the LIGHTEST peak in the envelope.
            // Its mass = realMasses[0] - DiffToMonoisotopic[i]
            //           = apex mass - (apex mass - monoisotopic mass)
            //           = monoisotopic mass.
            //
            // For each peak at mass m, its isotope index is:
            //   n = round((m - monoisotopicMass) / C13MinusC12)
            //
            // We then shift peak m to: m + n * perPeakOffset
            //
            // This means:
            //   monoisotopic peak (n=0): no shift — its mass is unchanged.
            //   +1 Da peak (n=1):        shifted by 1 * perPeakOffset.
            //   +2 Da peak (n=2):        shifted by 2 * perPeakOffset.
            //   etc.
            //
            // The intensity-descending order is preserved because the shift is
            // monotonically proportional to mass distance — peaks don't cross.

            double monoisotopicMass = realMasses[0] - realModel.GetDiffToMonoisotopic(i);

            double[] shiftedMasses = new double[realMasses.Length];
            for (int j = 0; j < realMasses.Length; j++)
            {
                // Isotope index: number of neutron-mass steps above monoisotopic.
                int n = (int)Math.Round(
                    (realMasses[j] - monoisotopicMass) / Constants.C13MinusC12);

                // Ensure n is non-negative (floating-point rounding guard).
                n = Math.Max(0, n);

                shiftedMasses[j] = realMasses[j] + n * perPeakOffset;
            }

            _allMasses[i] = shiftedMasses;
            // Apex is still at index 0 — same intensity ordering is preserved.
            _mostIntenseMasses[i] = shiftedMasses[0];
        }
    }
}