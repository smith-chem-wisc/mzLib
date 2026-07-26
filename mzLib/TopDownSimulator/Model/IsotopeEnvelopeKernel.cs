#nullable enable
using System;
using System.Collections.Generic;
using Chemistry;
using MassSpectrometry;

namespace TopDownSimulator.Model;

/// <summary>
/// φ(b; M, z, σ_m) — theoretical isotope envelope for a neutral monoisotopic mass M
/// at charge z, expressed as a sum of Gaussian peaks of width σ_m centered on the
/// predicted m/z of each isotopologue. Intensities are normalized so that the
/// sum over the isotopologue weights is 1.
/// </summary>
/// <remarks>
/// Uses mzLib's averagine precomputed tables. For a target mass M we look up the
/// averagine index with the closest most-intense mass, then shift every theoretical
/// isotopologue mass by (M − averagine_monoisotopic) so that the returned envelope
/// is anchored on the caller's monoisotopic mass.
/// <para>
/// Isotopologues are ordered by <b>ascending neutral mass</b>, and therefore by ascending m/z
/// for every charge: index 0 is the lowest-mass isotopologue, not the most intense one. Note
/// that this is the opposite of <see cref="Averagine"/>'s own tables, which are stored
/// most-intense-first; the reordering happens once in this constructor so that callers which
/// reason about isotopologue *position* — nearest-neighbour boundaries, m/z spans, binary
/// search — are correct by construction.
/// </para>
/// <para>
/// <b>Index 0 is not necessarily the monoisotopic isotopologue.</b> Averagine builds its tables
/// with an absolute probability cutoff of 1e-8, and the all-light composition falls below that
/// somewhere around 31 kDa, above which the monoisotopic peak is simply absent from the table.
/// The envelope is still correctly <i>anchored</i> — the shift below is derived from
/// <see cref="Averagine.DiffToMonoisotopic"/>, not from array position — but code that needs the
/// monoisotopic m/z must compute it from the mass rather than read <c>CentroidMzs(z)[0]</c>.
/// </para>
/// <para>
/// <b>Not thread-safe.</b> Centroid positions and per-width-model Gaussian constants are cached
/// lazily into plain dictionaries on first use, so a kernel — or a <see cref="ForwardModel"/>,
/// which holds kernels as instance state — must not be shared across threads. Callers that
/// parallelise over proteoforms build their own kernels and are unaffected.
/// </para>
/// </remarks>
public sealed class IsotopeEnvelopeKernel
{
    private const double EvaluationWindowInSigmas = 6.0;

    /// <summary>
    /// The per-isotopologue Gaussian constants for one charge under one width model. Each
    /// isotopologue is a Gaussian of the width appropriate to <b>its own</b> centre, which is what
    /// keeps every isotopologue unit-area under an m/z-dependent σ.
    /// </summary>
    private readonly record struct ChargeWidths(double[] Norms, double[] Inv2Sigma2, double MaxSigma);

    private readonly double[] _neutralMasses;
    private readonly double[] _normalizedIntensities;
    private readonly Dictionary<int, double[]> _centroidCacheByCharge = new();

    // Keyed on the width model by value, so a caller that rebuilds an equal model per call still
    // hits the cache. One model is used for the length of a simulation in practice.
    private IPeakWidthModel? _cachedWidthModel;
    private Dictionary<int, ChargeWidths>? _widthCacheByCharge;

    public IsotopeEnvelopeKernel(double monoisotopicMass)
    {
        int idx = new Averagine().GetMostIntenseMassIndex(monoisotopicMass);
        // Both of these are ordered most-intense-first by Averagine's static constructor, and
        // they are the shared static tables — never sort them, only the local copies below.
        double[] theorMasses = Averagine.AllMasses[idx];
        double[] theorIntensities = Averagine.AllIntensities[idx];

        // DiffToMonoisotopic is (most-intense mass − monoisotopic mass), so this must read the
        // most-intense entry — theorMasses[0] in Averagine's ordering. Computed before the local
        // arrays are reordered below, and from Averagine's array rather than ours, so it is
        // unaffected by that reordering.
        double averagineMonoisotopic = theorMasses[0] - Averagine.DiffToMonoisotopic[idx];
        double shift = monoisotopicMass - averagineMonoisotopic;

        _neutralMasses = new double[theorMasses.Length];
        _normalizedIntensities = new double[theorIntensities.Length];

        double sum = 0;
        for (int i = 0; i < theorIntensities.Length; i++)
        {
            _neutralMasses[i] = theorMasses[i] + shift;
            sum += theorIntensities[i];
        }
        for (int i = 0; i < theorIntensities.Length; i++)
        {
            _normalizedIntensities[i] = theorIntensities[i] / sum;
        }

        // Reorder both arrays together into ascending mass, establishing the ordering contract
        // documented above. Paired in one sort so mass i and intensity i stay the same isotopologue.
        Array.Sort(_neutralMasses, _normalizedIntensities);
    }

    public int IsotopologueCount => _neutralMasses.Length;
    public double NeutralMass(int i) => _neutralMasses[i];
    public double Intensity(int i) => _normalizedIntensities[i];

    /// <summary>
    /// Evaluate the Gaussian-convolved envelope at a single m/z point for charge z, with every
    /// isotopologue given the width <paramref name="widthModel"/> assigns to its own centroid.
    /// </summary>
    public double Evaluate(double mz, int charge, IPeakWidthModel widthModel)
    {
        if (widthModel is null)
            throw new ArgumentNullException(nameof(widthModel));

        // Ascending m/z, so the isotopologues within the evaluation window are one contiguous run.
        double[] centroids = CachedCentroidMzs(charge);
        var widths = GetChargeWidths(charge, widthModel);

        // Sized from the widest peak in the envelope rather than the local one: a narrower bound
        // would drop tail contributions from the wide end of the span.
        double window = EvaluationWindowInSigmas * widths.MaxSigma;
        int start = LowerBound(centroids, mz - window);
        int endExclusive = UpperBound(centroids, mz + window);

        if (start >= endExclusive)
            return 0;

        double[] intensities = _normalizedIntensities;
        double[] norms = widths.Norms;
        double[] inv2Sigma2 = widths.Inv2Sigma2;

        double sum = 0;
        for (int i = start; i < endExclusive; i++)
        {
            double d = mz - centroids[i];
            sum += intensities[i] * norms[i] * Math.Exp(-d * d * inv2Sigma2[i]);
        }
        return sum;
    }

    /// <summary>
    /// Convenience overload for a width that does not vary with m/z. Allocates a
    /// <see cref="ConstantPeakWidth"/> per call, so prefer the model overload on hot paths.
    /// </summary>
    public double Evaluate(double mz, int charge, double sigmaMz) =>
        Evaluate(mz, charge, new ConstantPeakWidth(sigmaMz));

    private ChargeWidths GetChargeWidths(int charge, IPeakWidthModel widthModel)
    {
        if (_widthCacheByCharge is null || !widthModel.Equals(_cachedWidthModel))
        {
            _cachedWidthModel = widthModel;
            _widthCacheByCharge = new Dictionary<int, ChargeWidths>();
        }
        else if (_widthCacheByCharge.TryGetValue(charge, out var cached))
        {
            return cached;
        }

        double[] centroids = CachedCentroidMzs(charge);
        var norms = new double[centroids.Length];
        var inv2Sigma2 = new double[centroids.Length];
        double maxSigma = 0;

        for (int i = 0; i < centroids.Length; i++)
        {
            double sigma = widthModel.SigmaAt(centroids[i]);
            if (sigma <= 0 || !double.IsFinite(sigma))
                throw new ArgumentOutOfRangeException(nameof(widthModel),
                    $"σ_m must be finite and positive; got {sigma} at m/z {centroids[i]}.");

            norms[i] = 1.0 / (sigma * Math.Sqrt(2 * Math.PI));
            inv2Sigma2[i] = 1.0 / (2 * sigma * sigma);
            if (sigma > maxSigma)
                maxSigma = sigma;
        }

        var result = new ChargeWidths(norms, inv2Sigma2, maxSigma);
        _widthCacheByCharge[charge] = result;
        return result;
    }

    /// <summary>
    /// Centroid m/z positions for a given charge state (neutral mass → m/z), ascending.
    /// </summary>
    /// <remarks>
    /// A copy. <see cref="Evaluate(double, int, IPeakWidthModel)"/> binary-searches the cached array
    /// and indexes the per-charge width tables by position, so a caller that sorted or rescaled a
    /// returned array in place would silently corrupt this kernel's evaluation. Callers keep these
    /// arrays — <c>GroundTruthExtractor</c> stores them straight into a public
    /// <c>ProteoformGroundTruth.CentroidMzs</c> — so the boundary is worth the allocation.
    /// </remarks>
    public double[] CentroidMzs(int charge) => (double[])CachedCentroidMzs(charge).Clone();

    /// <summary>
    /// The cached array itself, for internal callers that only read it. Never hand this out.
    /// </summary>
    private double[] CachedCentroidMzs(int charge)
    {
        if (_centroidCacheByCharge.TryGetValue(charge, out var cached))
            return cached;

        var result = new double[_neutralMasses.Length];
        for (int i = 0; i < _neutralMasses.Length; i++)
            result[i] = _neutralMasses[i].ToMz(charge);

        _centroidCacheByCharge[charge] = result;
        return result;
    }

    public (double MinMz, double MaxMz) GetMzBounds(int charge)
    {
        var centroids = CachedCentroidMzs(charge);
        return centroids.Length == 0
            ? (double.PositiveInfinity, double.NegativeInfinity)
            : (centroids[0], centroids[^1]);
    }

    private static int LowerBound(double[] values, double target)
    {
        int lo = 0;
        int hi = values.Length;
        while (lo < hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            if (values[mid] < target)
                lo = mid + 1;
            else
                hi = mid;
        }

        return lo;
    }

    private static int UpperBound(double[] values, double target)
    {
        int lo = 0;
        int hi = values.Length;
        while (lo < hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            if (values[mid] <= target)
                lo = mid + 1;
            else
                hi = mid;
        }

        return lo;
    }
}
