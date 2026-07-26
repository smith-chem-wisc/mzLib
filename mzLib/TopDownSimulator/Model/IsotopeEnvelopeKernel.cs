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
/// for every charge: index 0 is the lowest-mass (monoisotopic) isotopologue, not the most
/// intense one. Note that this is the opposite of <see cref="Averagine"/>'s own tables, which
/// are stored most-intense-first; the reordering happens once in this constructor so that
/// callers which reason about isotopologue *position* — nearest-neighbour boundaries, m/z
/// spans, binary search — are correct by construction.
/// </para>
/// </remarks>
public sealed class IsotopeEnvelopeKernel
{
    private const double EvaluationWindowInSigmas = 6.0;

    private readonly double[] _neutralMasses;
    private readonly double[] _normalizedIntensities;
    private readonly Dictionary<int, double[]> _centroidCacheByCharge = new();

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
    /// Evaluate the Gaussian-convolved envelope at a single m/z point for charge z.
    /// </summary>
    public double Evaluate(double mz, int charge, double sigmaMz)
    {
        if (sigmaMz <= 0)
            throw new ArgumentOutOfRangeException(nameof(sigmaMz));

        double inv2Sigma2 = 1.0 / (2 * sigmaMz * sigmaMz);
        double norm = 1.0 / (sigmaMz * Math.Sqrt(2 * Math.PI));
        double window = EvaluationWindowInSigmas * sigmaMz;

        // Ascending m/z, so the isotopologues within the evaluation window are one contiguous run.
        double[] centroids = CentroidMzs(charge);
        double[] intensities = _normalizedIntensities;
        int start = LowerBound(centroids, mz - window);
        int endExclusive = UpperBound(centroids, mz + window);

        if (start >= endExclusive)
            return 0;

        double sum = 0;
        for (int i = start; i < endExclusive; i++)
        {
            double center = centroids[i];
            double d = mz - center;
            sum += intensities[i] * norm * Math.Exp(-d * d * inv2Sigma2);
        }
        return sum;
    }

    /// <summary>
    /// Centroid m/z positions for a given charge state (neutral mass → m/z).
    /// </summary>
    public double[] CentroidMzs(int charge)
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
        var centroids = CentroidMzs(charge);
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
