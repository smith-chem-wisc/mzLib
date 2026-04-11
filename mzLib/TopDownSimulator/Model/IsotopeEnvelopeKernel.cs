using System;
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
/// </remarks>
public sealed class IsotopeEnvelopeKernel
{
    private readonly double[] _neutralMasses;
    private readonly double[] _normalizedIntensities;

    public IsotopeEnvelopeKernel(double monoisotopicMass)
    {
        int idx = new Averagine().GetMostIntenseMassIndex(monoisotopicMass);
        double[] theorMasses = Averagine.AllMasses[idx];
        double[] theorIntensities = Averagine.AllIntensities[idx];

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
    }

    public int IsotopologueCount => _neutralMasses.Length;
    public double NeutralMass(int i) => _neutralMasses[i];
    public double Intensity(int i) => _normalizedIntensities[i];

    /// <summary>
    /// Evaluate the Gaussian-convolved envelope at a single m/z point for charge z.
    /// </summary>
    public double Evaluate(double mz, int charge, double sigmaMz)
    {
        double inv2Sigma2 = 1.0 / (2 * sigmaMz * sigmaMz);
        double norm = 1.0 / (sigmaMz * Math.Sqrt(2 * Math.PI));
        double sum = 0;
        for (int i = 0; i < _neutralMasses.Length; i++)
        {
            double center = _neutralMasses[i].ToMz(charge);
            double d = mz - center;
            sum += _normalizedIntensities[i] * norm * Math.Exp(-d * d * inv2Sigma2);
        }
        return sum;
    }

    /// <summary>
    /// Centroid m/z positions for a given charge state (neutral mass → m/z).
    /// </summary>
    public double[] CentroidMzs(int charge)
    {
        var result = new double[_neutralMasses.Length];
        for (int i = 0; i < _neutralMasses.Length; i++)
            result[i] = _neutralMasses[i].ToMz(charge);
        return result;
    }
}
