using System;
using System.Collections.Generic;
using System.Linq;

namespace TopDownSimulator.Model;

/// <summary>
/// Evaluates the top-down MS1 forward model
/// I(t, mz) = Σ_p A_p · g_p(t) · Σ_z f_p(z) · φ(b; M_p, z, σ_m)
/// for a set of proteoforms on a caller-supplied (scanTime × m/z) grid.
/// </summary>
public sealed class ForwardModel
{
    private readonly IReadOnlyList<(ProteoformModel Model, IsotopeEnvelopeKernel Kernel)> _species;
    private readonly int _minCharge;
    private readonly int _maxCharge;
    private readonly double _sigmaMz;

    public ForwardModel(IEnumerable<ProteoformModel> proteoforms, int minCharge, int maxCharge, double sigmaMz)
    {
        if (minCharge < 1 || maxCharge < minCharge)
            throw new ArgumentException("Charge range must satisfy 1 ≤ minCharge ≤ maxCharge.");

        _species = proteoforms
            .Select(p => (p, new IsotopeEnvelopeKernel(p.MonoisotopicMass)))
            .ToList();
        _minCharge = minCharge;
        _maxCharge = maxCharge;
        _sigmaMz = sigmaMz;
    }

    /// <summary>
    /// Evaluate the forward model at a single (t, mz) point.
    /// </summary>
    public double Evaluate(double t, double mz)
    {
        double total = 0;
        foreach (var (model, kernel) in _species)
        {
            double rt = model.RtProfile.Evaluate(t);
            if (rt <= 0) continue;

            double chargeSum = 0;
            for (int z = _minCharge; z <= _maxCharge; z++)
            {
                double fz = model.ChargeDistribution.Evaluate(z);
                if (fz <= 0) continue;
                chargeSum += fz * kernel.Evaluate(mz, z, _sigmaMz);
            }
            total += model.Abundance * rt * chargeSum;
        }
        return total;
    }

    /// <summary>
    /// Rasterize the forward model onto a rectangular grid.
    /// Returns a [scan, mz] intensity array.
    /// </summary>
    public double[,] Rasterize(double[] scanTimes, double[] mzGrid)
    {
        var result = new double[scanTimes.Length, mzGrid.Length];
        for (int s = 0; s < scanTimes.Length; s++)
        {
            double t = scanTimes[s];
            for (int b = 0; b < mzGrid.Length; b++)
            {
                result[s, b] = Evaluate(t, mzGrid[b]);
            }
        }
        return result;
    }
}
