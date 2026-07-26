#nullable enable
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
    private const double EvaluationWindowInSigmas = 6.0;
    private const double MinimumContributionWeight = 1e-18;

    private readonly IReadOnlyList<(ProteoformModel Model, IsotopeEnvelopeKernel Kernel)> _species;
    private readonly int _minCharge;
    private readonly int _maxCharge;
    private readonly IPeakWidthModel _widthModel;

    public ForwardModel(IEnumerable<ProteoformModel> proteoforms, int minCharge, int maxCharge, IPeakWidthModel widthModel)
    {
        if (minCharge < 1 || maxCharge < minCharge)
            throw new ArgumentException("Charge range must satisfy 1 ≤ minCharge ≤ maxCharge.");

        _species = proteoforms
            .Select(p => (p, new IsotopeEnvelopeKernel(p.MonoisotopicMass)))
            .ToList();
        _minCharge = minCharge;
        _maxCharge = maxCharge;
        _widthModel = widthModel ?? throw new ArgumentNullException(nameof(widthModel));
    }

    public ForwardModel(IEnumerable<ProteoformModel> proteoforms, int minCharge, int maxCharge, double sigmaMz)
        : this(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz))
    {
    }

    /// <summary>The width law every peak in this model is rendered under.</summary>
    public IPeakWidthModel WidthModel => _widthModel;

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
                chargeSum += fz * kernel.Evaluate(mz, z, _widthModel);
            }
            total += model.Abundance * rt * chargeSum;
        }
        return total;
    }

    /// <summary>
    /// Rasterize the forward model onto a rectangular grid.
    /// Returns a [scan, mz] intensity array.
    /// </summary>
    /// <remarks>
    /// The envelope shape does not depend on scan time, so the kernel is evaluated once per
    /// (species, charge) across that envelope's m/z span and then scaled per scan, rather than
    /// re-evaluated for every scan. For any given output cell the (species, charge) terms are still
    /// accumulated in the same order as a scan-outer loop would, so the result is unchanged.
    /// </remarks>
    public double[,] Rasterize(double[] scanTimes, double[] mzGrid)
    {
        var result = new double[scanTimes.Length, mzGrid.Length];

        var rtScaledAbundances = new double[scanTimes.Length];
        double[] kernelValues = Array.Empty<double>();

        foreach (var (model, kernel) in _species)
        {
            bool anyScanContributes = false;
            for (int s = 0; s < scanTimes.Length; s++)
            {
                double rt = model.RtProfile.Evaluate(scanTimes[s]);
                double scaled = rt > 0 ? model.Abundance * rt : 0;
                rtScaledAbundances[s] = scaled;
                anyScanContributes |= scaled > MinimumContributionWeight;
            }

            if (!anyScanContributes)
                continue;

            for (int z = _minCharge; z <= _maxCharge; z++)
            {
                double fz = model.ChargeDistribution.Evaluate(z);
                if (fz <= 0)
                    continue;

                var (minMz, maxMz) = kernel.GetMzBounds(z);
                if (!double.IsFinite(minMz) || !double.IsFinite(maxMz))
                    continue;

                // Padding is per (species, charge) rather than hoisted, because under an
                // m/z-dependent width law the envelope at charge 6 and the one at charge 40 sit at
                // very different m/z and need very different tails. Taking the widest σ over the
                // envelope's own span keeps the window conservative at both ends.
                double mzPadding = EvaluationWindowInSigmas * _widthModel.MaxSigmaOver(minMz, maxMz);
                int start = LowerBound(mzGrid, minMz - mzPadding);
                int endExclusive = UpperBound(mzGrid, maxMz + mzPadding);
                int span = endExclusive - start;
                if (span <= 0)
                    continue;

                if (kernelValues.Length < span)
                    kernelValues = new double[span];

                for (int b = 0; b < span; b++)
                    kernelValues[b] = kernel.Evaluate(mzGrid[start + b], z, _widthModel);

                for (int s = 0; s < scanTimes.Length; s++)
                {
                    double weight = rtScaledAbundances[s] * fz;
                    if (rtScaledAbundances[s] <= MinimumContributionWeight || weight <= MinimumContributionWeight)
                        continue;

                    for (int b = 0; b < span; b++)
                        result[s, start + b] += weight * kernelValues[b];
                }
            }
        }
        return result;
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
