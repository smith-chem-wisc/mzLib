using System;
using System.Collections.Generic;
using System.Linq;
using TopDownSimulator.Model;

namespace TopDownSimulator.Simulation;

public sealed record RasterizedScanGrid(double[] ScanTimes, double[] MzGrid, double[,] Intensities);

/// <summary>
/// Builds an m/z grid for a set of proteoforms and evaluates the forward model on it.
/// </summary>
public sealed class GridRasterizer
{
    public RasterizedScanGrid Rasterize(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        var mzGrid = BuildMzGrid(proteoforms, minCharge, maxCharge, sigmaMz, pointsPerSigma, mzPaddingInSigmas);
        var model = new ForwardModel(proteoforms, minCharge, maxCharge, sigmaMz);
        return new RasterizedScanGrid((double[])scanTimes.Clone(), mzGrid, model.Rasterize(scanTimes, mzGrid));
    }

    /// <summary>
    /// Evaluates the forward model directly at the theoretical isotopologue m/z of every charge
    /// state of every proteoform, producing centroid spectra without ever building a profile grid.
    /// </summary>
    /// <remarks>
    /// This is both cheaper and more accurate than rasterizing a dense profile grid and then
    /// peak-picking it: there is no grid spacing to quantize peak positions, and the m/z axis is
    /// the number of isotopologues rather than the number of profile bins. Intensities still come
    /// from <see cref="ForwardModel.Rasterize"/>, so overlap between neighbouring isotopologues
    /// and between co-eluting proteoforms is accounted for exactly as in the profile path.
    /// </remarks>
    public RasterizedScanGrid RasterizeAtCentroids(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes)
    {
        var mzAxis = BuildCentroidMzAxis(proteoforms, minCharge, maxCharge, sigmaMz);
        var model = new ForwardModel(proteoforms, minCharge, maxCharge, sigmaMz);
        return new RasterizedScanGrid((double[])scanTimes.Clone(), mzAxis, model.Rasterize(scanTimes, mzAxis));
    }

    /// <summary>
    /// The sorted union of every isotopologue m/z across the supplied proteoforms and charges.
    /// Positions closer together than a small fraction of σ_m are collapsed, since they are one
    /// peak as far as the instrument is concerned and the forward model sums into either of them
    /// identically.
    /// </summary>
    public double[] BuildCentroidMzAxis(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz)
    {
        if (proteoforms.Count == 0)
            throw new ArgumentException("At least one proteoform is required.", nameof(proteoforms));
        if (sigmaMz <= 0)
            throw new ArgumentOutOfRangeException(nameof(sigmaMz));
        if (minCharge < 1 || maxCharge < minCharge)
            throw new ArgumentException("Charge range must satisfy 1 ≤ minCharge ≤ maxCharge.");

        var all = new List<double>();
        foreach (var proteoform in proteoforms)
        {
            var kernel = new IsotopeEnvelopeKernel(proteoform.MonoisotopicMass);
            for (int z = minCharge; z <= maxCharge; z++)
                all.AddRange(kernel.CentroidMzs(z));
        }

        all.Sort();

        double mergeTolerance = sigmaMz / 100.0;
        var axis = new List<double>(all.Count);
        foreach (double mz in all)
        {
            if (axis.Count == 0 || mz - axis[^1] > mergeTolerance)
                axis.Add(mz);
        }

        return axis.ToArray();
    }

    public double[] BuildMzGrid(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        if (proteoforms.Count == 0)
            throw new ArgumentException("At least one proteoform is required.", nameof(proteoforms));
        if (sigmaMz <= 0)
            throw new ArgumentOutOfRangeException(nameof(sigmaMz));
        if (pointsPerSigma < 1)
            throw new ArgumentOutOfRangeException(nameof(pointsPerSigma));

        double minMz = double.PositiveInfinity;
        double maxMz = double.NegativeInfinity;
        foreach (var proteoform in proteoforms)
        {
            var kernel = new IsotopeEnvelopeKernel(proteoform.MonoisotopicMass);
            for (int z = minCharge; z <= maxCharge; z++)
            {
                var centroids = kernel.CentroidMzs(z);
                if (centroids.Length == 0)
                    continue;

                minMz = Math.Min(minMz, centroids[0]);
                maxMz = Math.Max(maxMz, centroids[^1]);
            }
        }

        if (!double.IsFinite(minMz) || !double.IsFinite(maxMz))
            throw new InvalidOperationException("Could not determine an m/z range for the supplied proteoforms.");

        double padding = mzPaddingInSigmas * sigmaMz;
        double start = minMz - padding;
        double end = maxMz + padding;
        double step = sigmaMz / pointsPerSigma;
        int pointCount = Math.Max(2, (int)Math.Ceiling((end - start) / step) + 1);

        var mzGrid = new double[pointCount];
        for (int i = 0; i < pointCount; i++)
            mzGrid[i] = start + i * step;

        return mzGrid;
    }
}
