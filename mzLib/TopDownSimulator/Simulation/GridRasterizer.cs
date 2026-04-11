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
