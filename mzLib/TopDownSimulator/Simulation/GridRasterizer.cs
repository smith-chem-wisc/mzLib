#nullable enable
using System;
using System.Collections.Generic;
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
        IPeakWidthModel widthModel,
        double[] scanTimes,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        var mzGrid = BuildMzGrid(proteoforms, minCharge, maxCharge, widthModel, pointsPerSigma, mzPaddingInSigmas);
        var model = new ForwardModel(proteoforms, minCharge, maxCharge, widthModel);
        return new RasterizedScanGrid((double[])scanTimes.Clone(), mzGrid, model.Rasterize(scanTimes, mzGrid));
    }

    public RasterizedScanGrid Rasterize(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0) =>
        Rasterize(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes, pointsPerSigma, mzPaddingInSigmas);

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
        IPeakWidthModel widthModel,
        double[] scanTimes)
    {
        var mzAxis = BuildCentroidMzAxis(proteoforms, minCharge, maxCharge, widthModel);
        var model = new ForwardModel(proteoforms, minCharge, maxCharge, widthModel);
        return new RasterizedScanGrid((double[])scanTimes.Clone(), mzAxis, model.Rasterize(scanTimes, mzAxis));
    }

    public RasterizedScanGrid RasterizeAtCentroids(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes) =>
        RasterizeAtCentroids(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes);

    /// <summary>
    /// The sorted union of every isotopologue m/z across the supplied proteoforms and charges.
    /// Positions closer together than a small fraction of the local σ_m are collapsed, since they
    /// are one peak as far as the instrument is concerned and the forward model sums into either of
    /// them identically.
    /// </summary>
    public double[] BuildCentroidMzAxis(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel)
    {
        if (proteoforms.Count == 0)
            throw new ArgumentException("At least one proteoform is required.", nameof(proteoforms));
        if (widthModel is null)
            throw new ArgumentNullException(nameof(widthModel));
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

        var axis = new List<double>(all.Count);
        foreach (double mz in all)
        {
            // Evaluated at the candidate position rather than once for the whole axis: under a
            // width law the merge tolerance at m/z 1500 is several times the one at m/z 600.
            double mergeTolerance = widthModel.SigmaAt(mz) / 100.0;
            if (axis.Count == 0 || mz - axis[^1] > mergeTolerance)
                axis.Add(mz);
        }

        return axis.ToArray();
    }

    public double[] BuildCentroidMzAxis(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz) =>
        BuildCentroidMzAxis(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz));

    /// <summary>
    /// A profile m/z grid whose spacing tracks the local σ_m, so every peak is sampled
    /// <paramref name="pointsPerSigma"/> times regardless of where in the range it lands.
    /// </summary>
    public double[] BuildMzGrid(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        if (proteoforms.Count == 0)
            throw new ArgumentException("At least one proteoform is required.", nameof(proteoforms));
        if (widthModel is null)
            throw new ArgumentNullException(nameof(widthModel));
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

        double start = minMz - mzPaddingInSigmas * widthModel.SigmaAt(minMz);
        double end = maxMz + mzPaddingInSigmas * widthModel.SigmaAt(maxMz);

        var mzGrid = new List<double>();
        for (double mz = start; mz <= end; )
        {
            mzGrid.Add(mz);
            double step = widthModel.SigmaAt(mz) / pointsPerSigma;
            if (!(step > 0) || !double.IsFinite(step))
                throw new InvalidOperationException($"Width model produced a non-positive σ_m at m/z {mz}.");

            mz += step;
        }

        // The accumulating loop stops at the last position not past `end`, which would leave the
        // top of the 6σ padding uncovered. Under a constant width the old closed-form grid always
        // spanned it, so extend rather than truncate.
        if (mzGrid.Count < 2 || mzGrid[^1] < end)
            mzGrid.Add(end);

        return mzGrid.ToArray();
    }

    public double[] BuildMzGrid(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0) =>
        BuildMzGrid(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), pointsPerSigma, mzPaddingInSigmas);
}
