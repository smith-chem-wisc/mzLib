using System.Collections.Generic;
using MassSpectrometry;
using Readers;
using TopDownSimulator.Model;

namespace TopDownSimulator.Simulation;

public sealed record SimulationResult(RasterizedScanGrid Grid, MsDataScan[] Scans, GenericMsDataFile DataFile);

/// <summary>
/// High-level Phase 5 entrypoint for simulating synthetic MS1 data.
/// </summary>
public sealed class Simulator
{
    private readonly GridRasterizer _rasterizer;
    private readonly ScanBuilder _scanBuilder;

    public Simulator(GridRasterizer? rasterizer = null, ScanBuilder? scanBuilder = null)
    {
        _rasterizer = rasterizer ?? new GridRasterizer();
        _scanBuilder = scanBuilder ?? new ScanBuilder();
    }

    public SimulationResult Simulate(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        bool isCentroid = false,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        var grid = _rasterizer.Rasterize(proteoforms, minCharge, maxCharge, sigmaMz, scanTimes, pointsPerSigma, mzPaddingInSigmas);
        var scans = _scanBuilder.BuildMs1Scans(grid, isCentroid: isCentroid);
        var file = _scanBuilder.BuildMsDataFile(scans);
        return new SimulationResult(grid, scans, file);
    }

    public void WriteMzml(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        string outputPath,
        bool isCentroid = false,
        bool writeIndexed = true,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        var result = Simulate(proteoforms, minCharge, maxCharge, sigmaMz, scanTimes, isCentroid, pointsPerSigma, mzPaddingInSigmas);
        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(result.DataFile, outputPath, writeIndexed);
    }
}
