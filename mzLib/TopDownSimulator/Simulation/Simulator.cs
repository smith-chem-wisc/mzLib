#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using MassSpectrometry;
using Readers;
using TopDownSimulator.Model;

namespace TopDownSimulator.Simulation;

public sealed record SimulationResult(RasterizedScanGrid Grid, MsDataScan[] Scans, GenericMsDataFile DataFile);

/// <summary>
/// Outcome of writing a simulation to disk. <see cref="PeakCount"/> is counted after reduction,
/// so it reflects what actually landed in the file.
/// </summary>
public sealed record SimulationExportResult(
    string MzmlPath,
    string? GroundTruthPath,
    int ScanCount,
    int PeakCount);

/// <summary>
/// High-level Phase 5 entrypoint for simulating synthetic MS1 data and writing it as mzML.
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

    /// <summary>
    /// Simulates centroid spectra by evaluating the forward model at the theoretical isotopologue
    /// m/z of every charge state, with no intervening profile grid.
    /// </summary>
    public SimulationResult SimulateCentroided(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes)
    {
        var grid = _rasterizer.RasterizeAtCentroids(proteoforms, minCharge, maxCharge, sigmaMz, scanTimes);
        var scans = _scanBuilder.BuildMs1Scans(grid, isCentroid: true);
        var file = _scanBuilder.BuildMsDataFile(scans);
        return new SimulationResult(grid, scans, file);
    }

    /// <summary>
    /// Simulates centroid spectra and writes them as mzML. A companion ground-truth TSV holding
    /// the θ_p that generated the file is written alongside it unless
    /// <paramref name="writeGroundTruthSidecar"/> is false.
    /// </summary>
    /// <remarks>
    /// Output is centroided, which is not merely a preference: mzLib's own mzML reader rejects
    /// profile-mode files, so a profile simulation could be written but never read back by mzLib
    /// or MetaMorpheus.
    /// </remarks>
    public SimulationExportResult WriteMzml(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true)
    {
        if (string.IsNullOrWhiteSpace(outputPath))
            throw new ArgumentException("An output path is required.", nameof(outputPath));

        var simulation = SimulateCentroided(proteoforms, minCharge, maxCharge, sigmaMz, scanTimes);

        reduction ??= new ScanReductionOptions();
        var scans = SimulatedScanReducer.Reduce(simulation.Scans, reduction);
        var file = _scanBuilder.BuildMsDataFile(scans);

        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, outputPath, writeIndexed);

        int peakCount = 0;
        foreach (var scan in scans)
            peakCount += scan.MassSpectrum.XArray.Length;

        string? groundTruthPath = null;
        if (writeGroundTruthSidecar)
        {
            groundTruthPath = Path.ChangeExtension(outputPath, ".groundtruth.tsv");
            WriteGroundTruth(proteoforms, minCharge, maxCharge, sigmaMz, groundTruthPath);
        }

        return new SimulationExportResult(outputPath, groundTruthPath, scans.Length, peakCount);
    }

    /// <summary>
    /// Writes the parameter vector for every simulated proteoform so a consumer of the mzML can
    /// score itself against the values that generated it.
    /// </summary>
    public static void WriteGroundTruth(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        string outputPath)
    {
        var sb = new StringBuilder();
        sb.AppendLine(string.Join('\t',
            "Identifier", "MonoisotopicMass", "Abundance",
            "RtMu", "RtSigma", "RtTau",
            "ChargeMu", "ChargeSigma", "MinCharge", "MaxCharge", "SigmaMz"));

        foreach (var p in proteoforms)
        {
            var gaussian = p.ChargeDistribution as GaussianChargeDistribution;

            sb.AppendLine(string.Join('\t',
                p.Identifier ?? string.Empty,
                Num(p.MonoisotopicMass),
                Num(p.Abundance),
                Num(p.RtProfile.Mu),
                Num(p.RtProfile.Sigma),
                Num(p.RtProfile.Tau),
                gaussian is null ? string.Empty : Num(gaussian.MuZ),
                gaussian is null ? string.Empty : Num(gaussian.SigmaZ),
                minCharge.ToString(CultureInfo.InvariantCulture),
                maxCharge.ToString(CultureInfo.InvariantCulture),
                Num(sigmaMz)));
        }

        File.WriteAllText(outputPath, sb.ToString());
    }

    private static string Num(double value) => value.ToString("R", CultureInfo.InvariantCulture);
}
