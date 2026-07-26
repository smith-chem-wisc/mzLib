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
    int PeakCount,
    string? FeatureGroundTruthPath = null,
    int FeatureCount = 0);

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
        IPeakWidthModel widthModel,
        double[] scanTimes,
        bool isCentroid = false,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        var grid = _rasterizer.Rasterize(proteoforms, minCharge, maxCharge, widthModel, scanTimes, pointsPerSigma, mzPaddingInSigmas);
        var scans = _scanBuilder.BuildMs1Scans(grid, isCentroid: isCentroid);
        var file = _scanBuilder.BuildMsDataFile(scans);
        return new SimulationResult(grid, scans, file);
    }

    public SimulationResult Simulate(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        bool isCentroid = false,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0) =>
        Simulate(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes, isCentroid, pointsPerSigma, mzPaddingInSigmas);

    /// <summary>
    /// Simulates centroid spectra by evaluating the forward model at the theoretical isotopologue
    /// m/z of every charge state, with no intervening profile grid.
    /// </summary>
    public SimulationResult SimulateCentroided(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes)
    {
        var grid = _rasterizer.RasterizeAtCentroids(proteoforms, minCharge, maxCharge, widthModel, scanTimes);
        var scans = _scanBuilder.BuildMs1Scans(grid, isCentroid: true);
        var file = _scanBuilder.BuildMsDataFile(scans);
        return new SimulationResult(grid, scans, file);
    }

    public SimulationResult SimulateCentroided(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes) =>
        SimulateCentroided(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes);

    /// <summary>
    /// Simulates centroid spectra and writes them as mzML, together with two sidecars: the
    /// parameter vector θ_p that generated the run, and the feature-level truth a benchmark scores
    /// against. Both are suppressed by <paramref name="writeGroundTruthSidecar"/> = false.
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
        IPeakWidthModel widthModel,
        double[] scanTimes,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true)
    {
        if (string.IsNullOrWhiteSpace(outputPath))
            throw new ArgumentException("An output path is required.", nameof(outputPath));

        var simulation = SimulateCentroided(proteoforms, minCharge, maxCharge, widthModel, scanTimes);

        reduction ??= new ScanReductionOptions();
        double floor = SimulatedScanReducer.ComputeFloor(simulation.Scans, reduction);
        var scans = SimulatedScanReducer.Reduce(simulation.Scans, reduction);
        var file = _scanBuilder.BuildMsDataFile(scans);

        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, outputPath, writeIndexed);

        int peakCount = 0;
        foreach (var scan in scans)
            peakCount += scan.MassSpectrum.XArray.Length;

        string? groundTruthPath = null;
        string? featurePath = null;
        int featureCount = 0;
        if (writeGroundTruthSidecar)
        {
            groundTruthPath = Path.ChangeExtension(outputPath, ".groundtruth.tsv");
            WriteGroundTruth(proteoforms, minCharge, maxCharge, widthModel, groundTruthPath);

            var scanNumbers = new int[scans.Length];
            for (int s = 0; s < scans.Length; s++)
                scanNumbers[s] = scans[s].OneBasedScanNumber;

            var features = FeatureGroundTruth.Build(
                proteoforms, minCharge, maxCharge, widthModel,
                simulation.Grid.ScanTimes, scanNumbers, simulation.Grid.MzGrid, floor);

            featurePath = Path.ChangeExtension(outputPath, ".features.tsv");
            FeatureGroundTruth.Write(features, featurePath);
            featureCount = features.Length;
        }

        return new SimulationExportResult(outputPath, groundTruthPath, scans.Length, peakCount, featurePath, featureCount);
    }

    public SimulationExportResult WriteMzml(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true) =>
        WriteMzml(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes, outputPath,
            reduction, writeIndexed, writeGroundTruthSidecar);

    /// <summary>
    /// Writes the parameter vector for every simulated proteoform so a run can be reproduced.
    /// </summary>
    /// <remarks>
    /// This file is not scorable on its own — see <see cref="FeatureGroundTruth"/> for the
    /// feature-level truth. <c>SigmaMz</c> is the width at that proteoform's most abundant charge,
    /// which under a non-constant width law differs between rows; <c>PeakWidthModel</c> records the
    /// law itself.
    /// </remarks>
    public static void WriteGroundTruth(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        string outputPath)
    {
        var sb = new StringBuilder();
        sb.AppendLine(string.Join('\t',
            "Identifier", "MonoisotopicMass", "Abundance",
            "RtMu", "RtSigma", "RtTau",
            "ChargeMu", "ChargeSigma", "MinCharge", "MaxCharge", "SigmaMz", "PeakWidthModel"));

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
                Num(RepresentativeSigma(p, minCharge, maxCharge, widthModel)),
                widthModel.ToString()));
        }

        File.WriteAllText(outputPath, sb.ToString());
    }

    public static void WriteGroundTruth(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        string outputPath) =>
        WriteGroundTruth(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), outputPath);

    /// <summary>σ at the m/z of the proteoform's most abundant charge state.</summary>
    private static double RepresentativeSigma(
        ProteoformModel proteoform, int minCharge, int maxCharge, IPeakWidthModel widthModel)
    {
        int bestCharge = minCharge;
        double bestWeight = double.NegativeInfinity;
        for (int z = minCharge; z <= maxCharge; z++)
        {
            double fz = proteoform.ChargeDistribution.Evaluate(z);
            if (fz > bestWeight)
            {
                bestWeight = fz;
                bestCharge = z;
            }
        }

        var centroids = new IsotopeEnvelopeKernel(proteoform.MonoisotopicMass).CentroidMzs(bestCharge);
        return centroids.Length == 0 ? 0 : widthModel.SigmaAt(centroids[0]);
    }

    private static string Num(double value) => value.ToString("R", CultureInfo.InvariantCulture);
}
