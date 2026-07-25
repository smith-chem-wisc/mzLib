#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using Chemistry;
using MassSpectrometry;
using Readers;
using TopDownSimulator.Simulation;

namespace TopDownSimulator.Ms2;

/// <summary>Outcome of simulating a set of MS2 acquisitions.</summary>
public sealed record Ms2SimulationResult(
    MsDataScan[] Scans,
    GenericMsDataFile DataFile,
    IReadOnlyList<Ms2FragmentationResult> Fragmentations);

/// <summary>Outcome of writing simulated MS2 scans to disk.</summary>
public sealed record Ms2ExportResult(
    string MzmlPath,
    string? GroundTruthPath,
    int ScanCount,
    int PeakCount);

/// <summary>
/// High-level entrypoint for simulating synthetic centroided MS2 spectra and writing them as mzML.
/// </summary>
/// <remarks>
/// Depends only on <see cref="IMs2FragmentationModel"/>, so the whole spectrum model is replaceable
/// by passing a different implementation to the constructor.
/// </remarks>
public sealed class Ms2Simulator
{
    private readonly Ms2ScanBuilder _scanBuilder;
    private readonly ScanBuilder _fileBuilder;

    public Ms2Simulator(
        IMs2FragmentationModel? fragmentationModel = null,
        Ms2ScanBuilder? scanBuilder = null,
        ScanBuilder? fileBuilder = null)
    {
        FragmentationModel = fragmentationModel ?? new PropensityCascadeFragmentationModel();
        _scanBuilder = scanBuilder ?? new Ms2ScanBuilder();
        _fileBuilder = fileBuilder ?? new ScanBuilder();
    }

    public IMs2FragmentationModel FragmentationModel { get; }

    /// <summary>
    /// Dissociation type stamped on the written scans. Read from the model's options when it exposes
    /// them, otherwise <see cref="DissociationType.ETD"/>.
    /// </summary>
    public DissociationType DissociationType => FragmentationModel is PropensityCascadeFragmentationModel cascade
        ? cascade.Options.DissociationType
        : DissociationType.ETD;

    /// <summary>Simulates one centroided MS2 scan.</summary>
    public Ms2SimulationResult SimulateCentroided(
        Ms2ScanRequest request,
        ScanReductionOptions? reduction = null,
        int firstOneBasedScanNumber = 1) =>
        SimulateCentroided(new[] { request }, reduction, firstOneBasedScanNumber);

    /// <summary>
    /// Simulates centroided MS2 scans, one per request, numbered consecutively from
    /// <paramref name="firstOneBasedScanNumber"/>.
    /// </summary>
    public Ms2SimulationResult SimulateCentroided(
        IReadOnlyList<Ms2ScanRequest> requests,
        ScanReductionOptions? reduction = null,
        int firstOneBasedScanNumber = 1)
    {
        if (requests is null)
            throw new ArgumentNullException(nameof(requests));
        if (requests.Count == 0)
            throw new ArgumentException("At least one MS2 request is required.", nameof(requests));
        if (firstOneBasedScanNumber < 1)
            throw new ArgumentOutOfRangeException(nameof(firstOneBasedScanNumber), "Scan numbers are one-based.");

        var fragmentations = new Ms2FragmentationResult[requests.Count];
        for (int i = 0; i < requests.Count; i++)
            fragmentations[i] = FragmentationModel.Fragment(requests[i].Precursor);

        double floor = Ms2ScanBuilder.ComputeIntensityFloor(fragmentations, reduction);
        var dissociationType = DissociationType;

        var scans = new MsDataScan[requests.Count];
        for (int i = 0; i < requests.Count; i++)
        {
            scans[i] = _scanBuilder.BuildMs2Scan(
                requests[i],
                fragmentations[i],
                dissociationType,
                oneBasedScanNumber: firstOneBasedScanNumber + i,
                oneBasedPrecursorScanNumber: null,
                intensityFloor: floor);
        }

        var file = _fileBuilder.BuildMsDataFile(scans);
        return new Ms2SimulationResult(scans, file, fragmentations);
    }

    /// <summary>
    /// Simulates centroided MS2 scans and writes them as mzML, with a companion ground-truth TSV
    /// listing every simulated fragment ion unless <paramref name="writeGroundTruthSidecar"/> is
    /// false.
    /// </summary>
    public Ms2ExportResult WriteMzml(
        IReadOnlyList<Ms2ScanRequest> requests,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true,
        int firstOneBasedScanNumber = 1)
    {
        if (string.IsNullOrWhiteSpace(outputPath))
            throw new ArgumentException("An output path is required.", nameof(outputPath));

        var simulation = SimulateCentroided(requests, reduction, firstOneBasedScanNumber);

        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(simulation.DataFile, outputPath, writeIndexed);

        int peakCount = 0;
        foreach (var scan in simulation.Scans)
            peakCount += scan.MassSpectrum.XArray.Length;

        string? groundTruthPath = null;
        if (writeGroundTruthSidecar)
        {
            groundTruthPath = Path.ChangeExtension(outputPath, ".ms2.groundtruth.tsv");
            WriteGroundTruth(requests, simulation, groundTruthPath);
        }

        return new Ms2ExportResult(outputPath, groundTruthPath, simulation.Scans.Length, peakCount);
    }

    /// <summary>
    /// Writes one row per simulated fragment ion, so a consumer of the mzML can score itself against
    /// the assignment that generated it. Scan-level columns are repeated on every row.
    /// </summary>
    public static void WriteGroundTruth(
        IReadOnlyList<Ms2ScanRequest> requests,
        Ms2SimulationResult simulation,
        string outputPath)
    {
        var sb = new StringBuilder();
        sb.AppendLine(string.Join('\t',
            "ScanNumber", "RetentionTime", "BaseSequence", "FullSequence",
            "PrecursorCharge", "PrecursorMz", "StartingIonCount", "SurvivingPrecursorIntensity",
            "Annotation", "ProductType", "Terminus", "OneBasedBondNumber", "FragmentNumber",
            "FragmentNeutralMass", "FragmentCharge", "FragmentMz", "IsotopologueIndex", "Intensity"));

        for (int i = 0; i < requests.Count; i++)
        {
            var request = requests[i];
            var fragmentation = simulation.Fragmentations[i];
            var scan = simulation.Scans[i];
            double precursorMz = request.Precursor.MonoisotopicMass.ToMz(request.PrecursorCharge);

            foreach (var fragment in fragmentation.Fragments)
            {
                sb.AppendLine(string.Join('\t',
                    scan.OneBasedScanNumber.ToString(CultureInfo.InvariantCulture),
                    Num(request.RetentionTime),
                    fragmentation.BaseSequence,
                    request.Precursor.FullSequence ?? string.Empty,
                    request.PrecursorCharge.ToString(CultureInfo.InvariantCulture),
                    Num(precursorMz),
                    Num(fragmentation.StartingIonCount),
                    Num(fragmentation.SurvivingPrecursorIntensity),
                    fragment.Product.Annotation,
                    fragment.Product.ProductType.ToString(),
                    fragment.Product.Terminus.ToString(),
                    fragment.OneBasedBondNumber.ToString(CultureInfo.InvariantCulture),
                    fragment.Product.FragmentNumber.ToString(CultureInfo.InvariantCulture),
                    Num(fragment.Product.NeutralMass),
                    fragment.Charge.ToString(CultureInfo.InvariantCulture),
                    Num(fragment.Mz),
                    fragment.IsotopologueIndex.ToString(CultureInfo.InvariantCulture),
                    Num(fragment.Intensity)));
            }
        }

        File.WriteAllText(outputPath, sb.ToString());
    }

    private static string Num(double value) => value.ToString("R", CultureInfo.InvariantCulture);
}
