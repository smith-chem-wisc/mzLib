#nullable enable
using System;
using System.Collections.Generic;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using TopDownSimulator.Simulation;

namespace TopDownSimulator.Ms2;

/// <summary>
/// One MS2 acquisition to simulate.
/// </summary>
/// <param name="Precursor">The proteoform being fragmented.</param>
/// <param name="PrecursorCharge">Charge of the isolated precursor.</param>
/// <param name="RetentionTime">Retention time to stamp on the scan, in minutes.</param>
/// <param name="PrecursorIntensity">
/// MS1 intensity of the isolated precursor, written as the selected-ion peak intensity. When null
/// the model's starting ion count is used, which is the same unit the fragment intensities are in.
/// It is never written as absent: mzLib's mzML writer allocates a fixed three-element cvParam array
/// for the selected ion and leaves the peak-intensity slot null if this is missing.
/// </param>
/// <param name="IsolationWidth">Isolation window width in m/z. Defaults to 4.</param>
public sealed record Ms2ScanRequest(
    IBioPolymerWithSetMods Precursor,
    int PrecursorCharge,
    double RetentionTime,
    double? PrecursorIntensity = null,
    double IsolationWidth = 4.0);

/// <summary>
/// Converts simulated fragment ions into centroided MS2 <see cref="MsDataScan"/> objects.
/// </summary>
/// <remarks>
/// <para>
/// Output is centroid-only, matching the MS1 path: mzLib's own mzML reader throws on profile-mode
/// files, so a profile MS2 simulation could be written but never read back.
/// </para>
/// <para>
/// Thresholding uses <see cref="ScanReductionOptions"/> and
/// <see cref="SimulatedScanReducer.GlobalMaxIntensity"/>, but is applied to the peak lists
/// <em>before</em> the scans are constructed rather than by calling
/// <see cref="SimulatedScanReducer.Reduce"/> afterwards. <c>Reduce</c> rebuilds each scan through
/// the precursor-free <see cref="MsDataScan"/> overload, which would silently drop the selected-ion
/// m/z, charge and dissociation type that an MSn scan needs — and that the mzML writer dereferences
/// unconditionally.
/// </para>
/// </remarks>
public sealed class Ms2ScanBuilder
{
    /// <summary>
    /// Peaks closer together than this (in m/z) are merged into one centroid, intensity-summed.
    /// Distinct fragments essentially never collide; complementary series occasionally produce
    /// numerically identical m/z, and two peaks at the same position are not a thing a centroided
    /// spectrum can represent.
    /// </summary>
    public const double CentroidMergeTolerance = 1e-9;

    public MsDataScan BuildMs2Scan(
        Ms2ScanRequest request,
        Ms2FragmentationResult fragmentation,
        DissociationType dissociationType,
        int oneBasedScanNumber,
        int? oneBasedPrecursorScanNumber = null,
        double intensityFloor = 0,
        Polarity polarity = Polarity.Positive,
        MZAnalyzerType analyzer = MZAnalyzerType.Orbitrap,
        string scanFilter = "synthetic ms2",
        MzRange? scanWindowRange = null)
    {
        if (request is null)
            throw new ArgumentNullException(nameof(request));
        if (fragmentation is null)
            throw new ArgumentNullException(nameof(fragmentation));
        if (request.PrecursorCharge < 1)
            throw new ArgumentOutOfRangeException(nameof(request), "Precursor charge must be >= 1.");

        var (mz, intensities) = BuildPeakList(fragmentation.Fragments, intensityFloor);

        double tic = 0;
        for (int i = 0; i < intensities.Length; i++)
            tic += intensities[i];

        double precursorMz = request.Precursor.MonoisotopicMass.ToMz(request.PrecursorCharge);
        var window = scanWindowRange ?? DeriveScanWindow(mz, precursorMz);

        return new MsDataScan(
            massSpectrum: new MzSpectrum(mz, intensities, false),
            oneBasedScanNumber: oneBasedScanNumber,
            msnOrder: 2,
            isCentroid: true,
            polarity: polarity,
            retentionTime: request.RetentionTime,
            scanWindowRange: window,
            scanFilter: scanFilter,
            mzAnalyzer: analyzer,
            totalIonCurrent: tic,
            injectionTime: 1.0,
            noiseData: null,
            nativeId: $"scan={oneBasedScanNumber}",
            selectedIonMz: precursorMz,
            selectedIonChargeStateGuess: request.PrecursorCharge,
            selectedIonIntensity: request.PrecursorIntensity ?? fragmentation.StartingIonCount,
            isolationMZ: precursorMz,
            isolationWidth: request.IsolationWidth,
            dissociationType: dissociationType,
            oneBasedPrecursorScanNumber: oneBasedPrecursorScanNumber,
            selectedIonMonoisotopicGuessMz: precursorMz);
    }

    /// <summary>
    /// Absolute intensity floor implied by <paramref name="options"/> across a whole set of
    /// fragmentation results, so a faint spectrum is not given its own private floor.
    /// </summary>
    public static double ComputeIntensityFloor(
        IReadOnlyList<Ms2FragmentationResult> fragmentations,
        ScanReductionOptions? options)
    {
        options ??= new ScanReductionOptions();
        if (options.RelativeIntensityThreshold < 0)
            throw new ArgumentOutOfRangeException(nameof(options), "Relative intensity threshold must be non-negative.");
        if (options.MinimumIntensity < 0)
            throw new ArgumentOutOfRangeException(nameof(options), "Minimum intensity must be non-negative.");

        double max = 0;
        foreach (var fragmentation in fragmentations)
        {
            foreach (var fragment in fragmentation.Fragments)
            {
                if (fragment.Intensity > max)
                    max = fragment.Intensity;
            }
        }

        return Math.Max(options.MinimumIntensity, max * options.RelativeIntensityThreshold);
    }

    /// <summary>
    /// Sorts fragments by m/z, drops anything at or below <paramref name="intensityFloor"/>, and
    /// merges coincident peaks.
    /// </summary>
    public static (double[] Mz, double[] Intensities) BuildPeakList(
        IReadOnlyList<SimulatedFragmentIon> fragments, double intensityFloor)
    {
        var kept = new List<SimulatedFragmentIon>(fragments.Count);
        foreach (var fragment in fragments)
        {
            if (fragment.Intensity <= 0 || fragment.Intensity < intensityFloor)
                continue;

            kept.Add(fragment);
        }

        kept.Sort(static (a, b) => a.Mz.CompareTo(b.Mz));

        var mz = new List<double>(kept.Count);
        var intensities = new List<double>(kept.Count);
        foreach (var fragment in kept)
        {
            if (mz.Count > 0 && fragment.Mz - mz[^1] <= CentroidMergeTolerance)
            {
                intensities[^1] += fragment.Intensity;
                continue;
            }

            mz.Add(fragment.Mz);
            intensities.Add(fragment.Intensity);
        }

        return (mz.ToArray(), intensities.ToArray());
    }

    private static MzRange DeriveScanWindow(double[] mz, double precursorMz)
    {
        if (mz.Length == 0)
            return new MzRange(Math.Max(precursorMz - 1, 0), precursorMz + 1);

        double min = Math.Min(mz[0], precursorMz);
        double max = Math.Max(mz[^1], precursorMz);
        return new MzRange(Math.Max(min - 1, 0), max + 1);
    }
}
