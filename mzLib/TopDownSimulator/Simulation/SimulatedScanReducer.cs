#nullable enable
using System;
using System.Collections.Generic;
using MassSpectrometry;
using MzLibUtil;

namespace TopDownSimulator.Simulation;

/// <summary>
/// Controls how a simulated centroid spectrum is thinned before being written to file.
/// </summary>
public sealed record ScanReductionOptions
{
    /// <summary>
    /// Discard peaks below this fraction of the brightest peak in the entire simulation.
    /// The reference is global rather than per-scan so that a near-empty scan does not have
    /// its own faint signal promoted above the floor.
    /// </summary>
    public double RelativeIntensityThreshold { get; init; } = 1e-4;

    /// <summary>
    /// Absolute intensity floor, applied together with <see cref="RelativeIntensityThreshold"/>.
    /// The effective floor is the larger of the two.
    /// </summary>
    public double MinimumIntensity { get; init; }
}

/// <summary>
/// Drops negligible peaks from simulated centroid scans.
/// </summary>
/// <remarks>
/// The forward model is evaluated at every isotopologue of every charge state of every
/// proteoform, so each scan initially carries the union of all of them. Away from a
/// proteoform's elution apex nearly all of those peaks evaluate to effectively zero, and real
/// instruments do not record them. Thresholding is what keeps a full-run simulation to a
/// realistic peak count.
/// </remarks>
public static class SimulatedScanReducer
{
    public static MsDataScan[] Reduce(MsDataScan[] scans, ScanReductionOptions? options = null)
    {
        if (scans is null)
            throw new ArgumentNullException(nameof(scans));

        options ??= new ScanReductionOptions();
        if (options.RelativeIntensityThreshold < 0)
            throw new ArgumentOutOfRangeException(nameof(options), "Relative intensity threshold must be non-negative.");
        if (options.MinimumIntensity < 0)
            throw new ArgumentOutOfRangeException(nameof(options), "Minimum intensity must be non-negative.");

        double floor = Math.Max(options.MinimumIntensity, GlobalMaxIntensity(scans) * options.RelativeIntensityThreshold);

        var reduced = new MsDataScan[scans.Length];
        for (int s = 0; s < scans.Length; s++)
        {
            var scan = scans[s];
            var (mz, intensities) = Threshold(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, floor);
            reduced[s] = CloneWithSpectrum(scan, mz, intensities);
        }

        return reduced;
    }

    public static double GlobalMaxIntensity(MsDataScan[] scans)
    {
        double max = 0;
        foreach (var scan in scans)
        {
            var y = scan.MassSpectrum.YArray;
            for (int i = 0; i < y.Length; i++)
            {
                if (y[i] > max)
                    max = y[i];
            }
        }

        return max;
    }

    private static (double[] Mz, double[] Intensities) Threshold(double[] x, double[] y, double floor)
    {
        var mz = new List<double>();
        var intensities = new List<double>();

        for (int i = 0; i < y.Length; i++)
        {
            if (y[i] < floor)
                continue;

            mz.Add(x[i]);
            intensities.Add(y[i]);
        }

        return (mz.ToArray(), intensities.ToArray());
    }

    private static MsDataScan CloneWithSpectrum(MsDataScan source, double[] mz, double[] intensities)
    {
        double tic = 0;
        for (int i = 0; i < intensities.Length; i++)
            tic += intensities[i];

        var scanWindowRange = mz.Length > 0
            ? new MzRange(mz[0], mz[^1])
            : source.ScanWindowRange;

        return new MsDataScan(
            massSpectrum: new MzSpectrum(mz, intensities, false),
            oneBasedScanNumber: source.OneBasedScanNumber,
            msnOrder: source.MsnOrder,
            isCentroid: true,
            polarity: source.Polarity,
            retentionTime: source.RetentionTime,
            scanWindowRange: scanWindowRange,
            scanFilter: source.ScanFilter,
            mzAnalyzer: source.MzAnalyzer,
            totalIonCurrent: tic,
            injectionTime: source.InjectionTime,
            noiseData: null,
            nativeId: source.NativeId);
    }
}
