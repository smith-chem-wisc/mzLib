#nullable enable
using System;
using System.Collections.Generic;
using TopDownSimulator.Extraction;

namespace TopDownSimulator.Fitting;

public enum WidthFitMode
{
    /// <summary>Width was recovered from a multi-point profile peak shape.</summary>
    Profile,
    /// <summary>Input was centroided (≤1 peak per isotopologue window); returned the fallback.</summary>
    CentroidedFallback,
}

/// <summary>
/// One σ_m estimate, from one isotopologue's peak shape at one (charge, scan) cell.
/// </summary>
/// <param name="Mz">Theoretical centroid the shape was measured around.</param>
/// <param name="SigmaMz">Second central moment of the peak, as a standard deviation.</param>
/// <param name="Intensity">Summed intensity in the window — how bright this peak was.</param>
/// <param name="PeakCount">Profile samples the moment was computed from.</param>
public readonly record struct PeakWidthMeasurement(double Mz, double SigmaMz, double Intensity, int PeakCount);

public sealed record FittedEnvelopeWidth(
    double SigmaMz,
    WidthFitMode Mode,
    int PeaksUsed,
    IReadOnlyList<PeakWidthMeasurement>? Measurements = null,
    int WindowsTooCoarselySampled = 0)
{
    /// <summary>
    /// The individual (m/z, σ) observations behind <see cref="SigmaMz"/>. Pooling these across
    /// every record is what makes an m/z-dependent width law fittable; a per-record scalar rests on
    /// a few dozen peaks and is dominated by estimator variance.
    /// </summary>
    public IReadOnlyList<PeakWidthMeasurement> Measurements { get; init; } =
        Measurements ?? Array.Empty<PeakWidthMeasurement>();
}

/// <summary>
/// Estimates per-peak m/z Gaussian widths σ_m for a proteoform from the peak shapes at each
/// charge's apex scan. When the input is profile data, σ_m is recovered as the weighted second
/// central moment of the peak samples around each theoretical centroid. When the input is
/// centroided (only one peak falls in each window), a fallback σ_m is returned.
/// </summary>
public sealed class EnvelopeWidthFitter
{
    /// <summary>
    /// Windows with fewer samples than this are discarded. A second moment from two points carries
    /// essentially no information about width and would only add scatter to a pooled fit.
    /// </summary>
    public const int DefaultMinimumPeaksPerWindow = 3;

    /// <summary>
    /// How finely a window must be sampled, in samples per σ, for its second moment to be a peak
    /// width at all.
    /// </summary>
    /// <remarks>
    /// This is the guard against reading a width off centroided input. A few isolated centroids
    /// inside an isotopologue's boundaries still have a second moment, but that number describes
    /// how far apart the centroids are, not how wide a peak is — and because the boundaries are
    /// half the isotopologue spacing, 0.5/z, while m/z ≈ M/z, that spread grows in direct
    /// proportion to m/z. Pooled across records it therefore fits σ ∝ (m/z)^1 with high apparent
    /// precision, which is an artifact of the extraction window rather than a property of the
    /// instrument. A genuinely sampled peak has many points across one σ; scattered centroids are
    /// spaced comparably to the σ they appear to imply.
    /// <para>
    /// The threshold is 3 because the two populations separate there. For k roughly evenly spread
    /// centroids the implied σ is about gap·√((k²−1)/12), so samples-per-σ is ~1.2 at k = 3 and
    /// falls as k grows — always around 1 or below. Profile data from an Orbitrap carries several
    /// points across one σ by construction. Setting the bar too low does not merely admit noise: it
    /// admits the <i>wide</i> windows preferentially, because the test scales with σ, and that
    /// selection bias inflates the fitted width rather than just scattering it.
    /// </para>
    /// </remarks>
    public const double DefaultMinimumSamplesPerSigma = 3.0;

    private readonly double _fallbackSigmaMz;
    private readonly int _minPeaksPerWindow;
    private readonly double _minSamplesPerSigma;

    public EnvelopeWidthFitter(
        double fallbackSigmaMz = 0.01,
        int minPeaksPerWindow = DefaultMinimumPeaksPerWindow,
        double minSamplesPerSigma = DefaultMinimumSamplesPerSigma)
    {
        if (fallbackSigmaMz <= 0) throw new ArgumentOutOfRangeException(nameof(fallbackSigmaMz));
        if (minPeaksPerWindow < 2) throw new ArgumentOutOfRangeException(nameof(minPeaksPerWindow));
        if (minSamplesPerSigma < 0) throw new ArgumentOutOfRangeException(nameof(minSamplesPerSigma));
        _fallbackSigmaMz = fallbackSigmaMz;
        _minPeaksPerWindow = minPeaksPerWindow;
        _minSamplesPerSigma = minSamplesPerSigma;
    }

    public FittedEnvelopeWidth Fit(ProteoformGroundTruth truth)
    {
        var measurements = new List<PeakWidthMeasurement>();
        var gaps = new List<double>();
        int peaksUsed = 0;
        int windowsTooCoarselySampled = 0;

        for (int c = 0; c < truth.ChargeCount; c++)
        {
            int apexScan = ApexScanForCharge(truth, c);
            if (apexScan < 0)
                continue;

            int nIso = truth.CentroidMzs[c].Length;
            for (int i = 0; i < nIso; i++)
            {
                double centroid = truth.CentroidMzs[c][i];

                // Midpoint boundaries keep a neighbouring isotopologue that shares the extraction
                // window out of this peak's moment. They are only the true neighbours because
                // CentroidMzs is ordered by ascending m/z.
                double leftBoundary = i == 0
                    ? double.NegativeInfinity
                    : 0.5 * (truth.CentroidMzs[c][i - 1] + centroid);
                double rightBoundary = i == nIso - 1
                    ? double.PositiveInfinity
                    : 0.5 * (centroid + truth.CentroidMzs[c][i + 1]);

                var samples = truth.IsotopologuePeakWindows[c][i][apexScan];
                double sw = 0, swx = 0, swxx = 0;
                int n = 0;
                double previousMz = double.NaN;
                for (int k = 0; k < samples.Length; k++)
                {
                    var p = samples[k];
                    if (p.Intensity <= 0) continue;
                    if (p.Mz < leftBoundary || p.Mz > rightBoundary) continue;
                    sw += p.Intensity;
                    swx += p.Intensity * p.Mz;
                    swxx += p.Intensity * p.Mz * p.Mz;
                    if (n > 0)
                        gaps.Add(Math.Abs(p.Mz - previousMz));
                    previousMz = p.Mz;
                    n++;
                }

                if (n < _minPeaksPerWindow || sw <= 0)
                {
                    gaps.Clear();
                    continue;
                }

                double mean = swx / sw;
                double variance = swxx / sw - mean * mean;
                if (variance <= 0)
                {
                    gaps.Clear();
                    continue;
                }

                double sigma = Math.Sqrt(variance);
                double medianGap = MedianOf(gaps);
                gaps.Clear();

                // Reject anything too coarsely sampled to be a peak shape. See
                // DefaultMinimumSamplesPerSigma: a scatter of centroids has a second moment too,
                // and pooling those produces a confident fit to the wrong thing.
                if (_minSamplesPerSigma > 0 && medianGap * _minSamplesPerSigma > sigma)
                {
                    windowsTooCoarselySampled++;
                    continue;
                }

                measurements.Add(new PeakWidthMeasurement(centroid, sigma, sw, n));
                peaksUsed += n;
            }
        }

        if (measurements.Count == 0)
            return new FittedEnvelopeWidth(
                _fallbackSigmaMz, WidthFitMode.CentroidedFallback, 0,
                WindowsTooCoarselySampled: windowsTooCoarselySampled);

        double totalWeightedVar = 0;
        double totalWeight = 0;
        foreach (var m in measurements)
        {
            totalWeightedVar += m.Intensity * m.SigmaMz * m.SigmaMz;
            totalWeight += m.Intensity;
        }

        double pooled = totalWeight > 0 ? Math.Sqrt(totalWeightedVar / totalWeight) : 0;
        return pooled > 0
            ? new FittedEnvelopeWidth(pooled, WidthFitMode.Profile, peaksUsed, measurements, windowsTooCoarselySampled)
            : new FittedEnvelopeWidth(_fallbackSigmaMz, WidthFitMode.CentroidedFallback, peaksUsed, measurements, windowsTooCoarselySampled);
    }

    private static double MedianOf(List<double> values)
    {
        if (values.Count == 0)
            return 0;

        values.Sort();
        int n = values.Count;
        return n % 2 == 1 ? values[n / 2] : 0.5 * (values[n / 2 - 1] + values[n / 2]);
    }

    /// <summary>
    /// The scan whose peak windows the extractor populated for this charge, falling back to a
    /// recomputed XIC apex for truths built by hand.
    /// </summary>
    private static int ApexScanForCharge(ProteoformGroundTruth truth, int chargeIndex)
    {
        if (truth.ApexScanIndexByCharge.Length == truth.ChargeCount)
            return truth.ApexScanIndexByCharge[chargeIndex];

        int best = -1;
        double bestValue = 0;
        for (int s = 0; s < truth.ScanCount; s++)
        {
            if (truth.ChargeXics[chargeIndex][s] > bestValue)
            {
                bestValue = truth.ChargeXics[chargeIndex][s];
                best = s;
            }
        }

        return best;
    }
}
