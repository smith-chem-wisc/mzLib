#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

/// <param name="ReferenceMz">m/z at which the fitted model is clamped to a plausible σ.</param>
/// <param name="OutlierMadMultiple">
/// Measurements further than this many robust deviations from the fixed-slope law are dropped, so
/// one bad record cannot distort k. Non-positive disables the pass.
/// </param>
/// <param name="MinimumMzRatio">
/// Required ratio of the highest to the lowest measured m/z. A power law fitted over a narrow band
/// is an extrapolation dressed as a fit: the slope is barely constrained, so k inherits whatever
/// the scatter happened to do. Measured on a centroided run, the surviving windows spanned m/z
/// 601–725 — a ratio of 1.2 — and produced a k roughly 4× an Orbitrap's true width.
/// </param>
public sealed record PeakWidthModelFitOptions(
    double ReferenceMz = 800.0,
    double MinSigmaAtReferenceMz = 1e-4,
    double MaxSigmaAtReferenceMz = 0.2,
    double OutlierMadMultiple = 3.0,
    int MinimumMeasurements = 30,
    double MinimumMzRatio = 1.5);

/// <summary>
/// The outcome of fitting σ_m = k · (m/z)^1.5 to pooled width measurements.
/// </summary>
/// <param name="Model">The shipped model: slope pinned at 1.5, only k fitted.</param>
/// <param name="FreeSlope">
/// Slope of an unconstrained regression of log σ on log(m/z), reported rather than assumed away. If
/// it is not near 1.5 the data is saying something about the instrument or about the estimator, and
/// that is worth seeing before trusting <paramref name="Model"/>.
/// </param>
/// <param name="FreeSlopeStandardError">
/// Standard error of <paramref name="FreeSlope"/>. Without it the slope cannot be judged: an
/// estimate of 1.33 means something very different depending on whether it is ±0.03 or ±0.2, and
/// the m/z span the measurements cover (<paramref name="MinMz"/>..<paramref name="MaxMz"/>) is what
/// mostly decides which.
/// </param>
/// <param name="MeasurementsUsed">Measurements surviving the outlier pass.</param>
/// <param name="Clamped">Whether <paramref name="Model"/> hit the plausibility bounds.</param>
public sealed record PeakWidthModelFit(
    OrbitrapPeakWidth Model,
    double K,
    double FreeSlope,
    double FreeIntercept,
    double FreeSlopeStandardError,
    double MedianSigmaMz,
    double MedianMz,
    double MinMz,
    double MaxMz,
    int MeasurementCount,
    int MeasurementsUsed,
    bool Clamped)
{
    /// <summary>How far the free slope sits from the fixed exponent, in slope units.</summary>
    public double SlopeDeviation => Math.Abs(FreeSlope - OrbitrapPeakWidth.Exponent);

    /// <summary>
    /// How many standard errors the free slope sits from the fixed exponent.
    /// </summary>
    public double SlopeDeviationInStandardErrors =>
        FreeSlopeStandardError > 0
            ? SlopeDeviation / FreeSlopeStandardError
            : double.NaN;

    /// <summary>
    /// Whether the data actively contradicts σ ∝ (m/z)^1.5, as opposed to merely failing to confirm
    /// it precisely. Requires the slope to be off by a materially wrong amount <i>and</i> for that
    /// to be more than <paramref name="maxStandardErrors"/> of estimator noise.
    /// </summary>
    /// <remarks>
    /// The absolute term is not redundant. On data that lies exactly on the law the residuals are
    /// rounding error, so the standard error collapses toward zero and the ratio blows up — a
    /// perfect fit would be rejected by a standard-error test alone. The absolute term is what
    /// keeps the decision about whether the exponent is wrong enough to matter: at a slope of 1.25
    /// σ is off by 17% at m/z 1500, at 1.0 by 37%.
    /// </remarks>
    public bool ContradictsWidthLaw(double maxStandardErrors = 3.0, double maxSlopeDeviation = 0.25) =>
        SlopeDeviation > maxSlopeDeviation
        && SlopeDeviationInStandardErrors > maxStandardErrors;
}

/// <summary>
/// Fits the one free parameter of <see cref="OrbitrapPeakWidth"/> by pooled weighted regression
/// over width measurements from every record at once.
/// </summary>
/// <remarks>
/// Pooling rather than taking a median of per-record scalars is the point. A per-record σ rests on
/// the few dozen peaks of one charge state, and across a real run those scalars scatter by 3× with
/// no clean m/z trend — that spread is estimator variance, not physics, and a median of it throws
/// away the m/z dependence entirely.
/// <para>
/// Measurements are weighted by (samples − 1), which is what sets the precision of a variance
/// estimate. Weighting by intensity instead would let a handful of very bright peaks decide k.
/// </para>
/// </remarks>
public sealed class PeakWidthModelFitter
{
    private readonly PeakWidthModelFitOptions _options;

    public PeakWidthModelFitter(PeakWidthModelFitOptions? options = null)
    {
        _options = options ?? new PeakWidthModelFitOptions();
        if (_options.ReferenceMz <= 0)
            throw new ArgumentOutOfRangeException(nameof(options), "ReferenceMz must be positive.");
        if (_options.MinSigmaAtReferenceMz <= 0 || _options.MaxSigmaAtReferenceMz <= _options.MinSigmaAtReferenceMz)
            throw new ArgumentOutOfRangeException(nameof(options), "Sigma clamp bounds must satisfy 0 < min < max.");
    }

    /// <summary>
    /// Returns null when there are too few usable measurements to fit anything — the caller should
    /// fall back to a constant width rather than ship a one-point extrapolation of a power law.
    /// </summary>
    public PeakWidthModelFit? Fit(IEnumerable<PeakWidthMeasurement> measurements)
    {
        if (measurements is null) throw new ArgumentNullException(nameof(measurements));

        var xs = new List<double>();
        var ys = new List<double>();
        var ws = new List<double>();
        var sigmas = new List<double>();
        var mzs = new List<double>();

        foreach (var m in measurements)
        {
            if (!(m.Mz > 0) || !(m.SigmaMz > 0) || !double.IsFinite(m.Mz) || !double.IsFinite(m.SigmaMz))
                continue;

            xs.Add(Math.Log(m.Mz));
            ys.Add(Math.Log(m.SigmaMz));
            ws.Add(Math.Max(1.0, m.PeakCount - 1.0));
            sigmas.Add(m.SigmaMz);
            mzs.Add(m.Mz);
        }

        int total = xs.Count;
        if (total < _options.MinimumMeasurements)
            return null;

        // Residual of each measurement against the fixed-slope law. log σ − 1.5 log(m/z) = log k,
        // so these are direct estimates of log k and their spread is the scatter in k itself.
        var logK = new double[total];
        for (int i = 0; i < total; i++)
            logK[i] = ys[i] - OrbitrapPeakWidth.Exponent * xs[i];

        var keep = new bool[total];
        int used = 0;
        if (_options.OutlierMadMultiple > 0)
        {
            double median = Median((double[])logK.Clone());
            var deviations = new double[total];
            for (int i = 0; i < total; i++)
                deviations[i] = Math.Abs(logK[i] - median);

            // 1.4826 rescales the MAD to a standard deviation for normally distributed residuals.
            double mad = 1.4826 * Median(deviations);
            double cutoff = _options.OutlierMadMultiple * mad;

            for (int i = 0; i < total; i++)
            {
                keep[i] = mad <= 0 || Math.Abs(logK[i] - median) <= cutoff;
                if (keep[i]) used++;
            }
        }
        else
        {
            for (int i = 0; i < total; i++) keep[i] = true;
            used = total;
        }

        // Too few survivors is a refusal, not a licence to readmit the measurements the outlier
        // pass just identified as bad. Readmitting them would let precisely the worst data set k,
        // and "35/35 used" in the log would be indistinguishable from "no outliers found".
        if (used < _options.MinimumMeasurements)
            return null;

        // Evaluated over the surviving measurements, because those are the ones the fit actually
        // rests on. Checking the span before outlier rejection would let a handful of far-out
        // measurements satisfy the guard and then be discarded, leaving the narrow-band
        // extrapolation this exists to refuse.
        double minMz = double.PositiveInfinity, maxMz = double.NegativeInfinity;
        for (int i = 0; i < total; i++)
        {
            if (!keep[i]) continue;
            minMz = Math.Min(minMz, mzs[i]);
            maxMz = Math.Max(maxMz, mzs[i]);
        }

        if (minMz <= 0 || maxMz / minMz < _options.MinimumMzRatio)
            return null;

        // Fixed slope: the weighted mean of the per-measurement log k estimates.
        double sumW = 0, sumWLogK = 0;
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < total; i++)
        {
            if (!keep[i]) continue;
            double w = ws[i];
            sumW += w;
            sumWLogK += w * logK[i];
            sx += w * xs[i];
            sy += w * ys[i];
            sxx += w * xs[i] * xs[i];
            sxy += w * xs[i] * ys[i];
        }

        if (sumW <= 0)
            return null;

        double k = Math.Exp(sumWLogK / sumW);

        // Free slope and intercept, for diagnosis only.
        double denominator = sumW * sxx - sx * sx;
        double freeSlope = double.NaN, freeIntercept = double.NaN, slopeStandardError = double.NaN;
        if (Math.Abs(denominator) > 0)
        {
            freeSlope = (sumW * sxy - sx * sy) / denominator;
            freeIntercept = (sy - freeSlope * sx) / sumW;

            double weightedSquaredResiduals = 0;
            for (int i = 0; i < total; i++)
            {
                if (!keep[i]) continue;
                double r = ys[i] - (freeIntercept + freeSlope * xs[i]);
                weightedSquaredResiduals += ws[i] * r * r;
            }

            // Degrees of freedom is the number of MEASUREMENTS minus two, not the sum of the
            // weights. A window contributes one σ observation however many profile points it was
            // sampled at; treating its (samples − 1) weight as a replication count would shrink the
            // standard error by about √(samples), which is enough to make any decision built on it
            // depend on the detector's sampling rate rather than on the data.
            if (used > 2)
                slopeStandardError = Math.Sqrt(weightedSquaredResiduals / (used - 2) * (sumW / denominator));
        }

        bool clamped = false;
        double sigmaAtReference = k * Math.Pow(_options.ReferenceMz, OrbitrapPeakWidth.Exponent);
        if (sigmaAtReference < _options.MinSigmaAtReferenceMz)
        {
            k = _options.MinSigmaAtReferenceMz / Math.Pow(_options.ReferenceMz, OrbitrapPeakWidth.Exponent);
            clamped = true;
        }
        else if (sigmaAtReference > _options.MaxSigmaAtReferenceMz)
        {
            k = _options.MaxSigmaAtReferenceMz / Math.Pow(_options.ReferenceMz, OrbitrapPeakWidth.Exponent);
            clamped = true;
        }

        return new PeakWidthModelFit(
            Model: new OrbitrapPeakWidth(k),
            K: k,
            FreeSlope: freeSlope,
            FreeIntercept: freeIntercept,
            FreeSlopeStandardError: slopeStandardError,
            MedianSigmaMz: Median(sigmas.ToArray()),
            MedianMz: Median(mzs.ToArray()),
            MinMz: minMz,
            MaxMz: maxMz,
            MeasurementCount: total,
            MeasurementsUsed: used,
            Clamped: clamped);
    }

    private static double Median(double[] values)
    {
        Array.Sort(values);
        int n = values.Length;
        if (n == 0) return 0;
        return n % 2 == 1
            ? values[n / 2]
            : 0.5 * (values[n / 2 - 1] + values[n / 2]);
    }
}
