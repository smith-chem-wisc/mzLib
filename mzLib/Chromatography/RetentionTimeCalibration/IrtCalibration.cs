using MathNet.Numerics.Statistics;
using MzLibUtil;

namespace Chromatography.RetentionTimeCalibration;

/// <summary>The shape of a run-minutes to iRT calibration.</summary>
public enum IrtCalibrationKind
{
    /// <summary>A straight line, fitted robustly.</summary>
    Linear,

    /// <summary>A robust locally weighted regression (Cleveland 1979), for runs whose gradient bends.</summary>
    Lowess,
}

/// <param name="Kind">Straight line or LOWESS curve.</param>
/// <param name="MinimumAnchors">Fewer anchors than this are refused rather than fitted.</param>
/// <param name="Bandwidth">LOWESS only: the fraction of anchors in each local fit, in (0, 1].</param>
/// <param name="Knots">LOWESS only: how many points along the run the curve is evaluated at.</param>
/// <param name="RobustnessIterations">Rounds of downweighting anchors with large residuals.</param>
public sealed record IrtCalibrationOptions(
    IrtCalibrationKind Kind = IrtCalibrationKind.Lowess,
    int MinimumAnchors = 20,
    double Bandwidth = 0.3,
    int Knots = 100,
    int RobustnessIterations = 3);

/// <summary>
/// Fits one run's retention times onto a library's iRT scale from anchors: the observed apex RT, in minutes, of
/// confidently identified precursors paired with their library iRT. Anchors from a first pass include wrong
/// identifications, so the fit is robust. It starts from per-bin medians and then applies Tukey bisquare reweighting.
/// </summary>
public static class IrtCalibration
{
    private const int StartingBins = 20;

    /// <exception cref="ArgumentNullException"><paramref name="anchors"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">An option is out of range.</exception>
    /// <exception cref="ArgumentException">
    /// Too few anchors, a non-finite value, anchors that span no time, or a fit that does not increase with time.
    /// </exception>
    public static IrtCalibrationModel Fit(IReadOnlyList<(RtMinutes Rt, Irt Irt)> anchors, IrtCalibrationOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(anchors);
        options ??= new IrtCalibrationOptions();
        ValidateOptions(options);
        if (anchors.Count < options.MinimumAnchors)
            throw new ArgumentException($"Calibration needs at least {options.MinimumAnchors} anchors; got {anchors.Count}.", nameof(anchors));
        if (anchors.Any(a => !double.IsFinite(a.Rt.Value) || !double.IsFinite(a.Irt.Value)))
            throw new ArgumentException("Every anchor's retention time and iRT must be finite.", nameof(anchors));

        var sorted = anchors.OrderBy(a => a.Rt.Value).ThenBy(a => a.Irt.Value).ToArray();
        double[] x = sorted.Select(a => a.Rt.Value).ToArray();
        double[] y = sorted.Select(a => a.Irt.Value).ToArray();
        if (x[^1] <= x[0])
            throw new ArgumentException("The anchors span no time, so no calibration can be fitted.", nameof(anchors));

        // Robust start: residuals from per-bin medians decide the first weights
        double[] weights = BisquareWeights(Residuals(x, y, BinnedMedianCurve(x, y)));

        IrtCalibrationModel model = options.Kind == IrtCalibrationKind.Linear
            ? FitLinear(x, y, weights, options.RobustnessIterations)
            : FitLowess(x, y, weights, options);

        return model with { ResidualSd = InlierSd(Residuals(x, y, model)), AnchorCount = anchors.Count };
    }

    private static IrtCalibrationModel FitLinear(double[] x, double[] y, double[] weights, int iterations)
    {
        var (intercept, slope) = WeightedLine(x, y, weights, 0, x.Length);
        for (int iteration = 0; iteration < iterations; iteration++)
        {
            var line = LineModel(x, intercept, slope);
            weights = BisquareWeights(Residuals(x, y, line));
            (intercept, slope) = WeightedLine(x, y, weights, 0, x.Length);
        }
        return LineModel(x, intercept, slope);
    }

    private static IrtCalibrationModel FitLowess(double[] x, double[] y, double[] weights, IrtCalibrationOptions options)
    {
        int neighbours = Math.Clamp((int)Math.Ceiling(options.Bandwidth * x.Length), 3, x.Length);
        double[] knotRt = Enumerable.Range(0, options.Knots)
            .Select(k => x[0] + (x[^1] - x[0]) * k / (options.Knots - 1)).ToArray();

        IrtCalibrationModel model = LocalFits(x, y, weights, knotRt, neighbours);
        for (int iteration = 0; iteration < options.RobustnessIterations; iteration++)
        {
            weights = BisquareWeights(Residuals(x, y, model));
            model = LocalFits(x, y, weights, knotRt, neighbours);
        }
        return model;
    }

    /// <summary>A local weighted line at each knot, with tricube distance weights, made strictly increasing.</summary>
    private static IrtCalibrationModel LocalFits(double[] x, double[] y, double[] robustness, double[] knotRt, int neighbours)
    {
        var knotIrt = new double[knotRt.Length];
        double firstSlope = 0, lastSlope = 0;
        var local = new double[x.Length];
        for (int k = 0; k < knotRt.Length; k++)
        {
            // The neighbours nearest this knot form a contiguous run of the sorted anchors
            int lo = LowerBound(x, knotRt[k]), hi = lo;
            while (hi - lo < neighbours)
            {
                if (lo > 0 && (hi >= x.Length || knotRt[k] - x[lo - 1] <= x[hi] - knotRt[k])) lo--;
                else hi++;
            }
            double reach = Math.Max(knotRt[k] - x[lo], x[hi - 1] - knotRt[k]) * 1.0001 + 1e-12;
            for (int i = lo; i < hi; i++)
            {
                double u = Math.Abs(x[i] - knotRt[k]) / reach;
                double tricube = Math.Pow(1 - u * u * u, 3);
                local[i] = tricube * robustness[i];
            }
            var (intercept, slope) = WeightedLine(x, y, local, lo, hi);
            knotIrt[k] = intercept + slope * knotRt[k];
            if (k == 0) firstSlope = slope;
            if (k == knotRt.Length - 1) lastSlope = slope;
            Array.Clear(local, lo, hi - lo);
        }

        MakeStrictlyIncreasing(knotIrt);
        double fallback = (knotIrt[^1] - knotIrt[0]) / (knotRt[^1] - knotRt[0]);
        return new IrtCalibrationModel(IrtCalibrationKind.Lowess, knotRt, knotIrt,
            firstSlope > 0 ? firstSlope : fallback, lastSlope > 0 ? lastSlope : fallback);
    }

    private static IrtCalibrationModel LineModel(double[] x, double intercept, double slope)
    {
        if (!(slope > 0))
            throw new ArgumentException("The anchors' iRT does not increase with retention time, so the map cannot be inverted.");
        double[] knotRt = [x[0], x[^1]];
        return new IrtCalibrationModel(IrtCalibrationKind.Linear, knotRt, knotRt.Select(t => intercept + slope * t).ToArray(), slope, slope);
    }

    /// <summary>Weighted least-squares line over anchors [from, to).</summary>
    private static (double Intercept, double Slope) WeightedLine(double[] x, double[] y, double[] w, int from, int to)
    {
        double sw = 0, sx = 0, sy = 0;
        for (int i = from; i < to; i++) { sw += w[i]; sx += w[i] * x[i]; sy += w[i] * y[i]; }
        if (sw <= 0)
            return (y.Skip(from).Take(to - from).Average(), 0);
        double mx = sx / sw, my = sy / sw, sxx = 0, sxy = 0;
        for (int i = from; i < to; i++) { sxx += w[i] * (x[i] - mx) * (x[i] - mx); sxy += w[i] * (x[i] - mx) * (y[i] - my); }
        double slope = sxx > 0 ? sxy / sxx : 0;
        return (my - slope * mx, slope);
    }

    /// <summary>Medians of iRT within equal-count RT bins, as a piecewise-linear curve: a start that ignores wrong anchors.</summary>
    private static IrtCalibrationModel BinnedMedianCurve(double[] x, double[] y)
    {
        int bins = Math.Min(StartingBins, x.Length / 3);
        var rt = new List<double>();
        var irt = new List<double>();
        for (int b = 0; b < bins; b++)
        {
            int from = b * x.Length / bins, to = (b + 1) * x.Length / bins;
            double binRt = Statistics.Median(x[from..to]);
            if (rt.Count > 0 && binRt <= rt[^1])
                continue;
            rt.Add(binRt);
            irt.Add(Statistics.Median(y[from..to]));
        }
        if (rt.Count < 2)
        {
            rt = [x[0], x[^1]];
            irt = [Statistics.Median(y), Statistics.Median(y) + 1e-9];
        }
        double[] knotIrt = irt.ToArray();
        MakeStrictlyIncreasing(knotIrt);
        double slope = (knotIrt[^1] - knotIrt[0]) / (rt[^1] - rt[0]);
        return new IrtCalibrationModel(IrtCalibrationKind.Linear, rt.ToArray(), knotIrt, slope, slope);
    }

    private static double[] Residuals(double[] x, double[] y, IrtCalibrationModel model) =>
        x.Select((rt, i) => y[i] - model.ToIrt(new RtMinutes(rt)).Value).ToArray();

    /// <summary>Tukey bisquare weights with scale 6 × the median absolute residual (Cleveland's choice for LOWESS).</summary>
    private static double[] BisquareWeights(double[] residuals)
    {
        double scale = 6 * Statistics.Median(residuals.Select(Math.Abs));
        if (!(scale > 0))
            return residuals.Select(r => r == 0 ? 1.0 : 0.0).ToArray();
        return residuals.Select(r => { double u = r / scale; return Math.Abs(u) < 1 ? (1 - u * u) * (1 - u * u) : 0; }).ToArray();
    }

    /// <summary>
    /// The residual SD among anchors that fit: residuals beyond 3 robust SDs are set aside and the SD re-estimated, twice.
    /// A plain MAD would be inflated by wrong identifications.
    /// </summary>
    private static double InlierSd(double[] residuals)
    {
        double sd = 1.4826 * Statistics.Median(residuals.Select(Math.Abs));
        for (int round = 0; round < 2 && sd > 0; round++)
        {
            var inliers = residuals.Where(r => Math.Abs(r) <= 3 * sd).ToArray();
            if (inliers.Length < 2)
                break;
            sd = Math.Sqrt(inliers.Average(r => r * r));
        }
        return sd;
    }

    /// <summary>Isotonic (pool-adjacent-violators) then a tiny tilt on any flat run, so the map is invertible.</summary>
    private static void MakeStrictlyIncreasing(double[] values)
    {
        var blocks = new List<(double Mean, int Count)>();
        foreach (double v in values)
        {
            blocks.Add((v, 1));
            while (blocks.Count > 1 && blocks[^2].Mean >= blocks[^1].Mean)
            {
                var (m2, c2) = blocks[^1];
                var (m1, c1) = blocks[^2];
                blocks.RemoveRange(blocks.Count - 2, 2);
                blocks.Add(((m1 * c1 + m2 * c2) / (c1 + c2), c1 + c2));
            }
        }
        int i = 0;
        foreach (var (mean, count) in blocks)
            for (int j = 0; j < count; j++)
                values[i++] = mean;
        double span = Math.Max(1e-9, Math.Abs(values[^1] - values[0]));
        for (int k = 1; k < values.Length; k++)
            if (values[k] <= values[k - 1])
                values[k] = values[k - 1] + span * 1e-9;
    }

    private static int LowerBound(double[] sorted, double value)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi) { int mid = (lo + hi) >>> 1; if (sorted[mid] < value) lo = mid + 1; else hi = mid; }
        return lo;
    }

    private static void ValidateOptions(IrtCalibrationOptions options)
    {
        if (options.MinimumAnchors < 3)
            throw new ArgumentOutOfRangeException(nameof(options), options.MinimumAnchors, "MinimumAnchors must be at least 3.");
        if (!(options.Bandwidth > 0 && options.Bandwidth <= 1))
            throw new ArgumentOutOfRangeException(nameof(options), options.Bandwidth, "Bandwidth must be in (0, 1].");
        if (options.Knots < 3)
            throw new ArgumentOutOfRangeException(nameof(options), options.Knots, "Knots must be at least 3.");
        if (options.RobustnessIterations < 0)
            throw new ArgumentOutOfRangeException(nameof(options), options.RobustnessIterations, "RobustnessIterations cannot be negative.");
    }
}

/// <summary>
/// One run's calibrated map between retention time in minutes and library iRT. It is strictly increasing, piecewise
/// linear through its knots, and continues linearly beyond the anchors (see <see cref="IsExtrapolated"/>).
/// </summary>
public sealed record IrtCalibrationModel
{
    private readonly double[] _knotRt;
    private readonly double[] _knotIrt;
    private readonly double _slopeBefore;
    private readonly double _slopeAfter;

    internal IrtCalibrationModel(IrtCalibrationKind kind, double[] knotRt, double[] knotIrt, double slopeBefore, double slopeAfter)
    {
        Kind = kind;
        _knotRt = knotRt;
        _knotIrt = knotIrt;
        _slopeBefore = slopeBefore;
        _slopeAfter = slopeAfter;
    }

    public IrtCalibrationKind Kind { get; }

    /// <summary>Robust residual SD of the anchors about the fit, in iRT units (wrong identifications set aside).</summary>
    public double ResidualSd { get; init; }

    public int AnchorCount { get; init; }

    /// <summary>The earliest anchor's retention time; the map is fitted from here.</summary>
    public RtMinutes FirstAnchor => new(_knotRt[0]);

    /// <summary>The latest anchor's retention time; the map is fitted up to here.</summary>
    public RtMinutes LastAnchor => new(_knotRt[^1]);

    public Irt ToIrt(RtMinutes retentionTime) => new(Evaluate(_knotRt, _knotIrt, retentionTime.Value, _slopeBefore, _slopeAfter));

    public RtMinutes ToRtMinutes(Irt indexedRetentionTime) =>
        new(Evaluate(_knotIrt, _knotRt, indexedRetentionTime.Value, 1 / _slopeBefore, 1 / _slopeAfter));

    /// <summary>True outside the anchors' span, where the map is a straight-line continuation rather than fitted.</summary>
    public bool IsExtrapolated(RtMinutes retentionTime) => retentionTime.Value < _knotRt[0] || retentionTime.Value > _knotRt[^1];

    private static double Evaluate(double[] from, double[] to, double value, double slopeBefore, double slopeAfter)
    {
        if (value <= from[0])
            return to[0] + slopeBefore * (value - from[0]);
        if (value >= from[^1])
            return to[^1] + slopeAfter * (value - from[^1]);
        int hi = Array.BinarySearch(from, value);
        if (hi >= 0)
            return to[hi];
        hi = ~hi;
        int lo = hi - 1;
        double t = (value - from[lo]) / (from[hi] - from[lo]);
        return to[lo] + t * (to[hi] - to[lo]);
    }
}
