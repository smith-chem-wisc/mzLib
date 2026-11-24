using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Optimization;
using MathNet.Numerics.Statistics;

namespace Quantification.Differential;

/// <summary>
/// The normalization check (QuantProject NORMCHECK, <c>DEF-DIFF-NORMCHECK</c>): does a two-group comparison look shifted
/// by a normalization error? Each fitted protein's score is Efron's (2004) z-value of its test, sign(effect) Φ⁻¹(1 − p/2),
/// clamped to [−8, 8]. Efron's empirical null is fitted to the scores; if the normalization is right, the proteins that did
/// not change centre on 0, and the comparison warns when that centre δ0 is further from 0 than <see cref="Threshold"/>.
/// </summary>
/// <remarks>
/// Calibrated on label-free two-group comparisons at this runner's defaults (QuantProject NORMCHECK third attempt,
/// 2026-11-24): on 19 test deposits it warned on 4.7% of null comparisons and, with 2,000 or more fitted proteins, caught
/// 73% to 79% of normalization errors of 0.2 log2 or more. Of the comparisons it did not flag, 18.7% still had such an
/// error, so no comparison is ever called fine.
/// </remarks>
public static class NormalizationCheck
{
    /// <summary>The fewest fitted proteins the check runs on.</summary>
    public const int MinimumFittedProteins = 200;

    /// <summary>
    /// The |δ0| above which a comparison whose smaller group has <paramref name="groupSize"/> samples warns: the 95th
    /// percentile of |δ0| over null comparisons of 23 development deposits, deposits weighted equally. Groups of 5 or more
    /// take the value for 4; none larger was measured.
    /// </summary>
    public static double Threshold(int groupSize) => groupSize switch
    {
        < 2 => throw new ArgumentOutOfRangeException(nameof(groupSize), groupSize, "A group needs 2 samples or more."),
        2 => 0.410702939399882,
        3 => 0.363279703570296,
        _ => 0.429775562501974,
    };

    /// <summary>δ0 of the empirical null fitted to these proteins' scores, or NaN when it cannot be fitted.</summary>
    /// <param name="fitted">Each fitted protein's log2 effect and two-sided p-value.</param>
    public static double Centre(IEnumerable<(double Log2Effect, double PValue)> fitted) =>
        FitNull(fitted.Select(f => Math.Clamp(Math.Sign(f.Log2Effect) * -Normal.InvCDF(0, 1, Math.Max(f.PValue / 2, double.Epsilon)), -8, 8))
            .ToArray())?.Delta ?? double.NaN;

    /// <summary>The check of one comparison, as the phrase the metadata and the report carry.</summary>
    internal static string Verdict(IReadOnlyList<(double Log2Effect, double PValue)> fitted, int groupSize)
    {
        if (groupSize < 2) return "not checked: a group of one sample";
        if (fitted.Count < MinimumFittedProteins)
            return $"not checked: {fitted.Count} fitted proteins, fewer than {MinimumFittedProteins}";
        double delta0 = Centre(fitted);
        if (double.IsNaN(delta0)) return "not checked: no empirical null could be fitted";
        string d = delta0.ToString("0.###", CultureInfo.InvariantCulture);
        string t = Threshold(groupSize).ToString("0.###", CultureInfo.InvariantCulture);
        return Math.Abs(delta0) > Threshold(groupSize)
            ? $"warning: the unchanged proteins' scores centre at {d}, further from 0 than {t} for its group size, so its log2 changes may all be shifted"
            : $"no warning (centre {d}, threshold {t})";
    }

    /// <summary>
    /// Efron's empirical null N(δ0, σ0²) by maximum likelihood on the scores in a central window (Efron 2007, Ann. Stat.
    /// 35:1351, eq. 4.12: the normal truncated to the window). The window is R locfdr 1.1.8's default, which the papers do
    /// not state: the median ± b (Q3 − Q1) / (2 Φ⁻¹(0.75)), b = 4.3 exp(−0.26 log10 N), then refitted on δ0 ± b σ0 of
    /// the first fit. Written from that description, not from locfdr's source; null when either pass cannot be fitted.
    /// </summary>
    internal static (double Delta, double Sigma)? FitNull(double[] z)
    {
        if (z.Length < 3) return null;
        Array.Sort(z);
        double b = 4.3 * Math.Exp(-0.26 * Math.Log10(z.Length));
        double Q(double p) => SortedArrayStatistics.QuantileCustom(z, p, QuantileDefinition.R7);
        var first = FitWindow(z, Q(0.5), b * (Q(0.75) - Q(0.25)) / (2 * Normal.InvCDF(0, 1, 0.75)));
        return first is { } f ? FitWindow(z, f.Delta, b * f.Sigma) : null;
    }

    /// <summary>The normal truncated to [centre − half, centre + half] that maximises the likelihood of the scores in it,
    /// found by Nelder–Mead on the scores rescaled to [−1, 1], where it depends on them only through their mean and mean
    /// square.</summary>
    private static (double Delta, double Sigma)? FitWindow(double[] sorted, double centre, double half)
    {
        if (!(half > 0)) return null;
        double[] u = sorted.Where(v => v >= centre - half && v <= centre + half).Select(v => (v - centre) / half).ToArray();
        double m1 = u.Length < 3 ? 0 : u.Average(), m2 = u.Length < 3 ? 0 : u.Average(x => x * x);
        if (!(m2 - m1 * m1 > 0)) return null;
        double NegativeLogLikelihood(Vector<double> p)
        {
            double mu = p[0], s = Math.Exp(p[1]);
            double inside = Normal.CDF(0, 1, (1 - mu) / s) - Normal.CDF(0, 1, (-1 - mu) / s);
            if (!(inside > 0)) return double.PositiveInfinity;
            return Math.Log(s) + Math.Log(inside) + (m2 - 2 * mu * m1 + mu * mu) / (2 * s * s);
        }
        try
        {
            var start = Vector<double>.Build.DenseOfArray(new[] { m1, 0.5 * Math.Log(m2 - m1 * m1) });
            var best = NelderMeadSimplex.Minimum(ObjectiveFunction.Value(NegativeLogLikelihood), start,
                Vector<double>.Build.DenseOfArray(new[] { 0.1, 0.1 }), 1e-15, 20_000).MinimizingPoint;
            double sigma = Math.Exp(best[1]);
            return double.IsFinite(best[0]) && sigma < 1e4 ? (centre + half * best[0], half * sigma) : null;
        }
        catch (MaximumIterationsException)
        {
            return null;
        }
    }
}
