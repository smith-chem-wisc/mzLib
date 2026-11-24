using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using StatisticalModels;

namespace Quantification.Differential;

/// <summary>
/// The normalization check (QuantProject NORMCHECK, <c>DEF-DIFF-NORMCHECK</c>; <see cref="DifferentialNormalizationCheck"/>):
/// after fitting, whether a two-sided comparison's results look shifted by a normalization error. It changes nothing in
/// the results; a comparison it flags carries a warning, and one it does not flag is never called fine.
/// </summary>
/// <remarks>
/// <para>
/// <b>The score.</b> Each fitted protein of the comparison (status <c>fitted</c> or <c>fitted_single_peptide</c>, with a
/// finite effect, standard error above 0 and finite df) gets z = Φ⁻¹(T_df(t)), t = log2 effect / standard error, the
/// z-value of Efron (2004), clamped to [−8, 8]. It is taken as sign(t) × −Φ⁻¹(T_df(−|t|)), which keeps its digits in the
/// tails.
/// </para>
/// <para>
/// <b>The fit.</b> <see cref="EmpiricalNull.Fit"/> fits Efron's empirical null to the scores. If the normalization is
/// right, the proteins that did not change centre on 0; an error of e log2 moves them by about e over their standard
/// errors, and δ0 with them. A comparison whose |δ0| is above <see cref="Threshold"/> for its group size is a warning.
/// </para>
/// <para>
/// <b>Where it was calibrated</b> (QuantProject NORMCHECK third attempt, 2026-11-24): label-free comparisons of two groups
/// run by this runner at mzLib 1.0.600's defaults (<c>shared_peptide_median</c>, robust weighting on), on MS/MS values
/// only. The threshold is the 95th percentile of |δ0| over null comparisons (two groups drawn from one condition) of 23
/// development deposits, deposits weighted equally, for groups of 2, 3 and 4 samples. On 19 test deposits it warned on
/// 4.7% of null comparisons, and in comparisons with 2,000 or more fitted proteins it caught 79% (groups of 3) and 73%
/// (groups of 4 or more) of planted normalization errors of 0.2 log2 or more; in smaller comparisons it caught 16% to 46%.
/// Of the comparisons it did not flag, 18.7% still had an error of 0.2 or more, so no "fine" verdict is given. Anywhere
/// else the check is not run, and the reason is recorded.
/// </para>
/// </remarks>
public static class NormalizationCheck
{
    /// <summary>A comparison the check flags: its |δ0| is above the threshold.</summary>
    public const string Warning = "warning";
    /// <summary>A comparison the check ran on and did not flag. It is not a verdict that the normalization is right.</summary>
    public const string NoWarning = "no_warning";
    /// <summary>A comparison the check did not run on; the reason says why.</summary>
    public const string NotChecked = "not_checked";

    /// <summary>Fewer than <see cref="MinimumFittedProteins"/> fitted proteins.</summary>
    public const string FewerThanMinimumFitted = "fewer_than_200_fitted";
    /// <summary>The empirical null could not be fitted to the scores.</summary>
    public const string FitFailed = "fit_failed";
    /// <summary>A slope (a covariate contrast) has no two groups.</summary>
    public const string Slope = "slope";
    /// <summary>The smaller group has one sample.</summary>
    public const string GroupOfOne = "group_of_one";
    /// <summary>Prefix of a reason naming a setting the check was not calibrated on, e.g. <c>not_calibrated:mbr_kept</c>.</summary>
    public const string NotCalibratedPrefix = "not_calibrated:";

    /// <summary>The fewest fitted proteins the check runs on.</summary>
    public const int MinimumFittedProteins = 200;
    /// <summary>The bound the scores are clamped to.</summary>
    public const double ScoreBound = 8;
    /// <summary>The threshold table's name, written with every check so a result names the table it was judged by.</summary>
    public const string Yardstick = "normcheck_yardstick_1";

    /// <summary>The calibration's numbers, for the report: development and test deposits.</summary>
    public const int DevelopmentDeposits = 23, TestDeposits = 19;
    /// <summary>The share of test null comparisons it warned on.</summary>
    public const double TestFalseWarningRate = 0.047;
    /// <summary>The detection rates for a true error of 0.2 log2 or more, in comparisons with
    /// <see cref="PoweredFittedProteins"/> or more fitted proteins: groups of 4 or more, and groups of 3.</summary>
    public const double TestDetectionLow = 0.73, TestDetectionHigh = 0.79;
    /// <summary>The fitted proteins from which the check caught most errors of 0.2 log2 or more.</summary>
    public const int PoweredFittedProteins = 2000;

    /// <summary>
    /// The threshold on |δ0| for a comparison whose smaller group has <paramref name="groupSize"/> samples: 0.4107 for 2,
    /// 0.3633 for 3, 0.4298 for 4 or more. Groups of 5 or more take the value measured on groups of 4; no larger group was
    /// measured.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="groupSize"/> is below 2.</exception>
    public static double Threshold(int groupSize) => groupSize switch
    {
        < 2 => throw new ArgumentOutOfRangeException(nameof(groupSize), groupSize, "A group of fewer than 2 samples has no threshold."),
        2 => 0.410702939399882,
        3 => 0.363279703570296,
        _ => 0.429775562501974,
    };

    /// <summary>Whether a comparison with this δ0 and group size warns: |δ0| above <see cref="Threshold"/>, so a δ0 on the
    /// threshold does not.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="groupSize"/> is below 2.</exception>
    public static bool Warns(double delta0, int groupSize) => Math.Abs(delta0) > Threshold(groupSize);

    /// <summary>The clamped z-value of one fitted protein: sign(t) × −Φ⁻¹(T_df(−|t|)), t = effect / standard error.</summary>
    public static double Score(double log2Effect, double standardError, double df)
    {
        double t = log2Effect / standardError;
        double tail = double.IsPositiveInfinity(df) ? Normal.CDF(0, 1, -Math.Abs(t)) : StudentT.CDF(0, 1, df, -Math.Abs(t));
        double z = tail <= 0 ? ScoreBound : -Normal.InvCDF(0, 1, tail);
        z = Math.Min(z, ScoreBound);
        return t < 0 ? -z : z;
    }

    /// <summary>
    /// Runs the check on one comparison's fitted proteins.
    /// </summary>
    /// <param name="stratum">The stratum.</param>
    /// <param name="quantBasis">The quant basis, or null for a style with none.</param>
    /// <param name="contrastId">The contrast.</param>
    /// <param name="fitted">Each fitted protein's log2 effect, standard error and df; only those with a finite effect, a
    /// finite standard error above 0 and a finite df are scored.</param>
    /// <param name="groupSize">The smaller group's sample count.</param>
    /// <exception cref="ArgumentNullException">A string or the list is null.</exception>
    public static DifferentialNormalizationCheck Check(string stratum, string? quantBasis, string contrastId,
        IEnumerable<(double Log2Effect, double StandardError, double Df)> fitted, int groupSize)
    {
        ArgumentNullException.ThrowIfNull(stratum);
        ArgumentNullException.ThrowIfNull(contrastId);
        ArgumentNullException.ThrowIfNull(fitted);
        var usable = fitted.Where(f => double.IsFinite(f.Log2Effect) && double.IsFinite(f.StandardError) && f.StandardError > 0
                                      && double.IsFinite(f.Df)).ToList();
        if (groupSize < 2)
            return NotRun(stratum, quantBasis, contrastId, GroupOfOne, usable.Count, groupSize);
        if (usable.Count < MinimumFittedProteins)
            return NotRun(stratum, quantBasis, contrastId, FewerThanMinimumFitted, usable.Count, groupSize);

        var fit = EmpiricalNull.Fit(usable.Select(f => Score(f.Log2Effect, f.StandardError, f.Df)).ToArray());
        if (fit.Status != EmpiricalNullStatus.Estimated)
            return NotRun(stratum, quantBasis, contrastId, FitFailed, usable.Count, groupSize);

        double threshold = Threshold(groupSize);
        double eHat = fit.Delta * usable.Select(f => f.StandardError).Median();
        return new DifferentialNormalizationCheck(stratum, quantBasis, contrastId,
            Warns(fit.Delta, groupSize) ? Warning : NoWarning, null, usable.Count, groupSize, fit.Delta, fit.Sigma, eHat,
            threshold, Yardstick);
    }

    /// <summary>A comparison the check did not run on, with the reason.</summary>
    /// <exception cref="ArgumentNullException">A string is null.</exception>
    public static DifferentialNormalizationCheck NotRun(string stratum, string? quantBasis, string contrastId, string reason,
        int fittedProteins, int? groupSize)
    {
        ArgumentNullException.ThrowIfNull(stratum);
        ArgumentNullException.ThrowIfNull(contrastId);
        ArgumentNullException.ThrowIfNull(reason);
        return new DifferentialNormalizationCheck(stratum, quantBasis, contrastId, NotChecked, reason, fittedProteins, groupSize,
            null, null, null, null, Yardstick);
    }
}
