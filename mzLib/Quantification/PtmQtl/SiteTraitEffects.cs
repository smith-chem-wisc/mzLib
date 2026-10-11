using System;
using System.Collections.Generic;
using System.Linq;
using StatisticalModels;

namespace Quantification.PtmQtl;

/// <summary>One site's association between occupancy and a trait.</summary>
/// <param name="Site">The site.</param>
/// <param name="Status">
/// <c>Fitted</c>; <c>BelowSupport</c> (fewer runs or replicates than the options require); with
/// <see cref="BoundCellModel.Sampling"/>, <c>MostlyOneBound</c> (fewer than <see cref="SiteTraitEffects.MinRunsOffExtremes"/>
/// runs that are not 0 of N or N of N), <c>NoDetectionLink</c> or <c>NoScatterToBorrow</c>; or the
/// <see cref="FeatureFitStatus"/> name of a fit that could not be computed. Every non-fitted value is NaN.
/// </param>
/// <param name="Effect">Change in the log odds of occupancy per unit of the trait, other covariates held fixed, on <see cref="Scale"/>.</param>
/// <param name="StandardError">Standard error of <paramref name="Effect"/>.</param>
/// <param name="DegreesOfFreedom">Denominator df of the t test (containment rule, see <see cref="MixedModel"/>).</param>
/// <param name="PValue">Two-sided p-value of the trait coefficient.</param>
/// <param name="Q">Benjamini-Hochberg q over the fitted sites of one call.</param>
/// <param name="Runs">Runs contributing a response (with <see cref="BoundCellModel.Sampling"/>, measured or counted).</param>
/// <param name="Replicates">Distinct biological replicates among those runs.</param>
/// <param name="MedianFraction">Median occupancy over the contributing runs with both forms quantified.</param>
/// <param name="CovariateEffects">Coefficient of each covariate, in the order the covariates were named.</param>
/// <param name="ReplicateVariance">Between-replicate variance of the log odds of occupancy (the random intercept).</param>
/// <param name="ResidualVariance">Within-replicate (run-to-run) variance of the log odds of occupancy.</param>
public sealed record SiteTraitEffect(ModificationSite Site, string Status, double Effect, double StandardError,
    double DegreesOfFreedom, double PValue, double Q, int Runs, int Replicates, double MedianFraction,
    IReadOnlyList<double> CovariateEffects, double ReplicateVariance, double ResidualVariance)
{
    /// <summary>The scale of every effect, standard error and variance above (<see cref="SiteTraitOptions.Scale"/>).</summary>
    public OccupancyScale Scale { get; init; }

    /// <summary>Of <see cref="Runs"/>, those that entered as counts ("k of N"); 0 unless <see cref="BoundCellModel.Sampling"/>.</summary>
    public int CountedRuns { get; init; }

    /// <summary>
    /// Whether <see cref="ResidualVariance"/> was borrowed (the call's median over sites with their own) because the site
    /// has fewer than <see cref="SiteTraitEffects.MinMeasuredRunsForScatter"/> measured runs.
    /// </summary>
    public bool ResidualVarianceBorrowed { get; init; }

    /// <summary>The detection link the counted runs were read through, on natural log odds; null when none were.</summary>
    public DetectionLink? Link { get; init; }
}

/// <summary>The scale of a site-trait fit's response: the log odds of occupancy f, f / (1 − f), in natural log or log2.</summary>
public enum OccupancyScale
{
    /// <summary>ln(f / (1 − f)), the logit.</summary>
    NaturalLogOdds,

    /// <summary>log2(f / (1 − f)): an effect of +1 doubles the odds. The base of QuantProject's occupancy effects.</summary>
    Log2Odds,
}

/// <summary>How a site-trait fit uses a run where the site lacks one of the two intensities.</summary>
public enum BoundCellModel
{
    /// <summary>Left out: only runs with both forms quantified are responses.</summary>
    Exclude,

    /// <summary>
    /// Every covered run without both intensities (a covered-zero, a ceiling, a floor, a count-only cell) enters as
    /// "k of the N covering peptidoforms carried the modification" (<see cref="SiteRunOccupancy.CoveringPeptidoforms"/>,
    /// <see cref="SiteRunOccupancy.ModifiedPeptidoforms"/>), read through a detection link (<see cref="LogitNormalMixedModel"/>).
    /// A run whose counts were not supplied is left out.
    /// </summary>
    Sampling,
}

/// <summary>Which occupancy cells a trait fit uses, and how much support a site needs.</summary>
public sealed record SiteTraitOptions
{
    /// <summary>
    /// Leave out a run where the site was seen only modified (occupancy 1 as a ceiling,
    /// <see cref="SiteRunOccupancy.UnmodifiedQuantified"/> false). Null (default): true, or false with
    /// <see cref="BoundCellModel.Sampling"/>, where a ceiling is a count; true together with Sampling is refused.
    /// </summary>
    public bool? ExcludeCeiling { get; init; }

    /// <summary>A fraction is clamped to [ε, 1 − ε] before the logit. Default 5e-5, half the 4-dp resolution of stored fractions.</summary>
    public double LogitEpsilon { get; init; } = 5e-5;

    /// <summary>Minimum runs with a response. Default 6.</summary>
    public int MinRuns { get; init; } = 6;

    /// <summary>Minimum distinct biological replicates with a response. Default 4.</summary>
    public int MinReplicates { get; init; } = 4;

    /// <summary>
    /// The response scale. Default <see cref="OccupancyScale.NaturalLogOdds"/>. Under <see cref="OccupancyScale.Log2Odds"/>
    /// the effects and standard errors are the natural-log fit's / ln 2 and the variances / (ln 2)²; df, p and q are the
    /// natural-log fit's, unchanged.
    /// </summary>
    public OccupancyScale Scale { get; init; } = OccupancyScale.NaturalLogOdds;

    /// <summary>How runs without both intensities enter. Default <see cref="BoundCellModel.Exclude"/>.</summary>
    public BoundCellModel BoundCells { get; init; } = BoundCellModel.Exclude;

    /// <summary>
    /// With <see cref="BoundCellModel.Sampling"/>: a depth per run (e.g. log2 of its target PSMs) for the link's depth
    /// term, centred on its median over the call's runs. Every run of the trait needs one. Null (default): no depth term.
    /// </summary>
    public IReadOnlyDictionary<string, double>? RunDepth { get; init; }

    /// <summary>With <see cref="BoundCellModel.Sampling"/>: a link to use instead of fitting one from the call's sites.</summary>
    public DetectionLink? DetectionLink { get; init; }

    /// <summary>
    /// With <see cref="BoundCellModel.Sampling"/>: the link to use when the call has fewer than
    /// <see cref="SiteTraitEffects.MinLinkCells"/> cells to fit one from (e.g. a link fitted over many datasets). Null:
    /// sites with counted runs are then not fitted (<c>NoDetectionLink</c>).
    /// </summary>
    public DetectionLink? FallbackLink { get; init; }
}

/// <summary>
/// Per-site association between intensity occupancy and a consumer-defined trait (D6's feature_trait_effect),
/// within one dataset.
/// </summary>
/// <remarks>
/// <para>
/// Model, per site: log odds(occupancy) = β0 + β_trait·trait + Σ γ_k·covariate_k + u_replicate + ε (natural log, or
/// log2 by <see cref="SiteTraitOptions.Scale"/>), a random-intercept
/// linear mixed model (<see cref="MixedModel"/>, REML), so repeated runs of one biological replicate (technical
/// injections) are not counted as independent. The trait and covariates are per run; the trait is normally
/// constant within a replicate and is then tested at the between-replicate df.
/// </para>
/// <para>
/// By default only runs with both forms quantified are responses; Floor, CountOnly and NotDetected are left out,
/// never entered as 0. With <see cref="BoundCellModel.Sampling"/> (ptmQtl D37, D39, D44-D49) every other covered run
/// enters as a count, "k of N covering peptidoforms carried the modification", whose chance per peptidoform is
/// linked to the same latent log odds (<see cref="LogitNormalMixedModel"/>, approximate REML). The link is fitted
/// once per call, on whether the modified form was seen, from the sites with at least
/// <see cref="MinMeasuredRunsForScatter"/> measured runs (their mean and SD of log odds as each cell's latent).
/// A site needs at least <see cref="MinRunsOffExtremes"/> runs that are not 0 of N or N of N (else the trait could
/// separate them and its estimate run to infinity); one with fewer than <see cref="MinMeasuredRunsForScatter"/>
/// measured runs borrows the call's median run-to-run variance. A site with no counted run is fitted by
/// <see cref="MixedModel"/> exactly as under <see cref="BoundCellModel.Exclude"/>.
/// </para>
/// </remarks>
public static class SiteTraitEffects
{
    /// <summary>With <see cref="BoundCellModel.Sampling"/>: the fewest runs a site needs that are measured or counted with 0 &lt; k &lt; N.</summary>
    public const int MinRunsOffExtremes = 3;

    /// <summary>With <see cref="BoundCellModel.Sampling"/>: the fewest measured runs for a site's own run-to-run variance, and for it to inform the link.</summary>
    public const int MinMeasuredRunsForScatter = 3;

    /// <summary>With <see cref="BoundCellModel.Sampling"/>: the fewest cells to fit the link from; below this the fallback link is used.</summary>
    public const int MinLinkCells = 100;

    /// <summary>Fits the trait against every site's occupancy.</summary>
    /// <param name="occupancy">Per-run occupancy cells, any states; each site is fitted separately.</param>
    /// <param name="trait">Trait value per run. Runs absent here are not used.</param>
    /// <param name="replicate">Biological replicate (random-intercept group) per run. Every run in <paramref name="trait"/> needs one.</param>
    /// <param name="covariates">Optional covariate values per run, one per name in <paramref name="covariateNames"/>, all finite.</param>
    /// <param name="covariateNames">Names of the covariates, in order.</param>
    /// <param name="options">Cell selection and support thresholds; defaults when null.</param>
    public static IReadOnlyList<SiteTraitEffect> Fit(IEnumerable<SiteRunOccupancy> occupancy,
        IReadOnlyDictionary<string, double> trait, IReadOnlyDictionary<string, string> replicate,
        IReadOnlyDictionary<string, double[]>? covariates = null, IReadOnlyList<string>? covariateNames = null,
        SiteTraitOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(occupancy);
        ArgumentNullException.ThrowIfNull(trait);
        ArgumentNullException.ThrowIfNull(replicate);
        options ??= new SiteTraitOptions();
        covariateNames ??= Array.Empty<string>();
        int k = covariateNames.Count;
        if (k > 0 && covariates is null)
            throw new ArgumentException("Covariate names were given without covariate values.", nameof(covariates));
        bool sampling = options.BoundCells switch
        {
            BoundCellModel.Exclude => false,
            BoundCellModel.Sampling => true,
            _ => throw new ArgumentOutOfRangeException(nameof(options), options.BoundCells, "Unknown bound-cell model."),
        };
        if (sampling && options.ExcludeCeiling == true)
            throw new ArgumentException("With BoundCells = Sampling a ceiling is a count; leave ExcludeCeiling unset or false.", nameof(options));
        bool excludeCeiling = options.ExcludeCeiling ?? !sampling;

        var runs = trait.Keys.OrderBy(r => r, StringComparer.Ordinal).ToList();
        var runIndex = runs.Select((r, i) => (r, i)).ToDictionary(t => t.r, t => t.i, StringComparer.Ordinal);
        var design = new double[runs.Count, 2 + k];
        var groups = new string[runs.Count];
        for (int i = 0; i < runs.Count; i++)
        {
            string run = runs[i];
            if (!double.IsFinite(trait[run]))
                throw new ArgumentException($"The trait is not finite for run '{run}'.", nameof(trait));
            if (!replicate.TryGetValue(run, out var g) || string.IsNullOrEmpty(g))
                throw new ArgumentException($"Run '{run}' has no biological replicate.", nameof(replicate));
            groups[i] = g;
            design[i, 0] = 1;
            design[i, 1] = trait[run];
            if (k == 0) continue;
            if (!covariates!.TryGetValue(run, out var c) || c.Length != k || c.Any(v => !double.IsFinite(v)))
                throw new ArgumentException($"Run '{run}' needs {k} finite covariate values.", nameof(covariates));
            for (int j = 0; j < k; j++) design[i, 2 + j] = c[j];
        }
        double[]? depth = null;
        if (sampling && options.RunDepth != null)
        {
            var raw = runs.Select(r => options.RunDepth.TryGetValue(r, out double d) && double.IsFinite(d) ? d
                : throw new ArgumentException($"Run '{r}' needs a finite depth.", nameof(options))).ToArray();
            double median = Median(raw);
            depth = raw.Select(d => d - median).ToArray();
        }
        if (sampling && options.RunDepth == null && (options.DetectionLink?.DepthSlope ?? 0) != 0)
            throw new ArgumentException("The given detection link has a depth term; supply RunDepth.", nameof(options));

        double eps = options.LogitEpsilon;
        // The fit is always on the natural-log scale, and a Log2Odds result is divided by ln 2 afterwards. Refitting on
        // log2 would move the variances, and with them p, by the REML optimizer's tolerance (1.5e-7 relative in a test).
        double unit = options.Scale switch
        {
            OccupancyScale.NaturalLogOdds => 1,
            OccupancyScale.Log2Odds => Math.Log(2),
            _ => throw new ArgumentOutOfRangeException(nameof(options), options.Scale, "Unknown occupancy scale."),
        };

        // Per site and run: a measured log odds, or (with Sampling) a count, or nothing.
        var sites = occupancy
            .Where(o => runIndex.ContainsKey(o.Run))
            .GroupBy(o => o.Site)
            .OrderBy(g => g.Key.Key, StringComparer.Ordinal)
            .Select(g => (site: g.Key, cells: g.GroupBy(o => o.Run).Select(c => c.First())
                .Select(o => (cell: o, measured: Measured(o, excludeCeiling, sampling, eps), count: sampling ? Count(o, excludeCeiling, eps) : null))
                .Where(t => double.IsFinite(t.measured) || t.count != null).ToList()))
            .Where(s => s.cells.Count > 0)
            .ToList();

        var results = new SiteTraitEffect[sites.Count];
        var none = Enumerable.Repeat(double.NaN, k).ToArray();
        var measuredOnly = new List<int>();
        var withCounts = new List<int>();
        for (int s = 0; s < sites.Count; s++)
        {
            var cells = sites[s].cells;
            int reps = cells.Select(t => replicate[t.cell.Run]).Distinct(StringComparer.Ordinal).Count();
            int counted = cells.Count(t => t.count != null);
            int offExtremes = cells.Count(t => t.count is not BinomialCount c || (c.Successes > 0 && c.Successes < c.Trials));
            double median = Median(cells.Where(t => double.IsFinite(t.measured)).Select(t => t.cell.Fraction).ToArray());
            results[s] = new SiteTraitEffect(sites[s].site, "BelowSupport", double.NaN, double.NaN, double.NaN, double.NaN,
                double.NaN, cells.Count, reps, median, none, double.NaN, double.NaN) { Scale = options.Scale, CountedRuns = counted };
            if (cells.Count < options.MinRuns || reps < options.MinReplicates) continue;
            if (counted == 0) { measuredOnly.Add(s); continue; }
            if (offExtremes < MinRunsOffExtremes) { results[s] = results[s] with { Status = "MostlyOneBound" }; continue; }
            withCounts.Add(s);
        }

        var names = new[] { "intercept", "trait" }.Concat(covariateNames).ToList();
        var ownScatter = new List<double>();
        if (measuredOnly.Count > 0)
        {
            var matrix = new double[measuredOnly.Count, runs.Count];
            for (int f = 0; f < measuredOnly.Count; f++)
            {
                for (int r = 0; r < runs.Count; r++) matrix[f, r] = double.NaN;
                foreach (var t in sites[measuredOnly[f]].cells) matrix[f, runIndex[t.cell.Run]] = t.measured;
            }
            var fit = MixedModel.Fit(matrix, design, groups, names);
            for (int f = 0; f < measuredOnly.Count; f++)
            {
                int s = measuredOnly[f];
                results[s] = Read(results[s], fit.Status[f], fit.Coefficient, fit.StandardError, fit.DegreesOfFreedom, fit.PValue,
                    f, k, unit, fit.GroupVariance[f], fit.ResidualVariance[f]);
                if (fit.Status[f] == FeatureFitStatus.Fitted && sites[s].cells.Count(t => double.IsFinite(t.measured)) >= MinMeasuredRunsForScatter)
                    ownScatter.Add(fit.ResidualVariance[f]);
            }
        }

        if (withCounts.Count > 0)
        {
            var link = options.DetectionLink ?? FitLink(sites, runIndex, depth);
            link ??= options.FallbackLink;
            if (link is null)
                foreach (int s in withCounts) results[s] = results[s] with { Status = "NoDetectionLink" };
            else
            {
                if (link.DepthSlope != 0 && depth == null)
                    throw new ArgumentException("The fallback detection link has a depth term; supply RunDepth.", nameof(options));
                var own = withCounts.Where(s => sites[s].cells.Count(t => double.IsFinite(t.measured)) >= MinMeasuredRunsForScatter).ToList();
                var borrowing = withCounts.Except(own).ToList();
                foreach (var (batch, borrowed) in new[] { (own, false), (borrowing, true) })
                {
                    if (batch.Count == 0) continue;
                    double scatter = double.NaN;
                    if (borrowed)
                    {
                        if (ownScatter.Count == 0)
                        {
                            foreach (int s in batch) results[s] = results[s] with { Status = "NoScatterToBorrow" };
                            continue;
                        }
                        scatter = Median(ownScatter.ToArray());
                    }
                    var y = new double[batch.Count, runs.Count];
                    var counts = new BinomialCount?[batch.Count, runs.Count];
                    for (int f = 0; f < batch.Count; f++)
                    {
                        for (int r = 0; r < runs.Count; r++) y[f, r] = double.NaN;
                        foreach (var t in sites[batch[f]].cells)
                        {
                            if (double.IsFinite(t.measured)) y[f, runIndex[t.cell.Run]] = t.measured;
                            else counts[f, runIndex[t.cell.Run]] = t.count;
                        }
                    }
                    var fit = LogitNormalMixedModel.Fit(y, counts, design, groups, link, depth,
                        Enumerable.Repeat(scatter, batch.Count).ToArray(), names);
                    for (int f = 0; f < batch.Count; f++)
                    {
                        int s = batch[f];
                        results[s] = Read(results[s], fit.Status[f], fit.Coefficient, fit.StandardError, fit.DegreesOfFreedom, fit.PValue,
                            f, k, unit, fit.GroupVariance[f], fit.ResidualVariance[f]) with { ResidualVarianceBorrowed = borrowed, Link = link };
                        if (!borrowed && fit.Status[f] == FeatureFitStatus.Fitted) ownScatter.Add(fit.ResidualVariance[f]);
                    }
                }
            }
        }

        var tested = Enumerable.Range(0, results.Length).Where(i => double.IsFinite(results[i].PValue)).ToList();
        var q = MultipleTesting.BenjaminiHochberg(tested.Select(i => results[i].PValue).ToList());
        for (int t = 0; t < tested.Count; t++) results[tested[t]] = results[tested[t]] with { Q = q[t] };
        return results;
    }

    /// <summary>
    /// A Quantified run's log odds; else NaN. A ceiling is NaN when left out, and under Sampling, where it is a count;
    /// otherwise (ExcludeCeiling false, no Sampling) it reads 1, clamped, as before.
    /// </summary>
    private static double Measured(SiteRunOccupancy o, bool excludeCeiling, bool sampling, double eps)
    {
        if (o.State != OccupancyState.Quantified || !double.IsFinite(o.Fraction)) return double.NaN;
        if (!o.UnmodifiedQuantified && (excludeCeiling || sampling)) return double.NaN;
        double f = Math.Clamp(o.Fraction, eps, 1 - eps);
        return Math.Log(f / (1 - f));
    }

    /// <summary>With Sampling, a covered run that is not measured: its count of modified among covering peptidoforms; else null.</summary>
    private static BinomialCount? Count(SiteRunOccupancy o, bool excludeCeiling, double eps)
    {
        if (double.IsFinite(Measured(o, excludeCeiling, true, eps))) return null;
        if (!Enum.IsDefined(o.State))
            throw new ArgumentException($"{o.Site.Key} in {o.Run} has an unknown state {o.State}.", "occupancy");
        if (o.CoveringPeptidoforms is not int n || n < 1 || o.ModifiedPeptidoforms is not int m || m < 0 || m > n) return null;
        return new BinomialCount(m, n);
    }

    /// <summary>
    /// The call's link, fitted on whether the modified form was seen, from every counted cell of the sites with at least
    /// <see cref="MinMeasuredRunsForScatter"/> measured runs; null below <see cref="MinLinkCells"/> cells or when the fit fails.
    /// </summary>
    private static DetectionLink? FitLink(
        List<(ModificationSite site, List<(SiteRunOccupancy cell, double measured, BinomialCount? count)> cells)> sites,
        Dictionary<string, int> runIndex, double[]? depth)
    {
        var linkCells = new List<LinkCell>();
        foreach (var (_, cells) in sites)
        {
            var values = cells.Where(t => double.IsFinite(t.measured)).Select(t => t.measured).ToArray();
            if (values.Length < MinMeasuredRunsForScatter) continue;
            double mean = values.Average();
            double sd = Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1));
            foreach (var (cell, _, _) in cells)
            {
                if (cell.CoveringPeptidoforms is not int n || n < 1 || cell.ModifiedPeptidoforms is not int m || m < 0 || m > n) continue;
                linkCells.Add(new LinkCell(mean, sd, m, n, depth?[runIndex[cell.Run]] ?? 0));
            }
        }
        return linkCells.Count < MinLinkCells ? null : LogitNormalMixedModel.FitLink(linkCells, withDepth: depth != null);
    }

    private static SiteTraitEffect Read(SiteTraitEffect prior, FeatureFitStatus status, Func<int, int, double> coefficient,
        Func<int, int, double> standardError, Func<int, int, double> df, Func<int, int, double> p, int f, int k, double unit,
        double groupVariance, double residualVariance) =>
        status != FeatureFitStatus.Fitted
            ? prior with { Status = status.ToString() }
            : prior with
            {
                Status = "Fitted",
                Effect = coefficient(f, 1) / unit,
                StandardError = standardError(f, 1) / unit,
                DegreesOfFreedom = df(f, 1),
                PValue = p(f, 1),
                CovariateEffects = Enumerable.Range(0, k).Select(j => coefficient(f, 2 + j) / unit).ToArray(),
                ReplicateVariance = groupVariance / (unit * unit),
                ResidualVariance = residualVariance / (unit * unit),
            };

    /// <summary>
    /// Pooled occupancy per run over a set of sites: Σ fraction × covering intensity / Σ covering intensity over the
    /// run's Quantified cells, i.e. the modified share of all covering signal. A run-level index, e.g. a sample-handling
    /// covariate. Uses <see cref="SiteRunOccupancy.Fraction"/> (the reported fraction when there is one) for the numerator.
    /// </summary>
    /// <param name="occupancy">Per-run occupancy cells.</param>
    /// <param name="include">Which cells count (e.g. Gln deamidation with both forms seen); all Quantified cells when null.</param>
    public static IReadOnlyDictionary<string, double> PooledOccupancy(IEnumerable<SiteRunOccupancy> occupancy,
        Func<SiteRunOccupancy, bool>? include = null)
    {
        ArgumentNullException.ThrowIfNull(occupancy);
        return occupancy
            .Where(o => o.State == OccupancyState.Quantified && double.IsFinite(o.Fraction)
                        && double.IsFinite(o.CoveringIntensity) && o.CoveringIntensity > 0 && (include?.Invoke(o) ?? true))
            .GroupBy(o => o.Run, StringComparer.Ordinal)
            .ToDictionary(g => g.Key, g => g.Sum(o => o.Fraction * o.CoveringIntensity) / g.Sum(o => o.CoveringIntensity),
                StringComparer.Ordinal);
    }

    private static double Median(IReadOnlyList<double> values)
    {
        var v = values.OrderBy(x => x).ToArray();
        return v.Length == 0 ? double.NaN : v.Length % 2 == 1 ? v[v.Length / 2] : 0.5 * (v[v.Length / 2 - 1] + v[v.Length / 2]);
    }
}
