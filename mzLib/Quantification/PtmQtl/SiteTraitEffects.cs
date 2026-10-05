using System;
using System.Collections.Generic;
using System.Linq;
using StatisticalModels;

namespace Quantification.PtmQtl;

/// <summary>One site's association between occupancy and a trait.</summary>
/// <param name="Site">The site.</param>
/// <param name="Status">
/// <c>Fitted</c>, <c>BelowSupport</c> (fewer runs or replicates than the options require), or the
/// <see cref="FeatureFitStatus"/> name of a fit that could not be computed. Every non-fitted value is NaN.
/// </param>
/// <param name="Effect">Change in logit occupancy per unit of the trait, other covariates held fixed.</param>
/// <param name="StandardError">Standard error of <paramref name="Effect"/>.</param>
/// <param name="DegreesOfFreedom">Denominator df of the t test (containment rule, see <see cref="MixedModel"/>).</param>
/// <param name="PValue">Two-sided p-value of the trait coefficient.</param>
/// <param name="Q">Benjamini-Hochberg q over the fitted sites of one call.</param>
/// <param name="Runs">Runs contributing a response.</param>
/// <param name="Replicates">Distinct biological replicates among those runs.</param>
/// <param name="MedianFraction">Median occupancy over the contributing runs.</param>
/// <param name="CovariateEffects">Coefficient of each covariate, in the order the covariates were named.</param>
/// <param name="ReplicateVariance">Between-replicate variance of logit occupancy (the random intercept).</param>
/// <param name="ResidualVariance">Within-replicate (run-to-run) variance of logit occupancy.</param>
public sealed record SiteTraitEffect(ModificationSite Site, string Status, double Effect, double StandardError,
    double DegreesOfFreedom, double PValue, double Q, int Runs, int Replicates, double MedianFraction,
    IReadOnlyList<double> CovariateEffects, double ReplicateVariance, double ResidualVariance);

/// <summary>Which occupancy cells a trait fit uses, and how much support a site needs.</summary>
public sealed record SiteTraitOptions
{
    /// <summary>
    /// Leave out a run where the site was seen only modified (occupancy 1 as a ceiling,
    /// <see cref="SiteRunOccupancy.UnmodifiedQuantified"/> false). Default true.
    /// </summary>
    public bool ExcludeCeiling { get; init; } = true;

    /// <summary>A fraction is clamped to [ε, 1 − ε] before the logit. Default 5e-5, half the 4-dp resolution of stored fractions.</summary>
    public double LogitEpsilon { get; init; } = 5e-5;

    /// <summary>Minimum runs with a response. Default 6.</summary>
    public int MinRuns { get; init; } = 6;

    /// <summary>Minimum distinct biological replicates with a response. Default 4.</summary>
    public int MinReplicates { get; init; } = 4;
}

/// <summary>
/// Per-site association between intensity occupancy and a consumer-defined trait (D6's feature_trait_effect,
/// intensity part), within one dataset.
/// </summary>
/// <remarks>
/// <para>
/// Model, per site: logit(occupancy) = β0 + β_trait·trait + Σ γ_k·covariate_k + u_replicate + ε, a random-intercept
/// linear mixed model (<see cref="MixedModel"/>, REML), so repeated runs of one biological replicate (technical
/// injections) are not counted as independent. The trait and covariates are per run; the trait is normally
/// constant within a replicate and is then tested at the between-replicate df.
/// </para>
/// <para>
/// Only <see cref="OccupancyState.Quantified"/> cells are responses. Floor, CountOnly and NotDetected are left
/// out, never entered as 0: this is the intensity part of the hurdle only, and it says nothing about whether the
/// modification's detection changes with the trait.
/// </para>
/// </remarks>
public static class SiteTraitEffects
{
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

        double eps = options.LogitEpsilon;
        var sites = occupancy
            .Where(o => o.State == OccupancyState.Quantified && runIndex.ContainsKey(o.Run)
                        && (!options.ExcludeCeiling || o.UnmodifiedQuantified) && double.IsFinite(o.Fraction))
            .GroupBy(o => o.Site)
            .OrderBy(g => g.Key.Key, StringComparer.Ordinal)
            .ToList();

        var fitted = new List<int>();
        var responses = new List<double[]>();
        var results = new SiteTraitEffect[sites.Count];
        var none = Enumerable.Repeat(double.NaN, k).ToArray();
        for (int s = 0; s < sites.Count; s++)
        {
            var cells = sites[s].GroupBy(o => o.Run).Select(g => g.First()).ToList();
            int reps = cells.Select(o => replicate[o.Run]).Distinct(StringComparer.Ordinal).Count();
            double median = Median(cells.Select(o => o.Fraction));
            results[s] = new SiteTraitEffect(sites[s].Key, "BelowSupport", double.NaN, double.NaN, double.NaN, double.NaN,
                double.NaN, cells.Count, reps, median, none, double.NaN, double.NaN);
            if (cells.Count < options.MinRuns || reps < options.MinReplicates) continue;
            var y = Enumerable.Repeat(double.NaN, runs.Count).ToArray();
            foreach (var o in cells)
            {
                double f = Math.Clamp(o.Fraction, eps, 1 - eps);
                y[runIndex[o.Run]] = Math.Log(f / (1 - f));
            }
            fitted.Add(s);
            responses.Add(y);
        }

        if (fitted.Count > 0)
        {
            var matrix = new double[fitted.Count, runs.Count];
            for (int f = 0; f < fitted.Count; f++)
                for (int r = 0; r < runs.Count; r++) matrix[f, r] = responses[f][r];
            var names = new[] { "intercept", "trait" }.Concat(covariateNames).ToList();
            var fit = MixedModel.Fit(matrix, design, groups, names);
            for (int f = 0; f < fitted.Count; f++)
            {
                var status = fit.Status[f];
                var prior = results[fitted[f]];
                results[fitted[f]] = status != FeatureFitStatus.Fitted
                    ? prior with { Status = status.ToString() }
                    : prior with
                    {
                        Status = "Fitted",
                        Effect = fit.Coefficient(f, 1),
                        StandardError = fit.StandardError(f, 1),
                        DegreesOfFreedom = fit.DegreesOfFreedom(f, 1),
                        PValue = fit.PValue(f, 1),
                        CovariateEffects = Enumerable.Range(0, k).Select(j => fit.Coefficient(f, 2 + j)).ToArray(),
                        ReplicateVariance = fit.GroupVariance[f],
                        ResidualVariance = fit.ResidualVariance[f],
                    };
            }
        }

        var tested = Enumerable.Range(0, results.Length).Where(i => double.IsFinite(results[i].PValue)).ToList();
        var q = MultipleTesting.BenjaminiHochberg(tested.Select(i => results[i].PValue).ToList());
        for (int t = 0; t < tested.Count; t++) results[tested[t]] = results[tested[t]] with { Q = q[t] };
        return results;
    }

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

    private static double Median(IEnumerable<double> values)
    {
        var v = values.OrderBy(x => x).ToArray();
        return v.Length == 0 ? double.NaN : v.Length % 2 == 1 ? v[v.Length / 2] : 0.5 * (v[v.Length / 2 - 1] + v[v.Length / 2]);
    }
}
