using System;
using System.Collections.Generic;

namespace Quantification.Differential;

/// <summary>
/// One biological sample as the statistics see it: which stratum it is in, its level of each condition factor, the
/// individual it came from, its batch and its numeric covariates. A technical replicate or a fraction is never a sample.
/// </summary>
/// <param name="SampleId">The sample's key; unique within a dataset.</param>
public sealed record AnalysisSample(string SampleId)
{
    /// <summary>The sample's level of each condition factor (age group, treatment, ...), keyed by factor name.</summary>
    public IReadOnlyDictionary<string, string> Factors { get; init; } = new Dictionary<string, string>();

    /// <summary>What the sample is (organism part, cell type, preparation, ...), keyed by name; empty when one stratum.</summary>
    public IReadOnlyDictionary<string, string> Strata { get; init; } = new Dictionary<string, string>();

    /// <summary>The animal, donor or culture the sample came from (the random intercept, GR-6); null when not given.</summary>
    public string? Individual { get; init; }

    /// <summary>The sample's batch (TMT plex, acquisition batch); null when not given.</summary>
    public string? Batch { get; init; }

    /// <summary>Numeric covariates in their stored unit (e.g. age in years), keyed by name.</summary>
    public IReadOnlyDictionary<string, double> Covariates { get; init; } = new Dictionary<string, double>();
}

/// <summary>A condition factor of the model.</summary>
/// <param name="Name">The factor's name, as in <see cref="AnalysisSample.Factors"/>.</param>
/// <param name="ReferenceLevel">
/// The level every other level is compared against. Null when none was declared: every pair of levels is then
/// compared and the design warns (GR-14).
/// </param>
public sealed record FactorSpec(string Name, string? ReferenceLevel = null);

/// <summary>A numeric covariate of the model, entered as <c>(value − centre) / scale</c>.</summary>
/// <param name="Name">The covariate's name, as in <see cref="AnalysisSample.Covariates"/>.</param>
/// <param name="Unit">The stored unit, e.g. <c>years</c>.</param>
/// <param name="Scale">Stored units per scaled unit, e.g. 10 years per decade. Default 1.</param>
/// <param name="ScaledUnit">The unit a slope is per, e.g. <c>decade</c>; <paramref name="Unit"/> when null.</param>
/// <param name="Centre">The value subtracted first, in the stored unit; the stratum's mean when null.</param>
public sealed record CovariateSpec(string Name, string Unit, double Scale = 1, string? ScaledUnit = null, double? Centre = null);

/// <summary>What goes into the fixed-effects model.</summary>
public sealed record DesignSpec
{
    /// <summary>Condition factors, in column order. All enter one additive model (GR-1).</summary>
    public IReadOnlyList<FactorSpec> Factors { get; init; } = Array.Empty<FactorSpec>();

    /// <summary>Numeric covariates, in column order.</summary>
    public IReadOnlyList<CovariateSpec> Covariates { get; init; } = Array.Empty<CovariateSpec>();

    /// <summary>Two-factor interactions, only when asked for; each names two factors of <see cref="Factors"/>.</summary>
    public IReadOnlyList<(string First, string Second)> Interactions { get; init; } = Array.Empty<(string, string)>();

    /// <summary>Fit batch as a fixed block when the samples carry one. Default true.</summary>
    public bool FitBatch { get; init; } = true;
}

/// <summary>A comparison to estimate: a weight per design column (<c>DEF-DIFF-COLUMNS</c>' Contrast group).</summary>
/// <param name="Id">Unique within an analysis: <c>c1</c>, <c>c2</c>, … for the default set.</param>
/// <param name="Label">Readable, e.g. <c>age_group=old vs age_group=young</c> or <c>age (per decade)</c>.</param>
/// <param name="Weights">One weight per column of <see cref="AnalysisDesign.ColumnNames"/>.</param>
/// <param name="Numerator">The numerator level as <c>factor=level</c>; null for a slope or a custom contrast.</param>
/// <param name="Denominator">The reference level as <c>factor=level</c>; null for a slope or a custom contrast.</param>
/// <param name="Covariate">For a slope, the covariate.</param>
/// <param name="CovariateUnit">For a slope, the unit the effect is per.</param>
/// <param name="CovariateScale">For a slope, how the stored covariate became that unit, e.g. <c>years/10, centred at 50 years</c>.</param>
public sealed record Contrast(string Id, string Label, IReadOnlyList<double> Weights, string? Numerator, string? Denominator,
    string? Covariate = null, string? CovariateUnit = null, string? CovariateScale = null);

/// <summary>A contrast the specification implies that a stratum cannot run, and why (e.g. <c>not_in_stratum</c>).</summary>
public sealed record NotRunContrast(string Label, string Reason);
