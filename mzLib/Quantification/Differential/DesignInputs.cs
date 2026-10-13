using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using MassSpectrometry;
using Omics.BioPolymerGroup;
using StatisticalModels;

namespace Quantification.Differential;

/// <summary>
/// One biological sample as the statistics see it: which stratum it is in, its level of each condition factor, the
/// individual it came from, its batch and its numeric covariates. A technical replicate or a fraction is never a sample.
/// </summary>
/// <param name="SampleId">The sample's key; unique within a dataset.</param>
public sealed record AnalysisSample(string SampleId)
{
    /// <summary>The factor <see cref="FromExperimentalDesign"/> puts each sample's condition under.</summary>
    public const string ConditionFactor = "condition";

    /// <summary>
    /// The biological sample an experimental-design entry belongs to, as a <see cref="SampleId"/>.
    /// </summary>
    /// <remarks>
    /// A label-free run gives <c>{condition}_{biological replicate + 1}</c> (<see cref="SampleGroupLabels.ForSample"/>):
    /// the label MetaMorpheus's quantified tables use and the one the observation table's runs carry, so a sample's
    /// fractions and technical replicates share it. An isobaric channel gives its design's sample name, or
    /// <c>plex{plex}_{channel}</c> when the design names none, because one plex's channel is one sample in every fraction
    /// file of that plex. (<see cref="SampleGroupLabels.ForSample"/> labels each file's channel separately, for columns,
    /// so it is not used for channels.)
    /// </remarks>
    public static string IdFor(ISampleInfo sample)
    {
        ArgumentNullException.ThrowIfNull(sample);
        if (sample is IsobaricQuantSampleInfo channel)
            return string.IsNullOrWhiteSpace(channel.SampleName)
                ? string.Create(CultureInfo.InvariantCulture, $"plex{channel.PlexId}_{channel.ChannelLabel}")
                : channel.SampleName;
        return SampleGroupLabels.ForSample(sample);
    }

    /// <summary>
    /// The biological samples of a pipeline's experimental design, each with its condition under
    /// <see cref="ConditionFactor"/> (STATS-FRAMEWORK section 3.2: grouping comes only from the design). The runs of one
    /// sample (its fractions and technical replicates) become one sample, keyed by <see cref="IdFor"/>. An isobaric
    /// sample's batch is its plex, <c>plex{n}</c>; a pooled reference channel is not a biological sample and is left out.
    /// An experimental design has no individuals, strata or covariates; add them with <c>with</c>.
    /// </summary>
    /// <returns>The samples, ordered by id (ordinal).</returns>
    /// <exception cref="ArgumentException">
    /// The design has a null entry or mixes label-free files with isobaric channels; one sample's runs disagree on its
    /// condition; or, for isobaric channels, one sample name is on two channels or one channel carries two names. A sample
    /// name is display text, so channels are never merged by it.
    /// </exception>
    public static IReadOnlyList<AnalysisSample> FromExperimentalDesign(IExperimentalDesign design)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(design.FileNameSampleInfoDictionary, nameof(design));
        var entries = design.FileNameSampleInfoDictionary.Values.SelectMany(s => s ?? Array.Empty<ISampleInfo>()).ToList();
        if (entries.Any(e => e is null))
            throw new ArgumentException("The experimental design has a null sample.", nameof(design));
        if (entries.Any(e => e is IsobaricQuantSampleInfo) && entries.Any(e => e is not IsobaricQuantSampleInfo))
            throw new ArgumentException(
                "The experimental design mixes label-free files with isobaric channels; a design is one or the other.",
                nameof(design));

        var channels = entries.OfType<IsobaricQuantSampleInfo>()
            .GroupBy(c => (c.PlexId, c.ChannelLabel))
            .FirstOrDefault(g => g.Select(IdFor).Distinct(StringComparer.Ordinal).Count() > 1);
        if (channels is not null)
            throw new ArgumentException(
                $"Channel {channels.Key.ChannelLabel} of plex {channels.Key.PlexId} is named " +
                $"{string.Join(" and ", channels.Select(IdFor).Distinct(StringComparer.Ordinal).Select(n => $"'{n}'"))} " +
                "in different files; one channel of one plex is one sample.", nameof(design));

        var samples = new List<AnalysisSample>();
        foreach (var runs in entries
                     .Where(e => e is not IsobaricQuantSampleInfo { IsReferenceChannel: true })
                     .GroupBy(IdFor, StringComparer.Ordinal)
                     .OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            ISampleInfo first = runs.First();
            ISampleInfo? other = runs.FirstOrDefault(r => !string.Equals(r.Condition, first.Condition, StringComparison.Ordinal));
            if (other is not null)
                throw new ArgumentException(
                    $"Sample '{runs.Key}' is condition '{first.Condition}' in '{Path.GetFileName(first.FullFilePathWithExtension)}' " +
                    $"but '{other.Condition}' in '{Path.GetFileName(other.FullFilePathWithExtension)}'.", nameof(design));

            string? batch = null;
            if (first is IsobaricQuantSampleInfo)
            {
                var places = runs.Cast<IsobaricQuantSampleInfo>().Select(c => (c.PlexId, c.ChannelLabel)).Distinct()
                    .OrderBy(c => c.PlexId).ThenBy(c => c.ChannelLabel, StringComparer.Ordinal).ToList();
                if (places.Count > 1)
                    throw new ArgumentException(
                        $"Sample name '{runs.Key}' is on {string.Join(" and ", places.Select(c => $"channel {c.ChannelLabel} of plex {c.PlexId}"))}; " +
                        "a sample name is display text, so channels are not merged by it. Name each channel's sample once.",
                        nameof(design));
                batch = string.Create(CultureInfo.InvariantCulture, $"plex{places[0].PlexId}");
            }

            samples.Add(new AnalysisSample(runs.Key)
            {
                Factors = new Dictionary<string, string> { [ConditionFactor] = first.Condition },
                Batch = batch,
            });
        }
        return samples;
    }

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
/// <param name="ScaledUnit">
/// The unit a slope is per, e.g. <c>decade</c>. When null, a slope is per <paramref name="Scale"/> stored units:
/// <c>years</c> at scale 1, <c>10 years</c> at scale 10.
/// </param>
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
    string? Covariate = null, string? CovariateUnit = null, string? CovariateScale = null)
{
    /// <summary>
    /// The contrast as the model fits take it (<see cref="EmpiricalBayes.ModerateContrasts"/> and the mixed models), named
    /// by <see cref="Id"/>: ids are unique within an analysis, so each moderated result maps back to its contrast by name.
    /// </summary>
    public ContrastWeights ToWeights() => new(Id, Weights);
}

/// <summary>A contrast the specification implies that a stratum cannot run, and why (e.g. <c>not_in_stratum</c>).</summary>
public sealed record NotRunContrast(string Label, string Reason);
