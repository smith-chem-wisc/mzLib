using System.Collections.Generic;

namespace Quantification.Differential;

/// <summary>
/// One row of a differential analysis, as QuantProject <c>DEF-DIFF-COLUMNS</c> (DATA-DEFINITIONS v3.7) defines it: one
/// (analysis, stratum, feature, contrast, quant basis, method). Each property is one column; its machine and readable
/// header names live in <see cref="DifferentialColumns"/>. A null (or non-finite) number is written as an empty cell,
/// never 0 and never <c>NaN</c>; a 0 is a real value.
/// </summary>
public sealed record DifferentialResult
{
    // ── Identity ──
    /// <summary><c>QuantProject:DEF-DIFF-&lt;QUANTITY&gt; v&lt;n&gt;</c> for the row's quantity.</summary>
    public required string DefinitionId { get; init; }
    /// <summary>The analysis this row belongs to (<see cref="DifferentialMetadata.AnalysisId"/>).</summary>
    public required string AnalysisId { get; init; }
    /// <summary><c>lfq</c>, <c>tmt</c>, <c>itraq</c>, <c>dileu</c>, <c>silac</c>, <c>pulse_silac</c>.</summary>
    public string? QuantStyle { get; init; }
    /// <summary><c>ms2</c> or <c>ms3_sps</c> for isobaric styles; null otherwise.</summary>
    public string? ReporterAcquisition { get; init; }
    /// <summary>The stratum key (<c>DEF-DIFF-STRATUM</c>), or <c>all</c>.</summary>
    public required string Stratum { get; init; }
    /// <summary><c>protein_group</c>, <c>peptidoform</c>, <c>site</c>.</summary>
    public required string Grain { get; init; }
    /// <summary>The feature within the dataset (protein group name, full sequence, or site key).</summary>
    public required string FeatureId { get; init; }
    /// <summary>Every member accession; written ordinal-sorted, joined by <c>;</c>.</summary>
    public IReadOnlyList<string> FeatureAccessions { get; init; } = new List<string>();
    /// <summary>The members' gene names, aligned with <see cref="FeatureAccessions"/> as given; null or empty where none.</summary>
    public IReadOnlyList<string?> Genes { get; init; } = new List<string?>();
    /// <summary><c>abundance</c>, <c>occupancy</c>, …</summary>
    public required string Quantity { get; init; }
    /// <summary>The effect type (<c>DEF-DIFF-EFFECT</c>), e.g. <c>abundance_log2_ratio</c>.</summary>
    public required string EffectType { get; init; }
    /// <summary><c>msms_only</c> or <c>mbr_kept</c>; null for styles with no MBR.</summary>
    public string? QuantBasis { get; init; }

    // ── Contrast ──
    /// <summary>Unique within the analysis, e.g. <c>c1</c>.</summary>
    public required string ContrastId { get; init; }
    /// <summary>Readable, e.g. <c>age=old vs age=young</c>.</summary>
    public string? ContrastLabel { get; init; }
    /// <summary>The numerator level as <c>factor=level</c>; null for a slope.</summary>
    public string? Numerator { get; init; }
    /// <summary>The reference level as <c>factor=level</c>; null for a slope.</summary>
    public string? Denominator { get; init; }
    /// <summary>For a slope, the covariate.</summary>
    public string? Covariate { get; init; }
    /// <summary>For a slope, the unit the effect is per.</summary>
    public string? CovariateUnit { get; init; }
    /// <summary>For a slope, how the stored covariate became that unit.</summary>
    public string? CovariateScale { get; init; }

    // ── Effect ──
    /// <summary>The effect in log2 units (<c>DEF-DIFF-EFFECT</c>).</summary>
    public double? Log2Effect { get; init; }
    /// <summary>2^<see cref="Log2Effect"/>.</summary>
    public double? NaturalEffect { get; init; }
    /// <summary><c>ratio</c>, <c>odds_ratio</c>, <c>rate_ratio</c>, <c>ratio_per_&lt;unit&gt;</c>.</summary>
    public string? NaturalUnit { get; init; }
    /// <summary>Occupancy only: 100 × (numerator mean − denominator mean).</summary>
    public double? DeltaPercentagePoints { get; init; }
    /// <summary>The model's estimated mean for the numerator level.</summary>
    public double? MeanNumerator { get; init; }
    /// <summary>The model's estimated mean for the reference level.</summary>
    public double? MeanDenominator { get; init; }

    // ── Confidence ──
    /// <summary>Standard error of <see cref="Log2Effect"/>.</summary>
    public double? StandardError { get; init; }
    /// <summary>Lower confidence bound.</summary>
    public double? CiLow { get; init; }
    /// <summary>Upper confidence bound.</summary>
    public double? CiHigh { get; init; }
    /// <summary>The interval's level, e.g. 0.95.</summary>
    public double? CiLevel { get; init; }
    /// <summary>The test statistic (moderated t).</summary>
    public double? Statistic { get; init; }
    /// <summary>Degrees of freedom of the test.</summary>
    public double? Df { get; init; }
    /// <summary><c>moderated_residual</c> or <c>satterthwaite_moderated</c>.</summary>
    public string? DfMethod { get; init; }
    /// <summary>Two-sided p-value.</summary>
    public double? PValue { get; init; }
    /// <summary>Benjamini–Hochberg adjusted p within the row's family.</summary>
    public double? PAdjusted { get; init; }
    /// <summary><c>benjamini_hochberg</c>.</summary>
    public string? AdjustmentMethod { get; init; }
    /// <summary>Fitted rows in the row's family; must equal the count the writer finds.</summary>
    public int? FamilySize { get; init; }
    /// <summary>Bayesian method only: posterior error probability.</summary>
    public double? Pep { get; init; }
    /// <summary>Bayesian method only: Bayesian FDR (never a Benjamini–Hochberg value).</summary>
    public double? BayesianFdr { get; init; }
    /// <summary>Bayesian method only: Bayes factor.</summary>
    public double? BayesFactor { get; init; }
    /// <summary>Bayesian method only: half-width of the no-change interval, log2 units.</summary>
    public double? NullWidth { get; init; }

    // ── Evidence ──
    /// <summary>Numerator-side biological samples with at least one value.</summary>
    public int? NSamplesNumerator { get; init; }
    /// <summary>Reference-side biological samples with at least one value.</summary>
    public int? NSamplesDenominator { get; init; }
    /// <summary>Biological samples in the contrast with at least one value.</summary>
    public int? NSamplesWithValue { get; init; }
    /// <summary>Biological samples in the contrast, with or without a value.</summary>
    public int? NSamplesTotal { get; init; }
    /// <summary>Distinct peptides with at least one value that entered the model.</summary>
    public int? NPeptides { get; init; }
    /// <summary>(peptide, sample) values that entered the model.</summary>
    public int? NObservations { get; init; }
    /// <summary>Of <see cref="NObservations"/>, MBR transfers.</summary>
    public int? NMbrValues { get; init; }

    // ── Status and provenance ──
    /// <summary>The status (<see cref="DifferentialStatus"/>).</summary>
    public required string Status { get; init; }
    /// <summary>Free text: which rule failed and the counts it saw; null on a plain fitted row.</summary>
    public string? StatusDetail { get; init; }
    /// <summary><c>moderated</c> or <c>bayesian</c>.</summary>
    public required string Method { get; init; }
    /// <summary><c>peptide_mixed_model</c>, <c>moderated_t</c>, …</summary>
    public string? ModelUsed { get; init; }
    /// <summary>Whether outlying values were down-weighted.</summary>
    public bool? Robust { get; init; }
    /// <summary>The normalization that ran (<c>DEF-DIFF-NORM</c>).</summary>
    public string? Normalization { get; init; }
    /// <summary>Model terms besides the contrast's own factor, random ones marked, as given.</summary>
    public string? CovariatesFitted { get; init; }
    /// <summary>sha256 of the design the row was fitted on.</summary>
    public string? DesignSha256 { get; init; }
    /// <summary>The framework's version string.</summary>
    public string? MethodVersion { get; init; }
}
