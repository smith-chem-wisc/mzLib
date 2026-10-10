using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Omics.BioPolymerGroup;

namespace Quantification.Differential;

/// <summary>Which of a column's two names the header row uses (QuantProject GR-13).</summary>
public enum DifferentialHeaderStyle
{
    /// <summary>snake_case, e.g. <c>log2_effect</c>: for pipelines and anything that reads the file mechanically.</summary>
    Machine,
    /// <summary>Readable, e.g. <c>Log2 Effect</c>: for people running MetaMorpheus or the FlashLFQ app.</summary>
    Human,
}

/// <summary>
/// The status vocabulary of <c>DEF-DIFF-STATUS</c>. A fitted row (<see cref="Fitted"/>, <see cref="FittedSinglePeptide"/>)
/// carries every number; every other status is "evidence only": identity, contrast, evidence and status columns are written
/// and every effect and confidence cell is empty.
/// </summary>
public static class DifferentialStatus
{
    /// <summary>Fitted normally.</summary>
    public const string Fitted = "fitted";
    /// <summary>A protein with one usable peptide, fitted with the moderated t-test.</summary>
    public const string FittedSinglePeptide = "fitted_single_peptide";
    /// <summary>No value on the numerator side.</summary>
    public const string AbsentInNumerator = "absent_in_numerator";
    /// <summary>No value on the reference side.</summary>
    public const string AbsentInDenominator = "absent_in_denominator";
    /// <summary>Values in the stratum, but on neither side of this contrast.</summary>
    public const string AbsentInBoth = "absent_in_both";
    /// <summary>The fit did not converge.</summary>
    public const string NotConverged = "not_converged";

    private const string BelowSupportPrefix = "below_support:";
    private const string NotEstimablePrefix = "not_estimable:";

    /// <summary><c>below_support:&lt;rule&gt;</c>, e.g. <c>below_support:min_2_per_side</c>.</summary>
    public static string BelowSupport(string rule) => BelowSupportPrefix + Reason(rule, nameof(rule));

    /// <summary><c>not_estimable:&lt;reason&gt;</c>, e.g. <c>not_estimable:confounded</c>.</summary>
    public static string NotEstimable(string reason) => NotEstimablePrefix + Reason(reason, nameof(reason));

    /// <summary>Whether <paramref name="status"/> belongs to the vocabulary.</summary>
    public static bool IsValid(string status) =>
        status is Fitted or FittedSinglePeptide or AbsentInNumerator or AbsentInDenominator or AbsentInBoth or NotConverged
        || HasReason(status, BelowSupportPrefix) || HasReason(status, NotEstimablePrefix);

    /// <summary>Whether a row of this status carries an effect and its confidence.</summary>
    public static bool IsFitted(string status) => status is Fitted or FittedSinglePeptide;

    private static bool HasReason(string status, string prefix) =>
        status.StartsWith(prefix, StringComparison.Ordinal) && status.Length > prefix.Length && !status[prefix.Length..].Any(char.IsWhiteSpace);

    private static string Reason(string value, string name)
    {
        if (string.IsNullOrWhiteSpace(value) || value.Any(char.IsWhiteSpace))
            throw new ArgumentException("A status reason is one non-empty word, e.g. min_2_per_side.", name);
        return value;
    }
}

/// <summary>One column of the results table: both header names and the <c>DEF-DIFF-COLUMNS</c> group it belongs to.</summary>
/// <param name="MachineName">snake_case header.</param>
/// <param name="HumanName">Readable header.</param>
/// <param name="Group">Identity, Contrast, Effect, Confidence, Evidence or Status.</param>
public sealed record DifferentialColumn(string MachineName, string HumanName, string Group);

/// <summary>
/// The one registry of the results table's columns (<c>DEF-DIFF-COLUMNS</c> v1, 56 columns): order, both names, and how each
/// value is written. Numbers in invariant culture at round-trip precision; a null or non-finite number, and a null string,
/// as an empty cell; booleans <c>true</c>/<c>false</c>; accessions ordinal-sorted and joined by <c>;</c>, genes in that order.
/// </summary>
public static class DifferentialColumns
{
    private sealed record Entry(DifferentialColumn Column, Func<DifferentialResult, string?> Value);

    private static readonly Entry[] Entries =
    {
        E("definition_id", "Definition ID", "Identity", r => r.DefinitionId),
        E("analysis_id", "Analysis ID", "Identity", r => r.AnalysisId),
        E("quant_style", "Quantification Style", "Identity", r => r.QuantStyle),
        E("reporter_acquisition", "Reporter Ion Acquisition", "Identity", r => r.ReporterAcquisition),
        E("stratum", "Stratum", "Identity", r => r.Stratum),
        E("grain", "Grain", "Identity", r => r.Grain),
        E("feature_id", "Feature", "Identity", r => r.FeatureId),
        E("feature_accessions", "Accessions", "Identity", r => string.Join(";", SortedMembers(r).Select(m => m.Accession))),
        E("genes", "Genes", "Identity", Genes),
        E("quantity", "Quantity", "Identity", r => r.Quantity),
        E("effect_type", "Effect Type", "Identity", r => r.EffectType),
        E("quant_basis", "Quantification Basis", "Identity", r => r.QuantBasis),

        E("contrast_id", "Contrast ID", "Contrast", r => r.ContrastId),
        E("contrast_label", "Contrast", "Contrast", r => r.ContrastLabel),
        E("numerator", "Numerator", "Contrast", r => r.Numerator),
        E("denominator", "Denominator", "Contrast", r => r.Denominator),
        E("covariate", "Covariate", "Contrast", r => r.Covariate),
        E("covariate_unit", "Covariate Unit", "Contrast", r => r.CovariateUnit),
        E("covariate_scale", "Covariate Scale", "Contrast", r => r.CovariateScale),

        E("log2_effect", "Log2 Effect", "Effect", r => Number(r.Log2Effect)),
        E("natural_effect", "Effect (Natural Scale)", "Effect", r => Number(r.NaturalEffect)),
        E("natural_unit", "Natural Unit", "Effect", r => r.NaturalUnit),
        E("delta_percentage_points", "Change (Percentage Points)", "Effect", r => Number(r.DeltaPercentagePoints)),
        E("mean_numerator", "Numerator Mean", "Effect", r => Number(r.MeanNumerator)),
        E("mean_denominator", "Denominator Mean", "Effect", r => Number(r.MeanDenominator)),

        E("se", "Standard Error", "Confidence", r => Number(r.StandardError)),
        E("ci_low", "CI Lower", "Confidence", r => Number(r.CiLow)),
        E("ci_high", "CI Upper", "Confidence", r => Number(r.CiHigh)),
        E("ci_level", "CI Level", "Confidence", r => Number(r.CiLevel)),
        E("statistic", "Test Statistic", "Confidence", r => Number(r.Statistic)),
        E("df", "Degrees of Freedom", "Confidence", r => Number(r.Df)),
        E("df_method", "DF Method", "Confidence", r => r.DfMethod),
        E("p_value", "P-Value", "Confidence", r => Number(r.PValue)),
        E("p_adjusted", "Adjusted P-Value", "Confidence", r => Number(r.PAdjusted)),
        E("adjustment_method", "Adjustment Method", "Confidence", r => r.AdjustmentMethod),
        E("family_size", "Family Size", "Confidence", r => Integer(r.FamilySize)),
        E("pep", "Posterior Error Probability", "Confidence", r => Number(r.Pep)),
        E("bayesian_fdr", "Bayesian False Discovery Rate", "Confidence", r => Number(r.BayesianFdr)),
        E("bayes_factor", "Bayes Factor", "Confidence", r => Number(r.BayesFactor)),
        E("null_width", "Null Hypothesis Width", "Confidence", r => Number(r.NullWidth)),

        E("n_samples_numerator", "Numerator Samples With Value", "Evidence", r => Integer(r.NSamplesNumerator)),
        E("n_samples_denominator", "Denominator Samples With Value", "Evidence", r => Integer(r.NSamplesDenominator)),
        E("n_samples_with_value", "Samples With Value", "Evidence", r => Integer(r.NSamplesWithValue)),
        E("n_samples_total", "Samples In Contrast", "Evidence", r => Integer(r.NSamplesTotal)),
        E("n_peptides", "Peptides", "Evidence", r => Integer(r.NPeptides)),
        E("n_observations", "Peptide-Sample Values", "Evidence", r => Integer(r.NObservations)),
        E("n_mbr_values", "MBR Values", "Evidence", r => Integer(r.NMbrValues)),

        E("status", "Status", "Status", r => r.Status),
        E("status_detail", "Status Detail", "Status", r => r.StatusDetail),
        E("method", "Method", "Status", r => r.Method),
        E("model_used", "Model Used", "Status", r => r.ModelUsed),
        E("robust", "Robust Weighting", "Status", r => r.Robust is { } b ? (b ? "true" : "false") : null),
        E("normalization", "Normalization", "Status", r => r.Normalization),
        E("covariates_fitted", "Covariates Fitted", "Status", r => r.CovariatesFitted),
        E("design_sha256", "Design SHA-256", "Status", r => r.DesignSha256),
        E("method_version", "Method Version", "Status", r => r.MethodVersion),
    };

    /// <summary>Every column, in order.</summary>
    public static IReadOnlyList<DifferentialColumn> All { get; } = Entries.Select(e => e.Column).ToArray();

    /// <summary>The columns as a <see cref="TsvColumn{T}"/> schema with the chosen header names, for <see cref="TsvWriter"/>.</summary>
    public static IReadOnlyList<TsvColumn<DifferentialResult>> Schema(DifferentialHeaderStyle style) => Schema(style, null);

    internal static IReadOnlyList<TsvColumn<DifferentialResult>> Schema(DifferentialHeaderStyle style, Func<string?, string?>? clean) =>
        Entries.Select(e => new TsvColumn<DifferentialResult>(
                style == DifferentialHeaderStyle.Machine ? e.Column.MachineName : e.Column.HumanName,
                r => clean is null ? e.Value(r)! : clean(e.Value(r))!))
            .ToArray();

    private static Entry E(string machine, string human, string group, Func<DifferentialResult, string?> value) =>
        new(new DifferentialColumn(machine, human, group), value);

    /// <summary>Round-trip text in invariant culture; empty for null and for a non-finite value (never <c>NaN</c>).</summary>
    internal static string? Number(double? value) =>
        value is { } v && double.IsFinite(v) ? v.ToString("R", CultureInfo.InvariantCulture) : null;

    internal static string? Integer(int? value) => value?.ToString(CultureInfo.InvariantCulture);

    internal static IEnumerable<(string Accession, string? Gene)> SortedMembers(DifferentialResult r)
    {
        if (r.Genes.Count != 0 && r.Genes.Count != r.FeatureAccessions.Count)
            throw new ArgumentException(
                $"Feature '{r.FeatureId}' has {r.FeatureAccessions.Count} accessions but {r.Genes.Count} gene names; give one per accession, or none.");
        return r.FeatureAccessions
            .Select((a, i) => (Accession: a, Gene: r.Genes.Count == 0 ? null : r.Genes[i]))
            .OrderBy(m => m.Accession, StringComparer.Ordinal);
    }

    private static string? Genes(DifferentialResult r)
    {
        var members = SortedMembers(r).ToList();
        return members.All(m => string.IsNullOrEmpty(m.Gene)) ? null : string.Join(";", members.Select(m => m.Gene ?? ""));
    }
}

/// <summary>
/// Writes <c>DifferentialResults.tsv</c> (<c>DEF-DIFF-FILE</c>): UTF-8 without a byte-order mark, tab-separated, one header
/// row, <c>\n</c> line ends on every platform, rows in <c>DEF-DIFF-ROW</c>'s order whatever order they arrive in. Rows that
/// break the contract are refused before anything is written.
/// </summary>
/// <remarks>
/// Unlike <see cref="QuantificationWriter"/>, a missing number is an empty cell rather than <c>NaN</c>, doubles are written
/// at shortest round-trip precision rather than <c>G17</c>, and lines end in <c>\n</c> rather than the platform's newline,
/// so the same analysis gives the same bytes everywhere.
/// </remarks>
public static class DifferentialResultWriter
{
    /// <summary>The results table's file name.</summary>
    public const string ResultsFileName = "DifferentialResults.tsv";
    /// <summary>The metadata file's name (<see cref="DifferentialMetadataWriter"/>).</summary>
    public const string MetadataFileName = "DifferentialResults.metadata.json";

    /// <summary>
    /// <c>DEF-DIFF-ROW</c>'s order: stratum, grain, quantity, quant basis, contrast, method, then feature, each ordinal. All
    /// rows of one feature block are adjacent, so features never interleave.
    /// </summary>
    public static IReadOnlyList<DifferentialResult> Order(IEnumerable<DifferentialResult> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);
        return rows.OrderBy(r => r.Stratum, StringComparer.Ordinal)
            .ThenBy(r => r.Grain, StringComparer.Ordinal)
            .ThenBy(r => r.Quantity, StringComparer.Ordinal)
            .ThenBy(r => r.QuantBasis ?? "", StringComparer.Ordinal)
            .ThenBy(r => r.ContrastId, StringComparer.Ordinal)
            .ThenBy(r => r.Method, StringComparer.Ordinal)
            .ThenBy(r => r.FeatureId, StringComparer.Ordinal)
            .ToList();
    }

    /// <summary>
    /// Validates, orders and writes the table. A tab, CR or LF inside a text value is written as a space; the return value
    /// is how many cells that changed (the metadata's <c>cleaned_text_cells</c>).
    /// </summary>
    /// <exception cref="ArgumentException">A row breaks the contract: an unknown status; a number on an evidence-only row;
    /// a fitted row without an effect, SE or p-value; two rows with one key; a family size that disagrees with the rows;
    /// a missing required field.</exception>
    public static int Write(TextWriter output, IEnumerable<DifferentialResult> rows, DifferentialHeaderStyle style)
    {
        ArgumentNullException.ThrowIfNull(output);
        var ordered = Order(Validate(rows));
        int cleaned = 0;
        string? Clean(string? value)
        {
            if (value is null || value.IndexOfAny(Breaks) < 0) return value;
            cleaned++;
            return string.Concat(value.Select(c => Array.IndexOf(Breaks, c) >= 0 ? ' ' : c));
        }
        var schema = DifferentialColumns.Schema(style, Clean);
        output.Write(TsvWriter.HeaderLine(schema));
        output.Write('\n');
        foreach (var row in ordered)
        {
            output.Write(TsvWriter.RowLine(schema, row));
            output.Write('\n');
        }
        return cleaned;
    }

    /// <summary>
    /// Writes <see cref="ResultsFileName"/> (header style from the metadata's settings) and <see cref="MetadataFileName"/>
    /// into <paramref name="directory"/>, creating it, and returns the results file's path.
    /// </summary>
    public static string Write(string directory, IEnumerable<DifferentialResult> rows, DifferentialMetadata metadata)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(directory);
        ArgumentNullException.ThrowIfNull(metadata);
        var list = rows.ToList();
        var json = DifferentialMetadataWriter.ToBytes(metadata, list);
        using var buffer = new StringWriter();
        Write(buffer, list, metadata.Settings.HeaderStyle);
        Directory.CreateDirectory(directory);
        string path = Path.Combine(directory, ResultsFileName);
        File.WriteAllBytes(path, new UTF8Encoding(false).GetBytes(buffer.ToString()));
        File.WriteAllBytes(Path.Combine(directory, MetadataFileName), json);
        return path;
    }

    private static readonly char[] Breaks = { '\t', '\r', '\n' };

    internal static IReadOnlyList<DifferentialResult> Validate(IEnumerable<DifferentialResult> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);
        var list = rows.ToList();
        var keys = new HashSet<(string, string, string, string, string, string)>();
        foreach (var r in list)
        {
            ArgumentNullException.ThrowIfNull(r, nameof(rows));
            string at = $"Row for feature '{r.FeatureId}', contrast '{r.ContrastId}', stratum '{r.Stratum}'";
            foreach (var (name, value) in new[]
                     {
                         ("definition_id", r.DefinitionId), ("analysis_id", r.AnalysisId), ("stratum", r.Stratum),
                         ("grain", r.Grain), ("feature_id", r.FeatureId), ("quantity", r.Quantity),
                         ("effect_type", r.EffectType), ("contrast_id", r.ContrastId), ("status", r.Status), ("method", r.Method),
                     })
                if (string.IsNullOrWhiteSpace(value)) throw new ArgumentException($"{at} has no {name}.", nameof(rows));
            if (!DifferentialStatus.IsValid(r.Status))
                throw new ArgumentException($"{at} has status '{r.Status}', which is not in DEF-DIFF-STATUS.", nameof(rows));
            if (DifferentialStatus.IsFitted(r.Status))
            {
                if (!Finite(r.Log2Effect) || !Finite(r.StandardError) || !Finite(r.PValue))
                    throw new ArgumentException($"{at} is '{r.Status}' but lacks log2_effect, se or p_value.", nameof(rows));
            }
            else
            {
                var set = NumbersOnFittedRowsOnly(r).Where(c => c.Present).Select(c => c.Name).ToList();
                if (set.Count > 0)
                    throw new ArgumentException(
                        $"{at} is '{r.Status}', an evidence-only status, but sets {string.Join(", ", set)}; those cells must be empty.",
                        nameof(rows));
            }
            if (!keys.Add((r.AnalysisId, r.Stratum, r.FeatureId, r.ContrastId, r.QuantBasis ?? "", r.Method)))
                throw new ArgumentException($"{at} appears twice for quant basis '{r.QuantBasis}' and method '{r.Method}'.", nameof(rows));
            _ = DifferentialColumns.SortedMembers(r).ToList();
        }
        foreach (var family in list.Where(r => DifferentialStatus.IsFitted(r.Status)).GroupBy(FamilyKey))
        {
            int size = family.Count();
            var wrong = family.FirstOrDefault(r => r.FamilySize != size);
            if (wrong is not null)
                throw new ArgumentException(
                    $"Row for feature '{wrong.FeatureId}' says family_size {wrong.FamilySize?.ToString() ?? "(none)"}, but its family " +
                    $"({FamilyLabel(family.Key)}) has {size} fitted rows.", nameof(rows));
        }
        return list;
    }

    internal static (string Analysis, string Stratum, string Grain, string Quantity, string Basis, string Contrast, string Method)
        FamilyKey(DifferentialResult r) => (r.AnalysisId, r.Stratum, r.Grain, r.Quantity, r.QuantBasis ?? "", r.ContrastId, r.Method);

    private static string FamilyLabel((string, string, string, string, string, string, string) k) =>
        $"stratum {k.Item2}, grain {k.Item3}, quantity {k.Item4}, basis {k.Item5}, contrast {k.Item6}, method {k.Item7}";

    private static bool Finite(double? v) => v is { } x && double.IsFinite(x);

    private static IEnumerable<(string Name, bool Present)> NumbersOnFittedRowsOnly(DifferentialResult r)
    {
        yield return ("log2_effect", Finite(r.Log2Effect));
        yield return ("natural_effect", Finite(r.NaturalEffect));
        yield return ("delta_percentage_points", Finite(r.DeltaPercentagePoints));
        yield return ("mean_numerator", Finite(r.MeanNumerator));
        yield return ("mean_denominator", Finite(r.MeanDenominator));
        yield return ("se", Finite(r.StandardError));
        yield return ("ci_low", Finite(r.CiLow));
        yield return ("ci_high", Finite(r.CiHigh));
        yield return ("ci_level", Finite(r.CiLevel));
        yield return ("statistic", Finite(r.Statistic));
        yield return ("df", Finite(r.Df));
        yield return ("df_method", r.DfMethod is not null);
        yield return ("p_value", Finite(r.PValue));
        yield return ("p_adjusted", Finite(r.PAdjusted));
        yield return ("adjustment_method", r.AdjustmentMethod is not null);
        yield return ("family_size", r.FamilySize is not null);
        yield return ("pep", Finite(r.Pep));
        yield return ("bayesian_fdr", Finite(r.BayesianFdr));
        yield return ("bayes_factor", Finite(r.BayesFactor));
        yield return ("null_width", Finite(r.NullWidth));
    }
}

/// <summary>Software that produced an analysis.</summary>
public sealed record DifferentialSoftware(string Name, string Version);

/// <summary>An input of an analysis by role (<c>observation_table</c>, <c>design</c>, <c>contrast_spec</c>).</summary>
public sealed record DifferentialInput(string Role, string FileName, string Sha256);

/// <summary>The settings an analysis ran with; with the inputs, they determine <see cref="DifferentialMetadata.AnalysisId"/>.</summary>
public sealed record DifferentialSettings
{
    /// <summary>Header style of the results file.</summary>
    public DifferentialHeaderStyle HeaderStyle { get; init; } = DifferentialHeaderStyle.Machine;
    /// <summary>Two-sided interval level.</summary>
    public double CiLevel { get; init; } = 0.95;
    /// <summary><c>moderated</c> or <c>bayesian</c>.</summary>
    public string Method { get; init; } = "moderated";
    /// <summary>Whether outlying values were down-weighted.</summary>
    public bool Robust { get; init; } = true;
    /// <summary>The normalization setting (<c>DEF-DIFF-NORM</c>).</summary>
    public string Normalization { get; init; } = "shared_peptide_median";
    /// <summary>A curator's per-dataset override, if any.</summary>
    public string? NormalizationOverride { get; init; }
    /// <summary>A random seed, if any part of the analysis used one.</summary>
    public int? Seed { get; init; }
}

/// <summary>A contrast a stratum could not run, and why.</summary>
public sealed record DifferentialNotRun(string Label, string Reason);

/// <summary>
/// One stratum: its factors and values, which source chose them, samples per level, and contrasts not run there; and the
/// design facts the statistics report summarizes. A design fact left null is left out of the file.
/// </summary>
public sealed record DifferentialStratumInfo(string Stratum, IReadOnlyDictionary<string, string> Factors, string Source,
    IReadOnlyDictionary<string, int> SamplesPerLevel, IReadOnlyList<DifferentialNotRun> NotRun)
{
    /// <summary>Distinct individuals per level.</summary>
    public IReadOnlyDictionary<string, int>? IndividualsPerLevel { get; init; }
    /// <summary>The batch names in the stratum.</summary>
    public IReadOnlyList<string>? Batches { get; init; }
    /// <summary>How many distinct fraction numbers the stratum's runs carry.</summary>
    public int? Fractions { get; init; }
    /// <summary>How many distinct technical-replicate numbers the stratum's runs carry.</summary>
    public int? TechnicalReplicates { get; init; }
    /// <summary>Each factor's curated reference level; null when no reference is curated.</summary>
    public IReadOnlyDictionary<string, string>? ReferenceLevels { get; init; }
    /// <summary>The stratum's samples with no value at all.</summary>
    public IReadOnlyList<string>? SamplesWithoutValues { get; init; }
}

/// <summary>
/// A contrast's global shift in one stratum and quant basis: how far its two sides sit apart overall, before
/// normalization. Over the peptides with a value in every sample on both sides of the contrast in that stratum's table
/// for that basis, each peptide's mean log2 intensity on the numerator side minus its mean on the denominator side,
/// computed on intensities before normalization; <see cref="Value"/> is the median of those differences and
/// <see cref="Peptides"/> how many peptides it used. With no such peptide, <see cref="Value"/> is NaN (written
/// <c>null</c>) and <see cref="Peptides"/> is 0. <see cref="QuantBasis"/> is null for styles with no basis.
/// </summary>
public sealed record DifferentialGlobalShift(string Stratum, string? QuantBasis, double Value, int Peptides);

/// <summary>
/// One contrast, with its weight vector over model coefficients and, for a contrast with a numerator and a denominator,
/// its global shift in each stratum and quant basis it ran in. A covariate contrast has no two sides, so no global shift.
/// </summary>
public sealed record DifferentialContrastInfo(string Id, string Label, string? Numerator, string? Denominator,
    string? Covariate, string? CovariateUnit, string? CovariateScale, IReadOnlyDictionary<string, double> Weights,
    IReadOnlyList<DifferentialGlobalShift>? GlobalShift = null);

/// <summary>
/// The empirical-Bayes variance prior a fit used: the estimator, its degrees of freedom, whether it is trended, the prior
/// variance when it is one number, and how many features it was fitted on. A trended prior has a variance per feature,
/// so <see cref="Variance"/> is null and the values belong in the report's mean-variance table.
/// </summary>
public sealed record DifferentialPriorInfo(string Estimator, double Df, bool Trended, double? Variance, int Features)
{
    /// <summary>limma's legacy moment estimator.</summary>
    public const string MomentsLegacy = "moments_legacy";
    /// <summary>The marginal-likelihood estimator.</summary>
    public const string MarginalLikelihood = "marginal_likelihood";
}

/// <summary>
/// One model as fitted in one stratum and quant basis, for a method, the model it used and a grain: formula, REML or not,
/// df method, robust weighting rule, and its variance prior. <see cref="QuantBasis"/> is null for styles with no basis.
/// </summary>
public sealed record DifferentialModelInfo(string Stratum, string? QuantBasis, string Method, string ModelUsed, string Grain,
    string Formula, bool Reml, string DfMethod, string? RobustRule, DifferentialPriorInfo? Prior);

/// <summary>
/// One stratum's normalization in one quant basis: setting, reference set size, the floor that chose the setting, the
/// per-sample shifts and any warnings. The global shift is a property of a contrast, not of a stratum, so it is on
/// <see cref="DifferentialContrastInfo.GlobalShift"/>. <see cref="QuantBasis"/> is null for styles with no basis.
/// </summary>
public sealed record DifferentialNormalizationInfo(string Stratum, string? QuantBasis, string Setting, int ReferenceSetSize,
    int MinimumReferenceSetSize, IReadOnlyDictionary<string, double> PerSampleShift, IReadOnlyList<string> Warnings);

/// <summary>
/// Everything constant across an analysis's rows (<c>DEF-DIFF-META</c>), written by <see cref="DifferentialMetadataWriter"/>.
/// The columns, the Benjamini–Hochberg families, the row counts and the cleaned-cell count are not given here: the writer
/// derives them from the column registry and the rows, so they cannot disagree with the table.
/// </summary>
public sealed record DifferentialMetadata
{
    /// <summary>The definition version every row of this analysis follows.</summary>
    public const string Version = "DEF-DIFF v1";

    /// <summary>The calling program, mzLib, and the method version.</summary>
    public IReadOnlyList<DifferentialSoftware> Software { get; init; } = Array.Empty<DifferentialSoftware>();
    /// <summary>Each input by role, file name and sha256, in a fixed order.</summary>
    public IReadOnlyList<DifferentialInput> Inputs { get; init; } = Array.Empty<DifferentialInput>();
    /// <summary>The settings.</summary>
    public DifferentialSettings Settings { get; init; } = new();
    /// <summary>Per stratum.</summary>
    public IReadOnlyList<DifferentialStratumInfo> Strata { get; init; } = Array.Empty<DifferentialStratumInfo>();
    /// <summary>Per contrast.</summary>
    public IReadOnlyList<DifferentialContrastInfo> Contrasts { get; init; } = Array.Empty<DifferentialContrastInfo>();
    /// <summary>Per stratum, quant basis, method, model used and grain.</summary>
    public IReadOnlyList<DifferentialModelInfo> Models { get; init; } = Array.Empty<DifferentialModelInfo>();
    /// <summary>Per stratum and quant basis.</summary>
    public IReadOnlyList<DifferentialNormalizationInfo> Normalization { get; init; } = Array.Empty<DifferentialNormalizationInfo>();
    /// <summary>Every design warning (replicate, fraction or technical-replicate gaps).</summary>
    public IReadOnlyList<string> DesignWarnings { get; init; } = Array.Empty<string>();

    /// <summary>
    /// The first 16 hex characters of the sha256 of the inputs and settings, written canonically: the same inputs and
    /// settings always give the same id, and nothing else changes it.
    /// </summary>
    public string AnalysisId => ComputeAnalysisId(Inputs, Settings);

    /// <summary>See <see cref="AnalysisId"/>.</summary>
    public static string ComputeAnalysisId(IReadOnlyList<DifferentialInput> inputs, DifferentialSettings settings)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        ArgumentNullException.ThrowIfNull(settings);
        using var stream = new MemoryStream();
        using (var w = new Utf8JsonWriter(stream))
        {
            w.WriteStartObject();
            DifferentialMetadataWriter.WriteInputs(w, inputs);
            DifferentialMetadataWriter.WriteSettings(w, settings);
            w.WriteEndObject();
        }
        return Convert.ToHexString(SHA256.HashData(stream.ToArray()))[..16].ToLowerInvariant();
    }
}

/// <summary>
/// Writes <c>DifferentialResults.metadata.json</c> (<c>DEF-DIFF-META</c>): UTF-8 without a byte-order mark, indented, keys
/// in a fixed order, <c>\n</c> line ends and a final newline, no timestamps, so the same analysis gives the same bytes.
/// Optional values that are absent are left out; a non-finite number is written as <c>null</c>. Built key by key with
/// <see cref="Utf8JsonWriter"/>, as <c>OrthologySnapshotWriter</c> writes its manifest.
/// </summary>
public static class DifferentialMetadataWriter
{
    /// <summary>The metadata file's bytes for <paramref name="metadata"/> and the rows of its results table.</summary>
    /// <exception cref="ArgumentException">A row belongs to another analysis, or breaks the rows' contract; a global shift,
    /// model or normalization entry names a stratum the metadata does not list, or appears twice for its key; a global shift
    /// sits on a contrast without a numerator and a denominator; or a prior names an unknown estimator, or is trended with
    /// one variance, or constant without one.</exception>
    public static byte[] ToBytes(DifferentialMetadata metadata, IReadOnlyList<DifferentialResult> rows)
    {
        ArgumentNullException.ThrowIfNull(metadata);
        var list = DifferentialResultWriter.Validate(rows);
        string id = metadata.AnalysisId;
        var foreign = list.FirstOrDefault(r => r.AnalysisId != id);
        if (foreign is not null)
            throw new ArgumentException(
                $"Row for feature '{foreign.FeatureId}' belongs to analysis '{foreign.AnalysisId}', not '{id}'.", nameof(rows));
        ValidateMetadata(metadata);
        // The same cleaning the results table applies, so the count cannot disagree with that file.
        int cleaned = DifferentialResultWriter.Write(TextWriter.Null, list, DifferentialHeaderStyle.Machine);

        using var stream = new MemoryStream();
        using (var w = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = true, NewLine = "\n" }))
        {
            w.WriteStartObject();
            w.WriteString("definition_version", DifferentialMetadata.Version);
            w.WriteString("analysis_id", id);
            w.WriteStartArray("software");
            foreach (var s in metadata.Software)
            {
                w.WriteStartObject();
                w.WriteString("name", s.Name);
                w.WriteString("version", s.Version);
                w.WriteEndObject();
            }
            w.WriteEndArray();
            WriteInputs(w, metadata.Inputs);
            WriteSettings(w, metadata.Settings);

            w.WriteStartArray("columns");
            foreach (var c in DifferentialColumns.All)
            {
                w.WriteStartObject();
                w.WriteString("machine", c.MachineName);
                w.WriteString("human", c.HumanName);
                w.WriteString("group", c.Group);
                w.WriteEndObject();
            }
            w.WriteEndArray();

            w.WriteStartArray("strata");
            foreach (var s in metadata.Strata)
            {
                w.WriteStartObject();
                w.WriteString("stratum", s.Stratum);
                Map(w, "factors", s.Factors, (k, v) => w.WriteString(k, v));
                w.WriteString("source", s.Source);
                Map(w, "samples_per_level", s.SamplesPerLevel, (k, v) => w.WriteNumber(k, v));
                if (s.IndividualsPerLevel is { } individuals) Map(w, "individuals_per_level", individuals, (k, v) => w.WriteNumber(k, v));
                if (s.Batches is { } batches) Strings(w, "batches", batches);
                if (s.Fractions is { } fractions) w.WriteNumber("fractions", fractions);
                if (s.TechnicalReplicates is { } techReps) w.WriteNumber("technical_replicates", techReps);
                if (s.ReferenceLevels is { Count: > 0 } references) Map(w, "reference_levels", references, (k, v) => w.WriteString(k, v));
                if (s.SamplesWithoutValues is { } without) Strings(w, "samples_without_values", without);
                w.WriteStartArray("not_run");
                foreach (var n in s.NotRun)
                {
                    w.WriteStartObject();
                    w.WriteString("label", n.Label);
                    w.WriteString("reason", n.Reason);
                    w.WriteEndObject();
                }
                w.WriteEndArray();
                w.WriteEndObject();
            }
            w.WriteEndArray();

            w.WriteStartArray("contrasts");
            foreach (var c in metadata.Contrasts)
            {
                w.WriteStartObject();
                w.WriteString("id", c.Id);
                w.WriteString("label", c.Label);
                Optional(w, "numerator", c.Numerator);
                Optional(w, "denominator", c.Denominator);
                Optional(w, "covariate", c.Covariate);
                Optional(w, "covariate_unit", c.CovariateUnit);
                Optional(w, "covariate_scale", c.CovariateScale);
                Map(w, "weights", c.Weights, (k, v) => Number(w, k, v));
                if (c.GlobalShift is { Count: > 0 } shifts)
                {
                    w.WriteStartArray("global_shift");
                    foreach (var shift in shifts.OrderBy(s => s.Stratum, StringComparer.Ordinal)
                                 .ThenBy(s => s.QuantBasis ?? "", StringComparer.Ordinal))
                    {
                        w.WriteStartObject();
                        w.WriteString("stratum", shift.Stratum);
                        Optional(w, "quant_basis", shift.QuantBasis);
                        Number(w, "value", shift.Value);
                        w.WriteNumber("peptides", shift.Peptides);
                        w.WriteEndObject();
                    }
                    w.WriteEndArray();
                }
                w.WriteEndObject();
            }
            w.WriteEndArray();

            w.WriteStartArray("models");
            foreach (var m in metadata.Models.OrderBy(m => m.Stratum, StringComparer.Ordinal)
                         .ThenBy(m => m.QuantBasis ?? "", StringComparer.Ordinal).ThenBy(m => m.Method, StringComparer.Ordinal)
                         .ThenBy(m => m.ModelUsed, StringComparer.Ordinal).ThenBy(m => m.Grain, StringComparer.Ordinal))
            {
                w.WriteStartObject();
                w.WriteString("stratum", m.Stratum);
                Optional(w, "quant_basis", m.QuantBasis);
                w.WriteString("method", m.Method);
                w.WriteString("model_used", m.ModelUsed);
                w.WriteString("grain", m.Grain);
                w.WriteString("formula", m.Formula);
                w.WriteBoolean("reml", m.Reml);
                w.WriteString("df_method", m.DfMethod);
                Optional(w, "robust_rule", m.RobustRule);
                if (m.Prior is { } prior)
                {
                    w.WriteStartObject("prior");
                    w.WriteString("estimator", prior.Estimator);
                    Number(w, "df", prior.Df);
                    w.WriteBoolean("trended", prior.Trended);
                    if (prior.Variance is { } s0) Number(w, "variance", s0);
                    w.WriteNumber("features", prior.Features);
                    w.WriteEndObject();
                }
                w.WriteEndObject();
            }
            w.WriteEndArray();

            w.WriteStartArray("normalization");
            foreach (var n in metadata.Normalization.OrderBy(n => n.Stratum, StringComparer.Ordinal)
                         .ThenBy(n => n.QuantBasis ?? "", StringComparer.Ordinal))
            {
                w.WriteStartObject();
                w.WriteString("stratum", n.Stratum);
                Optional(w, "quant_basis", n.QuantBasis);
                w.WriteString("setting", n.Setting);
                w.WriteNumber("reference_set_size", n.ReferenceSetSize);
                w.WriteNumber("minimum_reference_set_size", n.MinimumReferenceSetSize);
                Map(w, "per_sample_shift", n.PerSampleShift, (k, v) => Number(w, k, v));
                w.WriteStartArray("warnings");
                foreach (var warning in n.Warnings) w.WriteStringValue(warning);
                w.WriteEndArray();
                w.WriteEndObject();
            }
            w.WriteEndArray();

            w.WriteStartArray("families");
            foreach (var f in list.Where(r => DifferentialStatus.IsFitted(r.Status))
                         .GroupBy(DifferentialResultWriter.FamilyKey)
                         .OrderBy(g => g.Key.Stratum, StringComparer.Ordinal).ThenBy(g => g.Key.Grain, StringComparer.Ordinal)
                         .ThenBy(g => g.Key.Quantity, StringComparer.Ordinal).ThenBy(g => g.Key.Basis, StringComparer.Ordinal)
                         .ThenBy(g => g.Key.Contrast, StringComparer.Ordinal).ThenBy(g => g.Key.Method, StringComparer.Ordinal))
            {
                w.WriteStartObject();
                w.WriteString("stratum", f.Key.Stratum);
                w.WriteString("grain", f.Key.Grain);
                w.WriteString("quantity", f.Key.Quantity);
                if (f.Key.Basis.Length > 0) w.WriteString("quant_basis", f.Key.Basis);
                w.WriteString("contrast_id", f.Key.Contrast);
                w.WriteString("method", f.Key.Method);
                w.WriteNumber("size", f.Count());
                w.WriteEndObject();
            }
            w.WriteEndArray();

            w.WriteStartArray("design_warnings");
            foreach (var warning in metadata.DesignWarnings) w.WriteStringValue(warning);
            w.WriteEndArray();

            w.WriteStartObject("row_count");
            w.WriteNumber("total", list.Count);
            w.WriteStartObject("by_status");
            foreach (var g in list.GroupBy(r => r.Status).OrderBy(g => g.Key, StringComparer.Ordinal))
                w.WriteNumber(g.Key, g.Count());
            w.WriteEndObject();
            w.WriteEndObject();

            w.WriteNumber("cleaned_text_cells", cleaned);

            w.WriteEndObject();
        }
        stream.WriteByte((byte)'\n');
        return stream.ToArray();
    }

    /// <summary>
    /// Every global shift, model and normalization entry names a listed stratum and appears once for its key; a global
    /// shift sits on a contrast with a numerator and a denominator; a prior names a known estimator and has one variance
    /// exactly when it is constant.
    /// </summary>
    private static void ValidateMetadata(DifferentialMetadata metadata)
    {
        var strata = metadata.Strata.Select(s => s.Stratum).ToHashSet(StringComparer.Ordinal);
        void Listed(string stratum, string what)
        {
            if (!strata.Contains(stratum))
                throw new ArgumentException($"{what} names stratum '{stratum}', which the metadata does not list.", nameof(metadata));
        }

        foreach (var c in metadata.Contrasts)
        {
            if (c.GlobalShift is not { Count: > 0 } shifts)
                continue;
            if (c.Numerator is null || c.Denominator is null)
                throw new ArgumentException(
                    $"Contrast '{c.Id}' has a global shift but no numerator and denominator to measure it between.", nameof(metadata));
            var seen = new HashSet<(string, string)>();
            foreach (var shift in shifts)
            {
                Listed(shift.Stratum, $"A global shift of contrast '{c.Id}'");
                if (!seen.Add((shift.Stratum, shift.QuantBasis ?? "")))
                    throw new ArgumentException($"Contrast '{c.Id}' has a global shift for stratum '{shift.Stratum}' and quant " +
                                                $"basis '{shift.QuantBasis}' twice.", nameof(metadata));
            }
        }

        var models = new HashSet<(string, string, string, string, string)>();
        foreach (var m in metadata.Models)
        {
            string at = $"The model for stratum '{m.Stratum}', quant basis '{m.QuantBasis}', method '{m.Method}', model used " +
                        $"'{m.ModelUsed}' and grain '{m.Grain}'";
            Listed(m.Stratum, at);
            if (!models.Add((m.Stratum, m.QuantBasis ?? "", m.Method, m.ModelUsed, m.Grain)))
                throw new ArgumentException($"{at} appears twice.", nameof(metadata));
            if (m.Prior is not { } prior)
                continue;
            if (prior.Estimator is not (DifferentialPriorInfo.MomentsLegacy or DifferentialPriorInfo.MarginalLikelihood))
                throw new ArgumentException($"{at} has prior estimator '{prior.Estimator}', which is neither " +
                                            $"'{DifferentialPriorInfo.MomentsLegacy}' nor '{DifferentialPriorInfo.MarginalLikelihood}'.",
                    nameof(metadata));
            if (prior.Trended && prior.Variance is not null)
                throw new ArgumentException($"{at} has a trended prior with one variance; a trended prior has a variance per feature.",
                    nameof(metadata));
            if (!prior.Trended && prior.Variance is null)
                throw new ArgumentException($"{at} has a constant prior without its variance.", nameof(metadata));
        }

        var normalization = new HashSet<(string, string)>();
        foreach (var n in metadata.Normalization)
        {
            string at = $"The normalization for stratum '{n.Stratum}' and quant basis '{n.QuantBasis}'";
            Listed(n.Stratum, at);
            if (!normalization.Add((n.Stratum, n.QuantBasis ?? "")))
                throw new ArgumentException($"{at} appears twice.", nameof(metadata));
        }
    }

    /// <summary>A list of strings as an array in ordinal order.</summary>
    private static void Strings(Utf8JsonWriter w, string name, IEnumerable<string> values)
    {
        w.WriteStartArray(name);
        foreach (var v in values.OrderBy(v => v, StringComparer.Ordinal)) w.WriteStringValue(v);
        w.WriteEndArray();
    }

    internal static void WriteInputs(Utf8JsonWriter w, IReadOnlyList<DifferentialInput> inputs)
    {
        w.WriteStartArray("inputs");
        foreach (var i in inputs)
        {
            w.WriteStartObject();
            w.WriteString("role", i.Role);
            w.WriteString("file", i.FileName);
            w.WriteString("sha256", i.Sha256);
            w.WriteEndObject();
        }
        w.WriteEndArray();
    }

    internal static void WriteSettings(Utf8JsonWriter w, DifferentialSettings s)
    {
        w.WriteStartObject("settings");
        w.WriteString("header_style", s.HeaderStyle == DifferentialHeaderStyle.Machine ? "machine" : "human");
        Number(w, "ci_level", s.CiLevel);
        w.WriteString("method", s.Method);
        w.WriteBoolean("robust", s.Robust);
        w.WriteString("normalization", s.Normalization);
        Optional(w, "normalization_override", s.NormalizationOverride);
        if (s.Seed is { } seed) w.WriteNumber("seed", seed);
        w.WriteEndObject();
    }

    private static void Optional(Utf8JsonWriter w, string name, string? value)
    {
        if (value is not null) w.WriteString(name, value);
    }

    /// <summary>A finite number at round-trip precision; a non-finite one as <c>null</c>.</summary>
    private static void Number(Utf8JsonWriter w, string name, double value)
    {
        if (double.IsFinite(value)) w.WriteNumber(name, value);
        else w.WriteNull(name);
    }

    /// <summary>A dictionary as an object with its keys in ordinal order.</summary>
    private static void Map<T>(Utf8JsonWriter w, string name, IReadOnlyDictionary<string, T> map, Action<string, T> write)
    {
        w.WriteStartObject(name);
        foreach (var (k, v) in map.OrderBy(kv => kv.Key, StringComparer.Ordinal)) write(k, v);
        w.WriteEndObject();
    }
}
