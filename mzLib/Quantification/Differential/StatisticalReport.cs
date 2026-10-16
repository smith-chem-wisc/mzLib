using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace Quantification.Differential;

/// <summary>
/// One fit's per-feature values for the mean-variance diagnostic, as limma's <c>plotSA</c> draws it: each feature's average
/// log2 intensity, its residual variance (s^2) and the empirical-Bayes prior variance it was moderated toward (s0^2; the
/// same value for every feature under a constant prior). The four lists are aligned on <see cref="FeatureIds"/>.
/// </summary>
public sealed record DiagnosticFit(string Stratum, string? QuantBasis, string Method, string ModelUsed, string Grain,
    IReadOnlyList<string> FeatureIds, IReadOnlyList<double> AverageLog2, IReadOnlyList<double> ResidualVariance,
    IReadOnlyList<double> PriorVariance);

/// <summary>
/// One stratum's per-sample log2 values in one quant basis, for the sample-correlation diagnostic. <see cref="Values"/> has
/// one list per sample, in <see cref="SampleIds"/> order, all aligned on the same features; NaN is a missing value.
/// </summary>
public sealed record SampleLog2(string Stratum, string? QuantBasis, IReadOnlyList<string> SampleIds,
    IReadOnlyList<IReadOnlyList<double>> Values);

/// <summary>The per-feature and per-sample values the report's diagnostics need beyond the metadata and the rows.</summary>
public sealed record DifferentialDiagnostics
{
    /// <summary>One entry per stratum, quant basis, method, model used and grain.</summary>
    public IReadOnlyList<DiagnosticFit> Fits { get; init; } = Array.Empty<DiagnosticFit>();
    /// <summary>One entry per stratum and quant basis.</summary>
    public IReadOnlyList<SampleLog2> Samples { get; init; } = Array.Empty<SampleLog2>();
}

/// <summary>A built report: the markdown and its figure files (SVG and TSV), in report order.</summary>
internal sealed record StatisticalReport(string Markdown, IReadOnlyList<ReportFile> Figures);

/// <summary>A file of the report's figures directory.</summary>
internal sealed record ReportFile(string Name, byte[] Bytes);

/// <summary>
/// Writes <c>StatisticalReport.md</c> (STAT1 milestone M6; QuantProject design/STATS-FRAMEWORK.md section 6), generated
/// from <see cref="DifferentialMetadata"/>, the <see cref="DifferentialResult"/> rows and <see cref="DifferentialDiagnostics"/>:
/// a methods paragraph ready to paste, a design summary, a results summary, diagnostic figures (SVG drawn by our own code,
/// each with the table of its data) and a journal checklist whose every item is filled from the run, a gap shown as N.
/// </summary>
/// <remarks>
/// Everything is built in memory and validated before any file is written. Like <see cref="DifferentialResultWriter"/>,
/// it writes UTF-8 without a byte-order mark, <c>\n</c> line ends, invariant-culture numbers and no timestamps, so the
/// same analysis gives the same bytes on every machine.
/// </remarks>
public static class StatisticalReportWriter
{
    /// <summary>The report's file name.</summary>
    public const string ReportFileName = "StatisticalReport.md";
    /// <summary>The directory, beside the report, that holds the figures and their tables.</summary>
    public const string FiguresDirectoryName = "figures";

    /// <summary>
    /// Writes <see cref="ReportFileName"/> and the <see cref="FiguresDirectoryName"/> directory into
    /// <paramref name="directory"/>, creating them, and returns the report's path.
    /// </summary>
    /// <exception cref="ArgumentException">The rows or metadata break their contract (see
    /// <see cref="DifferentialMetadataWriter.ToBytes"/>), or the diagnostics are inconsistent: a list of the wrong length, a
    /// stratum the metadata does not list, or an entry twice for its key.</exception>
    public static string Write(string directory, DifferentialMetadata metadata, IReadOnlyList<DifferentialResult> rows,
        DifferentialDiagnostics diagnostics)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(directory);
        var report = Build(metadata, rows, diagnostics);
        string figures = Path.Combine(directory, FiguresDirectoryName);
        Directory.CreateDirectory(figures);
        foreach (var file in report.Figures)
            File.WriteAllBytes(Path.Combine(figures, file.Name), file.Bytes);
        string path = Path.Combine(directory, ReportFileName);
        File.WriteAllBytes(path, new UTF8Encoding(false).GetBytes(report.Markdown));
        return path;
    }

    internal static StatisticalReport Build(DifferentialMetadata metadata, IReadOnlyList<DifferentialResult> rows,
        DifferentialDiagnostics diagnostics)
    {
        ArgumentNullException.ThrowIfNull(metadata);
        ArgumentNullException.ThrowIfNull(rows);
        ArgumentNullException.ThrowIfNull(diagnostics);
        _ = DifferentialMetadataWriter.ToBytes(metadata, rows);
        Validate(metadata, diagnostics);

        var ordered = DifferentialResultWriter.Order(rows);
        var md = new StringBuilder();
        md.Append("# Statistical report\n\n");
        md.Append($"Analysis `{metadata.AnalysisId}` (`{DifferentialMetadata.Version}`). Generated from `")
            .Append(DifferentialResultWriter.MetadataFileName).Append("` and `").Append(DifferentialResultWriter.ResultsFileName)
            .Append("`; nothing here is written by hand.\n\n");
        Methods(md, metadata, ordered);
        Design(md, metadata);
        Results(md, ordered);
        var figures = Diagnostics(md, metadata, ordered, diagnostics);
        Checklist(md, metadata, ordered);
        return new StatisticalReport(md.ToString(), figures);
    }

    private static void Validate(DifferentialMetadata metadata, DifferentialDiagnostics diagnostics)
    {
        var strata = metadata.Strata.Select(s => s.Stratum).ToHashSet(StringComparer.Ordinal);
        var fits = new HashSet<(string, string, string, string, string)>();
        foreach (var f in diagnostics.Fits)
        {
            string at = $"The diagnostic fit for stratum '{f.Stratum}', quant basis '{f.QuantBasis}', method '{f.Method}', model used " +
                        $"'{f.ModelUsed}' and grain '{f.Grain}'";
            if (!strata.Contains(f.Stratum))
                throw new ArgumentException($"{at} names stratum '{f.Stratum}', which the metadata does not list.", nameof(diagnostics));
            if (!fits.Add((f.Stratum, f.QuantBasis ?? "", f.Method, f.ModelUsed, f.Grain)))
                throw new ArgumentException($"{at} appears twice.", nameof(diagnostics));
            int n = f.FeatureIds.Count;
            if (f.AverageLog2.Count != n || f.ResidualVariance.Count != n || f.PriorVariance.Count != n)
                throw new ArgumentException($"{at} has {n} feature ids but {f.AverageLog2.Count} average intensities, " +
                                            $"{f.ResidualVariance.Count} residual variances and {f.PriorVariance.Count} prior variances.",
                    nameof(diagnostics));
        }
        var samples = new HashSet<(string, string)>();
        foreach (var s in diagnostics.Samples)
        {
            string at = $"The sample matrix for stratum '{s.Stratum}' and quant basis '{s.QuantBasis}'";
            if (!strata.Contains(s.Stratum))
                throw new ArgumentException($"{at} names stratum '{s.Stratum}', which the metadata does not list.", nameof(diagnostics));
            if (!samples.Add((s.Stratum, s.QuantBasis ?? "")))
                throw new ArgumentException($"{at} appears twice.", nameof(diagnostics));
            if (s.Values.Count != s.SampleIds.Count)
                throw new ArgumentException($"{at} has {s.SampleIds.Count} samples but {s.Values.Count} value lists.", nameof(diagnostics));
            if (s.SampleIds.Distinct(StringComparer.Ordinal).Count() != s.SampleIds.Count)
                throw new ArgumentException($"{at} names a sample twice.", nameof(diagnostics));
            if (s.Values.Select(v => v.Count).Distinct().Count() > 1)
                throw new ArgumentException($"{at} has samples with different numbers of features.", nameof(diagnostics));
        }
    }

    private static void Methods(StringBuilder md, DifferentialMetadata m, IReadOnlyList<DifferentialResult> rows)
    {
        md.Append("## 1. Methods\n\n");
        var versions = rows.Select(r => r.MethodVersion).Where(v => !string.IsNullOrEmpty(v)).Distinct(StringComparer.Ordinal)
            .OrderBy(v => v, StringComparer.Ordinal).ToList();
        md.Append("**Software.** ")
            .Append(m.Software.Count == 0 ? "Not recorded." : string.Join("; ", m.Software.Select(s => $"{s.Name} {s.Version}")) + ".")
            .Append(versions.Count == 0 ? "" : $" Method version: {string.Join(", ", versions)}.").Append("\n\n");
        if (m.Inputs.Count > 0)
            md.Append("**Inputs.** ")
                .Append(string.Join("; ", m.Inputs.Select(i => $"{i.Role} `{i.FileName}` (sha256 `{Prefix(i.Sha256)}`)"))).Append(".\n\n");

        md.Append("**Transformation and normalization.** All quantities were analysed on the log2 scale.");
        if (m.Settings.NormalizationOverride is { } over)
            md.Append($" A curator overrode the default normalization: `{over}`.");
        if (m.Normalization.Count == 0)
            md.Append($" Normalization `{m.Settings.Normalization}`; no per-stratum details were recorded.");
        foreach (var n in m.Normalization.OrderBy(n => n.Stratum, StringComparer.Ordinal).ThenBy(n => n.QuantBasis ?? "", StringComparer.Ordinal))
        {
            md.Append($" Stratum `{n.Stratum}`{Basis(n.QuantBasis)}: ").Append(NormalizationSentence(n));
            if (n.Warnings.Count > 0) md.Append(" Warnings: ").Append(string.Join("; ", n.Warnings)).Append('.');
        }
        if (m.Normalization.Count > 0) md.Append(" The per-sample shifts are in the metadata file.");
        md.Append("\n\n");

        md.Append("**Missing values.** Missing values were not imputed: each model was fitted on the observed values only. A ")
            .Append("missing value is an empty cell in the results file, never 0.\n\n");

        md.Append("**Models.**");
        if (m.Models.Count == 0) md.Append(" No model was recorded.");
        md.Append('\n');
        foreach (var model in m.Models.OrderBy(x => x.Stratum, StringComparer.Ordinal).ThenBy(x => x.QuantBasis ?? "", StringComparer.Ordinal)
                     .ThenBy(x => x.Method, StringComparer.Ordinal).ThenBy(x => x.ModelUsed, StringComparer.Ordinal)
                     .ThenBy(x => x.Grain, StringComparer.Ordinal))
        {
            md.Append($"- Stratum `{model.Stratum}`{Basis(model.QuantBasis)}, `{model.Grain}`, `{model.ModelUsed}` (method ")
                .Append($"`{model.Method}`): `{model.Formula}`, {(model.Reml ? "fitted by REML" : "fitted without REML")}; ")
                .Append($"degrees of freedom by `{model.DfMethod}`");
            if (model.RobustRule is { } robust) md.Append($"; robust weighting: {robust}");
            if (model.Prior is { } p)
            {
                md.Append(p.Trended
                    ? $"; a trended empirical-Bayes variance prior (`{p.Estimator}`) with {PriorDf(p.Df)} degrees of freedom, its variance " +
                      "depending on the average intensity (see the mean-variance figure)"
                    : $"; an empirical-Bayes variance prior (`{p.Estimator}`) with {PriorDf(p.Df)} degrees of freedom and variance " +
                      $"{Num(p.Variance ?? double.NaN)}");
                md.Append($", fitted on {Int(p.Features)} features");
            }
            md.Append(".\n");
        }
        md.Append('\n');
        if (m.Settings.Robust) md.Append("Outlying values were down-weighted (robust weighting on).\n\n");

        md.Append($"**Tests.** Every comparison is tested two-sided. Each estimate is reported with its standard error, a ")
            .Append($"{Num(m.Settings.CiLevel * 100)}% confidence interval (t-based, on the stated degrees of freedom), its test ")
            .Append("statistic, its degrees of freedom and an exact p-value.\n\n");

        var fitted = rows.Where(r => DifferentialStatus.IsFitted(r.Status)).ToList();
        var adjustments = fitted.Select(r => r.AdjustmentMethod).Where(a => a is not null).Distinct(StringComparer.Ordinal).ToList();
        int families = fitted.GroupBy(DifferentialResultWriter.FamilyKey).Count();
        md.Append("**Multiple testing.** ").Append(adjustments.Count == 0
            ? "No multiple-testing adjustment was recorded."
            : $"p-values were adjusted with the {string.Join(", ", adjustments.Select(AdjustmentName))} method within each family of " +
              $"features sharing stratum, grain, quantity, quant basis, comparison and method: {Int(families)} " +
              $"{(families == 1 ? "family" : "families")}, their sizes in the metadata file.").Append("\n\n");

        md.Append("**Thresholds.** No significance threshold was applied and no feature is labelled significant; the counts at ")
            .Append("adjusted p < 0.05 in section 3 are descriptive only (ASA statements, 2016 and 2019).\n\n");
        if (m.Settings.Seed is { } seed)
            md.Append($"**Seed.** Random seed: {seed.ToString(CultureInfo.InvariantCulture)}.\n\n");
    }

    private static string NormalizationSentence(DifferentialNormalizationInfo n) => n.Setting switch
    {
        "shared_peptide_median" =>
            $"each sample was shifted by its median log2 difference from the across-sample reference, over the {Int(n.ReferenceSetSize)} " +
            $"peptides with a value in every sample that has values (at least {Int(n.MinimumReferenceSetSize)} required).",
        "shared_peptide_median_half" =>
            $"each sample was shifted by its median log2 difference from the across-sample reference, over the {Int(n.ReferenceSetSize)} " +
            $"peptides with a value in at least half the samples that have values, because fewer than {Int(n.MinimumReferenceSetSize)} " +
            "had a value in every one; the differences were measured again on the shifted values until the shifts stopped changing, " +
            "and the shifts were centred to average 0.",
        "none" => "not normalized.",
        _ when n.Setting.StartsWith("background_set:", StringComparison.Ordinal) =>
            $"each sample was shifted on the declared background set `{n.Setting["background_set:".Length..]}` ({Int(n.ReferenceSetSize)} peptides).",
        _ => $"normalization `{n.Setting}` over {Int(n.ReferenceSetSize)} reference peptides.",
    };

    private static void Design(StringBuilder md, DifferentialMetadata m)
    {
        md.Append("## 2. Design\n\n");
        foreach (var s in m.Strata)
        {
            md.Append($"### Stratum `{s.Stratum}`\n\n");
            md.Append("Factors: ").Append(s.Factors.Count == 0
                    ? "none (the whole dataset)"
                    : string.Join(", ", s.Factors.OrderBy(f => f.Key, StringComparer.Ordinal).Select(f => $"{f.Key} = {f.Value}")))
                .Append($". Chosen by: {s.Source}.\n\n");
            md.Append("| level | samples | individuals |\n|---|---|---|\n");
            foreach (var (level, count) in s.SamplesPerLevel.OrderBy(kv => kv.Key, StringComparer.Ordinal))
            {
                string individuals = s.IndividualsPerLevel is { } ind && ind.TryGetValue(level, out int k) ? Int(k) : "not recorded";
                md.Append($"| {Cell(level)} | {Int(count)} | {individuals} |\n");
            }
            md.Append('\n');
            md.Append("- Batches: ").Append(s.Batches is null ? "not recorded" : s.Batches.Count == 0 ? "none"
                : string.Join(", ", s.Batches.OrderBy(b => b, StringComparer.Ordinal))).Append(".\n");
            md.Append("- Fractions: ").Append(s.Fractions is { } fr ? Int(fr) : "not recorded")
                .Append(". Technical replicates: ").Append(s.TechnicalReplicates is { } tr ? Int(tr) : "not recorded").Append(".\n");
            md.Append(s.ReferenceLevels is { Count: > 0 } refs
                ? "- Reference level: " + string.Join("; ", refs.OrderBy(r => r.Key, StringComparer.Ordinal).Select(r => $"{r.Key} = {r.Value}")) + ".\n"
                : "- No curated reference level: every pair of levels is compared.\n");
            md.Append("- Samples without values: ").Append(s.SamplesWithoutValues is null ? "not recorded" : s.SamplesWithoutValues.Count == 0
                ? "none" : string.Join(", ", s.SamplesWithoutValues.OrderBy(x => x, StringComparer.Ordinal))).Append(".\n");
            md.Append("- Comparisons not run here: ").Append(s.NotRun.Count == 0 ? "none"
                : string.Join("; ", s.NotRun.Select(n => $"{n.Label} ({n.Reason})"))).Append(".\n\n");
        }
        md.Append("### Comparisons\n\n| id | comparison | numerator | denominator | covariate |\n|---|---|---|---|---|\n");
        foreach (var c in m.Contrasts)
        {
            string covariate = c.Covariate is null ? "" :
                c.Covariate + (c.CovariateUnit is null && c.CovariateScale is null ? "" :
                    $" ({string.Join("; ", new[] { c.CovariateUnit, c.CovariateScale }.Where(x => x is not null))})");
            md.Append($"| {Cell(c.Id)} | {Cell(c.Label)} | {Cell(c.Numerator ?? "")} | {Cell(c.Denominator ?? "")} | {Cell(covariate)} |\n");
        }
        md.Append("\n### Design warnings\n\n");
        md.Append(m.DesignWarnings.Count == 0 ? "None.\n" : string.Concat(m.DesignWarnings.Select(w => $"- {w}\n")));
        md.Append('\n');
    }

    private static void Results(StringBuilder md, IReadOnlyList<DifferentialResult> rows)
    {
        md.Append("## 3. Results\n\n");
        md.Append("| stratum | basis | grain | quantity | comparison | method | features | fitted | not fitted | adjusted p < 0.05 |\n");
        md.Append("|---|---|---|---|---|---|---|---|---|---|\n");
        var blocks = rows.GroupBy(DifferentialResultWriter.FamilyKey).ToList();
        foreach (var b in blocks)
        {
            int fitted = b.Count(r => DifferentialStatus.IsFitted(r.Status));
            int below = b.Count(r => DifferentialStatus.IsFitted(r.Status) && r.PAdjusted is { } q && q < 0.05);
            md.Append($"| {Cell(b.Key.Stratum)} | {BasisCell(b.Key.Basis)} | {Cell(b.Key.Grain)} | {Cell(b.Key.Quantity)} | ")
                .Append($"{Cell(b.Key.Contrast)} | {Cell(b.Key.Method)} | {Int(b.Count())} | {Int(fitted)} | {Int(b.Count() - fitted)} | {Int(below)} |\n");
        }
        md.Append("\nNot fitted, by reason:\n\n| stratum | basis | comparison | method | status | features |\n|---|---|---|---|---|---|\n");
        foreach (var b in blocks)
            foreach (var s in b.Where(r => !DifferentialStatus.IsFitted(r.Status)).GroupBy(r => r.Status).OrderBy(g => g.Key, StringComparer.Ordinal))
                md.Append($"| {Cell(b.Key.Stratum)} | {BasisCell(b.Key.Basis)} | {Cell(b.Key.Contrast)} | {Cell(b.Key.Method)} | {Cell(s.Key)} | {Int(s.Count())} |\n");
        md.Append("\nCounts at adjusted p < 0.05 are descriptive only; no feature is labelled significant.\n\n");
    }

    private static IReadOnlyList<ReportFile> Diagnostics(StringBuilder md, DifferentialMetadata m, IReadOnlyList<DifferentialResult> rows,
        DifferentialDiagnostics diagnostics)
    {
        md.Append("## 4. Diagnostics\n\n");
        var files = new List<ReportFile>();
        int number = 0;
        void Add(string kind, string heading, (byte[] Svg, byte[] Tsv) figure, string caption)
        {
            number++;
            string stem = $"figure_{number.ToString("00", CultureInfo.InvariantCulture)}_{kind}";
            files.Add(new ReportFile(stem + ".svg", figure.Svg));
            files.Add(new ReportFile(stem + ".tsv", figure.Tsv));
            md.Append($"### Figure {number.ToString(CultureInfo.InvariantCulture)}. {heading}\n\n")
                .Append($"![Figure {number.ToString(CultureInfo.InvariantCulture)}]({FiguresDirectoryName}/{stem}.svg)\n\n")
                .Append(caption).Append($" Data: [`{stem}.tsv`]({FiguresDirectoryName}/{stem}.tsv).\n\n");
        }

        foreach (var family in rows.Where(r => DifferentialStatus.IsFitted(r.Status)).GroupBy(DifferentialResultWriter.FamilyKey))
        {
            var k = family.Key;
            string where = $"stratum `{k.Stratum}`{Basis(k.Basis.Length == 0 ? null : k.Basis)}, `{k.Grain}` `{k.Quantity}`, comparison `{k.Contrast}`, method `{k.Method}`";
            string shortTitle = $"{k.Stratum}{(k.Basis.Length == 0 ? "" : ", " + k.Basis)}, {k.Contrast}";
            var list = family.ToList();
            Add("pvalue_histogram", $"p-value histogram: {where}",
                DiagnosticFigures.PValueHistogram($"p-values: {shortTitle}", list.Select(r => r.PValue ?? double.NaN).ToList()),
                $"p-values of the {Int(list.Count)} fitted features in 20 bins of width 0.05 (bin k holds 0.05k ≤ p < 0.05(k+1); p = 1 is " +
                "in the last bin). With no change anywhere the bars are flat; a spike near 0 over a flat floor is expected when some " +
                "features change; a U shape, or a slope rising toward 1, signals a problem with the model or the data.");
            int zero = list.Count(r => r.PValue == 0);
            Add("volcano", $"Volcano: {where}",
                DiagnosticFigures.Volcano($"Volcano: {shortTitle}", list.Select(r => (r.FeatureId, r.Log2Effect ?? double.NaN, r.PValue ?? double.NaN)).ToList()),
                $"log2 effect against -log10(p) for the {Int(list.Count)} fitted features; no threshold lines are drawn." +
                (zero == 0 ? "" : $" {Int(zero)} {(zero == 1 ? "feature with p = 0 is" : "features with p = 0 are")} drawn as " +
                                  $"{(zero == 1 ? "a triangle" : "triangles")} at the top edge."));
        }

        foreach (var fit in diagnostics.Fits.OrderBy(f => f.Stratum, StringComparer.Ordinal).ThenBy(f => f.QuantBasis ?? "", StringComparer.Ordinal)
                     .ThenBy(f => f.Method, StringComparer.Ordinal).ThenBy(f => f.ModelUsed, StringComparer.Ordinal).ThenBy(f => f.Grain, StringComparer.Ordinal))
        {
            var prior = m.Models.FirstOrDefault(x => x.Stratum == fit.Stratum && x.QuantBasis == fit.QuantBasis && x.Method == fit.Method &&
                                                     x.ModelUsed == fit.ModelUsed && x.Grain == fit.Grain)?.Prior;
            int points = DiagnosticFigures.MeanVariancePoints(fit).Count;
            Add("mean_variance", $"Mean-variance: stratum `{fit.Stratum}`{Basis(fit.QuantBasis)}, `{fit.Grain}`, `{fit.ModelUsed}`",
                DiagnosticFigures.MeanVariance($"Mean-variance: {fit.Stratum}{(fit.QuantBasis is null ? "" : ", " + fit.QuantBasis)}, {fit.ModelUsed}", fit),
                $"Each feature's sqrt(sigma), the square root of its residual standard deviation, against its average log2 intensity, for " +
                $"{Int(points)} features, as limma's plotSA draws it. The line is the empirical-Bayes prior at the same scale ((s0^2)^(1/4)): " +
                (prior is null ? "the metadata records no prior for this fit." : prior.Trended ? "a curve, because the prior is trended." : "flat, because the prior is constant."));
        }

        foreach (var s in diagnostics.Samples.OrderBy(x => x.Stratum, StringComparer.Ordinal).ThenBy(x => x.QuantBasis ?? "", StringComparer.Ordinal))
            Add("sample_correlation", $"Sample correlation: stratum `{s.Stratum}`{Basis(s.QuantBasis)}",
                DiagnosticFigures.SampleCorrelationHeatmap($"Sample correlation: {s.Stratum}{(s.QuantBasis is null ? "" : ", " + s.QuantBasis)}", s),
                $"Pearson correlation of log2 values between each pair of the {Int(s.SampleIds.Count)} samples, over the features both " +
                "samples have; r needs at least 3 such features, and fewer leave the cell empty. The table gives r and the number of features.");

        var shifts = m.Contrasts.Where(c => c.GlobalShift is { Count: > 0 })
            .SelectMany(c => c.GlobalShift!.OrderBy(g => g.Stratum, StringComparer.Ordinal).ThenBy(g => g.QuantBasis ?? "", StringComparer.Ordinal)
                .Select(g => (c.Id, g))).ToList();
        if (shifts.Count > 0)
            Add("global_shift", "Global shift per comparison, stratum and basis",
                DiagnosticFigures.GlobalShift("Global shift", shifts),
                "Each comparison's global shift: the median, over peptides with a value in every sample on both sides, of the difference " +
                "in mean log2 intensity, before normalization; 0 means the two sides sit level. The number of peptides each rests on is " +
                "beside it.");

        if (number == 0) md.Append("No diagnostic figure: the run has no fitted feature and no diagnostic values.\n\n");
        return files;
    }

    private static void Checklist(StringBuilder md, DifferentialMetadata m, IReadOnlyList<DifferentialResult> rows)
    {
        var fitted = rows.Where(r => DifferentialStatus.IsFitted(r.Status)).ToList();
        static bool Finite(double? v) => v is { } x && double.IsFinite(x);
        var items = new (string Item, string AskedBy, string Where, bool Present)[]
        {
            ("Software and versions", "JPR; PROTEOMICS; MIAPE-Quant 4.1", "1. Methods", m.Software.Count > 0),
            ("Data transformation", "JPR; PROTEOMICS; MIAPE-Quant 4.4", "1. Methods", true),
            ("Normalization", "JPR; PROTEOMICS", "1. Methods", m.Normalization.Count > 0 || m.Settings.Normalization == "none"),
            ("Missing values", "PROTEOMICS; MIAPE-Quant 4.4", "1. Methods", true),
            ("Statistical tests and sidedness", "JPR; PROTEOMICS; Nature reporting summary", "1. Methods", m.Models.Count > 0),
            ("Degrees of freedom", "PROTEOMICS; Nature reporting summary", "1. Methods; results file `df`",
                m.Models.Count > 0 && fitted.All(r => Finite(r.Df))),
            ("Multiple-testing correction", "Nature reporting summary; MIAPE-Quant 4.5", "1. Methods",
                fitted.Count > 0 && fitted.All(r => r.AdjustmentMethod is not null)),
            ("Exact n per group", "Nature reporting summary; MCP", "2. Design",
                m.Strata.Count > 0 && m.Strata.All(s => s.SamplesPerLevel.Count > 0)),
            ("Biological and technical replicates", "MCP; Nature reporting summary", "2. Design",
                m.Strata.Count > 0 && m.Strata.All(s => s.IndividualsPerLevel is not null && s.TechnicalReplicates is not null)),
            ("Effect sizes with confidence intervals", "Nature reporting summary; PROTEOMICS; JPR; MIAPE-Quant 4.5",
                "results file `log2_effect`, `ci_low`, `ci_high`",
                fitted.Count > 0 && fitted.All(r => Finite(r.Log2Effect) && Finite(r.CiLow) && Finite(r.CiHigh))),
            ("Exact p-values", "Nature reporting summary; ASA 2016", "results file `p_value`", fitted.Count > 0 && fitted.All(r => Finite(r.PValue))),
            ("Thresholds and their justification", "MIAPE-Quant 4.5; PROTEOMICS; ASA 2019", "1. Methods", true),
            ("Priors", "Nature reporting summary", "1. Methods", m.Models.Count > 0 && m.Models.All(x => x.Prior is not null)),
            ("Level of each test in nested designs", "Nature reporting summary", "1. Methods (model formulas)",
                m.Models.Count > 0 && m.Models.All(x => !string.IsNullOrWhiteSpace(x.Formula))),
            ("Power analysis", "JPR (where appropriate)", "not computed", false),
        };
        md.Append("## 5. Reporting checklist\n\n| item | asked for by | answered in | this run |\n|---|---|---|---|\n");
        foreach (var (item, askedBy, where, present) in items)
            md.Append($"| {item} | {askedBy} | {where} | {(present ? "Y" : "N")} |\n");
        md.Append("\nAn N is a gap in this run, shown rather than left out.\n");
    }

    private static string Basis(string? basis) => basis is null ? "" : $", basis `{basis}`";

    private static string BasisCell(string basis) => basis.Length == 0 ? "none" : Cell(basis);

    private static string Cell(string text) => text.Replace("|", "\\|").Replace('\n', ' ').Replace('\r', ' ').Replace('\t', ' ');

    private static string Prefix(string sha256) => sha256.Length > 16 ? sha256[..16] + "..." : sha256;

    private static string AdjustmentName(string? method) => method switch
    {
        "benjamini_hochberg" => "Benjamini–Hochberg",
        null => "",
        _ => $"`{method}`",
    };

    /// <summary>
    /// A prior's degrees of freedom. Infinite is a real answer (the variances spread no more than sampling explains, so every
    /// feature takes the prior variance), not a missing one.
    /// </summary>
    private static string PriorDf(double df) => double.IsPositiveInfinity(df) ? "infinite" : Num(df);

    private static string Int(int value) => value.ToString("N0", CultureInfo.InvariantCulture);

    private static string Num(double value)
    {
        if (!double.IsFinite(value)) return "not recorded";
        double rounded = Math.Round(value, 3, MidpointRounding.AwayFromZero);
        if (rounded == 0) rounded = 0;
        return rounded.ToString("0.###", CultureInfo.InvariantCulture);
    }
}
