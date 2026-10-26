using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;

namespace Quantification.Differential;

/// <summary>What a frozen hypothesis predicts for its members (QuantProject <c>DEF-DIFF-SET-VERDICT</c>).</summary>
public enum SetPrediction
{
    /// <summary><c>up</c>: the set's average log2 change is above 0.</summary>
    Up,
    /// <summary><c>down</c>: below 0.</summary>
    Down,
    /// <summary><c>unchanged</c>: within ±δ (an equivalence test).</summary>
    Unchanged,
    /// <summary><c>not_up</c>: below +δ.</summary>
    NotUp,
    /// <summary><c>not_down</c>: above −δ.</summary>
    NotDown,
    /// <summary><c>attenuated_down</c>: falls, by less than f times the members' mRNA fall.</summary>
    AttenuatedDown,
    /// <summary><c>attenuated_up</c>: rises, by less than f times the members' mRNA rise.</summary>
    AttenuatedUp,
}

/// <summary>
/// How a dataset contrast may be pooled (<c>DEF-DIFF-SET-POOL</c>'s reading of <c>independence</c>, RULED 2026-10-19).
/// </summary>
public enum SetIndependence
{
    /// <summary>A value beginning <c>independent</c>. Only these datasets are pooled.</summary>
    Independent,
    /// <summary>A value containing <c>reproduces</c>: reported, never pooled.</summary>
    Reproduction,
    /// <summary>Any other value (<c>weak (not a test)</c>, ...): reported, never pooled.</summary>
    NotPooled,
}

/// <summary>One level set of a <c>contexts</c> cell: a factor and the levels it is restricted to.</summary>
/// <param name="Factor">The factor, as the design names it.</param>
/// <param name="Levels">Its listed levels, in the order written; at least one.</param>
public sealed record SetContext(string Factor, IReadOnlyList<string> Levels);

/// <summary>One row of a frozen <c>tests.tsv</c>: a hypothesis, one member gene, one dataset contrast.</summary>
public sealed record SetTestRow
{
    /// <summary>The row's line in the file (the header is line 1).</summary>
    public required int Line { get; init; }
    /// <summary><c>hypothesis_id</c>, e.g. <c>H3</c>.</summary>
    public required string HypothesisId { get; init; }
    /// <summary><c>kind</c>: <c>single</c>, <c>set</c>, or a label such as <c>concordance</c> (AGING-D9).</summary>
    public required string Kind { get; init; }
    /// <summary><c>prediction</c>.</summary>
    public required SetPrediction Prediction { get; init; }
    /// <summary><c>member</c>: the gene symbol, as frozen.</summary>
    public required string Member { get; init; }
    /// <summary><c>species</c>. Not used to find genes: the dataset's own databases are (<c>DEF-DIFF-SET-MEMBER</c>).</summary>
    public required string Species { get; init; }
    /// <summary><c>dataset</c>, e.g. <c>PXD057407</c>.</summary>
    public required string Dataset { get; init; }
    /// <summary><c>factor</c>: the factor the contrast compares.</summary>
    public required string Factor { get; init; }
    /// <summary><c>test_level</c>: the numerator.</summary>
    public required string TestLevel { get; init; }
    /// <summary><c>reference_level</c>: the denominator.</summary>
    public required string ReferenceLevel { get; init; }
    /// <summary><c>contexts</c> exactly as written; empty means every context of the stratum.</summary>
    public required string ContextsText { get; init; }
    /// <summary><c>contexts</c> parsed (<c>DEF-DIFF-SET-CONTRAST</c>).</summary>
    public required IReadOnlyList<SetContext> Contexts { get; init; }
    /// <summary><c>mrna_log2</c>, or null when empty.</summary>
    public double? MrnaLog2 { get; init; }
    /// <summary><c>independence</c> exactly as written; the reading is <see cref="Independence"/>.</summary>
    public required string IndependenceText { get; init; }
    /// <summary>The row's own reading of <see cref="IndependenceText"/>.</summary>
    public required SetIndependence Independence { get; init; }
    /// <summary><c>post_hoc</c> is <c>yes</c>: this member's result in this dataset was seen before the freeze.</summary>
    public required bool PostHoc { get; init; }
    /// <summary><c>primary</c>: true for <c>yes</c>, false for <c>no</c>, null when the column or the cell is empty.</summary>
    public bool? Primary { get; init; }
    /// <summary><c>margin_log2</c> (δ), or null for the default.</summary>
    public double? MarginLog2 { get; init; }
    /// <summary><c>attenuation_fraction</c> (f), or null for the default.</summary>
    public double? AttenuationFraction { get; init; }
    /// <summary><c>margin_reason</c>, or null.</summary>
    public string? MarginReason { get; init; }
    /// <summary><c>member_accessions</c>, split on <c>;</c>; empty when the member is found by gene name (AGING-D8).</summary>
    public IReadOnlyList<string> MemberAccessions { get; init; } = Array.Empty<string>();
    /// <summary><c>measure</c>: <c>protein</c> when the column is absent or empty (<c>DEF-DIFF-SHARE-INPUT</c>).</summary>
    public string Measure { get; init; } = SetTestFreeze.ProteinMeasure;
    /// <summary><c>isoform_accession</c>, or null.</summary>
    public string? IsoformAccession { get; init; }
    /// <summary><c>marker_region</c>, or null.</summary>
    public string? MarkerRegion { get; init; }
    /// <summary>The row's <c>note</c> begins <see cref="SetTestFreeze.MaturationPrefix"/> (<c>DEF-DIFF-SET-MATURATION</c>).</summary>
    public required bool Maturation { get; init; }
}

/// <summary>
/// One contrast of one dataset for one hypothesis part, with the member rows that name it, and the flags read over all of
/// them: a contrast is independent only if every member row is, post hoc or maturation if any member row is
/// (<c>DEF-DIFF-SET-POOL</c>, <c>-POSTHOC</c>, <c>-MATURATION</c>).
/// </summary>
public sealed record SetDatasetContrast
{
    /// <summary><c>dataset</c>.</summary>
    public required string Dataset { get; init; }
    /// <summary><c>factor</c>.</summary>
    public required string Factor { get; init; }
    /// <summary><c>test_level</c>: the numerator.</summary>
    public required string TestLevel { get; init; }
    /// <summary><c>reference_level</c>: the denominator.</summary>
    public required string ReferenceLevel { get; init; }
    /// <summary><c>contexts</c> as written.</summary>
    public required string ContextsText { get; init; }
    /// <summary><c>contexts</c> parsed.</summary>
    public required IReadOnlyList<SetContext> Contexts { get; init; }
    /// <summary>The member rows, in file order.</summary>
    public required IReadOnlyList<SetTestRow> Rows { get; init; }
    /// <summary>Independent only if every member row is; otherwise reproduction if any row reproduces, else not pooled.</summary>
    public required SetIndependence Independence { get; init; }
    /// <summary>Post hoc if any member row is (v3.13).</summary>
    public required bool PostHoc { get; init; }
    /// <summary>A maturation contrast if any member row is (v3.16).</summary>
    public required bool Maturation { get; init; }
    /// <summary>The rows' shared <c>primary</c> mark; null when none is given.</summary>
    public bool? Primary { get; init; }
    /// <summary>The members named, in file order, each once.</summary>
    public IReadOnlyList<string> Members => Rows.Select(r => r.Member).Distinct(StringComparer.Ordinal).ToList();
}

/// <summary>
/// One part of a hypothesis: its id and one prediction, with its own members, contrasts, verdict and place in the family
/// (<c>DEF-DIFF-SET-CONTRAST</c>).
/// </summary>
public sealed record SetHypothesisPart
{
    /// <summary><c>hypothesis_id</c>.</summary>
    public required string HypothesisId { get; init; }
    /// <summary><c>prediction</c>.</summary>
    public required SetPrediction Prediction { get; init; }
    /// <summary><c>kind</c>, the same on every row of the part.</summary>
    public required string Kind { get; init; }
    /// <summary>The members named on any row, in file order, each once.</summary>
    public required IReadOnlyList<string> Members { get; init; }
    /// <summary>The dataset contrasts, in file order.</summary>
    public required IReadOnlyList<SetDatasetContrast> Contrasts { get; init; }

    /// <summary>
    /// Tested as a set: kind <c>set</c>; or a label kind (AGING-D9) listing several members. Kind <c>single</c>, or a label
    /// kind with one member, is tested as a single.
    /// </summary>
    public bool IsSet => Kind switch
    {
        SetTestFreeze.SetKind => true,
        SetTestFreeze.SingleKind => false,
        _ => Members.Count > 1,
    };
}

/// <summary>
/// A frozen <c>tests.tsv</c> (QuantProject <c>DEF-DIFF-SET-INPUT</c>): its rows, identified by the sha256 of the bytes
/// read and by the freeze's tag, grouped into hypothesis parts and dataset contrasts.
/// </summary>
/// <remarks>
/// <para>
/// <b>The file:</b> UTF-8 (a byte-order mark is skipped), tab-separated, one header line, no quoting: a cell is the text
/// between tabs. Columns are matched by header name. The columns read are <see cref="RequiredColumns"/> and
/// <see cref="OptionalColumns"/>; <c>rank</c>, <c>grade</c>, <c>low_abundance</c> and <c>registered_in</c> are not read,
/// and <c>note</c> only for <see cref="MaturationPrefix"/>.
/// </para>
/// <para>
/// <b><c>contexts</c>:</b> factor groups separated by <c>;</c>; within a group, levels joined by <c>|</c>, each written
/// <c>factor=level</c> or, after the first, <c>level</c> alone (<c>treatment=UBX0101|treatment=vehicle</c> and
/// <c>treatment=UBX0101|vehicle</c> are the same).
/// </para>
/// <para>
/// <b>Refusals</b> (<see cref="FormatException"/>, naming the line): a missing required column or a header named twice; a
/// row with a different number of cells; an empty required cell other than <c>contexts</c> and <c>mrna_log2</c>; an
/// unknown <c>prediction</c>; <c>post_hoc</c> or <c>primary</c> other than <c>yes</c> or <c>no</c>; a number that does not
/// parse; a malformed <c>contexts</c> cell; rows of one part with different kinds; rows of one dataset contrast with
/// different <c>primary</c> marks.
/// </para>
/// </remarks>
public sealed class SetTestFreeze
{
    /// <summary>The <c>kind</c> tested as a set.</summary>
    public const string SetKind = "set";
    /// <summary>The <c>kind</c> tested as a single.</summary>
    public const string SingleKind = "single";
    /// <summary>The default <c>measure</c>.</summary>
    public const string ProteinMeasure = "protein";
    /// <summary>A row whose <c>note</c> begins with this is a maturation row (aging's protocol, hypotheses-splicing-v2).</summary>
    public const string MaturationPrefix = "MATURATION CONTRAST";

    /// <summary>Columns every freeze has.</summary>
    public static readonly IReadOnlyList<string> RequiredColumns = new[]
    {
        "hypothesis_id", "kind", "prediction", "member", "species", "dataset", "factor", "test_level", "reference_level",
        "contexts", "mrna_log2", "independence", "post_hoc",
    };

    /// <summary>Columns read when present.</summary>
    public static readonly IReadOnlyList<string> OptionalColumns = new[]
    {
        "member_accessions", "margin_log2", "attenuation_fraction", "margin_reason", "primary", "measure",
        "isoform_accession", "marker_region", "note",
    };

    private static readonly IReadOnlyDictionary<string, SetPrediction> Predictions = new Dictionary<string, SetPrediction>(StringComparer.Ordinal)
    {
        ["up"] = SetPrediction.Up,
        ["down"] = SetPrediction.Down,
        ["unchanged"] = SetPrediction.Unchanged,
        ["not_up"] = SetPrediction.NotUp,
        ["not_down"] = SetPrediction.NotDown,
        ["attenuated_down"] = SetPrediction.AttenuatedDown,
        ["attenuated_up"] = SetPrediction.AttenuatedUp,
    };

    private SetTestFreeze(string freeze, string sha256, IReadOnlyList<SetTestRow> rows, IReadOnlyList<SetHypothesisPart> parts)
    {
        Freeze = freeze;
        Sha256 = sha256;
        Rows = rows;
        Parts = parts;
    }

    /// <summary>The freeze's tag, as given, e.g. <c>hypotheses-v4</c>.</summary>
    public string Freeze { get; }

    /// <summary>Lower-case hex sha256 of the bytes read.</summary>
    public string Sha256 { get; }

    /// <summary>Every row, in file order.</summary>
    public IReadOnlyList<SetTestRow> Rows { get; }

    /// <summary>The hypothesis parts: one per (<c>hypothesis_id</c>, <c>prediction</c>), in file order.</summary>
    public IReadOnlyList<SetHypothesisPart> Parts { get; }

    /// <summary>The prediction's name as <c>tests.tsv</c> writes it.</summary>
    public static string PredictionName(SetPrediction prediction) =>
        Predictions.Single(p => p.Value == prediction).Key;

    /// <summary>Reads a freeze from disk.</summary>
    public static SetTestFreeze ReadFile(string path, string freeze)
    {
        ArgumentNullException.ThrowIfNull(path);
        return Read(File.ReadAllBytes(path), freeze);
    }

    /// <summary>Reads a freeze from its bytes, exactly as stored; <see cref="Sha256"/> is taken over them.</summary>
    /// <param name="content">The file's bytes.</param>
    /// <param name="freeze">The freeze's tag.</param>
    public static SetTestFreeze Read(byte[] content, string freeze)
    {
        ArgumentNullException.ThrowIfNull(content);
        if (string.IsNullOrWhiteSpace(freeze))
            throw new ArgumentException("Name the freeze (its tag).", nameof(freeze));

        string sha = Convert.ToHexString(SHA256.HashData(content)).ToLowerInvariant();
        string text = new UTF8Encoding(false, true).GetString(content);
        if (text.Length > 0 && text[0] == '\uFEFF')
            text = text[1..];
        var lines = text.Split('\n').Select(l => l.EndsWith('\r') ? l[..^1] : l).ToList();
        while (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);
        if (lines.Count == 0)
            throw new FormatException("The file is empty; a freeze has a header line.");

        var header = lines[0].Split('\t');
        var index = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int i = 0; i < header.Length; i++)
            if (!index.TryAdd(header[i], i))
                throw new FormatException($"Line 1: the column '{header[i]}' is named twice.");
        var missing = RequiredColumns.Where(c => !index.ContainsKey(c)).ToList();
        if (missing.Count > 0)
            throw new FormatException($"Line 1: missing column(s) {string.Join(", ", missing.Select(c => $"'{c}'"))}.");

        var rows = new List<SetTestRow>();
        for (int l = 1; l < lines.Count; l++)
        {
            int line = l + 1;
            var cells = lines[l].Split('\t');
            if (cells.Length != header.Length)
                throw new FormatException($"Line {line}: {cells.Length} cells, but the header has {header.Length}.");
            rows.Add(ReadRow(line, cells, index));
        }

        return new SetTestFreeze(freeze, sha, rows, BuildParts(rows));
    }

    private static SetTestRow ReadRow(int line, string[] cells, IReadOnlyDictionary<string, int> index)
    {
        string? Cell(string column) => index.TryGetValue(column, out int i) ? cells[i] : null;
        string Required(string column)
        {
            string value = Cell(column)!;
            if (value.Length == 0)
                throw new FormatException($"Line {line}: '{column}' is empty.");
            return value;
        }
        double? Number(string column)
        {
            string? value = Cell(column);
            if (string.IsNullOrEmpty(value))
                return null;
            if (!double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out double x) || !double.IsFinite(x))
                throw new FormatException($"Line {line}: '{column}' is '{value}', not a number.");
            return x;
        }
        bool? YesNo(string column, bool required)
        {
            string? value = Cell(column);
            if (string.IsNullOrEmpty(value))
            {
                if (required)
                    throw new FormatException($"Line {line}: '{column}' is empty.");
                return null;
            }
            return value switch
            {
                "yes" => true,
                "no" => false,
                _ => throw new FormatException($"Line {line}: '{column}' is '{value}'; write yes or no."),
            };
        }
        string? Optional(string column) => string.IsNullOrEmpty(Cell(column)) ? null : Cell(column);

        string predictionText = Required("prediction");
        if (!Predictions.TryGetValue(predictionText, out var prediction))
            throw new FormatException(
                $"Line {line}: prediction '{predictionText}' is not one of {string.Join(", ", Predictions.Keys)}.");

        string independence = Required("independence");
        string contexts = Cell("contexts")!;
        return new SetTestRow
        {
            Line = line,
            HypothesisId = Required("hypothesis_id"),
            Kind = Required("kind"),
            Prediction = prediction,
            Member = Required("member"),
            Species = Required("species"),
            Dataset = Required("dataset"),
            Factor = Required("factor"),
            TestLevel = Required("test_level"),
            ReferenceLevel = Required("reference_level"),
            ContextsText = contexts,
            Contexts = ParseContexts(contexts, line),
            MrnaLog2 = Number("mrna_log2"),
            IndependenceText = independence,
            Independence = ReadIndependence(independence),
            PostHoc = YesNo("post_hoc", required: true)!.Value,
            Primary = YesNo("primary", required: false),
            MarginLog2 = Number("margin_log2"),
            AttenuationFraction = Number("attenuation_fraction"),
            MarginReason = Optional("margin_reason"),
            MemberAccessions = (Optional("member_accessions") ?? "")
                .Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries).ToList(),
            Measure = Optional("measure") ?? ProteinMeasure,
            IsoformAccession = Optional("isoform_accession"),
            MarkerRegion = Optional("marker_region"),
            Maturation = (Cell("note") ?? "").StartsWith(MaturationPrefix, StringComparison.Ordinal),
        };
    }

    /// <summary>
    /// One row's reading of <c>independence</c> (RULED 2026-10-19): a value beginning <c>independent</c> is independent;
    /// one containing <c>reproduces</c> is a reproduction; anything else is reported and never pooled.
    /// </summary>
    public static SetIndependence ReadIndependence(string value)
    {
        ArgumentNullException.ThrowIfNull(value);
        if (value.StartsWith("independent", StringComparison.Ordinal))
            return SetIndependence.Independent;
        return value.Contains("reproduces", StringComparison.Ordinal) ? SetIndependence.Reproduction : SetIndependence.NotPooled;
    }

    /// <summary>
    /// Parses a <c>contexts</c> cell: groups separated by <c>;</c>, levels within a group joined by <c>|</c>, each
    /// <c>factor=level</c> or, after the first, <c>level</c> alone.
    /// </summary>
    /// <exception cref="FormatException">An empty group or level, a group's first level without its factor, two factors
    /// in one group, or one factor in two groups.</exception>
    public static IReadOnlyList<SetContext> ParseContexts(string text, int line = 0)
    {
        ArgumentNullException.ThrowIfNull(text);
        string where = line > 0 ? $"Line {line}: " : "";
        var result = new List<SetContext>();
        if (text.Trim().Length == 0)
            return result;
        foreach (string group in text.Split(';').Select(g => g.Trim()))
        {
            if (group.Length == 0)
                throw new FormatException($"{where}contexts '{text}' has an empty group.");
            string? factor = null;
            var levels = new List<string>();
            foreach (string item in group.Split('|').Select(i => i.Trim()))
            {
                int eq = item.IndexOf('=');
                string itemFactor = eq < 0 ? factor ?? throw new FormatException(
                        $"{where}contexts '{text}': '{item}' names no factor; write factor=level.")
                    : item[..eq].Trim();
                string level = eq < 0 ? item : item[(eq + 1)..].Trim();
                if (itemFactor.Length == 0 || level.Length == 0)
                    throw new FormatException($"{where}contexts '{text}': '{item}' has an empty factor or level.");
                if (factor is not null && !string.Equals(factor, itemFactor, StringComparison.Ordinal))
                    throw new FormatException(
                        $"{where}contexts '{text}': '{factor}' and '{itemFactor}' are joined by |, which joins levels of one factor; separate factors with ;.");
                factor = itemFactor;
                levels.Add(level);
            }
            if (result.Any(c => string.Equals(c.Factor, factor, StringComparison.Ordinal)))
                throw new FormatException($"{where}contexts '{text}' names '{factor}' twice.");
            result.Add(new SetContext(factor!, levels));
        }
        return result;
    }

    private static IReadOnlyList<SetHypothesisPart> BuildParts(IReadOnlyList<SetTestRow> rows)
    {
        var parts = new List<SetHypothesisPart>();
        foreach (var part in rows.GroupBy(r => (r.HypothesisId, r.Prediction)))
        {
            var kinds = part.Select(r => r.Kind).Distinct(StringComparer.Ordinal).ToList();
            if (kinds.Count > 1)
                throw new FormatException(
                    $"Lines {string.Join(", ", part.Select(r => r.Line))}: {part.Key.HypothesisId} {PredictionName(part.Key.Prediction)} " +
                    $"has kinds {string.Join(" and ", kinds.Select(k => $"'{k}'"))}; one part has one kind.");

            var contrasts = new List<SetDatasetContrast>();
            foreach (var contrast in part.GroupBy(r => (r.Dataset, r.Factor, r.TestLevel, r.ReferenceLevel, r.ContextsText)))
            {
                var members = contrast.ToList();
                var primaries = members.Select(r => r.Primary).Distinct().ToList();
                if (primaries.Count > 1)
                    throw new FormatException(
                        $"Lines {string.Join(", ", members.Select(r => r.Line))}: one dataset contrast of {part.Key.HypothesisId} " +
                        $"{PredictionName(part.Key.Prediction)} ({contrast.Key.Dataset}) carries different 'primary' marks.");
                contrasts.Add(new SetDatasetContrast
                {
                    Dataset = contrast.Key.Dataset,
                    Factor = contrast.Key.Factor,
                    TestLevel = contrast.Key.TestLevel,
                    ReferenceLevel = contrast.Key.ReferenceLevel,
                    ContextsText = contrast.Key.ContextsText,
                    Contexts = members[0].Contexts,
                    Rows = members,
                    Independence = members.All(r => r.Independence == SetIndependence.Independent) ? SetIndependence.Independent
                        : members.Any(r => r.Independence == SetIndependence.Reproduction) ? SetIndependence.Reproduction
                        : SetIndependence.NotPooled,
                    PostHoc = members.Any(r => r.PostHoc),
                    Maturation = members.Any(r => r.Maturation),
                    Primary = primaries[0],
                });
            }

            parts.Add(new SetHypothesisPart
            {
                HypothesisId = part.Key.HypothesisId,
                Prediction = part.Key.Prediction,
                Kind = kinds[0],
                Members = part.Select(r => r.Member).Distinct(StringComparer.Ordinal).ToList(),
                Contrasts = contrasts,
            });
        }
        return parts;
    }
}
