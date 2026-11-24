using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Omics.BioPolymerGroup;

namespace Quantification.Differential;

/// <summary>
/// Writes <c>DifferentialDecisions.tsv</c> (QuantProject <c>DEF-DIFF-DECISIONS</c>, user requirement R-DECISIONS): one line
/// for every place a rule of the STAT1 runner changed what happened to the analysis, a stratum, a feature, or a feature in
/// one comparison, with the ruling it follows, the facts that triggered it and what it did. Settings that apply to
/// everything are one analysis-scope line. The file joins to <c>DifferentialResults.tsv</c> on its keys (stratum, quant
/// basis, feature, contrast). UTF-8 without a byte-order mark, <c>\n</c> line ends, a fixed line order and no timestamps,
/// so the same analysis gives the same bytes; a tab, CR or LF in a text cell becomes a space, as in the results table.
/// </summary>
public static class DifferentialDecisionWriter
{
    /// <summary>The file's name, beside <see cref="DifferentialResultWriter.ResultsFileName"/>.</summary>
    public const string FileName = "DifferentialDecisions.tsv";

    /// <summary>A line about the whole analysis.</summary>
    public const string ScopeAnalysis = "analysis";
    /// <summary>A line about one stratum (and quant basis).</summary>
    public const string ScopeStratum = "stratum";
    /// <summary>A line about one comparison in one stratum, for every feature.</summary>
    public const string ScopeContrast = "contrast";
    /// <summary>A line about one feature in one stratum, every comparison.</summary>
    public const string ScopeFeature = "feature";
    /// <summary>A line about one feature in one comparison: one results row.</summary>
    public const string ScopeFeatureContrast = "feature_contrast";

    /// <summary>The stratum of an analysis-scope line.</summary>
    public const string AnalysisStratum = "";

    private static readonly IReadOnlyList<TsvColumn<(string Analysis, DifferentialDecision Decision)>> Schema = new[]
    {
        Column("analysis_id", d => d.Analysis),
        Column("scope", d => Scope(d.Decision)),
        Column("ruling", d => DifferentialDecisionRules.Ruling(d.Decision.Rule)),
        Column("rule", d => d.Decision.Rule),
        Column("stratum", d => d.Decision.Stratum),
        Column("quant_basis", d => d.Decision.QuantBasis),
        Column("feature_id", d => d.Decision.FeatureId),
        Column("contrast_id", d => d.Decision.ContrastId),
        Column("facts", d => d.Decision.Facts),
        Column("effect", d => d.Decision.Effect),
    };

    private static TsvColumn<(string Analysis, DifferentialDecision Decision)> Column(string header,
        Func<(string Analysis, DifferentialDecision Decision), string?> value) =>
        new(header, d => DifferentialResultWriter.CleanText(value(d)) ?? "");

    /// <summary>A line's scope, from which of its keys are set.</summary>
    public static string Scope(DifferentialDecision decision)
    {
        ArgumentNullException.ThrowIfNull(decision);
        if (decision.Stratum == AnalysisStratum) return ScopeAnalysis;
        return (decision.FeatureId, decision.ContrastId) switch
        {
            (null, null) => ScopeStratum,
            (null, _) => ScopeContrast,
            (_, null) => ScopeFeature,
            _ => ScopeFeatureContrast,
        };
    }

    /// <summary>
    /// The lines in file order: the analysis-scope lines, then by stratum, quant basis, feature and contrast (a missing key
    /// before any value, ordinally), then by the rules' order in <see cref="DifferentialDecisionRules.All"/>, then by facts.
    /// </summary>
    public static IReadOnlyList<DifferentialDecision> Order(IEnumerable<DifferentialDecision> decisions)
    {
        ArgumentNullException.ThrowIfNull(decisions);
        return decisions
            .OrderBy(d => d.Stratum == AnalysisStratum ? 0 : 1)
            .ThenBy(d => d.Stratum, StringComparer.Ordinal)
            .ThenBy(d => d.QuantBasis ?? "", StringComparer.Ordinal)
            .ThenBy(d => d.FeatureId is null ? 0 : 1).ThenBy(d => d.FeatureId ?? "", StringComparer.Ordinal)
            .ThenBy(d => d.ContrastId is null ? 0 : 1).ThenBy(d => d.ContrastId ?? "", StringComparer.Ordinal)
            .ThenBy(d => DifferentialDecisionRules.Index(d.Rule))
            .ThenBy(d => d.Facts, StringComparer.Ordinal)
            .ToList();
    }

    /// <summary>
    /// The file's bytes. Every line is checked first: its rule is one of <see cref="DifferentialDecisionRules"/>; an
    /// analysis-scope line has no other key; any other line names a stratum of <paramref name="metadata"/>, and a line
    /// naming a feature (and a comparison) matches a results row with those keys.
    /// </summary>
    /// <exception cref="ArgumentException">A line breaks one of those rules.</exception>
    public static byte[] ToBytes(DifferentialMetadata metadata, IReadOnlyList<DifferentialResult> rows,
        IEnumerable<DifferentialDecision> decisions)
    {
        ArgumentNullException.ThrowIfNull(metadata);
        ArgumentNullException.ThrowIfNull(rows);
        var ordered = Order(decisions);
        Validate(metadata, rows, ordered);
        string analysis = metadata.AnalysisId;
        var text = new StringBuilder(TsvWriter.HeaderLine(Schema)).Append('\n');
        foreach (var d in ordered) text.Append(TsvWriter.RowLine(Schema, (analysis, d))).Append('\n');
        return new UTF8Encoding(false).GetBytes(text.ToString());
    }

    private static void Validate(DifferentialMetadata metadata, IReadOnlyList<DifferentialResult> rows,
        IReadOnlyList<DifferentialDecision> decisions)
    {
        var strata = new HashSet<string>(metadata.Strata.Select(s => s.Stratum), StringComparer.Ordinal);
        var features = new HashSet<(string, string, string)>(rows.Select(r => (r.Stratum, r.QuantBasis ?? "", r.FeatureId)));
        var keys = new HashSet<(string, string, string, string)>(rows.Select(r => (r.Stratum, r.QuantBasis ?? "", r.FeatureId, r.ContrastId)));
        foreach (var d in decisions)
        {
            ArgumentNullException.ThrowIfNull(d, nameof(decisions));
            string at = $"Decision '{d.Rule}' (stratum '{d.Stratum}', feature '{d.FeatureId}', contrast '{d.ContrastId}')";
            if (!DifferentialDecisionRules.IsValid(d.Rule))
                throw new ArgumentException($"{at}: '{d.Rule}' is not a decision rule.", nameof(decisions));
            if (string.IsNullOrWhiteSpace(d.Facts) || string.IsNullOrWhiteSpace(d.Effect))
                throw new ArgumentException($"{at} has no facts or no effect.", nameof(decisions));
            if (d.Stratum == AnalysisStratum)
            {
                if (d.QuantBasis is not null || d.FeatureId is not null || d.ContrastId is not null)
                    throw new ArgumentException($"{at}: an analysis-scope line has no other key.", nameof(decisions));
                continue;
            }
            if (!strata.Contains(d.Stratum))
                throw new ArgumentException($"{at} names a stratum the metadata does not list.", nameof(decisions));
            if (d.FeatureId is null) continue;
            bool matched = d.ContrastId is null
                ? features.Contains((d.Stratum, d.QuantBasis ?? "", d.FeatureId))
                : keys.Contains((d.Stratum, d.QuantBasis ?? "", d.FeatureId, d.ContrastId));
            if (!matched)
                throw new ArgumentException($"{at} matches no results row.", nameof(decisions));
        }
    }

    /// <summary>
    /// How many distinct features, and feature x comparison rows, each rule touched (<c>StatisticalReport.md</c>'s
    /// "Decisions applied"), in rule order; rules with no line are left out.
    /// </summary>
    public static IReadOnlyList<(string Rule, int Lines, int Features, int Rows)> Counts(IEnumerable<DifferentialDecision> decisions)
    {
        ArgumentNullException.ThrowIfNull(decisions);
        return decisions.GroupBy(d => d.Rule, StringComparer.Ordinal)
            .OrderBy(g => DifferentialDecisionRules.Index(g.Key))
            .Select(g => (g.Key, g.Count(),
                g.Where(d => d.FeatureId is not null).Select(d => (d.Stratum, d.QuantBasis, d.FeatureId)).Distinct().Count(),
                g.Where(d => d.FeatureId is not null && d.ContrastId is not null)
                    .Select(d => (d.Stratum, d.QuantBasis, d.FeatureId, d.ContrastId)).Distinct().Count()))
            .ToList();
    }

    internal static string Count(int n) => n.ToString(CultureInfo.InvariantCulture);
}
