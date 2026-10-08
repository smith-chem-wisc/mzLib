using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace Quantification.Differential;

/// <summary>
/// One stratum's samples, the fixed-effects design matrix built from them, and the contrasts that can be tested on it
/// (STAT1 milestone M1; <c>DEF-DIFF-STRATUM</c>, <c>DEF-DIFF-STATUS</c>).
/// </summary>
/// <remarks>
/// <para>
/// <b>Model.</b> One additive model over every condition factor (GR-1): an intercept; each factor treatment-coded
/// against its reference level; batch, treatment-coded against its first level in ordinal order; each numeric covariate
/// as <c>(value − centre) / scale</c>; and any interactions named in <see cref="DesignSpec.Interactions"/>. Columns
/// appear in that order, and within a factor the reference comes first and the other levels follow in ordinal order.
/// The matrix equals R's <c>stats::model.matrix</c> with the same level order and formula
/// (<c>~ a + b + batch + I((x − centre)/scale) + a:b</c>).
/// </para>
/// <para>
/// <b>What is not in the matrix.</b> The individual (subject) is carried on each sample for the random intercept a
/// fit adds (GR-6), never as a fixed column. Strata are not columns either: <see cref="ByStratum"/> builds one design
/// per stratum (GR-5).
/// </para>
/// <para>
/// <b>Refusals and non-fits.</b> Input a design cannot be built from (a duplicate sample, a missing factor value, a
/// non-finite covariate, batch on some samples only, an interaction naming an unknown factor) throws
/// <see cref="ArgumentException"/> naming the sample and the column. A design that can be built but not fully
/// estimated never throws: <see cref="Estimability"/> says which contrasts cannot be estimated and why, and
/// <see cref="NotRun"/> lists the contrasts a stratum cannot run at all.
/// </para>
/// <para>
/// <b>Determinism.</b> Samples are ordered by <see cref="AnalysisSample.SampleId"/> (ordinal), levels ordinally, so the
/// result does not depend on input order.
/// </para>
/// </remarks>
public sealed class AnalysisDesign
{
    /// <summary>The stratum of a dataset with no stratum factor.</summary>
    public const string AllSamples = "all";

    /// <summary>Contrast estimability: estimable.</summary>
    public const string Estimable = "estimable";

    /// <summary>
    /// Contrast estimability: the contrast lies along a direction two or more terms share (e.g. a condition that only
    /// ever occurs in one batch).
    /// </summary>
    public const string NotEstimableConfounded = "not_estimable:confounded";

    /// <summary>Contrast estimability: the contrast needs a column no sample informs (e.g. an empty interaction cell).</summary>
    public const string NotEstimableRankDeficient = "not_estimable:rank_deficient";

    /// <summary>Why a contrast was not run in a stratum: one side has no sample there.</summary>
    public const string NotInStratum = "not_in_stratum";

    private const string Intercept = "intercept";
    private const string BatchTerm = "batch";

    /// <summary>A singular value below this fraction of the largest counts as zero.</summary>
    private const double RankTolerance = 1e-10;

    /// <summary>A null-space loading (or a contrast's projection on it) below this counts as zero.</summary>
    private const double LoadingTolerance = 1e-8;

    private readonly double[,] _matrix;
    private readonly List<string> _columnTerms;
    private readonly List<FactorCoding> _factors;
    private readonly List<CovariateCoding> _covariates;
    private readonly Matrix<double>? _nullSpace;

    private AnalysisDesign(string stratum, IReadOnlyList<AnalysisSample> samples, double[,] matrix,
        List<string> columnNames, List<string> columnTerms, List<FactorCoding> factors, List<CovariateCoding> covariates,
        List<string> warnings, List<NotRunContrast> notRun)
    {
        Stratum = stratum;
        Samples = samples;
        _matrix = matrix;
        ColumnNames = columnNames;
        _columnTerms = columnTerms;
        _factors = factors;
        _covariates = covariates;
        Warnings = warnings;
        NotRun = notRun;

        var x = Matrix<double>.Build.DenseOfArray(matrix);
        var svd = x.Svd(computeVectors: true);
        double largest = svd.S.Count == 0 ? 0 : svd.S.Maximum();
        Rank = svd.S.Count(s => s > RankTolerance * largest);
        int p = columnNames.Count;
        _nullSpace = Rank < p ? svd.VT.SubMatrix(Rank, p - Rank, 0, p).Transpose() : null;
        AliasedTerms = _nullSpace is null
            ? Array.Empty<string>()
            : Enumerable.Range(0, _nullSpace.ColumnCount)
                .Select(TermsOf)
                .Where(t => t.Count > 1)
                .SelectMany(t => t)
                .Distinct(StringComparer.Ordinal)
                .OrderBy(t => _columnTerms.IndexOf(t))
                .ToArray();
    }

    /// <summary>The stratum this design covers (<see cref="StratumKey"/>), or <see cref="AllSamples"/>.</summary>
    public string Stratum { get; }

    /// <summary>The samples, in row order (ordinal by <see cref="AnalysisSample.SampleId"/>).</summary>
    public IReadOnlyList<AnalysisSample> Samples { get; }

    /// <summary>
    /// Column names, in column order: <c>intercept</c>; <c>factor=level</c> for each non-reference level;
    /// <c>batch=level</c>; <c>name (per unit)</c> for a covariate; <c>a=level:b=level</c> for an interaction.
    /// </summary>
    public IReadOnlyList<string> ColumnNames { get; }

    /// <summary>A copy of the design matrix, <c>[sample, column]</c>, in the layout <c>LinearModel.Fit</c> takes.</summary>
    public double[,] Matrix => (double[,])_matrix.Clone();

    /// <summary>The numerical rank of <see cref="Matrix"/>.</summary>
    public int Rank { get; }

    /// <summary>
    /// The terms (factor names, <c>batch</c>, covariate names, <c>intercept</c>) that share a direction with another
    /// term, so cannot be separated from it; empty when the design has full column rank or its only redundancy lies
    /// within one term.
    /// </summary>
    public IReadOnlyList<string> AliasedTerms { get; }

    /// <summary>What was left out of the model or changed from the specification, in plain words.</summary>
    public IReadOnlyList<string> Warnings { get; }

    /// <summary>Contrasts the specification asks for that cannot run in this stratum, each with its reason.</summary>
    public IReadOnlyList<NotRunContrast> NotRun { get; }

    /// <summary>Builds one design over every sample, as one stratum.</summary>
    /// <param name="samples">The samples; each needs a value for every factor and covariate in <paramref name="spec"/>.</param>
    /// <param name="spec">Factors, covariates, interactions and whether batch is fitted.</param>
    /// <param name="stratum">The stratum's key, written on every result row.</param>
    public static AnalysisDesign Create(IEnumerable<AnalysisSample> samples, DesignSpec spec, string stratum = AllSamples)
    {
        ArgumentNullException.ThrowIfNull(samples);
        ArgumentNullException.ThrowIfNull(spec);
        ArgumentException.ThrowIfNullOrWhiteSpace(stratum);

        var rows = samples.OrderBy(s => s.SampleId, StringComparer.Ordinal).ToList();
        if (rows.Count == 0) throw new ArgumentException("A design needs at least one sample.", nameof(samples));
        for (int i = 1; i < rows.Count; i++)
            if (string.Equals(rows[i].SampleId, rows[i - 1].SampleId, StringComparison.Ordinal))
                throw new ArgumentException($"Sample '{rows[i].SampleId}' is given twice.", nameof(samples));

        var warnings = new List<string>();
        var notRun = new List<NotRunContrast>();
        var factors = new List<FactorCoding>();
        foreach (var f in spec.Factors)
        {
            var coding = FactorCoding.Build(f, rows, stratum, warnings, notRun);
            if (coding is not null) factors.Add(coding);
        }

        var known = spec.Factors.Select(f => f.Name).ToHashSet(StringComparer.Ordinal);
        foreach (var (a, b) in spec.Interactions)
            if (!known.Contains(a) || !known.Contains(b))
                throw new ArgumentException(
                    $"The interaction {a}:{b} names a factor that is not in the design ({string.Join(", ", known)}).",
                    nameof(spec));

        FactorCoding? batch = spec.FitBatch ? BuildBatch(rows, stratum, warnings) : null;
        var covariates = spec.Covariates.Select(c => CovariateCoding.Build(c, rows)).ToList();

        var names = new List<string> { Intercept };
        var terms = new List<string> { Intercept };
        var columns = new List<Func<AnalysisSample, double>> { _ => 1 };
        foreach (var f in factors) f.AddColumns(names, terms, columns, f.Name);
        batch?.AddColumns(names, terms, columns, BatchTerm);
        foreach (var c in covariates)
        {
            names.Add(c.ColumnName);
            terms.Add(c.Spec.Name);
            columns.Add(c.Value);
        }
        foreach (var (a, b) in spec.Interactions)
        {
            var fa = factors.FirstOrDefault(f => f.Name == a);
            var fb = factors.FirstOrDefault(f => f.Name == b);
            if (fa is null || fb is null)
            {
                warnings.Add($"The interaction {a}:{b} is not in the model in stratum '{stratum}', because one of its " +
                             "factors has a single level there.");
                continue;
            }
            // R's order for a:b: the first factor's level varies fastest.
            foreach (var lb in fb.Levels.Skip(1))
                foreach (var la in fa.Levels.Skip(1))
                {
                    names.Add($"{a}={la}:{b}={lb}");
                    terms.Add($"{a}:{b}");
                    columns.Add(s => s.Factors[a] == la && s.Factors[b] == lb ? 1 : 0);
                }
        }

        var matrix = new double[rows.Count, columns.Count];
        for (int i = 0; i < rows.Count; i++)
            for (int j = 0; j < columns.Count; j++)
                matrix[i, j] = columns[j](rows[i]);

        return new AnalysisDesign(stratum, rows, matrix, names, terms, factors, covariates, warnings, notRun);
    }

    /// <summary>
    /// Builds one design per stratum (GR-5), keyed by <see cref="StratumKey"/> of each sample's
    /// <see cref="AnalysisSample.Strata"/>, in ordinal order of the key. Nothing is shared between strata.
    /// </summary>
    public static IReadOnlyList<AnalysisDesign> ByStratum(IEnumerable<AnalysisSample> samples, DesignSpec spec)
    {
        ArgumentNullException.ThrowIfNull(samples);
        ArgumentNullException.ThrowIfNull(spec);
        return samples
            .GroupBy(s => StratumKey(s.Strata), StringComparer.Ordinal)
            .OrderBy(g => g.Key, StringComparer.Ordinal)
            .Select(g => Create(g, spec, g.Key))
            .ToList();
    }

    /// <summary>
    /// The stratum's key as <c>DEF-DIFF-STRATUM</c> writes it: the stratum factors sorted by name, each
    /// <c>name=value</c>, joined by <c>;</c>, with <c>%</c>, <c>;</c> and <c>=</c> in a value written <c>%25</c>,
    /// <c>%3B</c> and <c>%3D</c>. No stratum factor gives <see cref="AllSamples"/>.
    /// </summary>
    public static string StratumKey(IReadOnlyDictionary<string, string> values)
    {
        ArgumentNullException.ThrowIfNull(values);
        if (values.Count == 0) return AllSamples;
        return string.Join(";", values
            .OrderBy(kv => kv.Key, StringComparer.Ordinal)
            .Select(kv => kv.Key + "=" + kv.Value.Replace("%", "%25").Replace(";", "%3B").Replace("=", "%3D")));
    }

    /// <summary>A factor's levels in coding order: the reference first, then the others in ordinal order.</summary>
    public IReadOnlyList<string> Levels(string factor) => Factor(factor).Levels;

    /// <summary>The centre subtracted from a covariate: the one stated, else the mean over this stratum's samples.</summary>
    public double CovariateCentre(string covariate) => Covariate(covariate).Centre;

    /// <summary>
    /// The contrasts run when none is declared (ST-4, GR-14): each factor's non-reference levels against its reference,
    /// adjusted for every other term; every pair of levels, later against earlier in level order, for a factor with no
    /// declared reference; and the slope of each covariate. Ids are <c>c1</c>, <c>c2</c>, … in that order. With an
    /// interaction in the model, a factor's contrast is its effect at the reference level of the other factor.
    /// </summary>
    public IReadOnlyList<Contrast> DefaultContrasts()
    {
        var list = new List<Contrast>();
        foreach (var f in _factors)
        {
            if (f.HasReference)
                foreach (var level in f.Levels.Skip(1)) list.Add(Between(f.Name, level, f.Levels[0]));
            else if (f.Declared is null)
                for (int i = 0; i < f.Levels.Count; i++)
                    for (int j = i + 1; j < f.Levels.Count; j++)
                        list.Add(Between(f.Name, f.Levels[j], f.Levels[i]));
        }
        foreach (var c in _covariates) list.Add(Slope(c.Spec.Name));
        return list.Select((c, i) => c with { Id = "c" + (i + 1).ToString(CultureInfo.InvariantCulture) }).ToList();
    }

    /// <summary>
    /// <c>numerator</c> against <c>denominator</c> for one factor: the difference of their coefficients (a level at the
    /// coding baseline has none).
    /// </summary>
    public Contrast Between(string factor, string numerator, string denominator)
    {
        var f = Factor(factor);
        var weights = new double[ColumnNames.Count];
        Add(weights, f, numerator, +1);
        Add(weights, f, denominator, -1);
        string num = $"{factor}={numerator}", den = $"{factor}={denominator}";
        string label = $"{num} vs {den}";
        return new Contrast(label, label, weights, num, den);
    }

    /// <summary>The slope of a covariate, per its scaled unit, other terms held fixed.</summary>
    public Contrast Slope(string covariate)
    {
        var c = Covariate(covariate);
        var weights = new double[ColumnNames.Count];
        weights[ColumnNames.ToList().IndexOf(c.ColumnName)] = 1;
        string unit = c.Spec.ScaledUnit ?? c.Spec.Unit;
        string scale = string.Create(CultureInfo.InvariantCulture,
            $"{c.Spec.Unit}/{c.Spec.Scale:R}, centred at {c.Centre:R} {c.Spec.Unit}");
        return new Contrast(c.ColumnName, c.ColumnName, weights, null, null, c.Spec.Name, unit, scale);
    }

    /// <summary>An explicit contrast: a weight per named column, zero elsewhere.</summary>
    public Contrast Custom(string id, string label, IReadOnlyDictionary<string, double> weights)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(id);
        ArgumentException.ThrowIfNullOrWhiteSpace(label);
        ArgumentNullException.ThrowIfNull(weights);
        var w = new double[ColumnNames.Count];
        foreach (var (column, weight) in weights)
        {
            int j = ColumnNames.ToList().IndexOf(column);
            if (j < 0)
                throw new ArgumentException(
                    $"'{column}' is not a column of the design ({string.Join(", ", ColumnNames)}).", nameof(weights));
            if (!double.IsFinite(weight))
                throw new ArgumentException($"The weight on '{column}' is not finite.", nameof(weights));
            w[j] = weight;
        }
        return new Contrast(id, label, w, null, null);
    }

    /// <summary>
    /// Whether a contrast can be estimated from this design: <see cref="Estimable"/> when its weights are orthogonal to
    /// the matrix's null space; otherwise <see cref="NotEstimableConfounded"/> when the direction it needs is shared by
    /// two or more terms, else <see cref="NotEstimableRankDeficient"/>.
    /// </summary>
    public string Estimability(Contrast contrast)
    {
        ArgumentNullException.ThrowIfNull(contrast);
        if (contrast.Weights.Count != ColumnNames.Count)
            throw new ArgumentException(
                $"The contrast has {contrast.Weights.Count} weights but the design has {ColumnNames.Count} columns.",
                nameof(contrast));
        if (_nullSpace is null) return Estimable;
        var c = Vector<double>.Build.DenseOfEnumerable(contrast.Weights);
        double scale = Math.Max(1, c.L2Norm());
        var projection = _nullSpace.TransposeThisAndMultiply(c);
        var touched = Enumerable.Range(0, projection.Count)
            .Where(k => Math.Abs(projection[k]) > LoadingTolerance * scale).ToList();
        if (touched.Count == 0) return Estimable;
        return touched.Any(k => TermsOf(k).Count > 1) ? NotEstimableConfounded : NotEstimableRankDeficient;
    }

    private List<string> TermsOf(int nullVector) =>
        Enumerable.Range(0, _nullSpace!.RowCount)
            .Where(j => Math.Abs(_nullSpace[j, nullVector]) > LoadingTolerance)
            .Select(j => _columnTerms[j])
            .Distinct(StringComparer.Ordinal)
            .ToList();

    private void Add(double[] weights, FactorCoding f, string level, double sign)
    {
        if (!f.Levels.Contains(level, StringComparer.Ordinal))
            throw new ArgumentException(
                $"'{level}' is not a level of '{f.Name}' in stratum '{Stratum}' ({string.Join(", ", f.Levels)}).",
                nameof(level));
        if (level == f.Levels[0]) return;
        weights[ColumnNames.ToList().IndexOf($"{f.Name}={level}")] += sign;
    }

    private FactorCoding Factor(string name) =>
        _factors.FirstOrDefault(f => f.Name == name)
        ?? throw new ArgumentException(
            $"'{name}' is not a factor of the model in stratum '{Stratum}' ({string.Join(", ", _factors.Select(f => f.Name))}).",
            nameof(name));

    private CovariateCoding Covariate(string name) =>
        _covariates.FirstOrDefault(c => c.Spec.Name == name)
        ?? throw new ArgumentException(
            $"'{name}' is not a covariate of the model ({string.Join(", ", _covariates.Select(c => c.Spec.Name))}).",
            nameof(name));

    private static FactorCoding? BuildBatch(List<AnalysisSample> rows, string stratum, List<string> warnings)
    {
        var withBatch = rows.Where(r => !string.IsNullOrWhiteSpace(r.Batch)).ToList();
        if (withBatch.Count == 0) return null;
        var missing = rows.FirstOrDefault(r => string.IsNullOrWhiteSpace(r.Batch));
        if (missing is not null)
            throw new ArgumentException(
                $"Sample '{missing.SampleId}' has no batch, but other samples do; a batch is needed on every sample or none.",
                nameof(rows));
        var levels = withBatch.Select(r => r.Batch!).Distinct(StringComparer.Ordinal).OrderBy(l => l, StringComparer.Ordinal).ToList();
        if (levels.Count < 2)
        {
            warnings.Add($"Every sample in stratum '{stratum}' is in batch '{levels[0]}', so batch is not in the model.");
            return null;
        }
        return new FactorCoding(BatchTerm, null, levels, HasReference: true, BatchOf: true);
    }

    private sealed record FactorCoding(string Name, string? Declared, List<string> Levels, bool HasReference, bool BatchOf = false)
    {
        public static FactorCoding? Build(FactorSpec spec, List<AnalysisSample> rows, string stratum,
            List<string> warnings, List<NotRunContrast> notRun)
        {
            ArgumentException.ThrowIfNullOrWhiteSpace(spec.Name);
            foreach (var r in rows)
                if (!r.Factors.TryGetValue(spec.Name, out var v) || string.IsNullOrWhiteSpace(v))
                    throw new ArgumentException($"Sample '{r.SampleId}' has no value for factor '{spec.Name}'.", nameof(rows));

            var levels = rows.Select(r => r.Factors[spec.Name]).Distinct(StringComparer.Ordinal)
                .OrderBy(l => l, StringComparer.Ordinal).ToList();
            if (levels.Count < 2)
            {
                warnings.Add($"Factor '{spec.Name}' has one level ('{levels[0]}') in stratum '{stratum}', so it is not in " +
                             "the model and none of its contrasts is run there.");
                notRun.Add(new NotRunContrast($"{spec.Name} (one level: {levels[0]})", NotInStratum));
                return null;
            }

            if (spec.ReferenceLevel is null)
            {
                warnings.Add($"No reference level was declared for factor '{spec.Name}', so every pair of its levels is " +
                             $"compared (GR-14). Levels: {string.Join(", ", levels)}.");
                return new FactorCoding(spec.Name, null, levels, HasReference: false);
            }

            if (!levels.Contains(spec.ReferenceLevel, StringComparer.Ordinal))
            {
                warnings.Add($"The reference level '{spec.ReferenceLevel}' of factor '{spec.Name}' has no sample in stratum " +
                             $"'{stratum}', so no contrast against it is run there.");
                foreach (var l in levels)
                    notRun.Add(new NotRunContrast($"{spec.Name}={l} vs {spec.Name}={spec.ReferenceLevel}", NotInStratum));
                return new FactorCoding(spec.Name, spec.ReferenceLevel, levels, HasReference: false);
            }

            levels.Remove(spec.ReferenceLevel);
            levels.Insert(0, spec.ReferenceLevel);
            return new FactorCoding(spec.Name, spec.ReferenceLevel, levels, HasReference: true);
        }

        public void AddColumns(List<string> names, List<string> terms, List<Func<AnalysisSample, double>> columns, string term)
        {
            foreach (var level in Levels.Skip(1))
            {
                names.Add($"{Name}={level}");
                terms.Add(term);
                columns.Add(BatchOf ? s => s.Batch == level ? 1 : 0 : s => s.Factors[Name] == level ? 1 : 0);
            }
        }
    }

    private sealed record CovariateCoding(CovariateSpec Spec, double Centre)
    {
        public string ColumnName => $"{Spec.Name} (per {Spec.ScaledUnit ?? Spec.Unit})";

        public double Value(AnalysisSample s) => (s.Covariates[Spec.Name] - Centre) / Spec.Scale;

        public static CovariateCoding Build(CovariateSpec spec, List<AnalysisSample> rows)
        {
            ArgumentException.ThrowIfNullOrWhiteSpace(spec.Name);
            ArgumentException.ThrowIfNullOrWhiteSpace(spec.Unit);
            if (!double.IsFinite(spec.Scale) || spec.Scale <= 0)
                throw new ArgumentException($"The scale of covariate '{spec.Name}' must be finite and positive.", nameof(spec));
            if (spec.Centre is { } stated && !double.IsFinite(stated))
                throw new ArgumentException($"The centre of covariate '{spec.Name}' is not finite.", nameof(spec));
            foreach (var r in rows)
                if (!r.Covariates.TryGetValue(spec.Name, out var v) || !double.IsFinite(v))
                    throw new ArgumentException($"Sample '{r.SampleId}' needs a finite value for covariate '{spec.Name}'.",
                        nameof(rows));
            return new CovariateCoding(spec, spec.Centre ?? rows.Average(r => r.Covariates[spec.Name]));
        }
    }
}
