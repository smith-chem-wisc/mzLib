using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using StatisticalModels;

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
/// <see cref="NotRun"/> lists the contrasts a stratum cannot run at all. <c>LinearModel.Fit</c> refuses such a
/// <see cref="Matrix"/>, so a fit takes <see cref="FittedMatrix"/> (the later redundant columns left out) and each
/// estimable contrast's <see cref="FittedWeights"/>.
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

    /// <summary>
    /// A null-space loading (or a unit-length contrast's projection on it) below this counts as zero, both on the columns
    /// divided by their lengths.
    /// </summary>
    private const double LoadingTolerance = 1e-8;

    private readonly double[,] _matrix;
    private readonly List<string> _columnTerms;
    private readonly List<FactorCoding> _factors;
    private readonly List<CovariateCoding> _covariates;
    private readonly List<DefaultContrast> _defaultsRunHere;
    private readonly double[] _columnLengths;
    private readonly Matrix<double>? _nullProjector;

    private AnalysisDesign(string stratum, IReadOnlyList<AnalysisSample> samples, double[,] matrix,
        List<string> columnNames, List<string> columnTerms, List<FactorCoding> factors, List<CovariateCoding> covariates,
        List<DefaultContrast> defaultsRunHere, List<string> warnings, List<NotRunContrast> notRun)
    {
        Stratum = stratum;
        Samples = samples;
        _matrix = matrix;
        ColumnNames = columnNames;
        _columnTerms = columnTerms;
        _factors = factors;
        _covariates = covariates;
        _defaultsRunHere = defaultsRunHere;
        Warnings = warnings;
        NotRun = notRun;

        int p = columnNames.Count;
        // One rule decides everything: LinearModel's full-rank test, the one its Fit applies. The columns it keeps size the
        // rank, and the columns it leaves out give the null space, so a contrast is estimable here exactly when a fit of
        // FittedMatrix can give it (review round 4; an SVD tolerance beside it could disagree).
        FittedColumns = LinearModel.HasFullColumnRank(matrix) ? Enumerable.Range(0, p).ToArray() : KeptColumns(matrix);
        Rank = FittedColumns.Count;
        // What a redundancy involves is read from the projection onto the null space, which depends on the space alone,
        // not on the basis that spans it (review round 3). It is the null space of the columns divided by their lengths,
        // so neither it nor a verdict depends on the unit a column is stored in (review round 5). An all-zero column keeps
        // length 1. Scaling a column does not change which columns are redundant, so FittedColumns stand.
        _columnLengths = Matrix<double>.Build.DenseOfArray(matrix).ColumnNorms(2).Select(l => l > 0 ? l : 1).ToArray();
        var nullSpace = Rank < p ? NullSpace(matrix, _columnLengths, FittedColumns) : null;
        _nullProjector = nullSpace?.TransposeAndMultiply(nullSpace);
        AliasedTerms = _nullProjector is null
            ? Array.Empty<string>()
            : Enumerable.Range(0, p)
                .Select(j => TermsAlong(_nullProjector.Column(j)))
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

    /// <summary>
    /// The columns a fit takes: all of them when the design has full column rank. Otherwise <see cref="Matrix"/>'s columns
    /// left to right, each kept unless it is redundant with those already kept (by
    /// <see cref="LinearModel.HasFullColumnRank"/>), so the later of two aliased columns is left out.
    /// </summary>
    public IReadOnlyList<int> FittedColumns { get; }

    /// <summary>
    /// A copy of <see cref="Matrix"/> with only <see cref="FittedColumns"/>. It has full column rank, so
    /// <c>LinearModel.Fit</c> takes it; an estimable contrast is fitted on it with <see cref="FittedWeights"/>.
    /// </summary>
    public double[,] FittedMatrix
    {
        get
        {
            int n = _matrix.GetLength(0);
            var fitted = new double[n, FittedColumns.Count];
            for (int i = 0; i < n; i++)
                for (int k = 0; k < FittedColumns.Count; k++)
                    fitted[i, k] = _matrix[i, FittedColumns[k]];
            return fitted;
        }
    }

    /// <summary>
    /// An estimable contrast's weights on <see cref="FittedColumns"/>. Its estimate and standard error from a fit of
    /// <see cref="FittedMatrix"/> are those of every least-squares fit of <see cref="Matrix"/>: the left-out coefficients
    /// can be taken as 0, and an estimable contrast does not depend on which solution is taken. Hand them on as
    /// <c>new ContrastWeights(contrast.Id, weights)</c>.
    /// </summary>
    /// <exception cref="ArgumentException">The contrast is not estimable; the message gives its <see cref="Estimability"/>.</exception>
    public double[] FittedWeights(Contrast contrast)
    {
        string estimability = Estimability(contrast);
        if (estimability != Estimable)
            throw new ArgumentException($"Contrast '{contrast.Id}' is {estimability}, so no fit gives it an estimate.",
                nameof(contrast));
        return FittedColumns.Select(j => contrast.Weights[j]).ToArray();
    }

    /// <summary>
    /// The numerical rank of <see cref="Matrix"/>, by <see cref="LinearModel.HasFullColumnRank"/>: the number of
    /// <see cref="FittedColumns"/>, so the column count exactly when <c>LinearModel.Fit</c> will fit the matrix.
    /// </summary>
    public int Rank { get; }

    /// <summary>
    /// The terms (factor names, <c>batch</c>, covariate names, <c>intercept</c>) that share a direction with another
    /// term, so cannot be separated from it; empty when the design has full column rank or its only redundancy lies
    /// within one term.
    /// </summary>
    public IReadOnlyList<string> AliasedTerms { get; }

    /// <summary>What was left out of the model or changed from the specification, in plain words.</summary>
    public IReadOnlyList<string> Warnings { get; }

    /// <summary>
    /// Contrasts the specification asks for that cannot run in this stratum, each with its reason: every default contrast
    /// of the analysis (<see cref="DefaultContrasts"/>) with no sample on one side here, or none at the level of an
    /// interacting factor it is taken at, under its label, its <see cref="NotRunContrast.Id"/> and
    /// <see cref="NotInStratum"/>, in the order of their ids. A factor the analysis has no default contrast for is listed
    /// with no id: one with one level here as <c>factor (one level: level)</c>, one whose declared reference has no sample
    /// as each level here against it.
    /// </summary>
    public IReadOnlyList<NotRunContrast> NotRun { get; }

    /// <summary>Builds one design over every sample, as one stratum.</summary>
    /// <param name="samples">The samples; each needs a value for every factor and covariate in <paramref name="spec"/>.</param>
    /// <param name="spec">Factors, covariates, interactions and whether batch is fitted.</param>
    /// <param name="stratum">The stratum's key, written on every result row.</param>
    public static AnalysisDesign Create(IEnumerable<AnalysisSample> samples, DesignSpec spec, string stratum = AllSamples) =>
        CreateStratum(samples, spec, stratum, defaultSet: null);

    /// <param name="defaultSet">The analysis's default contrasts; null to take them from these samples alone.</param>
    private static AnalysisDesign CreateStratum(IEnumerable<AnalysisSample> samples, DesignSpec spec, string stratum,
        List<DefaultContrast>? defaultSet)
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
        var factors = new List<FactorCoding>();
        foreach (var f in spec.Factors)
        {
            var coding = FactorCoding.Build(f, rows, stratum, warnings);
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

        defaultSet ??= DefaultSet(rows, spec);
        var levelsHere = new Dictionary<string, List<string>>(StringComparer.Ordinal);
        foreach (var f in spec.Factors) levelsHere[f.Name] = LevelsOf(rows, f.Name);
        var notRun = NotRunHere(spec, defaultSet, levelsHere, stratum, warnings);
        var defaultsRunHere = defaultSet.Where(d => d.RunsIn(levelsHere)).ToList();

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

        return new AnalysisDesign(stratum, rows, matrix, names, terms, factors, covariates, defaultsRunHere, warnings, notRun);
    }

    /// <summary>
    /// Builds one design per stratum (GR-5), keyed by <see cref="StratumKey"/> of each sample's
    /// <see cref="AnalysisSample.Strata"/>, in ordinal order of the key. Each stratum has its own design, centres and
    /// fit; only the default contrasts' ids are shared: they are numbered once over all the samples
    /// (<see cref="DefaultContrasts"/>), so an id is the same comparison in every stratum.
    /// </summary>
    public static IReadOnlyList<AnalysisDesign> ByStratum(IEnumerable<AnalysisSample> samples, DesignSpec spec)
    {
        ArgumentNullException.ThrowIfNull(samples);
        ArgumentNullException.ThrowIfNull(spec);
        var all = samples.ToList();
        var defaultSet = DefaultSet(all, spec);
        return all
            .GroupBy(s => StratumKey(s.Strata), StringComparer.Ordinal)
            .OrderBy(g => g.Key, StringComparer.Ordinal)
            .Select(g => CreateStratum(g, spec, g.Key, defaultSet))
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
    /// declared reference; and the slope of each covariate. With an interaction in the model, a factor's contrast is its
    /// effect at the reference level of the other factor (its first level in ordinal order when it has no reference).
    /// </summary>
    /// <remarks>
    /// Ids are <c>c1</c>, <c>c2</c>, … in that order over the whole analysis: the levels of a design from
    /// <see cref="ByStratum"/> are taken over all its samples, so an id is the same comparison in every stratum. A
    /// stratum with no sample on one side of a contrast skips its id and lists the contrast in <see cref="NotRun"/>, as
    /// does one with no sample at the other factor's reference when the two interact: there the same weights would give
    /// the effect at another level.
    /// The weights are over this design's own columns.
    /// </remarks>
    public IReadOnlyList<Contrast> DefaultContrasts() =>
        _defaultsRunHere
            .Select(d => (d.Covariate is null ? Between(d.Factor!, d.Numerator!, d.Denominator!) : Slope(d.Covariate)) with { Id = d.Id })
            .ToList();

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
        string label = BetweenLabel(factor, numerator, denominator);
        return new Contrast(label, label, weights, $"{factor}={numerator}", $"{factor}={denominator}");
    }

    private static string BetweenLabel(string factor, string numerator, string denominator) =>
        $"{factor}={numerator} vs {factor}={denominator}";

    /// <summary>The slope of a covariate, per its scaled unit, other terms held fixed.</summary>
    public Contrast Slope(string covariate)
    {
        var c = Covariate(covariate);
        var weights = new double[ColumnNames.Count];
        weights[ColumnNames.ToList().IndexOf(c.ColumnName)] = 1;
        string scale = string.Create(CultureInfo.InvariantCulture,
            $"{c.Spec.Unit}/{c.Spec.Scale:R}, centred at {c.Centre:R} {c.Spec.Unit}");
        return new Contrast(c.ColumnName, c.ColumnName, weights, null, null, c.Spec.Name, c.PerUnit, scale);
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
    /// two or more terms, else <see cref="NotEstimableRankDeficient"/>. It is judged on the columns divided by their
    /// lengths, so it is the same in whatever unit a covariate is stored.
    /// </summary>
    public string Estimability(Contrast contrast)
    {
        ArgumentNullException.ThrowIfNull(contrast);
        if (contrast.Weights.Count != ColumnNames.Count)
            throw new ArgumentException(
                $"The contrast has {contrast.Weights.Count} weights but the design has {ColumnNames.Count} columns.",
                nameof(contrast));
        if (_nullProjector is null) return Estimable;
        // A column divided by its length has its coefficient multiplied by that length, so the contrast's weight on it
        // is divided by it.
        var c = Vector<double>.Build.Dense(ColumnNames.Count, j => contrast.Weights[j] / _columnLengths[j]);
        double length = c.L2Norm();
        // The part of the contrast that lies in the null space, which no fit can see.
        var unseen = _nullProjector * c;
        if (unseen.L2Norm() <= LoadingTolerance * length) return Estimable;
        return TermsAlong(unseen / length).Count > 1 ? NotEstimableConfounded : NotEstimableRankDeficient;
    }

    /// <summary>The terms a direction in the null space loads on: its entries above <see cref="LoadingTolerance"/> of its length.</summary>
    private List<string> TermsAlong(Vector<double> direction)
    {
        double length = direction.L2Norm();
        if (length <= LoadingTolerance) return new List<string>();
        return Enumerable.Range(0, direction.Count)
            .Where(j => Math.Abs(direction[j]) > LoadingTolerance * length)
            .Select(j => _columnTerms[j])
            .Distinct(StringComparer.Ordinal)
            .ToList();
    }

    /// <summary>Left to right, each column kept unless it is redundant with those already kept, by LinearModel's own test.</summary>
    private static int[] KeptColumns(double[,] matrix)
    {
        int n = matrix.GetLength(0), p = matrix.GetLength(1);
        var kept = new List<int>();
        for (int j = 0; j < p; j++)
        {
            var candidate = new double[n, kept.Count + 1];
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < kept.Count; k++)
                    candidate[i, k] = matrix[i, kept[k]];
                candidate[i, kept.Count] = matrix[i, j];
            }
            if (LinearModel.HasFullColumnRank(candidate))
                kept.Add(j);
        }
        return kept.ToArray();
    }

    /// <summary>
    /// An orthonormal basis of the null space of the columns divided by <paramref name="lengths"/>,
    /// <c>[column, direction]</c>: for each left-out column, its unit vector less its least-squares coefficients on the
    /// kept columns. Those directions are independent (each is 1 on its own left-out column and 0 on the others), and
    /// for an exact redundancy they span the same space an SVD gives.
    /// </summary>
    private static Matrix<double> NullSpace(double[,] matrix, double[] lengths, IReadOnlyList<int> kept)
    {
        int p = matrix.GetLength(1);
        var x = Matrix<double>.Build.DenseOfArray(matrix);
        x.MapIndexedInplace((i, j, v) => v / lengths[j]);
        var keptColumns = Matrix<double>.Build.DenseOfColumnVectors(kept.Select(j => x.Column(j)));
        var keptQr = keptColumns.QR();
        var directions = new List<Vector<double>>();
        for (int j = 0; j < p; j++)
        {
            if (kept.Contains(j)) continue;
            var coefficients = keptQr.Solve(x.Column(j));
            var direction = Vector<double>.Build.Dense(p);
            direction[j] = 1;
            for (int k = 0; k < kept.Count; k++)
                direction[kept[k]] = -coefficients[k];
            directions.Add(direction);
        }
        return Matrix<double>.Build.DenseOfColumnVectors(directions).QR(QRMethod.Thin).Q;
    }

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
        return new FactorCoding(BatchTerm, levels, BatchOf: true);
    }

    /// <summary>
    /// The default contrasts (<see cref="DefaultContrasts"/>) over all of an analysis's samples, numbered <c>c1</c>,
    /// <c>c2</c>, …: a factor's levels are those of every sample, so each stratum runs the ones it has samples for under
    /// these ids. A sample with no value for a factor adds no level; <see cref="Create"/> refuses it.
    /// </summary>
    private static List<DefaultContrast> DefaultSet(IReadOnlyList<AnalysisSample> samples, DesignSpec spec)
    {
        // The level each factor is coded against over the whole analysis: its declared reference when some sample has it,
        // else its first level in ordinal order, as FactorCoding.Build codes a stratum.
        var baseline = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (var f in spec.Factors)
        {
            if (string.IsNullOrWhiteSpace(f.Name)) continue;
            var levels = LevelsOf(samples, f.Name);
            if (levels.Count > 0)
                baseline[f.Name] = f.ReferenceLevel is not null && levels.Contains(f.ReferenceLevel, StringComparer.Ordinal)
                    ? f.ReferenceLevel
                    : levels[0];
        }

        var set = new List<DefaultContrast>();
        foreach (var f in spec.Factors)
        {
            if (string.IsNullOrWhiteSpace(f.Name)) continue;
            var levels = LevelsOf(samples, f.Name);
            // With an interaction, a factor's contrast is its effect at the baseline of each factor it interacts with.
            var at = spec.Interactions
                .Select(i => i.First == f.Name ? i.Second : i.Second == f.Name ? i.First : null)
                .Where(other => other is not null && other != f.Name && baseline.ContainsKey(other))
                .Distinct(StringComparer.Ordinal)
                .Select(other => (Factor: other!, Level: baseline[other!]))
                .ToList();
            if (f.ReferenceLevel is null)
            {
                for (int i = 0; i < levels.Count; i++)
                    for (int j = i + 1; j < levels.Count; j++)
                        set.Add(new DefaultContrast("", f.Name, levels[j], levels[i]) { At = at });
            }
            else if (levels.Contains(f.ReferenceLevel, StringComparer.Ordinal))
            {
                foreach (var level in levels.Where(l => l != f.ReferenceLevel))
                    set.Add(new DefaultContrast("", f.Name, level, f.ReferenceLevel) { At = at });
            }
        }
        foreach (var c in spec.Covariates) set.Add(new DefaultContrast("", Covariate: c.Name));
        return set.Select((d, i) => d with { Id = "c" + (i + 1).ToString(CultureInfo.InvariantCulture) }).ToList();
    }

    /// <summary>A factor's distinct levels over these samples, in ordinal order; a sample with no value adds none.</summary>
    private static List<string> LevelsOf(IEnumerable<AnalysisSample> samples, string factor) =>
        samples
            .Select(s => s.Factors.TryGetValue(factor, out var v) ? v : null)
            .Where(v => !string.IsNullOrWhiteSpace(v))
            .Select(v => v!)
            .Distinct(StringComparer.Ordinal).OrderBy(l => l, StringComparer.Ordinal).ToList();

    /// <summary>
    /// The one place a stratum's <see cref="NotRun"/> is decided. A factor with default contrasts lists each one this
    /// stratum cannot run, under its label and id (#1450 review); a factor the analysis has none for is listed by what it
    /// lacks here, with no id: one level, or no sample at its declared reference.
    /// </summary>
    private static List<NotRunContrast> NotRunHere(DesignSpec spec, List<DefaultContrast> defaultSet,
        Dictionary<string, List<string>> levelsHere, string stratum, List<string> warnings)
    {
        var notRun = new List<NotRunContrast>();
        foreach (var f in spec.Factors)
        {
            var own = defaultSet.Where(d => d.Factor == f.Name).ToList();
            var levels = levelsHere[f.Name];
            if (own.Count > 0)
            {
                foreach (var d in own.Where(d => !d.RunsIn(levelsHere)))
                {
                    notRun.Add(new NotRunContrast(d.Label, NotInStratum) { Id = d.Id });
                    var missing = d.At.Where(a => !levelsHere[a.Factor].Contains(a.Level, StringComparer.Ordinal)).ToList();
                    if (d.HasBothSidesIn(levelsHere) && missing.Count > 0)
                        warnings.Add($"Contrast {d.Id} ({d.Label}) is {f.Name}'s effect at " +
                                     $"{string.Join(" and ", d.At.Select(a => $"{a.Factor}={a.Level}"))}, where it interacts, and " +
                                     $"stratum '{stratum}' has no sample with {string.Join(" or ", missing.Select(a => $"{a.Factor}={a.Level}"))}, " +
                                     "so it is not run there.");
                }
            }
            else if (levels.Count < 2)
                notRun.Add(new NotRunContrast($"{f.Name} (one level: {levels[0]})", NotInStratum));
            else if (f.ReferenceLevel is not null && !levels.Contains(f.ReferenceLevel, StringComparer.Ordinal))
                foreach (var l in levels)
                    notRun.Add(new NotRunContrast(BetweenLabel(f.Name, l, f.ReferenceLevel), NotInStratum));
        }
        return notRun;
    }

    /// <summary>
    /// One contrast of an analysis's default set: <see cref="Numerator"/> against <see cref="Denominator"/> of
    /// <see cref="Factor"/>, or the slope of <see cref="Covariate"/>.
    /// </summary>
    private sealed record DefaultContrast(string Id, string? Factor = null, string? Numerator = null,
        string? Denominator = null, string? Covariate = null)
    {
        /// <summary>
        /// For a factor in an interaction, the level of each factor it interacts with that its effect is taken at: that
        /// factor's baseline over the whole analysis.
        /// </summary>
        public IReadOnlyList<(string Factor, string Level)> At { get; init; } = Array.Empty<(string, string)>();

        /// <summary>A factor contrast's label, as <see cref="Between"/> writes it.</summary>
        public string Label => BetweenLabel(Factor!, Numerator!, Denominator!);

        /// <summary>Whether a stratum with these levels has a sample at both of a factor contrast's levels.</summary>
        public bool HasBothSidesIn(Dictionary<string, List<string>> levels) =>
            levels[Factor!].Contains(Numerator!, StringComparer.Ordinal)
            && levels[Factor!].Contains(Denominator!, StringComparer.Ordinal);

        /// <summary>
        /// Whether a stratum with these levels runs the contrast: a sample at both of its levels, and at each level
        /// <see cref="At"/> names, so that it is the same comparison as in every other stratum. A stratum with only those
        /// levels of an interacting factor drops the interaction and still compares at them. Every stratum has every
        /// covariate.
        /// </summary>
        public bool RunsIn(Dictionary<string, List<string>> levels) =>
            Covariate is not null
            || HasBothSidesIn(levels) && At.All(a => levels[a.Factor].Contains(a.Level, StringComparer.Ordinal));
    }

    private sealed record FactorCoding(string Name, List<string> Levels, bool BatchOf = false)
    {
        public static FactorCoding? Build(FactorSpec spec, List<AnalysisSample> rows, string stratum,
            List<string> warnings)
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
                return null;
            }

            if (spec.ReferenceLevel is null)
            {
                warnings.Add($"No reference level was declared for factor '{spec.Name}', so every pair of its levels is " +
                             $"compared (GR-14). Levels: {string.Join(", ", levels)}.");
                return new FactorCoding(spec.Name, levels);
            }

            if (!levels.Contains(spec.ReferenceLevel, StringComparer.Ordinal))
            {
                warnings.Add($"The reference level '{spec.ReferenceLevel}' of factor '{spec.Name}' has no sample in stratum " +
                             $"'{stratum}', so no contrast against it is run there.");
                return new FactorCoding(spec.Name, levels);
            }

            levels.Remove(spec.ReferenceLevel);
            levels.Insert(0, spec.ReferenceLevel);
            return new FactorCoding(spec.Name, levels);
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
        /// <summary>The unit one step of the column is: the scaled unit, else <c>Scale</c> stored units ("10 years").</summary>
        public string PerUnit => Spec.ScaledUnit
            ?? (Spec.Scale == 1 ? Spec.Unit : string.Create(CultureInfo.InvariantCulture, $"{Spec.Scale:R} {Spec.Unit}"));

        public string ColumnName => $"{Spec.Name} (per {PerUnit})";

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
            // The mean, except that equal values are centred at that value: their floating-point mean can be a rounding
            // error off (ten 18.2s average 18.200000000000003), which would leave a column of 3.6e-15 instead of 0.
            double first = rows.Count == 0 ? double.NaN : rows[0].Covariates[spec.Name];
            double mean = rows.Count > 0 && rows.All(r => r.Covariates[spec.Name] == first)
                ? first
                : rows.Average(r => r.Covariates[spec.Name]);
            return new CovariateCoding(spec, spec.Centre ?? mean);
        }
    }
}
