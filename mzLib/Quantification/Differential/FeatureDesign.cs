using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using StatisticalModels;

namespace Quantification.Differential;

/// <summary>
/// One feature's own design (STAT1 M5b ruling 1): the stratum's design columns on the feature's observations, beside the
/// feature's nuisance columns (its peptides), reduced to the columns those observations can estimate, so a protein missing
/// a group or a batch is fitted on what it has instead of being refused as a whole.
/// </summary>
/// <remarks>
/// <para>
/// <b>One rule, the design's.</b> Columns are kept left to right, common before nuisance, each unless it is redundant
/// with those already kept (<see cref="AnalysisDesign.KeptColumns"/>, LinearModel's own full-rank test), and a contrast
/// is estimable exactly when its weights are orthogonal to the null space of the columns divided by their lengths, judged
/// WITH the nuisance columns: a peptide seen in one group only confounds that group with the peptide. So every verdict
/// agrees with <see cref="AnalysisDesign.Estimability"/> on the same matrix.
/// </para>
/// <para>
/// <b>Judged by what is estimable, never by which column survived.</b> With a factor's reference level absent no column
/// is zero: the other levels' columns add up to the intercept, the later one is left out, and the survivor no longer means
/// "against the reference". A contrast is therefore always given on the full stratum columns and mapped by
/// <see cref="KeptWeights"/>, which is exact for an estimable contrast (the left-out coefficients can be taken as 0).
/// </para>
/// </remarks>
internal sealed class FeatureDesign
{
    private readonly double[] _lengths;
    private readonly Matrix<double>? _nullProjector;
    private readonly IReadOnlyList<string> _terms;

    private FeatureDesign(int commonCount, int nuisanceCount, IReadOnlyList<int> kept, double[] lengths,
        Matrix<double>? nullProjector, IReadOnlyList<string> terms)
    {
        CommonCount = commonCount;
        NuisanceCount = nuisanceCount;
        KeptColumns = kept;
        KeptCommon = kept.Where(j => j < commonCount).ToArray();
        KeptNuisance = kept.Where(j => j >= commonCount).Select(j => j - commonCount).ToArray();
        _lengths = lengths;
        _nullProjector = nullProjector;
        _terms = terms;
    }

    /// <summary>The stratum's design columns.</summary>
    public int CommonCount { get; }

    /// <summary>The feature's nuisance columns.</summary>
    public int NuisanceCount { get; }

    /// <summary>Kept columns of [common | nuisance], ascending.</summary>
    public IReadOnlyList<int> KeptColumns { get; }

    /// <summary>Kept common columns, as indices into the stratum's columns, ascending.</summary>
    public IReadOnlyList<int> KeptCommon { get; }

    /// <summary>Kept nuisance columns, as indices into the feature's nuisance columns, ascending.</summary>
    public IReadOnlyList<int> KeptNuisance { get; }

    /// <summary>Whether every column was kept.</summary>
    public bool IsFull => KeptColumns.Count == CommonCount + NuisanceCount;

    /// <summary>A key for the kept common columns: features with the same key can share one fit.</summary>
    public string CommonKey => string.Join(",", KeptCommon);

    /// <summary>Builds a feature's design from its observations.</summary>
    /// <param name="common">[observation, stratum column]: each observation's row of the stratum's design.</param>
    /// <param name="nuisance">[observation, nuisance column], or null for none.</param>
    /// <param name="commonTerms">Each stratum column's term (as <see cref="AnalysisDesign"/> names them).</param>
    /// <param name="nuisanceTerm">The term every nuisance column belongs to.</param>
    public static FeatureDesign Build(double[,] common, double[,]? nuisance, IReadOnlyList<string> commonTerms,
        string nuisanceTerm = "peptide")
    {
        ArgumentNullException.ThrowIfNull(common);
        ArgumentNullException.ThrowIfNull(commonTerms);
        int n = common.GetLength(0), p = common.GetLength(1), q = nuisance?.GetLength(1) ?? 0;
        if (commonTerms.Count != p)
            throw new ArgumentException($"There are {commonTerms.Count} terms for {p} common columns.", nameof(commonTerms));
        if (nuisance is not null && nuisance.GetLength(0) != n)
            throw new ArgumentException($"The nuisance columns have {nuisance.GetLength(0)} rows but there are {n} observations.",
                nameof(nuisance));

        var matrix = new double[n, p + q];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < p; j++) matrix[i, j] = common[i, j];
            for (int j = 0; j < q; j++) matrix[i, p + j] = nuisance![i, j];
        }
        var terms = commonTerms.Concat(Enumerable.Repeat(nuisanceTerm, q)).ToArray();
        IReadOnlyList<int> kept = n > 0 && LinearModel.HasFullColumnRank(matrix)
            ? Enumerable.Range(0, p + q).ToArray()
            : AnalysisDesign.KeptColumns(matrix);
        var lengths = Matrix<double>.Build.DenseOfArray(matrix).ColumnNorms(2).Select(l => l > 0 ? l : 1).ToArray();
        Matrix<double>? projector = null;
        if (kept.Count < p + q)
        {
            if (kept.Count == 0)
                projector = Matrix<double>.Build.DenseIdentity(p + q);
            else
            {
                var nullSpace = AnalysisDesign.NullSpace(matrix, lengths, kept);
                projector = nullSpace.TransposeAndMultiply(nullSpace);
            }
        }
        return new FeatureDesign(p, q, kept, lengths, projector, terms);
    }

    /// <summary>
    /// Whether a contrast on the stratum's columns can be estimated from this feature's observations:
    /// <see cref="AnalysisDesign.Estimable"/>, <see cref="AnalysisDesign.NotEstimableConfounded"/> (the direction it needs
    /// is shared by two or more terms, a peptide's included) or <see cref="AnalysisDesign.NotEstimableRankDeficient"/>.
    /// </summary>
    public string Estimability(IReadOnlyList<double> commonWeights)
    {
        var unseen = Unseen(commonWeights, out double length);
        if (unseen is null || unseen.L2Norm() <= AnalysisDesign.LoadingTolerance * length) return AnalysisDesign.Estimable;
        return AnalysisDesign.TermsAlong(unseen / length, _terms).Count > 1
            ? AnalysisDesign.NotEstimableConfounded
            : AnalysisDesign.NotEstimableRankDeficient;
    }

    /// <summary>An estimable contrast's weights on <see cref="KeptCommon"/> (the left-out coefficients taken as 0).</summary>
    public double[] KeptWeights(IReadOnlyList<double> commonWeights)
    {
        Check(commonWeights);
        return KeptCommon.Select(j => commonWeights[j]).ToArray();
    }

    /// <summary>
    /// The variance-df rows of a feature (ruling 7): the unit vector of each of <paramref name="testedColumns"/> that was
    /// kept, as weights on <see cref="KeptCommon"/>, so the joint test is of the reduced fit's own coefficients for the
    /// tested terms, as lme4 tests them when it drops the same columns (it keeps columns left to right too). With every
    /// column kept they are the tested columns' unit vectors; empty when no tested column was kept. A kept coefficient of a
    /// reduced fit is the fit's own parameter (with a reference level absent, the present levels against the last present
    /// one), never read as a contrast: contrasts come only through <see cref="KeptWeights"/> after <see cref="Estimability"/>.
    /// </summary>
    public IReadOnlyList<double[]> VarianceDfRows(IReadOnlyList<int> testedColumns)
    {
        ArgumentNullException.ThrowIfNull(testedColumns);
        return testedColumns.Where(j => KeptCommon.Contains(j))
            .Select(j => KeptCommon.Select(k => k == j ? 1.0 : 0.0).ToArray())
            .ToArray();
    }

    private Vector<double>? Unseen(IReadOnlyList<double> commonWeights, out double length)
    {
        Check(commonWeights);
        int width = CommonCount + NuisanceCount;
        var c = Vector<double>.Build.Dense(width, j => j < CommonCount ? commonWeights[j] / _lengths[j] : 0);
        length = c.L2Norm();
        return _nullProjector is null ? null : _nullProjector * c;
    }

    private void Check(IReadOnlyList<double> commonWeights)
    {
        ArgumentNullException.ThrowIfNull(commonWeights);
        if (commonWeights.Count != CommonCount)
            throw new ArgumentException($"The contrast has {commonWeights.Count} weights but the stratum has {CommonCount} columns.",
                nameof(commonWeights));
    }
}
