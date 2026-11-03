using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MathNet.Numerics.Distributions;
using StatisticalModels;

namespace Quantification.Differential;

/// <summary>Reasons a dataset contrast of a hypothesis part is not run (<c>DEF-DIFF-SET-CONTRAST</c>, <c>-MODEL</c>).</summary>
public static class SetNotRunReason
{
    /// <summary>The <c>contexts</c> select samples of more than one stratum; STAT1 never pools strata (GR-5).</summary>
    public const string SpansStrata = "spans_strata";
    /// <summary>The stratum has more than one sample of one individual, which needs a third random level (not in v1).</summary>
    public const string RepeatedIndividuals = "repeated_individuals";
    /// <summary>The set model's fit failed.</summary>
    public const string NotConverged = "not_converged";
    /// <summary>
    /// The freeze names a factor or level the dataset's design does not have, followed by what was missing.
    /// Ruled by the user 2026-10-26 (DATA-DEFINITIONS v3.18).
    /// </summary>
    public const string NotInDesignPrefix = "not_in_design:";
    /// <summary>
    /// An <c>attenuated_*</c> prediction whose members' frozen mRNA change is missing or on the other side of 0.
    /// Ruled by the user 2026-10-26 (DATA-DEFINITIONS v3.18).
    /// </summary>
    public const string MrnaNotUsablePrefix = "mrna_not_usable:";
}

/// <summary>
/// One dataset contrast turned into a stratum, its samples and a fixed design (<c>DEF-DIFF-SET-CONTRAST</c>): stratum
/// factors in <c>contexts</c> select the stratum; other factors restrict the samples to their listed levels and enter the
/// model as an additive term (GR-1); the contrast is <c>factor=test_level</c> against <c>factor=reference_level</c>, run as
/// STAT1 runs a contrast in one stratum.
/// </summary>
public sealed record SetContrastSetup
{
    /// <summary>Why the contrast is not run, or null.</summary>
    public string? NotRun { get; init; }
    /// <summary>The stratum's key (<see cref="AnalysisDesign.StratumKey"/>); empty when not run before one was found.</summary>
    public string Stratum { get; init; } = "";
    /// <summary>The samples the model is fitted on: the stratum's, restricted by <c>contexts</c>; ordinal by id.</summary>
    public IReadOnlyList<AnalysisSample> Samples { get; init; } = Array.Empty<AnalysisSample>();
    /// <summary>The fixed design over <see cref="Samples"/>: STAT1's for the stratum, plus the contexts' factors.</summary>
    public AnalysisDesign? Design { get; init; }
    /// <summary>The contrast on <see cref="Design"/>.</summary>
    public Contrast? Contrast { get; init; }
    /// <summary>The samples at <c>test_level</c>.</summary>
    public IReadOnlyList<string> NumeratorSamples { get; init; } = Array.Empty<string>();
    /// <summary>The samples at <c>reference_level</c>.</summary>
    public IReadOnlyList<string> DenominatorSamples { get; init; } = Array.Empty<string>();

    /// <summary>The samples on either side of the contrast.</summary>
    public IReadOnlyList<string> ContrastSamples => NumeratorSamples.Concat(DenominatorSamples).Order(StringComparer.Ordinal).ToList();

    /// <summary>Builds the setup of one dataset contrast.</summary>
    /// <param name="contrast">The dataset contrast.</param>
    /// <param name="datasetSamples">Every sample of the dataset, with its strata, factors and individual.</param>
    /// <param name="spec">STAT1's fixed design for the dataset's analysis.</param>
    public static SetContrastSetup Build(SetDatasetContrast contrast, IReadOnlyList<AnalysisSample> datasetSamples, DesignSpec spec)
    {
        ArgumentNullException.ThrowIfNull(contrast);
        ArgumentNullException.ThrowIfNull(datasetSamples);
        ArgumentNullException.ThrowIfNull(spec);

        var stratumFactors = datasetSamples.SelectMany(s => s.Strata.Keys).ToHashSet(StringComparer.Ordinal);
        var modelFactors = datasetSamples.SelectMany(s => s.Factors.Keys).ToHashSet(StringComparer.Ordinal);
        var rows = datasetSamples.ToList();
        var factors = spec.Factors.ToList();
        foreach (var context in contrast.Contexts)
        {
            var levels = context.Levels.ToHashSet(StringComparer.Ordinal);
            if (stratumFactors.Contains(context.Factor))
                rows = rows.Where(s => s.Strata.TryGetValue(context.Factor, out var v) && levels.Contains(v)).ToList();
            else if (modelFactors.Contains(context.Factor))
            {
                rows = rows.Where(s => s.Factors.TryGetValue(context.Factor, out var v) && levels.Contains(v)).ToList();
                if (factors.All(f => f.Name != context.Factor))
                    factors.Add(new FactorSpec(context.Factor, context.Levels[0]));
            }
            else
                return new SetContrastSetup { NotRun = SetNotRunReason.NotInDesignPrefix + "context_factor=" + context.Factor };
        }

        var strata = rows.Select(s => AnalysisDesign.StratumKey(s.Strata)).Distinct(StringComparer.Ordinal).ToList();
        if (strata.Count == 0)
            return new SetContrastSetup { NotRun = SetNotRunReason.NotInDesignPrefix + "no_samples" };
        if (strata.Count > 1)
            return new SetContrastSetup { NotRun = SetNotRunReason.SpansStrata };
        string stratum = strata[0];

        bool repeated = datasetSamples.Where(s => AnalysisDesign.StratumKey(s.Strata) == stratum && s.Individual is not null)
            .GroupBy(s => s.Individual, StringComparer.Ordinal).Any(g => g.Count() > 1);
        if (repeated)
            return new SetContrastSetup { NotRun = SetNotRunReason.RepeatedIndividuals, Stratum = stratum };

        var numerator = rows.Where(s => s.Factors.TryGetValue(contrast.Factor, out var v) && v == contrast.TestLevel)
            .Select(s => s.SampleId).Order(StringComparer.Ordinal).ToList();
        var denominator = rows.Where(s => s.Factors.TryGetValue(contrast.Factor, out var v) && v == contrast.ReferenceLevel)
            .Select(s => s.SampleId).Order(StringComparer.Ordinal).ToList();
        if (!modelFactors.Contains(contrast.Factor))
            return new SetContrastSetup { NotRun = SetNotRunReason.NotInDesignPrefix + "factor=" + contrast.Factor, Stratum = stratum };
        if (numerator.Count == 0)
            return new SetContrastSetup { NotRun = SetNotRunReason.NotInDesignPrefix + $"{contrast.Factor}={contrast.TestLevel}", Stratum = stratum };
        if (denominator.Count == 0)
            return new SetContrastSetup { NotRun = SetNotRunReason.NotInDesignPrefix + $"{contrast.Factor}={contrast.ReferenceLevel}", Stratum = stratum };

        if (factors.All(f => f.Name != contrast.Factor))
            factors.Add(new FactorSpec(contrast.Factor, contrast.ReferenceLevel));
        var design = AnalysisDesign.Create(rows, spec with { Factors = factors }, stratum);
        var between = design.Between(contrast.Factor, contrast.TestLevel, contrast.ReferenceLevel);
        string estimability = design.Estimability(between);
        return new SetContrastSetup
        {
            NotRun = estimability == AnalysisDesign.Estimable ? null : estimability,
            Stratum = stratum,
            Samples = design.Samples,
            Design = design,
            Contrast = between,
            NumeratorSamples = numerator,
            DenominatorSamples = denominator,
        };
    }
}

/// <summary>A member's state after STAT1's row and the values are read (<c>DEF-DIFF-SET-MEMBER</c>), with its counts.</summary>
public sealed record SetMemberResult
{
    /// <summary>The member's resolution (<see cref="SetMemberResolver"/>), with its final <see cref="SetMemberAssignment.State"/>.</summary>
    public required SetMemberAssignment Assignment { get; init; }
    /// <summary>The member's STAT1 row, when it has a group.</summary>
    public DifferentialResult? Stat1Row { get; init; }
    /// <summary>Numerator-side samples in which one of the group's unique peptides has a value.</summary>
    public int SamplesWithValueNumerator { get; init; }
    /// <summary>Reference-side samples in which one of the group's unique peptides has a value.</summary>
    public int SamplesWithValueDenominator { get; init; }
    /// <summary>The member enters the set's estimate β̄ (or is the single's estimate).</summary>
    public bool InEstimate { get; init; }

    /// <summary>The member's final state.</summary>
    public string State => Assignment.State!;
}

/// <summary>The states read after member resolution (<c>DEF-DIFF-SET-MEMBER</c>).</summary>
public static class SetMemberStateFromRow
{
    /// <summary>STAT1's row for the group is <c>below_support:*</c>.</summary>
    public const string TooFew = "too_few";
    /// <summary>STAT1's row is <c>absent_in_denominator</c>: values only at the test level.</summary>
    public const string SeenNumeratorOnly = "seen_one_side:numerator";
    /// <summary>STAT1's row is <c>absent_in_numerator</c>: values only at the reference level.</summary>
    public const string SeenDenominatorOnly = "seen_one_side:denominator";
    /// <summary>No-change kinds only: values in fewer than half, or fewer than 2, of the samples on either side.</summary>
    public const string Incomplete = "incomplete";
    /// <summary>None of the above.</summary>
    public const string Tested = "tested";

    /// <summary>
    /// Reads each resolved member's final state, in the contract's order: <see cref="TooFew"/>, the two
    /// <c>seen_one_side</c> states, <see cref="Incomplete"/> (no-change predictions only), else <see cref="Tested"/>.
    /// </summary>
    /// <remarks>
    /// <b>Readings ruled by the user 2026-10-26</b> (DATA-DEFINITIONS v3.18; v3.12 mapped only <c>below_support:*</c>
    /// and <c>absent_in_*</c>):
    /// a STAT1 row <c>absent_in_both</c> is <see cref="TooFew"/>; a row <c>not_converged</c> or <c>not_estimable:*</c> does
    /// not stop the member, whose values still enter the set's own model (a single, which is STAT1's row, is then not run).
    /// </remarks>
    /// <param name="resolved">The members as <see cref="SetMemberResolver.Resolve"/> left them.</param>
    /// <param name="prediction">The part's prediction.</param>
    /// <param name="setup">The contrast's setup (sides).</param>
    /// <param name="normalized">The stratum's normalized observation table (<c>DEF-DIFF-NORM</c>).</param>
    /// <param name="stat1Rows">STAT1's rows for this setup and contrast, by protein group (<c>feature_id</c>).</param>
    /// <exception cref="ArgumentException">A member's group has no STAT1 row.</exception>
    public static IReadOnlyList<SetMemberResult> Apply(IReadOnlyList<SetMemberAssignment> resolved, SetPrediction prediction,
        SetContrastSetup setup, ObservationTable normalized, IReadOnlyDictionary<string, DifferentialResult> stat1Rows)
    {
        ArgumentNullException.ThrowIfNull(resolved);
        ArgumentNullException.ThrowIfNull(setup);
        ArgumentNullException.ThrowIfNull(normalized);
        ArgumentNullException.ThrowIfNull(stat1Rows);
        bool noChange = SetVerdict.IsNoChange(prediction);
        var result = new List<SetMemberResult>();
        foreach (var a in resolved)
        {
            if (a.ProteinGroup is null)
            {
                result.Add(new SetMemberResult { Assignment = a });
                continue;
            }
            if (!stat1Rows.TryGetValue(a.ProteinGroup, out var row))
                throw new ArgumentException($"Protein group '{a.ProteinGroup}' (member {a.Member}) has no STAT1 row.", nameof(stat1Rows));

            int num = SamplesWithValue(normalized, a.ProteinGroup, setup.NumeratorSamples);
            int den = SamplesWithValue(normalized, a.ProteinGroup, setup.DenominatorSamples);
            string state = row.Status.StartsWith("below_support:", StringComparison.Ordinal) || row.Status == DifferentialStatus.AbsentInBoth
                ? TooFew
                : row.Status == DifferentialStatus.AbsentInDenominator ? SeenNumeratorOnly
                : row.Status == DifferentialStatus.AbsentInNumerator ? SeenDenominatorOnly
                : noChange && (IncompleteSide(num, setup.NumeratorSamples.Count) || IncompleteSide(den, setup.DenominatorSamples.Count))
                    ? Incomplete
                    : Tested;
            result.Add(new SetMemberResult
            {
                Assignment = a with { State = state },
                Stat1Row = row,
                SamplesWithValueNumerator = num,
                SamplesWithValueDenominator = den,
                // Ruled by the user 2026-10-26: an incomplete member stays in the estimate (ruling 3 of v3.12 declined
                // leaving it out) and blocks supported.
                InEstimate = state is Tested or Incomplete,
            });
        }
        return result;
    }

    /// <summary>Fewer than half, or fewer than 2, of the side's samples have a value.</summary>
    internal static bool IncompleteSide(int withValue, int sideSize) => withValue < 2 || 2 * withValue < sideSize;

    private static int SamplesWithValue(ObservationTable table, string group, IReadOnlyList<string> samples)
    {
        var columns = samples.Select(s => table.Samples.ToList().IndexOf(s)).Where(i => i >= 0).ToList();
        var unique = Enumerable.Range(0, table.Peptides.Count)
            .Where(p => string.Equals(table.Peptides[p].UniqueProteinGroup, group, StringComparison.Ordinal)).ToList();
        return columns.Count(s => unique.Any(p => double.IsFinite(table.Log2Intensity(p, s))));
    }
}

/// <summary>A set's or a single's estimate: β̄, its standard error and degrees of freedom (<c>DEF-DIFF-SET-MODEL</c>).</summary>
public sealed record SetEstimate
{
    /// <summary>β̄: the equal-weight average of the members' log2 effects (a single: the member's).</summary>
    public required double Log2Effect { get; init; }
    /// <summary>Its standard error.</summary>
    public required double StandardError { get; init; }
    /// <summary>Its degrees of freedom.</summary>
    public required double Df { get; init; }
    /// <summary><c>satterthwaite</c> (set model), or the single's STAT1 <c>df_method</c>.</summary>
    public required string DfMethod { get; init; }
    /// <summary><c>set_mixed_model</c> or <c>stat1_row</c>.</summary>
    public required string ModelUsed { get; init; }
    /// <summary>The random terms fitted, e.g. <c>sample;sample:member</c>; empty for a single.</summary>
    public string RandomTerms { get; init; } = "";
    /// <summary>The protein groups that entered, each once.</summary>
    public int Terms { get; init; }
    /// <summary>Distinct peptides that entered.</summary>
    public int Peptides { get; init; }
    /// <summary>(peptide, sample) values that entered.</summary>
    public int Observations { get; init; }
    /// <summary>Residual variance of the set model; NaN for a single.</summary>
    public double ResidualVariance { get; init; } = double.NaN;
    /// <summary>Variance of the shared per-sample shift, <c>(1|sample)</c>; NaN for a single.</summary>
    public double SampleVariance { get; init; } = double.NaN;
    /// <summary>Variance of the per-protein sample term, <c>(1|sample:member)</c>; NaN for a single.</summary>
    public double SampleMemberVariance { get; init; } = double.NaN;
    /// <summary>
    /// Whether outlying values were down-weighted (GR-8): the set model's robust setting, or a single's STAT1 row's
    /// <see cref="DifferentialResult.Robust"/>.
    /// </summary>
    public bool? Robust { get; init; }
    /// <summary>
    /// Each value that entered the set model, with the weight it carried in the final fit, in the model's order (groups as
    /// given, then peptides and samples in ordinal order): 1 without robust weighting. Empty for a single.
    /// </summary>
    public IReadOnlyList<SetValueWeight> Weights { get; init; } = Array.Empty<SetValueWeight>();
}

/// <summary>One (peptide, sample) value of the set model and the weight it carried in the final fit (GR-8).</summary>
/// <param name="ProteinGroup">The member term the value belongs to.</param>
/// <param name="Peptide">The peptide's full sequence.</param>
/// <param name="Sample">The sample.</param>
/// <param name="Weight">1 without robust weighting; with it, the converged Huber weight, in (0, 1].</param>
public sealed record SetValueWeight(string ProteinGroup, string Peptide, string Sample, double Weight);

/// <summary>
/// The set model (<c>DEF-DIFF-SET-MODEL</c>): one mixed model over the tested members' unique peptides,
/// <c>log2 value ~ peptide + member:(the stratum's fixed design) + (1|sample) + (1|sample:member)</c>, REML with robust
/// weighting as STAT1's setting for the stratum (GR-8), no moderation.
/// It is <see cref="NestedMixedModel"/>'s <see cref="RandomStructure.IndividualAndSample"/> with each sample as the
/// "individual" and each member as the "sample" within it. β̄ is the equal-weight contrast over the members' contrasts;
/// its standard error is from the fit and its df Satterthwaite's.
/// </summary>
public static class SetModel
{
    /// <summary>The set model's name in <c>model_used</c>.</summary>
    public const string SetMixedModel = "set_mixed_model";
    /// <summary>A single's estimate, taken from its STAT1 row.</summary>
    public const string Stat1Row = "stat1_row";

    /// <summary>
    /// Fits the set model. Each protein group in <paramref name="groups"/> is one member term, so a group several members
    /// share through the all-members exception enters once.
    /// </summary>
    /// <remarks>
    /// A fit that fails gives no estimate and <paramref name="notRun"/> <c>not_converged</c>. Under robust weighting that
    /// includes reweighting that cannot finish (its round limit, a reweighted refit that fails, or residuals whose MAD is
    /// zero): the set is not run rather than refitted without weights, so every set of a family is fitted one way.
    /// </remarks>
    /// <param name="setup">The contrast's setup; must be run.</param>
    /// <param name="groups">The protein groups in the estimate, each once, in order.</param>
    /// <param name="normalized">The stratum's normalized observation table; only <paramref name="setup"/>'s samples are read.</param>
    /// <param name="notRun">Why no estimate came out, when it is null.</param>
    /// <param name="robust">
    /// STAT1's robust setting for the stratum (<see cref="DifferentialSettings.Robust"/>, GR-8), on by default as STAT1's:
    /// down-weight outlying values with <see cref="NestedMixedModel"/>'s robust option (msqrob2's Huber weights, iterated
    /// to convergence, GR-22), then REML and Satterthwaite on the weighted fit.
    /// </param>
    public static SetEstimate? FitSet(SetContrastSetup setup, IReadOnlyList<string> groups, ObservationTable normalized,
        out string? notRun, bool robust = true)
    {
        ArgumentNullException.ThrowIfNull(setup);
        ArgumentNullException.ThrowIfNull(groups);
        ArgumentNullException.ThrowIfNull(normalized);
        if (setup.NotRun is not null || setup.Design is null || setup.Contrast is null)
            throw new ArgumentException("The setup is not run, so there is no model to fit.", nameof(setup));
        if (groups.Distinct(StringComparer.Ordinal).Count() != groups.Count)
            throw new ArgumentException("A group enters a set at most once.", nameof(groups));

        var design = setup.Design;
        double[,] x = design.FittedMatrix;
        double[] w = design.FittedWeights(setup.Contrast);
        int p = design.FittedColumns.Count, m = groups.Count;
        var sampleRow = design.Samples.Select((s, i) => (s.SampleId, i)).ToDictionary(t => t.SampleId, t => t.i, StringComparer.Ordinal);
        var tableColumn = normalized.Samples.Select((s, i) => (s, i)).ToDictionary(t => t.s, t => t.i, StringComparer.Ordinal);

        var response = new List<double>();
        var commonRows = new List<double[]>();
        var nuisanceIndex = new List<int>(); // -1: the member's first peptide (no column)
        var sampleOf = new List<string>();
        var memberOf = new List<string>();
        var peptideOf = new List<string>();
        int nuisanceColumns = 0;
        int peptides = 0;
        for (int g = 0; g < m; g++)
        {
            var unique = Enumerable.Range(0, normalized.Peptides.Count)
                .Where(k => string.Equals(normalized.Peptides[k].UniqueProteinGroup, groups[g], StringComparison.Ordinal))
                .Where(k => design.Samples.Any(s => tableColumn.TryGetValue(s.SampleId, out int c) && double.IsFinite(normalized.Log2Intensity(k, c))))
                .OrderBy(k => normalized.Peptides[k].FullSequence, StringComparer.Ordinal)
                .ToList();
            for (int u = 0; u < unique.Count; u++)
            {
                int column = u == 0 ? -1 : nuisanceColumns++;
                peptides++;
                foreach (var sample in design.Samples)
                {
                    if (!tableColumn.TryGetValue(sample.SampleId, out int c))
                        throw new ArgumentException($"The table has no sample '{sample.SampleId}'.", nameof(normalized));
                    double value = normalized.Log2Intensity(unique[u], c);
                    if (!double.IsFinite(value)) continue;
                    var row = new double[m * p];
                    int r = sampleRow[sample.SampleId];
                    for (int j = 0; j < p; j++) row[g * p + j] = x[r, j];
                    response.Add(value);
                    commonRows.Add(row);
                    nuisanceIndex.Add(column);
                    sampleOf.Add(sample.SampleId);
                    memberOf.Add(groups[g]);
                    peptideOf.Add(normalized.Peptides[unique[u]].FullSequence);
                }
            }
        }

        int n = response.Count;
        var common = new double[n, m * p];
        var nuisance = nuisanceColumns > 0 ? new double[n, nuisanceColumns] : null;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m * p; j++) common[i, j] = commonRows[i][j];
            if (nuisanceIndex[i] >= 0) nuisance![i, nuisanceIndex[i]] = 1;
        }
        var names = groups.SelectMany(g => design.FittedColumns.Select(j => $"member={g}:{design.ColumnNames[j]}")).ToList();
        var problem = new NestedMixedModelProblem(response, common, nuisance, memberOf, sampleOf);
        var fit = NestedMixedModel.Fit(new[] { problem }, names, maxThreads: 1, robust: robust);
        if (fit.Status[0] != FeatureFitStatus.Fitted)
        {
            notRun = fit.Status[0] == FeatureFitStatus.NotConverged
                ? SetNotRunReason.NotConverged
                : SetNotRunReason.NotConverged + ":" + SnakeCase(fit.Status[0].ToString());
            return null;
        }

        // β̄: each member's contrast, weighted 1/m.
        var weights = new double[m * p];
        for (int g = 0; g < m; g++)
            for (int j = 0; j < p; j++)
                weights[g * p + j] = w[j] / m;
        double estimate = 0, quad = 0;
        for (int i = 0; i < m * p; i++)
        {
            estimate += weights[i] * fit.Coefficient(0, i);
            for (int j = 0; j < m * p; j++) quad += weights[i] * fit.UnscaledCovariance(0, i, j) * weights[j];
        }
        notRun = null;
        var structure = fit.Structure[0];
        var fitted = fit.Weights(0)!;
        return new SetEstimate
        {
            Log2Effect = estimate,
            StandardError = Math.Sqrt(quad * fit.ResidualVariance[0]),
            Df = fit.SatterthwaiteDf(0, weights),
            DfMethod = "satterthwaite",
            ModelUsed = SetMixedModel,
            RandomTerms = structure == RandomStructure.IndividualAndSample ? "sample;sample:member" : "sample",
            Terms = m,
            Peptides = peptides,
            Observations = n,
            ResidualVariance = fit.ResidualVariance[0],
            SampleVariance = fit.IndividualVariance[0],
            SampleMemberVariance = structure == RandomStructure.IndividualAndSample ? fit.SampleVariance[0] : double.NaN,
            Robust = robust,
            Weights = Enumerable.Range(0, n).Select(i => new SetValueWeight(memberOf[i], peptideOf[i], sampleOf[i], fitted[i])).ToList(),
        };
    }

    /// <summary>A single's estimate: its member's STAT1 row, not refitted (<c>DEF-DIFF-SET-MODEL</c>).</summary>
    /// <exception cref="ArgumentException">The row carries no effect.</exception>
    public static SetEstimate FromStat1Row(DifferentialResult row)
    {
        ArgumentNullException.ThrowIfNull(row);
        if (!DifferentialStatus.IsFitted(row.Status) || row.Log2Effect is null || row.StandardError is null || row.Df is null)
            throw new ArgumentException($"STAT1's row for '{row.FeatureId}' is {row.Status} and carries no effect.", nameof(row));
        return new SetEstimate
        {
            Log2Effect = row.Log2Effect.Value,
            StandardError = row.StandardError.Value,
            Df = row.Df.Value,
            DfMethod = row.DfMethod ?? "",
            ModelUsed = Stat1Row,
            Terms = 1,
            Peptides = row.NPeptides ?? 0,
            Observations = row.NObservations ?? 0,
            Robust = row.Robust,
        };
    }

    private static string SnakeCase(string name) =>
        string.Concat(name.Select((ch, i) => char.IsUpper(ch) && i > 0 ? "_" + char.ToLowerInvariant(ch) : char.ToLowerInvariant(ch).ToString()));
}

/// <summary>
/// One prediction's test on an estimate (<c>DEF-DIFF-SET-VERDICT</c>): the p (one-sided, or the larger of two one-sided
/// p's), the 90% interval, the uncorrected <c>contradicted</c> rule and <c>smallest_margin_passed</c>.
/// </summary>
public sealed record SetVerdict
{
    /// <summary>The default δ, log2.</summary>
    public const double DefaultMargin = 0.5;
    /// <summary>The default f.</summary>
    public const double DefaultAttenuation = 0.5;
    /// <summary>The level of the 90% interval's one-sided tails, 0.05.</summary>
    public const double Alpha = 0.05;

    /// <summary>The prediction tested.</summary>
    public required SetPrediction Prediction { get; init; }
    /// <summary>The estimate.</summary>
    public required SetEstimate Estimate { get; init; }
    /// <summary>δ as used, or null where the prediction has none.</summary>
    public double? Margin { get; init; }
    /// <summary>f as used, or null where the prediction has none.</summary>
    public double? Attenuation { get; init; }
    /// <summary>m̄ as used, or null where the prediction has none.</summary>
    public double? MrnaMean { get; init; }
    /// <summary>The null hypothesis's boundary or boundaries, e.g. <c>0</c>, <c>+0.5</c>, <c>±0.5</c>, <c>0;-1.1</c>.</summary>
    public required string Bound { get; init; }
    /// <summary>The p: one-sided, or the larger of two one-sided p's.</summary>
    public required double PValue { get; init; }
    /// <summary>Lower end of the 90% interval.</summary>
    public required double Ci90Low { get; init; }
    /// <summary>Upper end of the 90% interval.</summary>
    public required double Ci90High { get; init; }
    /// <summary>The 90% interval lies wholly where the prediction is refuted (not adjusted for multiplicity).</summary>
    public required bool Contradicted { get; init; }
    /// <summary>The smallest δ (or f) at which the unadjusted interval rule would pass; null for up, down, or none passing.</summary>
    public double? SmallestMarginPassed { get; init; }

    /// <summary>unchanged, not_up, not_down, attenuated_down and attenuated_up: the predictions of no change or of less change.</summary>
    public static bool IsNoChange(SetPrediction prediction) => prediction is not (SetPrediction.Up or SetPrediction.Down);

    /// <summary>Tests <paramref name="prediction"/> on <paramref name="estimate"/>.</summary>
    /// <param name="prediction">The prediction.</param>
    /// <param name="estimate">β̄ with its se and df.</param>
    /// <param name="margin">δ; <see cref="DefaultMargin"/> when null.</param>
    /// <param name="attenuation">f; <see cref="DefaultAttenuation"/> when null.</param>
    /// <param name="mrnaMean">m̄; required for the attenuated predictions, with the prediction's sign.</param>
    public static SetVerdict Test(SetPrediction prediction, SetEstimate estimate, double? margin = null, double? attenuation = null,
        double? mrnaMean = null)
    {
        ArgumentNullException.ThrowIfNull(estimate);
        double b = estimate.Log2Effect, se = estimate.StandardError, df = estimate.Df;
        if (!double.IsFinite(b) || !(se > 0) || !(df > 0))
            throw new ArgumentException("The estimate needs a finite effect, a positive standard error and positive df.", nameof(estimate));
        double d = margin ?? DefaultMargin, f = attenuation ?? DefaultAttenuation;
        if (!(d > 0) || !double.IsFinite(d)) throw new ArgumentOutOfRangeException(nameof(margin), "δ must be positive.");
        if (!(f > 0 && f < 1)) throw new ArgumentOutOfRangeException(nameof(attenuation), "f must lie strictly between 0 and 1.");
        bool attenuated = prediction is SetPrediction.AttenuatedDown or SetPrediction.AttenuatedUp;
        if (attenuated && (mrnaMean is null || !double.IsFinite(mrnaMean.Value)
                           || (prediction == SetPrediction.AttenuatedDown ? mrnaMean >= 0 : mrnaMean <= 0)))
            throw new ArgumentException($"{SetTestFreeze.PredictionName(prediction)} needs m̄ on its side of 0.", nameof(mrnaMean));

        double Cdf(double t) => double.IsPositiveInfinity(df) ? Normal.CDF(0, 1, t) : StudentT.CDF(0, 1, df, t);
        double Upper(double bound) => 1 - Cdf((b - bound) / se); // H0: β̄ ≤ bound
        double Lower(double bound) => Cdf((b - bound) / se);     // H0: β̄ ≥ bound
        double tq = double.IsPositiveInfinity(df) ? Normal.InvCDF(0, 1, 1 - Alpha) : StudentT.InvCDF(0, 1, df, 1 - Alpha);
        double low = b - tq * se, high = b + tq * se;
        double fm = attenuated ? f * mrnaMean!.Value : double.NaN;
        string R(double v) => v.ToString("R", CultureInfo.InvariantCulture);

        (double p, bool contradicted, string bound, double? smallest) = prediction switch
        {
            SetPrediction.Up => (Upper(0), high < 0, "0", (double?)null),
            SetPrediction.Down => (Lower(0), low > 0, "0", null),
            SetPrediction.Unchanged => (Math.Max(Upper(-d), Lower(d)), low > d || high < -d, "±" + R(d),
                Math.Max(Math.Abs(low), Math.Abs(high))),
            SetPrediction.NotUp => (Lower(d), low > d, "+" + R(d), Math.Max(high, 0)),
            SetPrediction.NotDown => (Upper(-d), high < -d, "-" + R(d), Math.Max(-low, 0)),
            SetPrediction.AttenuatedDown => (Math.Max(Lower(0), Upper(fm)), low >= 0 || high < fm, "0;" + R(fm),
                high < 0 ? low / mrnaMean!.Value : null),
            SetPrediction.AttenuatedUp => (Math.Max(Upper(0), Lower(fm)), high <= 0 || low > fm, "0;" + R(fm),
                low > 0 ? high / mrnaMean!.Value : null),
            _ => throw new ArgumentOutOfRangeException(nameof(prediction)),
        };
        return new SetVerdict
        {
            Prediction = prediction,
            Estimate = estimate,
            Margin = prediction is SetPrediction.Unchanged or SetPrediction.NotUp or SetPrediction.NotDown ? d : null,
            Attenuation = attenuated ? f : null,
            MrnaMean = attenuated ? mrnaMean : null,
            Bound = bound,
            PValue = p,
            Ci90Low = low,
            Ci90High = high,
            Contradicted = contradicted,
            SmallestMarginPassed = smallest,
        };
    }
}

/// <summary>The verdict words (<c>DEF-DIFF-SET-VERDICT</c>).</summary>
public static class SetVerdictName
{
    /// <summary>The adjusted p is at most 0.05.</summary>
    public const string Supported = "supported";
    /// <summary>The 90% interval lies wholly where the prediction is refuted; not adjusted for multiplicity.</summary>
    public const string Contradicted = "contradicted";
    /// <summary>Anything else, and any guard.</summary>
    public const string NotInformative = "not_informative";
    /// <summary>Fewer than 3 tested members (a set), or no tested member (a single).</summary>
    public const string NotTestable = "not_testable";
    /// <summary>A reason from the contrast or the model.</summary>
    public const string NotRun = "not_run";
}

/// <summary>
/// One hypothesis part on one dataset contrast, tested (<c>DEF-DIFF-SET-MODEL</c>, <c>-VERDICT</c>): the members'
/// states, the estimate, the test, and the guards. <see cref="Verdict"/> needs the family's adjusted p for
/// <c>supported</c>, which pooling and the family give.
/// </summary>
public sealed record SetContrastOutcome
{
    /// <summary>The part.</summary>
    public required SetHypothesisPart Part { get; init; }
    /// <summary>The dataset contrast.</summary>
    public required SetDatasetContrast Contrast { get; init; }
    /// <summary>The setup.</summary>
    public required SetContrastSetup Setup { get; init; }
    /// <summary>Every member of the contrast, with its final state; empty when the contrast is not run.</summary>
    public IReadOnlyList<SetMemberResult> Members { get; init; } = Array.Empty<SetMemberResult>();
    /// <summary><see cref="SetVerdictName.NotRun"/> or <see cref="SetVerdictName.NotTestable"/> before any test; null when tested.</summary>
    public string? Untested { get; init; }
    /// <summary>The reason for <see cref="Untested"/>.</summary>
    public string? UntestedReason { get; init; }
    /// <summary>The test, when there was one.</summary>
    public SetVerdict? Test { get; init; }
    /// <summary>The members whose state blocks <c>supported</c> (the guards), in member order.</summary>
    public IReadOnlyList<string> BlockedBy { get; init; } = Array.Empty<string>();
    /// <summary>For up and down: <c>seen_one_side</c> members seen on the predicted side.</summary>
    public int OneSideAgree { get; init; }
    /// <summary>For up and down: <c>seen_one_side</c> members seen on the other side.</summary>
    public int OneSideDisagree { get; init; }

    /// <summary>
    /// The verdict and its detail, given the family's Benjamini–Hochberg adjusted p (null outside the family): untested
    /// first, then <c>contradicted</c> (detail <c>uncorrected</c>), then a guard (<c>not_informative</c>, naming the members),
    /// then <c>supported</c> when the adjusted p is at most 0.05, else <c>not_informative</c>.
    /// </summary>
    public (string Verdict, string Detail) Verdict(double? pAdjusted)
    {
        if (Untested is not null) return (Untested, UntestedReason ?? "");
        if (Test!.Contradicted) return (SetVerdictName.Contradicted, "uncorrected");
        if (BlockedBy.Count > 0) return (SetVerdictName.NotInformative, "blocked_by:" + string.Join(";", BlockedBy));
        if (pAdjusted is <= SetVerdict.Alpha) return (SetVerdictName.Supported, "bh_adjusted_p<=0.05");
        return (SetVerdictName.NotInformative, pAdjusted is null ? "not_in_family" : "bh_adjusted_p>0.05");
    }
}

/// <summary>Runs one hypothesis part on one dataset contrast (<c>DEF-DIFF-SET-MODEL</c>, <c>-VERDICT</c>).</summary>
public static class SetTest
{
    /// <summary>The fewest tested members a set needs.</summary>
    public const int MinimumTestedMembers = 3;

    /// <summary>Tests one part on one dataset contrast.</summary>
    /// <param name="part">The hypothesis part.</param>
    /// <param name="contrast">One of its dataset contrasts.</param>
    /// <param name="setup">The contrast's setup (<see cref="SetContrastSetup.Build"/>).</param>
    /// <param name="resolved">Its members, resolved (<see cref="SetMemberResolver.Resolve"/>) on the same samples.</param>
    /// <param name="normalized">The stratum's normalized observation table.</param>
    /// <param name="stat1Rows">STAT1's rows for this setup and contrast, by protein group.</param>
    /// <param name="robust">
    /// STAT1's robust setting for the stratum (<see cref="DifferentialSettings.Robust"/>, GR-8), on by default as STAT1's;
    /// the set model is fitted with it (<see cref="SetModel.FitSet"/>). A single takes its STAT1 row as it is.
    /// </param>
    public static SetContrastOutcome Run(SetHypothesisPart part, SetDatasetContrast contrast, SetContrastSetup setup,
        IReadOnlyList<SetMemberAssignment> resolved, ObservationTable normalized, IReadOnlyDictionary<string, DifferentialResult> stat1Rows,
        bool robust = true)
    {
        ArgumentNullException.ThrowIfNull(part);
        ArgumentNullException.ThrowIfNull(contrast);
        ArgumentNullException.ThrowIfNull(setup);
        var outcome = new SetContrastOutcome { Part = part, Contrast = contrast, Setup = setup };
        if (setup.NotRun is not null)
            return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = setup.NotRun };

        var members = SetMemberStateFromRow.Apply(resolved, part.Prediction, setup, normalized, stat1Rows);
        outcome = outcome with { Members = members };
        var inEstimate = members.Where(x => x.InEstimate).ToList();

        // Margins are frozen per hypothesis: every member row of the contrast must agree.
        double? margin = One(contrast.Rows.Select(r => r.MarginLog2), "margin_log2", contrast);
        double? attenuation = One(contrast.Rows.Select(r => r.AttenuationFraction), "attenuation_fraction", contrast);
        double? mrnaMean = null;
        if (part.Prediction is SetPrediction.AttenuatedDown or SetPrediction.AttenuatedUp)
        {
            // Ruled by the user 2026-10-26: m̄ over the members in the estimate, so the bound matches what β̄ averages.
            var names = inEstimate.Select(x => x.Assignment.Member).ToHashSet(StringComparer.Ordinal);
            var mrna = contrast.Rows.Where(r => names.Contains(r.Member)).Select(r => r.MrnaLog2).ToList();
            if (mrna.Count > 0 && mrna.Any(v => v is null))
                return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = SetNotRunReason.MrnaNotUsablePrefix + "missing" };
            if (mrna.Count > 0)
            {
                mrnaMean = mrna.Average(v => v!.Value);
                if (part.Prediction == SetPrediction.AttenuatedDown ? mrnaMean >= 0 : mrnaMean <= 0)
                    return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = SetNotRunReason.MrnaNotUsablePrefix + "direction" };
            }
        }

        SetEstimate? estimate;
        if (!part.IsSet)
        {
            if (members.Count != 1)
                throw new ArgumentException(
                    $"{part.HypothesisId} is tested as a single but lists {members.Count} members in {contrast.Dataset}.", nameof(contrast));
            var only = members[0];
            if (!only.InEstimate)
                return outcome with { Untested = SetVerdictName.NotTestable, UntestedReason = "no_tested_member" };
            if (!DifferentialStatus.IsFitted(only.Stat1Row!.Status))
                return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = only.Stat1Row.Status };
            estimate = SetModel.FromStat1Row(only.Stat1Row);
        }
        else
        {
            // Ruled by the user 2026-10-26: the minimum counts protein groups, each once, so one family group shared through
            // the all-members exception counts once.
            var groups = inEstimate.Select(x => x.Assignment.ProteinGroup!).Distinct(StringComparer.Ordinal).ToList();
            if (groups.Count < MinimumTestedMembers)
                return outcome with { Untested = SetVerdictName.NotTestable, UntestedReason = "fewer_than_3_tested" };
            estimate = SetModel.FitSet(setup, groups, normalized, out string? notRun, robust);
            if (estimate is null)
                return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = notRun };
        }

        var test = SetVerdict.Test(part.Prediction, estimate, margin, attenuation, mrnaMean);
        var blocked = new List<string>();
        int agree = 0, disagree = 0;
        foreach (var x in members)
        {
            string? state = x.Assignment.State;
            bool numeratorOnly = state == SetMemberStateFromRow.SeenNumeratorOnly;
            bool denominatorOnly = state == SetMemberStateFromRow.SeenDenominatorOnly;
            switch (part.Prediction)
            {
                case SetPrediction.Up or SetPrediction.Down:
                    bool up = part.Prediction == SetPrediction.Up;
                    if (numeratorOnly) { if (up) agree++; else disagree++; }
                    if (denominatorOnly) { if (up) disagree++; else agree++; }
                    break;
                case SetPrediction.NotUp:
                    if (numeratorOnly || state == SetMemberStateFromRow.Incomplete) blocked.Add(x.Assignment.Member);
                    break;
                case SetPrediction.NotDown:
                    if (denominatorOnly || state == SetMemberStateFromRow.Incomplete) blocked.Add(x.Assignment.Member);
                    break;
                default:
                    if (numeratorOnly || denominatorOnly || state == SetMemberStateFromRow.Incomplete) blocked.Add(x.Assignment.Member);
                    break;
            }
        }
        return outcome with { Test = test, BlockedBy = blocked, OneSideAgree = agree, OneSideDisagree = disagree };
    }

    private static double? One(IEnumerable<double?> values, string column, SetDatasetContrast contrast)
    {
        var distinct = values.Distinct().ToList();
        if (distinct.Count > 1)
            throw new ArgumentException(
                $"The member rows of {contrast.Dataset} (lines {string.Join(", ", contrast.Rows.Select(r => r.Line))}) give " +
                $"different '{column}' values; margins are frozen per hypothesis.", nameof(contrast));
        return distinct[0];
    }
}
