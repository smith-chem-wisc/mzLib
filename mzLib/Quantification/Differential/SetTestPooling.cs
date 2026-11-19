using System;
using System.Collections.Generic;
using System.Linq;
using StatisticalModels;

namespace Quantification.Differential;

/// <summary>
/// The words of pooling and the families (<c>DEF-DIFF-SET-POOL</c>, <c>-FAMILY</c>, <c>-POSTHOC</c>, <c>-MATURATION</c>).
/// </summary>
public static class SetPoolingName
{
    /// <summary><c>contrast_kind</c> of every contrast that is not a maturation contrast, and of a pooled row.</summary>
    public const string Aging = "aging";
    /// <summary><c>contrast_kind</c> of a maturation contrast, and the name of its family (<c>DEF-DIFF-SET-MATURATION</c>).</summary>
    public const string Maturation = "maturation";
    /// <summary><c>model_used</c> of a pooled row.</summary>
    public const string PooledModel = "random_effects_hartung_knapp";
    /// <summary><c>df_method</c> of a pooled row.</summary>
    public const string PooledDfMethod = "hartung_knapp";

    /// <summary>Left out: <c>independence</c> says the contrast reproduces an earlier result.</summary>
    public const string Reproduction = "reproduction";
    /// <summary>Left out: <c>independence</c> is neither independent nor a reproduction.</summary>
    public const string NotPooled = "not_pooled";
    /// <summary>Left out: a member row says <c>post_hoc = yes</c>; followed by those members, <c>;</c>-joined (v3.13).</summary>
    public const string PostHocPrefix = "post_hoc:";
    /// <summary>Left out: the part names several contrasts of this dataset and does not mark exactly one <c>primary = yes</c>.</summary>
    public const string SeveralContrastsPerDataset = "several_contrasts_per_dataset";
    /// <summary>Left out: <c>primary = no</c>, or another contrast of the dataset is the one marked <c>primary = yes</c>.</summary>
    public const string NotPrimary = "not_primary";
    /// <summary>The pooled row is not run: every pooled dataset gives the same estimate, so the Hartung–Knapp standard error is 0.</summary>
    public const string PooledStandardErrorZero = "pooled_standard_error_zero";
    /// <summary>A dataset row of a selected hypothesis: its Benjamini–Bogomolov adjusted p is at most the level.</summary>
    public const string BogomolovPassed = "bb_adjusted_p<=level";
    /// <summary>A dataset row of a selected hypothesis: its Benjamini–Bogomolov adjusted p is above the level.</summary>
    public const string BogomolovFailed = "bb_adjusted_p>level";
}

/// <summary>One dataset contrast of a part after pooling and the families: its place, its adjusted p and its verdict.</summary>
public sealed record SetDatasetResult
{
    /// <summary>The tested contrast (<see cref="SetTest.Run"/>).</summary>
    public required SetContrastOutcome Outcome { get; init; }
    /// <summary><see cref="SetPoolingName.Aging"/> or <see cref="SetPoolingName.Maturation"/>.</summary>
    public string ContrastKind => Outcome.Contrast.Maturation ? SetPoolingName.Maturation : SetPoolingName.Aging;
    /// <summary>Why the contrast cannot build a primary result (a <see cref="SetPoolingName"/> word); null when it can.</summary>
    public string? LeftOut { get; init; }
    /// <summary>One of the k ≥ 2 datasets of its part's pooled row.</summary>
    public bool Pooled { get; init; }
    /// <summary>
    /// The contrast is a primary result itself: the one dataset of an aging part (k = 1), or a maturation contrast in the
    /// maturation family.
    /// </summary>
    public bool IsPrimary { get; init; }
    /// <summary>The family its p is in (<see cref="SetPoolingName.Aging"/> or <see cref="SetPoolingName.Maturation"/>); null outside.</summary>
    public string? Family { get; init; }
    /// <summary><c>in_family</c>.</summary>
    public bool InFamily => Family is not null;
    /// <summary>
    /// A primary's Benjamini–Hochberg adjusted p within its family; for a pooled dataset of a <c>supported</c> part, the
    /// Benjamini–Hochberg adjusted p within the part's pooled datasets, compared at <see cref="BogomolovLevel"/>; else null.
    /// </summary>
    public double? PAdjusted { get; init; }
    /// <summary>A primary's family size; null otherwise.</summary>
    public int? FamilySize { get; init; }
    /// <summary><c>bb_level</c>: 0.05·R/m for a pooled dataset of a <c>supported</c> part; null otherwise.</summary>
    public double? BogomolovLevel { get; init; }
    /// <summary>The verdict (<see cref="SetVerdictName"/>).</summary>
    public required string Verdict { get; init; }
    /// <summary>The verdict's detail: the rule, guard or reason.</summary>
    public required string VerdictDetail { get; init; }
}

/// <summary>One part's pooled row (<c>DEF-DIFF-SET-POOL</c>): k ≥ 2 independent datasets pooled by random effects.</summary>
public sealed record SetPooledResult
{
    /// <summary>The part.</summary>
    public required SetHypothesisPart Part { get; init; }
    /// <summary>The datasets pooled, in the part's contrast order.</summary>
    public required IReadOnlyList<SetContrastOutcome> Datasets { get; init; }
    /// <summary>Why the pooled row is not run, or null.</summary>
    public string? NotRun { get; init; }
    /// <summary>
    /// The pooled estimate: DerSimonian–Laird's estimate with the Hartung–Knapp standard error on k − 1 degrees of freedom;
    /// null when not run.
    /// </summary>
    public SetEstimate? Estimate { get; init; }
    /// <summary>The prediction's test on <see cref="Estimate"/> (<see cref="SetVerdict.Test"/>); null when not run.</summary>
    public SetVerdict? Test { get; init; }
    /// <summary>τ², the DerSimonian–Laird between-dataset variance.</summary>
    public double Tau2 { get; init; }
    /// <summary>
    /// The guards: the pooled datasets' members whose state blocks <c>supported</c>, each written <c>dataset:member</c>,
    /// in contrast then member order. Ruled by the user 2026-11-14: any of them blocks the pooled row.
    /// </summary>
    public IReadOnlyList<string> BlockedBy { get; init; } = Array.Empty<string>();
    /// <summary>For up and down: <c>seen_one_side</c> members on the predicted side, summed over the pooled datasets.</summary>
    public int OneSideAgree { get; init; }
    /// <summary>For up and down: <c>seen_one_side</c> members on the other side, summed over the pooled datasets.</summary>
    public int OneSideDisagree { get; init; }
    /// <summary><c>n_datasets_agree</c>: pooled datasets whose own 90% interval lies on the predicted side of the bound.</summary>
    public int DatasetsAgree { get; init; }
    /// <summary>
    /// <c>p_supported_in_2</c>: the partial conjunction p that at least 2 of the pooled datasets support the prediction,
    /// Fisher's combination of the k − 1 largest dataset p's (Benjamini and Heller 2008).
    /// </summary>
    public double PSupportedInTwo { get; init; } = double.NaN;
    /// <summary>The aging family when tested; null when not run.</summary>
    public string? Family { get; init; }
    /// <summary>Its Benjamini–Hochberg adjusted p within the aging family.</summary>
    public double? PAdjusted { get; init; }
    /// <summary>The aging family's size.</summary>
    public int? FamilySize { get; init; }
    /// <summary>The verdict (<see cref="SetVerdictName"/>).</summary>
    public required string Verdict { get; init; }
    /// <summary>The verdict's detail.</summary>
    public required string VerdictDetail { get; init; }
}

/// <summary>One hypothesis part after pooling: its dataset rows and, with k ≥ 2 pooled datasets, its pooled row.</summary>
public sealed record SetPartResult
{
    /// <summary>The part.</summary>
    public required SetHypothesisPart Part { get; init; }
    /// <summary>Every dataset contrast supplied for the part, in the part's contrast order.</summary>
    public required IReadOnlyList<SetDatasetResult> Datasets { get; init; }
    /// <summary>The pooled row, when k ≥ 2 datasets pooled; null otherwise.</summary>
    public SetPooledResult? Pooled { get; init; }
}

/// <summary>Every part of one freeze on one quantification basis, after pooling and the two families.</summary>
public sealed record SetPoolingResult
{
    /// <summary>The freeze's tag.</summary>
    public required string Freeze { get; init; }
    /// <summary>The sha256 of the freeze's <c>tests.tsv</c>.</summary>
    public required string TestsSha256 { get; init; }
    /// <summary>The quantification basis every outcome was computed on.</summary>
    public required QuantBasis Basis { get; init; }
    /// <summary>The parts with at least one outcome, in the freeze's order.</summary>
    public required IReadOnlyList<SetPartResult> Parts { get; init; }
    /// <summary>m: primary results in the aging family.</summary>
    public int AgingFamilySize { get; init; }
    /// <summary>R: aging primaries <c>supported</c>.</summary>
    public int AgingSupported { get; init; }
    /// <summary>0.05·R/m, the level of the dataset rows of a <c>supported</c> part; NaN when none is supported.</summary>
    public double BogomolovLevel { get; init; } = double.NaN;
    /// <summary>Maturation contrasts in the maturation family.</summary>
    public int MaturationFamilySize { get; init; }
}

/// <summary>
/// Pools each hypothesis part's datasets and forms the families (QuantProject <c>DEF-DIFF-SET-POOL</c>, <c>-FAMILY</c>,
/// <c>-POSTHOC</c>, <c>-MATURATION</c>).
/// </summary>
/// <remarks>
/// <para>
/// <b>Which contrasts build a primary result.</b> A contrast is left out when it is not independent (a reproduction, or
/// <c>not_pooled</c>), post hoc, or not its dataset's primary: a part with one contrast of a dataset uses it unless it says
/// <c>primary = no</c>; a part with several uses the one marked <c>primary = yes</c>, and none of them unless exactly one is.
/// Maturation contrasts are read the same way among themselves and never with aging contrasts. A contrast that was not
/// tested (<see cref="SetContrastOutcome.Untested"/>) is reported and builds nothing.
/// </para>
/// <para>
/// <b>An aging part's primary result:</b> none with k = 0 tested contrasts left in; the one dataset's result with k = 1; and
/// with k ≥ 2 a pooled row: random effects with the DerSimonian–Laird τ² and the Hartung–Knapp standard error on k − 1
/// degrees of freedom (IntHout et al. 2014), tested by the same rules as a dataset (<see cref="SetVerdict.Test"/>).
/// </para>
/// <para>
/// <b>Ruled by the user 2026-11-14 (DATA-DEFINITIONS v3.26):</b> a pooled dataset's guard blocks the pooled row's
/// <c>supported</c>, naming <c>dataset:member</c>; an <c>attenuated_*</c> pooled row tests against the average of the
/// datasets' m̄ with the weights that average their β̄, 1/(seᵢ² + τ²). <b>Taken by QuantProject, open to challenge:</b>
/// <c>n_datasets_agree</c>, the "supported in at least 2" p and the Benjamini–Bogomolov rows cover the pooled datasets
/// only; R counts the aging primaries <c>supported</c>; a pooled row whose datasets all give one estimate is not run.
/// </para>
/// <para>
/// <b>The families:</b> one Benjamini–Hochberg family of aging primaries and one of maturation primaries, each at 0.05.
/// A pooled dataset of a <c>supported</c> part is adjusted by Benjamini–Bogomolov (2014): Benjamini–Hochberg within the
/// part's pooled datasets, compared at 0.05·R/m. Every other row has no adjusted p.
/// </para>
/// </remarks>
public static class SetPooling
{
    /// <summary>The "supported in at least" count for <see cref="SetPooledResult.PSupportedInTwo"/>.</summary>
    public const int SupportedInAtLeast = 2;

    /// <summary>Pools and adjusts every part of one freeze on one quantification basis.</summary>
    /// <param name="freeze">The freeze the parts came from.</param>
    /// <param name="basis">The quantification basis of every outcome; one call per basis, one pair of families each.</param>
    /// <param name="outcomes">The tested contrasts (<see cref="SetTest.Run"/>), at most one per part and contrast.</param>
    /// <exception cref="ArgumentException">An outcome whose part is not in <paramref name="freeze"/>, or whose contrast is
    /// not one of its part's; one supplied twice; or pooled datasets with different margins.</exception>
    public static SetPoolingResult Run(SetTestFreeze freeze, QuantBasis basis, IReadOnlyList<SetContrastOutcome> outcomes)
    {
        ArgumentNullException.ThrowIfNull(freeze);
        ArgumentNullException.ThrowIfNull(outcomes);
        if (!Enum.IsDefined(basis)) throw new ArgumentOutOfRangeException(nameof(basis));

        var byPart = new Dictionary<SetHypothesisPart, List<SetContrastOutcome>>(ReferenceEqualityComparer.Instance);
        foreach (var o in outcomes)
        {
            ArgumentNullException.ThrowIfNull(o, nameof(outcomes));
            var part = freeze.Parts.FirstOrDefault(p => ReferenceEquals(p, o.Part))
                ?? throw new ArgumentException($"{o.Part.HypothesisId} {SetTestFreeze.PredictionName(o.Part.Prediction)} is not a part of freeze {freeze.Freeze}.", nameof(outcomes));
            if (!part.Contrasts.Any(c => ReferenceEquals(c, o.Contrast)))
                throw new ArgumentException($"{o.Contrast.Dataset} is not a contrast of {part.HypothesisId} {SetTestFreeze.PredictionName(part.Prediction)}.", nameof(outcomes));
            if (!byPart.TryGetValue(part, out var list)) byPart[part] = list = new List<SetContrastOutcome>();
            if (list.Any(x => ReferenceEquals(x.Contrast, o.Contrast)))
                throw new ArgumentException($"{part.HypothesisId} {SetTestFreeze.PredictionName(part.Prediction)} on {o.Contrast.Dataset} is supplied twice.", nameof(outcomes));
            list.Add(o);
        }

        // 1. Each part: what is left out, what pools, the pooled row's test.
        var parts = new List<PartWork>();
        foreach (var part in freeze.Parts.Where(byPart.ContainsKey))
        {
            var supplied = part.Contrasts.SelectMany(c => byPart[part].Where(o => ReferenceEquals(o.Contrast, c))).ToList();
            var work = new PartWork(part, supplied.Select(o => new DatasetWork(o, LeftOut(part, o.Contrast))).ToList());
            var pool = work.Datasets.Where(d => !d.Outcome.Contrast.Maturation && d.LeftOut is null && d.Outcome.Test is not null).ToList();
            if (pool.Count == 1) pool[0].IsPrimary = true;
            if (pool.Count >= 2)
            {
                foreach (var d in pool) d.Pooled = true;
                work.Pooled = Pool(part, pool.Select(d => d.Outcome).ToList());
            }
            foreach (var d in work.Datasets.Where(d => d.Outcome.Contrast.Maturation && d.LeftOut is null && d.Outcome.Test is not null))
                d.IsPrimary = true;
            parts.Add(work);
        }

        // 2. The two families.
        var aging = new List<(PartWork Part, DatasetWork? Dataset, double P)>();
        foreach (var w in parts)
        {
            if (w.Pooled?.Test is { } test) aging.Add((w, null, test.PValue));
            foreach (var d in w.Datasets.Where(d => d.IsPrimary && !d.Outcome.Contrast.Maturation)) aging.Add((w, d, d.Outcome.Test!.PValue));
        }
        var agingAdjusted = MultipleTesting.BenjaminiHochberg(aging.Select(a => a.P).ToArray());
        for (int i = 0; i < aging.Count; i++)
        {
            if (aging[i].Dataset is { } d) { d.PAdjusted = agingAdjusted[i]; d.FamilySize = aging.Count; d.Family = SetPoolingName.Aging; }
            else { aging[i].Part.PooledPAdjusted = agingAdjusted[i]; }
        }
        var maturation = parts.SelectMany(w => w.Datasets).Where(d => d.IsPrimary && d.Outcome.Contrast.Maturation).ToList();
        var maturationAdjusted = MultipleTesting.BenjaminiHochberg(maturation.Select(d => d.Outcome.Test!.PValue).ToArray());
        for (int i = 0; i < maturation.Count; i++)
        {
            maturation[i].PAdjusted = maturationAdjusted[i];
            maturation[i].FamilySize = maturation.Count;
            maturation[i].Family = SetPoolingName.Maturation;
        }

        // 3. The primaries' verdicts, then Benjamini–Bogomolov on the pooled datasets of the supported parts.
        var pooledVerdict = new Dictionary<PartWork, (string Verdict, string Detail)>();
        foreach (var w in parts.Where(w => w.Pooled is not null)) pooledVerdict[w] = PooledVerdict(w.Pooled!, w.PooledPAdjusted);
        var supported = aging.Select(a => a.Dataset is { } d
            ? d.Outcome.Verdict(d.PAdjusted).Verdict == SetVerdictName.Supported
            : pooledVerdict[a.Part].Verdict == SetVerdictName.Supported).ToArray();
        var withinFamily = aging.Select(a => a.Dataset is null
            ? (IReadOnlyList<double>)a.Part.Datasets.Where(d => d.Pooled).Select(d => d.Outcome.Test!.PValue).ToArray()
            : new[] { a.P }).ToArray();
        var bogomolov = MultipleTesting.BenjaminiBogomolov(withinFamily, supported, SetVerdict.Alpha);
        for (int i = 0; i < aging.Count; i++)
        {
            if (aging[i].Dataset is not null || !supported[i]) continue;
            var pooled = aging[i].Part.Datasets.Where(d => d.Pooled).ToList();
            for (int j = 0; j < pooled.Count; j++)
            {
                pooled[j].PAdjusted = bogomolov.Adjusted[i][j];
                pooled[j].BogomolovLevel = bogomolov.Level;
            }
        }

        var results = parts.Select(w => new SetPartResult
        {
            Part = w.Part,
            Datasets = w.Datasets.Select(DatasetResult).ToList(),
            Pooled = w.Pooled is null ? null : w.Pooled with
            {
                Family = w.Pooled.Test is null ? null : SetPoolingName.Aging,
                PAdjusted = w.PooledPAdjusted,
                FamilySize = w.Pooled.Test is null ? null : aging.Count,
                Verdict = pooledVerdict[w].Verdict,
                VerdictDetail = pooledVerdict[w].Detail,
            },
        }).ToList();

        return new SetPoolingResult
        {
            Freeze = freeze.Freeze,
            TestsSha256 = freeze.Sha256,
            Basis = basis,
            Parts = results,
            AgingFamilySize = aging.Count,
            AgingSupported = bogomolov.Selected,
            BogomolovLevel = bogomolov.Level,
            MaturationFamilySize = maturation.Count,
        };
    }

    /// <summary>Why a contrast cannot build a primary result, or null. Reads the freeze's flags only, never a result.</summary>
    private static string? LeftOut(SetHypothesisPart part, SetDatasetContrast contrast)
    {
        if (contrast.Independence == SetIndependence.Reproduction) return SetPoolingName.Reproduction;
        if (contrast.Independence == SetIndependence.NotPooled) return SetPoolingName.NotPooled;
        if (contrast.PostHoc)
            return SetPoolingName.PostHocPrefix + string.Join(";", contrast.Rows.Where(r => r.PostHoc).Select(r => r.Member).Distinct(StringComparer.Ordinal));
        var sameDataset = part.Contrasts
            .Where(c => c.Maturation == contrast.Maturation && string.Equals(c.Dataset, contrast.Dataset, StringComparison.Ordinal))
            .ToList();
        if (sameDataset.Count == 1) return contrast.Primary == false ? SetPoolingName.NotPrimary : null;
        var marked = sameDataset.Where(c => c.Primary == true).ToList();
        if (marked.Count != 1) return SetPoolingName.SeveralContrastsPerDataset;
        return ReferenceEquals(marked[0], contrast) ? null : SetPoolingName.NotPrimary;
    }

    private static SetPooledResult Pool(SetHypothesisPart part, IReadOnlyList<SetContrastOutcome> datasets)
    {
        var tests = datasets.Select(d => d.Test!).ToList();
        string Of(SetContrastOutcome d) => d.Contrast.Dataset;
        if (tests.Select(t => t.Margin).Distinct().Count() > 1 || tests.Select(t => t.Attenuation).Distinct().Count() > 1)
            throw new ArgumentException(
                $"The pooled datasets of {part.HypothesisId} {SetTestFreeze.PredictionName(part.Prediction)} " +
                $"({string.Join(", ", datasets.Select(Of))}) test different margins; margins are frozen per hypothesis.", nameof(datasets));

        var meta = RandomEffectsMeta.Pool(tests.Select(t => t.Estimate.Log2Effect).ToArray(),
            tests.Select(t => t.Estimate.StandardError).ToArray(), 1 - 2 * SetVerdict.Alpha, RandomEffectsInference.HartungKnapp);
        var blocked = datasets.SelectMany(d => d.BlockedBy.Select(m => Of(d) + ":" + m)).ToList();
        var pooled = new SetPooledResult
        {
            Part = part,
            Datasets = datasets,
            Tau2 = meta.Tau2,
            BlockedBy = blocked,
            OneSideAgree = datasets.Sum(d => d.OneSideAgree),
            OneSideDisagree = datasets.Sum(d => d.OneSideDisagree),
            DatasetsAgree = tests.Count(t => t.PValue <= SetVerdict.Alpha),
            PSupportedInTwo = PValueCombination.FisherPartialConjunction(tests.Select(t => t.PValue).ToArray(), SupportedInAtLeast).PValue,
            Verdict = "",
            VerdictDetail = "",
        };
        if (!(meta.StandardError > 0)) return pooled with { NotRun = SetPoolingName.PooledStandardErrorZero };

        double? mrnaMean = null;
        if (part.Prediction is SetPrediction.AttenuatedDown or SetPrediction.AttenuatedUp)
        {
            // Ruled by the user 2026-11-14: the bound averages the datasets' m̄ with the weights that average their β̄.
            double[] w = tests.Select(t => 1 / (t.Estimate.StandardError * t.Estimate.StandardError + meta.Tau2)).ToArray();
            mrnaMean = w.Zip(tests, (wi, t) => wi * t.MrnaMean!.Value).Sum() / w.Sum();
        }
        var estimate = new SetEstimate
        {
            Log2Effect = meta.Estimate,
            StandardError = meta.StandardError,
            Df = meta.DegreesOfFreedom,
            DfMethod = SetPoolingName.PooledDfMethod,
            ModelUsed = SetPoolingName.PooledModel,
        };
        return pooled with
        {
            Estimate = estimate,
            Test = SetVerdict.Test(part.Prediction, estimate, tests[0].Margin, tests[0].Attenuation, mrnaMean),
        };
    }

    /// <summary>
    /// The pooled row's verdict, in <see cref="SetContrastOutcome.Verdict"/>'s order. A tested pooled row is always in the
    /// aging family, so it always has an adjusted p.
    /// </summary>
    private static (string Verdict, string Detail) PooledVerdict(SetPooledResult pooled, double? pAdjusted)
    {
        if (pooled.NotRun is not null) return (SetVerdictName.NotRun, pooled.NotRun);
        if (pooled.Test!.Contradicted) return (SetVerdictName.Contradicted, "uncorrected");
        if (pooled.BlockedBy.Count > 0) return (SetVerdictName.NotInformative, "blocked_by:" + string.Join(";", pooled.BlockedBy));
        if (pAdjusted!.Value <= SetVerdict.Alpha) return (SetVerdictName.Supported, "bh_adjusted_p<=0.05");
        return (SetVerdictName.NotInformative, "bh_adjusted_p>0.05");
    }

    private static SetDatasetResult DatasetResult(DatasetWork d)
    {
        var o = d.Outcome;
        (string verdict, string detail) = o.Verdict(d.Family is null ? null : d.PAdjusted);
        bool reachedTheFamilyStep = o.Untested is null && !o.Test!.Contradicted && o.BlockedBy.Count == 0;
        if (reachedTheFamilyStep && d.Family is null)
        {
            if (d.BogomolovLevel is { } level)
                (verdict, detail) = d.PAdjusted <= level
                    ? (SetVerdictName.Supported, SetPoolingName.BogomolovPassed)
                    : (SetVerdictName.NotInformative, SetPoolingName.BogomolovFailed);
            else if (d.LeftOut is not null)
                detail = d.LeftOut;
        }
        return new SetDatasetResult
        {
            Outcome = o,
            LeftOut = d.LeftOut,
            Pooled = d.Pooled,
            IsPrimary = d.IsPrimary,
            Family = d.Family,
            PAdjusted = d.PAdjusted,
            FamilySize = d.FamilySize,
            BogomolovLevel = d.BogomolovLevel,
            Verdict = verdict,
            VerdictDetail = detail,
        };
    }

    private sealed class PartWork
    {
        public PartWork(SetHypothesisPart part, List<DatasetWork> datasets) { Part = part; Datasets = datasets; }
        public SetHypothesisPart Part { get; }
        public List<DatasetWork> Datasets { get; }
        public SetPooledResult? Pooled { get; set; }
        public double? PooledPAdjusted { get; set; }
    }

    private sealed class DatasetWork
    {
        public DatasetWork(SetContrastOutcome outcome, string? leftOut) { Outcome = outcome; LeftOut = leftOut; }
        public SetContrastOutcome Outcome { get; }
        public string? LeftOut { get; }
        public bool Pooled { get; set; }
        public bool IsPrimary { get; set; }
        public string? Family { get; set; }
        public double? PAdjusted { get; set; }
        public int? FamilySize { get; set; }
        public double? BogomolovLevel { get; set; }
    }
}
