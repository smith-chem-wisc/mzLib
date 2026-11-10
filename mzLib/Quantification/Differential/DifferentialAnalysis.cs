using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using StatisticalModels;

namespace Quantification.Differential;

/// <summary>A comparison the analysis asks for in every stratum: between two levels of a factor, or a covariate's slope.</summary>
public sealed record DifferentialContrastRequest
{
    /// <summary>The factor, for a comparison between two of its levels.</summary>
    public string? Factor { get; init; }
    /// <summary>The numerator level.</summary>
    public string? Numerator { get; init; }
    /// <summary>The reference level.</summary>
    public string? Denominator { get; init; }
    /// <summary>The covariate, for a slope.</summary>
    public string? Covariate { get; init; }

    /// <summary><paramref name="numerator"/> against <paramref name="denominator"/> for one factor.</summary>
    public static DifferentialContrastRequest Between(string factor, string numerator, string denominator) =>
        new() { Factor = factor, Numerator = numerator, Denominator = denominator };

    /// <summary>The slope of a covariate.</summary>
    public static DifferentialContrastRequest Slope(string covariate) => new() { Covariate = covariate };
}

/// <summary>
/// Everything one differential analysis needs (STAT1 M5b): exactly one observation table, BEFORE normalization (the global
/// shift is measured on raw values), the analysis samples, the design, the comparisons, the settings and the inputs that,
/// with the settings, determine the analysis id.
/// </summary>
public sealed record DifferentialRequest
{
    /// <summary>A label-free table (peptide x sample); null when <see cref="IsobaricTable"/> is given.</summary>
    public ObservationTable? Table { get; init; }
    /// <summary>An isobaric table (ion in its plex run x channel sample); null when <see cref="Table"/> is given.</summary>
    public IsobaricObservationTable? IsobaricTable { get; init; }
    /// <summary>The analysis samples, each with its factors, strata, individual, batch and covariates.</summary>
    public required IReadOnlyList<AnalysisSample> Samples { get; init; }
    /// <summary>The design: factors, covariates, interactions and whether batch is fitted.</summary>
    public required DesignSpec Spec { get; init; }
    /// <summary>The comparisons, in order (ids <c>c1</c>, <c>c2</c>, ...); null for each stratum's default contrasts.</summary>
    public IReadOnlyList<DifferentialContrastRequest>? Contrasts { get; init; }
    /// <summary>The settings.</summary>
    public DifferentialSettings Settings { get; init; } = new();
    /// <summary>The inputs by role, file name and sha256.</summary>
    public IReadOnlyList<DifferentialInput> Inputs { get; init; } = Array.Empty<DifferentialInput>();
    /// <summary><c>lfq</c>, <c>tmt</c>, ...; null takes <c>lfq</c> for a label-free table and <c>tmt</c> for an isobaric one.</summary>
    public string? QuantStyle { get; init; }
    /// <summary><c>ms2</c> or <c>ms3_sps</c> for an isobaric table; null otherwise.</summary>
    public string? ReporterAcquisition { get; init; }
    /// <summary>Which source chose the strata's factors (written on each stratum).</summary>
    public string StratumSource { get; init; } = "default_list";
}

/// <summary>
/// One place a rule changed what happened (R-DECISIONS): its rule, its keys (stratum, quant basis, feature, contrast; a
/// null key means the whole level above), the facts that triggered it and what it did.
/// </summary>
public sealed record DifferentialDecision(string Rule, string Stratum, string? QuantBasis, string? FeatureId,
    string? ContrastId, string Facts, string Effect);

/// <summary>The decision rule ids (R-DECISIONS), one per STAT1 M5b ruling that can change what happens.</summary>
public static class DifferentialDecisionRules
{
    /// <summary>Groups left out of the row set, by reason (ruling 12).</summary>
    public const string RowSet = "row_set";
    /// <summary>A contrast the stratum's design cannot estimate, or cannot run (ruling 2 step 1).</summary>
    public const string DesignContrast = "design_contrast";
    /// <summary>No usable unique peptide (ruling 2 step 2, ruling 6).</summary>
    public const string NoUniquePeptides = "no_unique_peptides";
    /// <summary>A side with no value (ruling 1, ruling 2 step 3).</summary>
    public const string AbsentSide = "absent_side";
    /// <summary>Fewer than 2 samples with a value on a side (ruling 2 step 4).</summary>
    public const string MinTwoPerSide = "min_2_per_side";
    /// <summary>The feature was fitted on a reduced design (ruling 1).</summary>
    public const string ReducedDesign = "reduced_design";
    /// <summary>The feature's own data cannot separate a contrast (ruling 2 step 5).</summary>
    public const string FeatureNotEstimable = "feature_not_estimable";
    /// <summary>One usable peptide: the moderated t-test (ruling 6).</summary>
    public const string SinglePeptide = "single_peptide";
    /// <summary>The fit's own outcome when it is not a fit (ruling 2 step 6).</summary>
    public const string FitOutcome = "fit_outcome";
    /// <summary>Robust reweighting fell back to one round (ruling 8).</summary>
    public const string RobustFellBack = "robust_fell_back";
    /// <summary>Fewer than 3 features fitted in a moderation: no prior (ruling 5).</summary>
    public const string NoPrior = "no_prior";

    private static readonly HashSet<string> All = new(StringComparer.Ordinal)
    {
        RowSet, DesignContrast, NoUniquePeptides, AbsentSide, MinTwoPerSide, ReducedDesign, FeatureNotEstimable,
        SinglePeptide, FitOutcome, RobustFellBack, NoPrior,
    };

    /// <summary>Whether <paramref name="rule"/> is one of these ids.</summary>
    public static bool IsValid(string rule) => All.Contains(rule);
}

/// <summary>The fixed strings the runner writes (DATA-DEFINITIONS <c>DEF-DIFF-*</c>; STAT1 M5b ruling 10).</summary>
public static class DifferentialConstants
{
    /// <summary>The definition every abundance row follows.</summary>
    public const string DefinitionId = "QuantProject:DEF-DIFF-ABUNDANCE v1";
    /// <summary>The framework's version string.</summary>
    public const string MethodVersion = "diff-1.0";
    /// <summary>A comparison between two levels.</summary>
    public const string EffectRatio = "abundance_log2_ratio";
    /// <summary>A covariate's slope.</summary>
    public const string EffectSlope = "abundance_log2_slope";
    /// <summary>The quantity.</summary>
    public const string Quantity = "abundance";
    /// <summary>The grain.</summary>
    public const string Grain = "protein_group";
    /// <summary>The method.</summary>
    public const string Method = "moderated";
    /// <summary>The peptide-level mixed model.</summary>
    public const string MixedModel = "peptide_mixed_model";
    /// <summary>The moderated t-test, for one usable peptide.</summary>
    public const string ModeratedT = "moderated_t";
    /// <summary>The mixed model's df.</summary>
    public const string DfSatterthwaite = "satterthwaite_moderated";
    /// <summary>The t-test's df.</summary>
    public const string DfResidual = "moderated_residual";
    /// <summary>The robust rule (DATA-DEFINITIONS v3.21).</summary>
    public const string RobustRule = "huber 1.345, MAD about 0, to convergence (100 rounds), else one round";
    /// <summary>The adjustment.</summary>
    public const string Adjustment = "benjamini_hochberg";
    /// <summary>Support rule: no usable unique peptide.</summary>
    public const string NoUniquePeptides = "no_unique_peptides";
    /// <summary>Support rule: fewer than 2 samples with a value on a side.</summary>
    public const string MinTwoPerSide = "min_2_per_side";
    /// <summary>Support rule: too little data for the variance (ruling 2 step 6).</summary>
    public const string TooFewObservations = "too_few_observations";
}

/// <summary>A differential analysis's rows, metadata and decisions.</summary>
public sealed class DifferentialAnalysisResult
{
    internal DifferentialAnalysisResult(IReadOnlyList<DifferentialResult> rows, DifferentialMetadata metadata,
        IReadOnlyList<DifferentialDecision> decisions)
    {
        Rows = rows;
        Metadata = metadata;
        Decisions = decisions;
    }

    /// <summary>The result rows, in DEF-DIFF-ROW order.</summary>
    public IReadOnlyList<DifferentialResult> Rows { get; }
    /// <summary>Everything constant across the rows.</summary>
    public DifferentialMetadata Metadata { get; }
    /// <summary>Every place a rule changed what happened (R-DECISIONS).</summary>
    public IReadOnlyList<DifferentialDecision> Decisions { get; }
}

/// <summary>
/// The STAT1 runner (milestone M5b): turns an observation table, its samples, a design and the comparisons into
/// <see cref="DifferentialResult"/> rows, their metadata and a decisions record, applying the 15 rulings of
/// QuantProject <c>design/stat1/M5B_CONTRACTS_2026-11-04.md</c>.
/// </summary>
/// <remarks>
/// <para>
/// <b>Per stratum</b> (<see cref="AnalysisDesign.ByStratum"/>): the stratum's table is normalized (label-free:
/// <see cref="SharedPeptideNormalization"/>; isobaric: <see cref="IsobaricNormalization"/>, once over the table); each
/// protein group in the row set gets one row per comparison and quant basis.
/// </para>
/// <para>
/// <b>Per protein group</b>, one status, first match wins (ruling 2): the stratum's design cannot estimate the comparison;
/// no usable unique peptide; a side with no value; fewer than 2 samples with a value on a side; the group's own data
/// cannot separate the comparison (its design reduced to what its observations estimate, ruling 1); the fit's outcome.
/// Groups with two or more usable peptides (label-free) or ions in plex runs (isobaric) are fitted by
/// <see cref="NestedMixedModel"/>; one by <see cref="LinearModel"/> and the moderated t-test (ruling 6). Each model's
/// fits share one variance prior over every reduced design (<c>ModerateContrasts</c> over parts; ruling 5: fewer than 3
/// fitted, no prior). Benjamini–Hochberg runs once per family over both models' rows (ruling 4).
/// </para>
/// </remarks>
public static class DifferentialAnalysis
{
    /// <summary>The usable-peptide minimum (ruling 6): a unique level with a value in at least this many samples.</summary>
    public const int UsableMinimumSamples = 2;

    /// <summary>Runs the analysis.</summary>
    /// <param name="request">Exactly one table, the samples, the design, the comparisons and the settings.</param>
    /// <param name="maxThreads">As <see cref="NestedMixedModel.Fit"/>; the output does not depend on it.</param>
    /// <exception cref="ArgumentException">A request that cannot be read (see the messages).</exception>
    public static DifferentialAnalysisResult Run(DifferentialRequest request, int maxThreads = -1)
    {
        ArgumentNullException.ThrowIfNull(request);
        ArgumentNullException.ThrowIfNull(request.Samples, nameof(request));
        ArgumentNullException.ThrowIfNull(request.Spec, nameof(request));
        if ((request.Table is null) == (request.IsobaricTable is null))
            throw new ArgumentException("Give exactly one observation table: label-free or isobaric.", nameof(request));
        var withIndividual = request.Samples.Where(s => s.Individual is not null).Select(s => s.SampleId).ToList();
        if (withIndividual.Count > 0 && withIndividual.Count < request.Samples.Count)
            throw new ArgumentException(
                $"Sample '{request.Samples.First(s => s.Individual is null).SampleId}' has no individual, but others do; an " +
                "individual is needed on every sample or none.", nameof(request));

        return new Analysis(request, maxThreads).Execute();
    }

    /// <summary>One analysis in progress.</summary>
    private sealed class Analysis
    {
        private readonly DifferentialRequest _request;
        private readonly int _maxThreads;
        private readonly bool _isobaric;
        private readonly string _analysisId;
        private readonly string? _basis;
        private readonly string _style;
        private readonly List<DifferentialResult> _rows = new();
        private readonly List<DifferentialDecision> _decisions = new();
        private readonly List<DifferentialStratumInfo> _strata = new();
        private readonly List<DifferentialModelInfo> _models = new();
        private readonly List<DifferentialNormalizationInfo> _normalization = new();
        private readonly Dictionary<string, DifferentialContrastInfo> _contrastInfo = new(StringComparer.Ordinal);
        private readonly List<string> _designWarnings = new();
        private IsobaricObservationTable? _tmt;
        private IsobaricNormalizationResult? _tmtNormalization;

        public Analysis(DifferentialRequest request, int maxThreads)
        {
            _request = request;
            _maxThreads = maxThreads;
            _isobaric = request.IsobaricTable is not null;
            _analysisId = DifferentialMetadata.ComputeAnalysisId(request.Inputs, request.Settings);
            _basis = _isobaric ? null : request.Table!.Basis.ToMachineName();
            _style = request.QuantStyle ?? (_isobaric ? "tmt" : "lfq");
        }

        public DifferentialAnalysisResult Execute()
        {
            if (_isobaric)
            {
                var table = _request.IsobaricTable!;
                if (_request.Settings.Normalization == SharedPeptideNormalization.None)
                    _tmt = table;
                else
                {
                    _tmtNormalization = IsobaricNormalization.Measure(table);
                    _tmt = IsobaricNormalization.Apply(table, _tmtNormalization);
                }
                _designWarnings.AddRange(table.Warnings);
            }
            else _designWarnings.AddRange(_request.Table!.Warnings);

            foreach (var design in AnalysisDesign.ByStratum(_request.Samples, _request.Spec))
                Stratum(design);

            var rows = Adjust(_rows);
            var metadata = new DifferentialMetadata
            {
                Software = new[] { new DifferentialSoftware("mzLib", DifferentialConstants.MethodVersion) },
                Inputs = _request.Inputs,
                Settings = _request.Settings,
                Strata = _strata,
                Contrasts = _contrastInfo.Values.OrderBy(c => c.Id, StringComparer.Ordinal).ToList(),
                Models = _models,
                Normalization = _normalization,
                DesignWarnings = _designWarnings.Distinct(StringComparer.Ordinal).ToList(),
            };
            return new DifferentialAnalysisResult(DifferentialResultWriter.Order(rows), metadata, _decisions);
        }

        // ---- one stratum ----

        private void Stratum(AnalysisDesign design)
        {
            string stratum = design.Stratum;
            var sampleMap = design.Samples.ToDictionary(s => s.SampleId, StringComparer.Ordinal);
            var contrasts = Contrasts(design);
            _designWarnings.AddRange(design.Warnings);

            IReadOnlyDictionary<string, FeatureObservations> features;
            IReadOnlyDictionary<string, ProteinGroupInfo> groupInfo;
            List<string> candidates;
            string normalization;
            if (_isobaric)
            {
                features = FeatureObservationBuilder.FromIsobaric(_tmt!, sampleMap, UsableMinimumSamples);
                groupInfo = _tmt!.ProteinGroups;
                candidates = IsobaricCandidates(_tmt!, sampleMap);
                normalization = _tmtNormalization?.Setting ?? SharedPeptideNormalization.None;
                if (_tmtNormalization is not null) _normalization.Add(IsobaricNormalizationInfo(stratum, sampleMap));
            }
            else
            {
                var table = _request.Table!.ForSamples(sampleMap.Keys);
                var normalized = table;
                normalization = SharedPeptideNormalization.None;
                if (_request.Settings.Normalization != SharedPeptideNormalization.None)
                {
                    var measured = SharedPeptideNormalization.Measure(table);
                    normalized = SharedPeptideNormalization.Apply(table, measured);
                    normalization = measured.Setting;
                    _normalization.Add(new DifferentialNormalizationInfo(stratum, _basis, measured.Setting,
                        measured.ReferenceSetSize, measured.MinimumReferenceSetSize, measured.PerSampleShift, measured.Warnings));
                }
                features = FeatureObservationBuilder.FromTable(normalized, sampleMap, UsableMinimumSamples);
                groupInfo = normalized.ProteinGroups;
                candidates = LabelFreeCandidates(normalized);
            }

            // The row set (ruling 12).
            var included = new List<string>();
            var left = new SortedDictionary<string, int>(StringComparer.Ordinal);
            foreach (string group in candidates)
            {
                string? why = FeatureObservationBuilder.RowSetExclusion(group, groupInfo);
                if (why is null) included.Add(group);
                else left[why] = left.GetValueOrDefault(why) + 1;
            }
            if (left.Count > 0)
                _decisions.Add(new DifferentialDecision(DifferentialDecisionRules.RowSet, stratum, _basis, null, null,
                    string.Join("; ", left.Select(kv => $"{kv.Key} {kv.Value}")), "left out of the row set: no rows"));

            var tested = TestedColumns(design, contrasts.Where(c => c.Estimability == AnalysisDesign.Estimable).Select(c => c.Contrast));
            var plans = new List<FeaturePlan>();
            foreach (string group in included)
                plans.Add(Plan(design, contrasts, group, features.GetValueOrDefault(group), groupInfo[group], tested));

            FitAndModerate(design, contrasts, plans, tested, normalization);
            _strata.Add(StratumInfo(design));
        }

        private List<StratumContrast> Contrasts(AnalysisDesign design)
        {
            var list = new List<StratumContrast>();
            if (_request.Contrasts is null)
            {
                foreach (var c in design.DefaultContrasts()) list.Add(new StratumContrast(c, design.Estimability(c)));
            }
            else
            {
                for (int i = 0; i < _request.Contrasts.Count; i++)
                {
                    var r = _request.Contrasts[i];
                    string id = string.Create(CultureInfo.InvariantCulture, $"c{i + 1}");
                    try
                    {
                        var c = (r.Covariate is not null ? design.Slope(r.Covariate) : design.Between(r.Factor!, r.Numerator!, r.Denominator!))
                            with { Id = id };
                        list.Add(new StratumContrast(c, design.Estimability(c)));
                    }
                    catch (ArgumentException e)
                    {
                        _decisions.Add(new DifferentialDecision(DifferentialDecisionRules.DesignContrast, design.Stratum, _basis,
                            null, id, e.Message, "not run in this stratum: no rows"));
                    }
                }
            }
            foreach (var c in list)
            {
                if (c.Estimability != AnalysisDesign.Estimable)
                    _decisions.Add(new DifferentialDecision(DifferentialDecisionRules.DesignContrast, design.Stratum, _basis,
                        null, c.Contrast.Id, $"the stratum's design: {c.Estimability}", $"every row of the contrast is {c.Estimability}"));
                if (!_contrastInfo.ContainsKey(c.Contrast.Id))
                    _contrastInfo[c.Contrast.Id] = new DifferentialContrastInfo(c.Contrast.Id, c.Contrast.Label, c.Contrast.Numerator,
                        c.Contrast.Denominator, c.Contrast.Covariate, c.Contrast.CovariateUnit, c.Contrast.CovariateScale,
                        Enumerable.Range(0, design.ColumnNames.Count).Where(j => c.Contrast.Weights[j] != 0)
                            .ToDictionary(j => design.ColumnNames[j], j => c.Contrast.Weights[j], StringComparer.Ordinal));
            }
            return list;
        }

        /// <summary>Every column of every term a contrast puts weight on (ruling 7).</summary>
        private static int[] TestedColumns(AnalysisDesign design, IEnumerable<Contrast> contrasts)
        {
            var terms = new HashSet<string>(StringComparer.Ordinal);
            foreach (var c in contrasts)
                for (int j = 0; j < c.Weights.Count; j++)
                    if (c.Weights[j] != 0) terms.Add(design.ColumnTerms[j]);
            return Enumerable.Range(0, design.ColumnNames.Count).Where(j => terms.Contains(design.ColumnTerms[j])).ToArray();
        }

        private static List<string> LabelFreeCandidates(ObservationTable table)
        {
            var groups = new SortedSet<string>(StringComparer.Ordinal);
            for (int p = 0; p < table.Peptides.Count; p++)
                if (Enumerable.Range(0, table.Samples.Count).Any(s => !double.IsNaN(table.Log2Intensity(p, s))))
                    foreach (string g in table.Peptides[p].ProteinGroups)
                        if (g != ObservedPeptide.UndefinedProteinGroup) groups.Add(g);
            return groups.ToList();
        }

        private static List<string> IsobaricCandidates(IsobaricObservationTable table, IReadOnlyDictionary<string, AnalysisSample> samples)
        {
            var groups = new SortedSet<string>(StringComparer.Ordinal);
            for (int r = 0; r < table.Rows.Count; r++)
                if (table.SamplesOfPlex(table.Rows[r].Plex).Any(s => samples.ContainsKey(table.Samples[s]) && !double.IsNaN(table.Log2Intensity(r, s))))
                    foreach (string g in table.Rows[r].Peptide.ProteinGroups)
                        if (g != ObservedPeptide.UndefinedProteinGroup) groups.Add(g);
            return groups.ToList();
        }

        // ---- one protein group: the status of each contrast before any fit ----

        private FeaturePlan Plan(AnalysisDesign design, List<StratumContrast> contrasts, string group,
            FeatureObservations? values, ProteinGroupInfo info, int[] tested)
        {
            var plan = new FeaturePlan(group, info, values);
            bool usable = values is not null && values.NuisanceLevels.Count > 0;
            FeatureDesign? featureDesign = null;
            if (usable)
            {
                var (common, nuisance) = Matrices(design, values!);
                featureDesign = FeatureDesign.Build(common, nuisance, design.ColumnTerms);
                plan.Design = featureDesign;
                plan.CommonRows = common;
                if (!featureDesign.IsFull)
                    Decide(DifferentialDecisionRules.ReducedDesign, design.Stratum, group, null,
                        $"kept {featureDesign.KeptCommon.Count} of {featureDesign.CommonCount} design columns and " +
                        $"{featureDesign.KeptNuisance.Count} of {featureDesign.NuisanceCount} nuisance columns",
                        "fitted on the columns its observations estimate");
            }

            foreach (var sc in contrasts)
            {
                var c = sc.Contrast;
                var (numerator, denominator) = Sides(design, c);
                var counts = new SideCounts(
                    numerator is null ? null : values?.SamplesWithValue(numerator) ?? 0,
                    denominator is null ? null : values?.SamplesWithValue(denominator) ?? 0,
                    values?.SamplesWithValue(Union(design, numerator, denominator)) ?? 0,
                    Union(design, numerator, denominator).Count);
                string? status = null, detail = null;
                if (sc.Estimability != AnalysisDesign.Estimable)
                {
                    status = sc.Estimability;
                    detail = "the stratum's design cannot estimate this contrast";
                }
                else if (!usable)
                {
                    status = DifferentialStatus.BelowSupport(DifferentialConstants.NoUniquePeptides);
                    detail = values is null || values.UniqueLevels == 0
                        ? "no peptide unique to the group"
                        : $"{values.UniqueLevels} peptide(s) unique to the group, none valued in {UsableMinimumSamples} or more of the stratum's samples";
                    Decide(DifferentialDecisionRules.NoUniquePeptides, design.Stratum, group, c.Id, detail, "no fit");
                }
                else if (counts.Numerator is int n && counts.Denominator is int d && (n == 0 || d == 0))
                {
                    status = n == 0 && d == 0 ? DifferentialStatus.AbsentInBoth
                        : n == 0 ? DifferentialStatus.AbsentInNumerator : DifferentialStatus.AbsentInDenominator;
                    detail = $"samples with a value: numerator {n}, denominator {d}";
                    Decide(DifferentialDecisionRules.AbsentSide, design.Stratum, group, c.Id, detail, $"reported {status}, no number");
                }
                else if (counts.Numerator is int n2 && counts.Denominator is int d2 && (n2 < 2 || d2 < 2))
                {
                    status = DifferentialStatus.BelowSupport(DifferentialConstants.MinTwoPerSide);
                    detail = $"samples with a value: numerator {n2}, denominator {d2}";
                    Decide(DifferentialDecisionRules.MinTwoPerSide, design.Stratum, group, c.Id, detail, "no number");
                }
                else
                {
                    string estimability = featureDesign!.Estimability(c.Weights);
                    if (estimability != AnalysisDesign.Estimable)
                    {
                        status = estimability;
                        detail = "the group's own observations cannot separate this contrast";
                        Decide(DifferentialDecisionRules.FeatureNotEstimable, design.Stratum, group, c.Id, detail, $"reported {status}");
                    }
                }
                plan.Contrasts.Add(new PlannedContrast(c, counts, status, detail));
            }
            return plan;
        }

        /// <summary>Each value's stratum design row, and the nuisance columns treatment-coded on the first level (ordinal).</summary>
        private static (double[,] Common, double[,]? Nuisance) Matrices(AnalysisDesign design, FeatureObservations values)
        {
            var matrix = design.Matrix;
            var rowOf = design.Samples.Select((s, i) => (s.SampleId, i)).ToDictionary(t => t.SampleId, t => t.i, StringComparer.Ordinal);
            int n = values.Values.Count, p = design.ColumnNames.Count, q = values.NuisanceLevels.Count - 1;
            var level = values.NuisanceLevels.Select((l, i) => (l, i)).ToDictionary(t => t.l, t => t.i, StringComparer.Ordinal);
            var common = new double[n, p];
            var nuisance = q > 0 ? new double[n, q] : null;
            for (int i = 0; i < n; i++)
            {
                var v = values.Values[i];
                int row = rowOf[v.SampleId];
                for (int j = 0; j < p; j++) common[i, j] = matrix[row, j];
                int k = level[v.Nuisance];
                if (k > 0) nuisance![i, k - 1] = 1;
            }
            return (common, nuisance);
        }

        /// <summary>The samples on each side of a between-levels contrast; nulls for a slope.</summary>
        private static (List<string>? Numerator, List<string>? Denominator) Sides(AnalysisDesign design, Contrast c)
        {
            if (c.Numerator is null || c.Denominator is null) return (null, null);
            List<string> Of(string side)
            {
                int eq = side.IndexOf('=');
                string factor = side[..eq], level = side[(eq + 1)..];
                return design.Samples.Where(s => s.Factors.TryGetValue(factor, out var v) && v == level).Select(s => s.SampleId).ToList();
            }
            return (Of(c.Numerator), Of(c.Denominator));
        }

        private static List<string> Union(AnalysisDesign design, List<string>? numerator, List<string>? denominator) =>
            numerator is null || denominator is null
                ? design.Samples.Select(s => s.SampleId).ToList()
                : numerator.Concat(denominator).Distinct(StringComparer.Ordinal).ToList();

        private void Decide(string rule, string stratum, string? feature, string? contrast, string facts, string effect) =>
            _decisions.Add(new DifferentialDecision(rule, stratum, _basis, feature, contrast, facts, effect));

        // ---- fits, moderation and rows ----

        private void FitAndModerate(AnalysisDesign design, List<StratumContrast> contrasts, List<FeaturePlan> plans,
            int[] tested, string normalization)
        {
            var toFit = plans.Where(p => p.Design is not null && p.Contrasts.Any(c => c.Status is null)).ToList();
            var mixed = toFit.Where(p => p.Values!.NuisanceLevels.Count > 1).GroupBy(p => p.Design!.CommonKey, StringComparer.Ordinal)
                .OrderBy(g => g.Key, StringComparer.Ordinal).ToList();
            var single = toFit.Where(p => p.Values!.NuisanceLevels.Count == 1).GroupBy(p => p.Design!.CommonKey, StringComparer.Ordinal)
                .OrderBy(g => g.Key, StringComparer.Ordinal).ToList();
            var estimable = contrasts.Where(c => c.Estimability == AnalysisDesign.Estimable).Select(c => c.Contrast).ToList();

            // Mixed model per kept-column group, one pooled moderation.
            var mixedParts = new List<MixedModelModerationPart>();
            var mixedMembers = new List<List<FeaturePlan>>();
            foreach (var g in mixed)
            {
                var members = g.ToList();
                var first = members[0].Design!;
                var names = first.KeptCommon.Select(j => design.ColumnNames[j]).ToArray();
                var problems = members.Select(p => Problem(p)).ToList();
                var fit = NestedMixedModel.Fit(problems, names, _maxThreads, _request.Settings.Robust);
                var weights = PartContrasts(first, estimable);
                var rows = first.VarianceDfRows(tested);
                for (int i = 0; i < members.Count; i++)
                {
                    members[i].MixedFit = fit;
                    members[i].Index = i;
                }
                if (rows.Count == 0)
                {
                    foreach (var p in members) p.NoVarianceDf = true;
                    continue;
                }
                mixedParts.Add(new MixedModelModerationPart(fit, weights, rows.Select(r => (IReadOnlyList<double>)r).ToList()));
                mixedMembers.Add(members);
            }
            if (mixedParts.Count > 0)
            {
                var moderated = NestedMixedModel.ModerateContrasts(mixedParts, _request.Settings.CiLevel);
                for (int part = 0; part < moderated.Count; part++)
                    foreach (var p in mixedMembers[part]) p.Moderated = moderated[part];
                Prior(design, moderated, DifferentialConstants.MixedModel, mixedParts.Sum(m => m.Fit.Status.Count(s => s == FeatureFitStatus.Fitted)));
            }

            // One usable level: the moderated t-test, its own pooled prior.
            var linearParts = new List<LinearModelModerationPart>();
            var linearMembers = new List<List<FeaturePlan>>();
            var sampleIds = design.Samples.Select(s => s.SampleId).ToList();
            var matrix = design.Matrix;
            foreach (var g in single)
            {
                var members = g.ToList();
                var first = members[0].Design!;
                var kept = first.KeptCommon;
                var x = new double[sampleIds.Count, kept.Count];
                for (int s = 0; s < sampleIds.Count; s++)
                    for (int k = 0; k < kept.Count; k++) x[s, k] = matrix[s, kept[k]];
                var y = new double[members.Count, sampleIds.Count];
                for (int f = 0; f < members.Count; f++)
                {
                    for (int s = 0; s < sampleIds.Count; s++) y[f, s] = double.NaN;
                    foreach (var v in members[f].Values!.Values) y[f, sampleIds.IndexOf(v.SampleId)] = v.Log2;
                    members[f].Index = f;
                }
                var fit = LinearModel.Fit(y, x, kept.Select(j => design.ColumnNames[j]).ToArray());
                foreach (var p in members)
                {
                    p.LinearFit = fit;
                    Decide(DifferentialDecisionRules.SinglePeptide, design.Stratum, p.Group, null,
                        "one usable peptide", "tested with the moderated t-test (fitted_single_peptide)");
                }
                linearParts.Add(new LinearModelModerationPart(fit, PartContrasts(first, estimable)));
                linearMembers.Add(members);
            }
            if (linearParts.Count > 0)
            {
                var moderated = EmpiricalBayes.ModerateContrasts(linearParts, _request.Settings.CiLevel);
                for (int part = 0; part < moderated.Count; part++)
                    foreach (var p in linearMembers[part]) p.Moderated = moderated[part];
                Prior(design, moderated, DifferentialConstants.ModeratedT, linearParts.Sum(l => l.Fit.Status.Count(s => s == FeatureFitStatus.Fitted)));
            }

            foreach (var plan in plans) Rows(design, plan, normalization);
        }

        private static NestedMixedModelProblem Problem(FeaturePlan plan)
        {
            var values = plan.Values!.Values;
            var d = plan.Design!;
            var (common, nuisance) = (new double[values.Count, d.KeptCommon.Count], d.KeptNuisance.Count > 0 ? new double[values.Count, d.KeptNuisance.Count] : null);
            var level = plan.Values.NuisanceLevels.Select((l, i) => (l, i)).ToDictionary(t => t.l, t => t.i, StringComparer.Ordinal);
            for (int i = 0; i < values.Count; i++)
            {
                for (int k = 0; k < d.KeptCommon.Count; k++) common[i, k] = plan.CommonRows![i, d.KeptCommon[k]];
                int lv = level[values[i].Nuisance];
                for (int k = 0; k < d.KeptNuisance.Count; k++) nuisance![i, k] = d.KeptNuisance[k] == lv - 1 ? 1 : 0;
            }
            bool individuals = values.Any(v => v.RandomIndividual is not null);
            return new NestedMixedModelProblem(values.Select(v => v.Log2).ToArray(), common, nuisance,
                values.Select(v => v.RandomSample).ToArray(), individuals ? values.Select(v => v.RandomIndividual).ToArray() : null);
        }

        /// <summary>Each estimable contrast's weights on the part's kept columns, unless they are all zero there.</summary>
        private static List<ContrastWeights> PartContrasts(FeatureDesign design, List<Contrast> contrasts) =>
            contrasts.Select(c => new ContrastWeights(c.Id, design.KeptWeights(c.Weights)))
                .Where(w => w.Weights.Any(x => x != 0)).ToList();

        private void Prior(AnalysisDesign design, IReadOnlyList<IReadOnlyList<ModeratedContrast>> moderated, string model, int fitted)
        {
            var any = moderated.SelectMany(m => m).FirstOrDefault();
            DifferentialPriorInfo? prior = null;
            if (any is not null)
            {
                var finite = any.Prior.Scale.Where(double.IsFinite).ToList();
                double? variance = finite.Count > 0 ? finite[0] : null;
                prior = new DifferentialPriorInfo(DifferentialPriorInfo.MomentsLegacy, any.Prior.Df, false,
                    any.Prior.Df == 0 ? null : variance, fitted);
                if (any.Prior.Df == 0)
                    Decide(DifferentialDecisionRules.NoPrior, design.Stratum, null, null, $"{fitted} feature(s) fitted ({model})",
                        "no variance prior: each keeps its own variance and df");
            }
            bool mixedModel = model == DifferentialConstants.MixedModel;
            _models.Add(new DifferentialModelInfo(design.Stratum, _basis, DifferentialConstants.Method, model, DifferentialConstants.Grain,
                Formula(design, mixedModel), true, mixedModel ? DifferentialConstants.DfSatterthwaite : DifferentialConstants.DfResidual,
                mixedModel && _request.Settings.Robust ? DifferentialConstants.RobustRule : null, prior));
        }

        private string Formula(AnalysisDesign design, bool mixedModel)
        {
            var terms = design.ColumnTerms.Where(t => t != "intercept").Distinct(StringComparer.Ordinal).ToList();
            string nuisance = _isobaric ? "ion_in_run" : "peptide";
            string random = _isobaric
                ? (_request.Samples.Any(s => s.Individual is not null) ? "(1|individual/sample)" : "(1|sample/run)")
                : (_request.Samples.Any(s => s.Individual is not null) ? "(1|individual/sample)" : "(1|sample)");
            return mixedModel
                ? $"log2 y ~ {string.Join(" + ", new[] { nuisance }.Concat(terms))} + {random}"
                : $"log2 y ~ {string.Join(" + ", terms.DefaultIfEmpty("1"))}";
        }

        private void Rows(AnalysisDesign design, FeaturePlan plan, string normalization)
        {
            foreach (var pc in plan.Contrasts)
            {
                var c = pc.Contrast;
                string? status = pc.Status, detail = pc.Detail;
                ModeratedContrast? result = null;
                bool single = plan.Values?.NuisanceLevels.Count == 1;
                if (status is null)
                {
                    result = plan.Moderated?.FirstOrDefault(m => m.Name == c.Id);
                    var fitStatus = single ? plan.LinearFit!.Status[plan.Index] : plan.MixedFit!.Status[plan.Index];
                    if (plan.NoVarianceDf)
                        (status, detail) = (DifferentialStatus.BelowSupport(DifferentialConstants.TooFewObservations), "no tested term kept for the variance df");
                    else if (fitStatus == FeatureFitStatus.NotConverged)
                        (status, detail) = (DifferentialStatus.NotConverged, "the fit did not converge");
                    else if (fitStatus is FeatureFitStatus.TooFewObservations or FeatureFitStatus.TooFewGroups)
                        (status, detail) = (DifferentialStatus.BelowSupport(DifferentialConstants.TooFewObservations), $"the model reported {fitStatus}");
                    else if (fitStatus == FeatureFitStatus.RankDeficient)
                        (status, detail) = (DifferentialStatus.NotEstimable("rank_deficient"), "the model reported RankDeficient");
                    else if (result is null || !double.IsFinite(result.Estimate[plan.Index]) || !double.IsFinite(result.PValue[plan.Index]))
                        (status, detail) = (DifferentialStatus.NotConverged, "no moderated result (variance or contrast df not finite)");
                    else
                        status = single ? DifferentialStatus.FittedSinglePeptide : DifferentialStatus.Fitted;
                    if (!DifferentialStatus.IsFitted(status))
                        Decide(DifferentialDecisionRules.FitOutcome, design.Stratum, plan.Group, c.Id, detail!, $"reported {status}");
                }

                bool fitted = DifferentialStatus.IsFitted(status);
                int f = plan.Index;
                var values = plan.Values;
                var row = new DifferentialResult
                {
                    DefinitionId = DifferentialConstants.DefinitionId,
                    AnalysisId = _analysisId,
                    QuantStyle = _style,
                    ReporterAcquisition = _isobaric ? _request.ReporterAcquisition : null,
                    Stratum = design.Stratum,
                    Grain = DifferentialConstants.Grain,
                    FeatureId = plan.Group,
                    FeatureAccessions = plan.Info.Accessions,
                    Quantity = DifferentialConstants.Quantity,
                    EffectType = c.Covariate is null ? DifferentialConstants.EffectRatio : DifferentialConstants.EffectSlope,
                    QuantBasis = _basis,
                    ContrastId = c.Id,
                    ContrastLabel = c.Label,
                    Numerator = c.Numerator,
                    Denominator = c.Denominator,
                    Covariate = c.Covariate,
                    CovariateUnit = c.CovariateUnit,
                    CovariateScale = c.CovariateScale,
                    Log2Effect = fitted ? result!.Estimate[f] : null,
                    NaturalEffect = fitted ? Math.Pow(2, result!.Estimate[f]) : null,
                    NaturalUnit = fitted ? (c.CovariateUnit is null ? "ratio" : "ratio_per_" + c.CovariateUnit.Replace(' ', '_')) : null,
                    StandardError = fitted ? result!.StandardError[f] : null,
                    CiLow = fitted ? result!.ConfidenceLow[f] : null,
                    CiHigh = fitted ? result!.ConfidenceHigh[f] : null,
                    CiLevel = fitted ? _request.Settings.CiLevel : null,
                    Statistic = fitted ? result!.T[f] : null,
                    Df = fitted ? result!.DfTotal[f] : null,
                    DfMethod = fitted ? (single ? DifferentialConstants.DfResidual : DifferentialConstants.DfSatterthwaite) : null,
                    PValue = fitted ? result!.PValue[f] : null,
                    AdjustmentMethod = fitted ? DifferentialConstants.Adjustment : null,
                    NSamplesNumerator = pc.Counts.Numerator,
                    NSamplesDenominator = pc.Counts.Denominator,
                    NSamplesWithValue = pc.Counts.WithValue,
                    NSamplesTotal = pc.Counts.Total,
                    NPeptides = values is null ? 0 : values.NuisanceLevels.Count,
                    NObservations = values?.Values.Count ?? 0,
                    NMbrValues = _isobaric ? null : values?.Values.Count(v => v.IsMbr) ?? 0,
                    Status = status!,
                    StatusDetail = fitted ? null : detail,
                    Method = DifferentialConstants.Method,
                    ModelUsed = values is null || values.NuisanceLevels.Count == 0 ? null
                        : single ? DifferentialConstants.ModeratedT : DifferentialConstants.MixedModel,
                    Robust = values is null || values.NuisanceLevels.Count == 0 ? null : !single && _request.Settings.Robust,
                    Normalization = normalization,
                    MethodVersion = DifferentialConstants.MethodVersion,
                };
                _rows.Add(row);
            }
        }

        /// <summary>Benjamini–Hochberg once per family over both models' fitted rows (ruling 4), and the family sizes.</summary>
        private static List<DifferentialResult> Adjust(List<DifferentialResult> rows)
        {
            var output = new List<DifferentialResult>(rows.Count);
            var adjusted = new Dictionary<DifferentialResult, (double, int)>(ReferenceEqualityComparer.Instance);
            foreach (var family in rows.Where(DifferentialResultWriter.InBenjaminiHochbergFamily).GroupBy(DifferentialResultWriter.FamilyKey))
            {
                var members = family.ToList();
                var q = MultipleTesting.BenjaminiHochberg(members.Select(r => r.PValue!.Value).ToArray());
                for (int i = 0; i < members.Count; i++) adjusted[members[i]] = (q[i], members.Count);
            }
            foreach (var r in rows)
                output.Add(adjusted.TryGetValue(r, out var a) ? r with { PAdjusted = a.Item1, FamilySize = a.Item2 } : r);
            return output;
        }

        private DifferentialStratumInfo StratumInfo(AnalysisDesign design)
        {
            var first = design.Samples[0];
            var perLevel = new SortedDictionary<string, int>(StringComparer.Ordinal);
            foreach (var factor in _request.Spec.Factors)
                foreach (var s in design.Samples)
                    if (s.Factors.TryGetValue(factor.Name, out var level))
                        perLevel[$"{factor.Name}={level}"] = perLevel.GetValueOrDefault($"{factor.Name}={level}") + 1;
            return new DifferentialStratumInfo(design.Stratum, first.Strata, _request.StratumSource, perLevel,
                design.NotRun.Select(n => new DifferentialNotRun(n.Label, n.Reason)).ToList());
        }

        private DifferentialNormalizationInfo IsobaricNormalizationInfo(string stratum, IReadOnlyDictionary<string, AnalysisSample> samples)
        {
            var shifts = _tmtNormalization!.Shifts.Where(s => samples.ContainsKey(s.SampleId) && double.IsFinite(s.Shift))
                .ToDictionary(s => string.Create(CultureInfo.InvariantCulture, $"{s.SampleId}@plex{s.Plex}.run{s.TechnicalReplicate}"),
                    s => s.Shift, StringComparer.Ordinal);
            return new DifferentialNormalizationInfo(stratum, null, _tmtNormalization.Setting, 0, 0, shifts, _tmtNormalization.Warnings);
        }
    }

    private sealed record StratumContrast(Contrast Contrast, string Estimability);

    private sealed record SideCounts(int? Numerator, int? Denominator, int WithValue, int Total);

    private sealed record PlannedContrast(Contrast Contrast, SideCounts Counts, string? Status, string? Detail);

    /// <summary>One protein group's plan, then its fit.</summary>
    private sealed class FeaturePlan
    {
        public FeaturePlan(string group, ProteinGroupInfo info, FeatureObservations? values)
        {
            Group = group;
            Info = info;
            Values = values;
        }

        public string Group { get; }
        public ProteinGroupInfo Info { get; }
        public FeatureObservations? Values { get; }
        public List<PlannedContrast> Contrasts { get; } = new();
        public FeatureDesign? Design { get; set; }
        public double[,]? CommonRows { get; set; }
        public NestedMixedModelFit? MixedFit { get; set; }
        public LinearModelFit? LinearFit { get; set; }
        public IReadOnlyList<ModeratedContrast>? Moderated { get; set; }
        public int Index { get; set; }
        public bool NoVarianceDf { get; set; }
    }
}
