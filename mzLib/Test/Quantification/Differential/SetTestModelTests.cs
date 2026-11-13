using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// SETTEST1 S2 (QuantProject <c>DEF-DIFF-SET-CONTRAST</c>, <c>-MEMBER</c>, <c>-MODEL</c>, <c>-VERDICT</c>): contexts to a
/// stratum and a design; the set model against lmerTest (<c>ReferenceData/make_settest_model_fixtures.R</c>, R never run
/// here), plain and (S2b, GR-8) under msqrob2-style robust weights run to convergence, or one round when they never settle
/// (GR-22 amended 2026-11-03); the member states read from STAT1's
/// rows; every prediction's p, 90% interval, <c>contradicted</c> rule and smallest margin against values from R's
/// <c>pt</c>/<c>qt</c>; the guards; and the order of the verdict words.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SetTestModelTests
{
    private const double Fit = 1e-6;

    /// <summary>
    /// Satterthwaite df: lmerTest differentiates numerically, so df are compared as <c>ReferenceComparisonNestedLmmTests</c>
    /// compares them (#1471), relative to max(1, |expected|).
    /// </summary>
    private const double Df = 1e-4;

    private static void Close(double actual, string expectedText, double tolerance, string what)
    {
        double expected = D(expectedText);
        TestContext.Out.WriteLine($"{what}: {actual:R} vs {expected:R} (|diff| {Math.Abs(actual - expected):E2})");
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Abs(expected))),
            $"{what}: {actual:R} vs {expected:R}");
    }

    private const string Header =
        "hypothesis_id\tkind\tprediction\tmember\tspecies\tdataset\tfactor\ttest_level\treference_level\tcontexts\tmrna_log2\tindependence\tpost_hoc";

    // ── fixtures ─────────────────────────────────────────────────────────────

    private static List<string[]> Tsv(string file) =>
        File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential",
            "ReferenceData", file)).Skip(1).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static SetTestFreeze Freeze(params string[] rows) =>
        SetTestFreeze.Read(Encoding.UTF8.GetBytes(string.Join('\n', new[] { Header }.Concat(rows)) + "\n"), "t");

    private static string Row(string member, string prediction = "up", string kind = "set", string contexts = "",
        string mrna = "", string test = "old", string reference = "young", string factor = "age group") =>
        string.Join('\t', "H1", kind, prediction, member, "mouse", "PXD1", factor, test, reference, contexts, mrna,
            "independent", "no");

    private static ObservationTable Table(IEnumerable<string> samples, IEnumerable<(string Member, string Peptide, string Sample, double Log2)> values,
        IEnumerable<string>? groups = null)
    {
        var runs = samples.Select(s => new ObservationRun(s + ".raw", s)).ToList();
        var peptides = values.GroupBy(v => (v.Member, v.Peptide))
            .Select(g => new PeptideRunValues(g.Key.Peptide, g.Key.Peptide, new[] { g.Key.Member },
                g.ToDictionary(v => v.Sample + ".raw", v => new RunValue(Math.Pow(2, v.Log2), ObservationState.Quantified))))
            .ToList();
        var infos = (groups ?? values.Select(v => v.Member)).Distinct().Select(m => new ProteinGroupInfo(m, 0.001, false, false));
        return ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, infos, QuantBasis.MsmsOnly);
    }

    private static readonly DesignSpec AgeSpec = new() { Factors = new[] { new FactorSpec("age group", "young") } };

    private static (SetContrastSetup Setup, ObservationTable Table, List<string> Groups, string[] R) Scenario(string name,
        Func<int, string>? contextsOfRow = null)
    {
        var samples = Tsv($"settest_model_{name}_samples.tsv")
            .Select(c => new AnalysisSample(c[0])
            {
                Factors = new Dictionary<string, string> { ["age group"] = c[1], ["treatment"] = c[2] },
            }).ToList();
        var values = Tsv($"settest_model_{name}_values.tsv").Select(c => (c[0], c[1], c[2], D(c[3]))).ToList();
        var table = Table(samples.Select(s => s.SampleId), values);
        var groups = values.Select(v => v.Item1).Distinct().Order(StringComparer.Ordinal).ToList();
        string contexts = name == "context" ? "treatment=a|b" : "";
        var contrast = Freeze(groups.Select((g, i) => Row(g, contexts: contextsOfRow?.Invoke(i) ?? contexts)).ToArray())
            .Parts.Single().Contrasts.Single();
        return (SetContrastSetup.Build(contrast, samples, AgeSpec), table, groups, Tsv($"settest_model_{name}_result.tsv").Single());
    }

    // ── the set model against lmerTest ───────────────────────────────────────

    [TestCase("balanced")]
    [TestCase("missing")]
    [TestCase("context")]
    public void TheSetModel_MatchesLmerTest(string name)
    {
        var (setup, table, groups, r) = Scenario(name);
        Assert.That(setup.NotRun, Is.Null);

        var estimate = SetModel.FitSet(setup, groups, table, out string? notRun, robust: false);

        Assert.That(notRun, Is.Null);
        Close(estimate!.Log2Effect, r[0], Fit, name + " estimate");
        Close(estimate.StandardError, r[1], Fit, name + " se");
        Close(estimate.Df, r[2], Df, name + " Satterthwaite df");
        Close(estimate.ResidualVariance, r[3], Fit, name + " sigma^2");
        Close(estimate.SampleVariance, r[4], Fit, name + " var (1|sample)");
        Close(estimate.SampleMemberVariance, r[5], Fit, name + " var (1|sample:member)");
        Assert.That(estimate.Terms, Is.EqualTo(int.Parse(r[7])));
        Assert.That(estimate.Peptides, Is.EqualTo(int.Parse(r[8])));
        Assert.That(estimate.Observations, Is.EqualTo(int.Parse(r[9])));
        Assert.That(estimate.RandomTerms, Is.EqualTo("sample;sample:member"));
        Assert.That(estimate.DfMethod, Is.EqualTo("satterthwaite"));
        Assert.That(estimate.ModelUsed, Is.EqualTo(SetModel.SetMixedModel));
    }

    [Test]
    public void TheContextFactor_EntersTheDesignAsAnAdditiveTerm()
    {
        var (setup, _, _, _) = Scenario("context");
        Assert.That(setup.Design!.ColumnNames, Is.EqualTo(new[] { "intercept", "age group=old", "treatment=b" }));
        Assert.That(setup.NumeratorSamples, Has.Count.EqualTo(6));
        Assert.That(setup.DenominatorSamples, Has.Count.EqualTo(6));
    }

    [Test]
    public void ContextLevelsInAnotherOrder_MoveTheTermsReferenceLevel_NotTheEstimate()
    {
        // The first row lists b first and the others a first. The freeze makes them one contrast with the first row's
        // order, so the added term's reference level is b.
        var (setup, table, groups, r) = Scenario("context", i => i == 0 ? "treatment=b|a" : "treatment=a|b");
        Assert.That(setup.Design!.ColumnNames, Is.EqualTo(new[] { "intercept", "age group=old", "treatment=a" }));

        var estimate = SetModel.FitSet(setup, groups, table, out string? notRun, robust: false);

        // Recoding a nuisance factor spans the same columns, so lmerTest's fit with reference a still holds.
        Assert.That(notRun, Is.Null);
        Close(estimate!.Log2Effect, r[0], Fit, "reference b estimate");
        Close(estimate.StandardError, r[1], Fit, "reference b se");
        Close(estimate.Df, r[2], Df, "reference b Satterthwaite df");
        Close(estimate.ResidualVariance, r[3], Fit, "reference b sigma^2");
        Close(estimate.SampleVariance, r[4], Fit, "reference b var (1|sample)");
        Close(estimate.SampleMemberVariance, r[5], Fit, "reference b var (1|sample:member)");
    }

    [Test]
    public void AGroupEntersTheSetAtMostOnce()
    {
        var (setup, table, groups, _) = Scenario("balanced");
        Assert.That(() => SetModel.FitSet(setup, groups.Append(groups[0]).ToList(), table, out _),
            Throws.ArgumentException.With.Message.Contains("at most once"));
    }

    // ── the set model under robust weighting (S2b, GR-8) ─────────────────────

    [TestCase("balanced")]
    [TestCase("missing")]
    [TestCase("context", Explicit = true, Reason = "6 s (25 rounds on 174 values); balanced and missing run the same code in CI")]
    public void TheRobustSetModel_MatchesMsqrob2StyleWeightsAndLmerTest(string name)
    {
        var (setup, table, groups, _) = Scenario(name);
        var r = Tsv($"settest_model_{name}_robust_result.tsv").Single();
        var expected = Tsv($"settest_model_{name}_robust_weights.tsv").ToDictionary(c => (c[0], c[1], c[2]), c => D(c[3]));

        var estimate = SetModel.FitSet(setup, groups, table, out string? notRun, robust: true);

        Assert.That(notRun, Is.Null);
        Assert.That(estimate!.Robust, Is.True);
        Assert.That(bool.Parse(r[8]), Is.True, "R settles");
        Assert.That(estimate.RobustFellBack, Is.False, "settled, so no fallback");
        Assert.That(estimate.Weights.Select(v => (v.ProteinGroup, v.Peptide, v.Sample)), Is.EquivalentTo(expected.Keys));
        double worst = estimate.Weights.Max(v => Math.Abs(v.Weight - expected[(v.ProteinGroup, v.Peptide, v.Sample)]));
        TestContext.Out.WriteLine($"{name} weights: largest |diff| {worst:E2}");
        foreach (var v in estimate.Weights)
            Assert.That(v.Weight, Is.EqualTo(expected[(v.ProteinGroup, v.Peptide, v.Sample)]).Within(Fit),
                $"{name} weight {v.Peptide} {v.Sample}");
        Assert.That(estimate.Weights.Count(v => v.Weight < 1), Is.EqualTo(int.Parse(r[9])), "values down-weighted, as R");
        Close(estimate.Log2Effect, r[0], Fit, name + " robust estimate");
        Close(estimate.StandardError, r[1], Fit, name + " robust se");
        Close(estimate.Df, r[2], Df, name + " robust Satterthwaite df");
        Close(estimate.ResidualVariance, r[3], Fit, name + " robust sigma^2");
        Close(estimate.SampleVariance, r[4], Fit, name + " robust var (1|sample)");
        Close(estimate.SampleMemberVariance, r[5], Fit, name + " robust var (1|sample:member)");
        Assert.That(estimate.Observations, Is.EqualTo(int.Parse(r[12])));
    }

    [Test]
    public void TheSetFit_IsRobustByDefault_AsStat1Is()
    {
        var (setup, table, groups) = SmallSet();

        var byDefault = SetModel.FitSet(setup, groups, table, out _);
        var plain = SetModel.FitSet(setup, groups, table, out _, robust: false);

        Assert.That(byDefault!.Robust, Is.True);
        Assert.That(byDefault.RobustFellBack, Is.False);
        Assert.That(byDefault.Weights.Count(v => v.Weight < 1), Is.GreaterThan(0));
        Assert.That(plain!.Robust, Is.False);
        Assert.That(plain.RobustFellBack, Is.Null, "no reweighting, so nothing to fall back from");
        Assert.That(plain.Weights, Has.Count.EqualTo(plain.Observations));
        Assert.That(plain.Weights.Select(v => v.Weight), Is.All.EqualTo(1.0));
        Assert.That(byDefault.Log2Effect, Is.Not.EqualTo(plain.Log2Effect));
    }

    [Test]
    public void FittingASetupThatIsNotRun_IsRefused()
    {
        var (setup, table, groups, _) = Scenario("balanced");
        foreach (var bad in new[] { setup with { NotRun = SetNotRunReason.SpansStrata }, setup with { Design = null }, setup with { Contrast = null } })
            Assert.That(() => SetModel.FitSet(bad, groups, table, out _), Throws.ArgumentException.With.Message.Contains("not run"));
    }

    [Test]
    public void ATableMissingASampleOfTheSetup_IsRefused()
    {
        var setup = SetContrastSetup.Build(Freeze(Row("G1")).Parts.Single().Contrasts.Single(), FourSamples, AgeSpec);
        var three = new[] { "o1", "o2", "y1" };
        var table = Table(three, three.Select((s, i) => ("G1", "G1PEPA", s, All4[i])));
        Assert.That(() => SetModel.FitSet(setup, new[] { "G1" }, table, out _),
            Throws.ArgumentException.With.Message.Contains("no sample 'y2'"));
    }

    [TestCase("not_converged", null, null, null)]
    [TestCase("fitted", null, 0.1, 12.0)]
    [TestCase("fitted", 0.4, null, 12.0)]
    [TestCase("fitted", 0.4, 0.1, null)]
    public void ASinglesStat1RowWithoutAnEffect_IsRefused(string status, double? effect, double? se, double? df) =>
        Assert.That(() => SetModel.FromStat1Row(Stat1("G1", status, effect, se, df)),
            Throws.ArgumentException.With.Message.Contains("carries no effect"));

    // ── contexts to a stratum and a design ───────────────────────────────────

    private static List<AnalysisSample> Samples(params (string Id, string Age, string Tissue, string Treatment, string? Individual)[] s) =>
        s.Select(x => new AnalysisSample(x.Id)
        {
            Factors = new Dictionary<string, string> { ["age group"] = x.Age, ["treatment"] = x.Treatment },
            Strata = x.Tissue.Length == 0 ? new Dictionary<string, string>() : new Dictionary<string, string> { ["organism part"] = x.Tissue },
            Individual = x.Individual,
        }).ToList();

    private static readonly List<AnalysisSample> TwoTissues = Samples(
        ("l1", "young", "liver", "a", null), ("l2", "young", "liver", "b", null), ("l3", "old", "liver", "a", null), ("l4", "old", "liver", "b", null),
        ("k1", "young", "kidney", "a", null), ("k2", "young", "kidney", "b", null), ("k3", "old", "kidney", "a", null), ("k4", "old", "kidney", "b", null));

    private static SetContrastSetup Setup(string contexts, List<AnalysisSample> samples, string test = "old", string factor = "age group",
        string reference = "young") =>
        SetContrastSetup.Build(Freeze(Row("GeneA", contexts: contexts, test: test, reference: reference, factor: factor)).Parts.Single()
            .Contrasts.Single(), samples, AgeSpec);

    [Test]
    public void AStratumFactorInContexts_SelectsTheStratum()
    {
        var setup = Setup("organism part=liver", TwoTissues);
        Assert.That(setup.NotRun, Is.Null);
        Assert.That(setup.Stratum, Is.EqualTo("organism part=liver"));
        Assert.That(setup.ContrastSamples, Is.EqualTo(new[] { "l1", "l2", "l3", "l4" }));
        Assert.That(setup.Design!.ColumnNames, Is.EqualTo(new[] { "intercept", "age group=old" }));
    }

    [Test]
    public void ContextsSpanningTwoStrata_AreNotRun()
    {
        Assert.That(Setup("", TwoTissues).NotRun, Is.EqualTo(SetNotRunReason.SpansStrata));
        Assert.That(Setup("organism part=liver|kidney", TwoTissues).NotRun, Is.EqualTo(SetNotRunReason.SpansStrata));
    }

    [Test]
    public void AModelFactorInContexts_RestrictsTheSamples()
    {
        var setup = Setup("organism part=liver; treatment=a", TwoTissues);
        Assert.That(setup.ContrastSamples, Is.EqualTo(new[] { "l1", "l3" }));
        // One level of treatment is left, so it adds no column.
        Assert.That(setup.Design!.ColumnNames, Is.EqualTo(new[] { "intercept", "age group=old" }));
    }

    [Test]
    public void AStratumWithTwoSamplesOfOneIndividual_IsNotRun()
    {
        var paired = Samples(("a1", "young", "", "a", "m1"), ("a2", "old", "", "a", "m1"), ("b1", "young", "", "a", "m2"), ("b2", "old", "", "a", "m2"));
        Assert.That(Setup("", paired).NotRun, Is.EqualTo(SetNotRunReason.RepeatedIndividuals));
        var oneEach = Samples(("a1", "young", "", "a", "m1"), ("b2", "old", "", "a", "m2"), ("c1", "young", "", "a", "m3"), ("d2", "old", "", "a", "m4"));
        Assert.That(Setup("", oneEach).NotRun, Is.Null);
    }

    [TestCase("diet=chow", "age group", "old", "not_in_design:context_factor=diet")]
    [TestCase("organism part=liver", "genotype", "old", "not_in_design:factor=genotype")]
    [TestCase("organism part=liver", "age group", "ancient", "not_in_design:age group=ancient")]
    [TestCase("organism part=heart", "age group", "old", "not_in_design:no_samples")]
    public void SomethingTheDesignDoesNotHave_IsNotRunAndNamed(string contexts, string factor, string test, string reason) =>
        Assert.That(Setup(contexts, TwoTissues, test, factor).NotRun, Is.EqualTo(reason));

    [Test]
    public void AContrastAliasedWithAContextFactor_IsNotRunAsNotEstimable()
    {
        // treatment follows age group exactly, and contexts add it to the model.
        var aliased = Samples(("a", "young", "", "x", null), ("b", "young", "", "x", null), ("c", "old", "", "y", null), ("d", "old", "", "y", null));
        Assert.That(Setup("treatment=x|y", aliased).NotRun, Is.EqualTo(AnalysisDesign.NotEstimableConfounded));
    }

    [Test]
    public void AReferenceLevelTheDesignDoesNotHave_IsNotRunAndNamed() =>
        Assert.That(Setup("organism part=liver", TwoTissues, reference: "middle").NotRun, Is.EqualTo("not_in_design:age group=middle"));

    [Test]
    public void AContrastOnAFactorTheSpecLacks_EntersTheDesignWithTheContrastsReferenceLevel()
    {
        // AgeSpec names only age group; a contrast on treatment adds treatment, with reference level a.
        var setup = Setup("organism part=liver", TwoTissues, test: "b", factor: "treatment", reference: "a");
        Assert.That(setup.NotRun, Is.Null);
        Assert.That(setup.Design!.ColumnNames, Is.EqualTo(new[] { "intercept", "age group=old", "treatment=b" }));
        Assert.That(setup.NumeratorSamples, Is.EqualTo(new[] { "l2", "l4" }));
        Assert.That(setup.DenominatorSamples, Is.EqualTo(new[] { "l1", "l3" }));
    }

    [Test]
    public void ASampleLackingAContextsFactor_IsLeftOut()
    {
        // x1 has no organism part and x2 no treatment, so "organism part=liver; treatment=a|b" keeps neither.
        var samples = TwoTissues.Take(4)
            .Append(new AnalysisSample("x1") { Factors = new Dictionary<string, string> { ["age group"] = "old", ["treatment"] = "a" } })
            .Append(new AnalysisSample("x2")
            {
                Factors = new Dictionary<string, string> { ["age group"] = "old" },
                Strata = new Dictionary<string, string> { ["organism part"] = "liver" },
            }).ToList();
        var setup = Setup("organism part=liver; treatment=a|b", samples);
        Assert.That(setup.NotRun, Is.Null);
        Assert.That(setup.ContrastSamples, Is.EqualTo(new[] { "l1", "l2", "l3", "l4" }));
    }

    // ── the member states read from STAT1's rows ─────────────────────────────

    private static readonly string[] Four = { "o1", "o2", "y1", "y2" };
    private static readonly List<AnalysisSample> FourSamples = Samples(
        ("o1", "old", "", "a", null), ("o2", "old", "", "a", null), ("y1", "young", "", "a", null), ("y2", "young", "", "a", null));

    private static DifferentialResult Stat1(string group, string status, double? effect = null, double? se = null, double? df = null) => new()
    {
        DefinitionId = "QuantProject:DEF-DIFF-ABUNDANCE v1", AnalysisId = "a", Stratum = "all", Grain = "protein_group",
        FeatureId = group, Quantity = "abundance", EffectType = "abundance_log2_ratio", ContrastId = "c1", Status = status,
        Method = "moderated",
        Log2Effect = effect, StandardError = se, Df = df, DfMethod = df is null ? null : "moderated_residual",
    };

    /// <summary>
    /// A member per group with two unique peptides: the given values in o1, o2, y1, y2 (NaN = none), and a second peptide
    /// 0.4 higher with a small fixed jitter, so <c>(1|sample:member)</c> is separable from the residual.
    /// </summary>
    private static (ObservationTable Table, IReadOnlyList<SetMemberAssignment> Resolved) Members(params (string Group, double[] Values)[] groups)
    {
        double[] jitter = { 0.03, -0.05, 0.02, -0.01, 0.04, -0.03 };
        var values = groups.SelectMany((g, k) => Four.SelectMany((s, i) => new[]
            {
                (g.Group, g.Group + "PEPA", s, g.Values[i]),
                (g.Group, g.Group + "PEPB", s, g.Values[i] + 0.4 + jitter[(i + k) % jitter.Length]),
            }))
            .Where(v => double.IsFinite(v.Item4)).ToList();
        var table = Table(Four, values, groups.Select(g => g.Group));
        var resolved = groups.Select(g => new SetMemberAssignment { Member = g.Group, ProteinGroup = g.Group }).ToList();
        return (table, resolved);
    }

    private static readonly double[] All4 = { 20, 21, 20.5, 20.7 };

    [TestCase("below_support:min_2_per_side", "too_few")]
    [TestCase("absent_in_both", "too_few")]
    [TestCase("absent_in_denominator", "seen_one_side:numerator")]
    [TestCase("absent_in_numerator", "seen_one_side:denominator")]
    [TestCase("fitted", "tested")]
    [TestCase("not_converged", "tested")]
    public void AMembersState_IsReadFromItsStat1Row(string status, string state)
    {
        var (table, resolved) = Members(("G1", All4));
        var setup = SetContrastSetup.Build(Freeze(Row("G1")).Parts.Single().Contrasts.Single(), FourSamples, AgeSpec);
        var result = SetMemberStateFromRow.Apply(resolved, SetPrediction.Up, setup, table,
            new Dictionary<string, DifferentialResult> { ["G1"] = Stat1("G1", status) }).Single();
        Assert.That(result.State, Is.EqualTo(state));
        Assert.That(result.InEstimate, Is.EqualTo(state == "tested"));
    }

    [Test]
    public void Incomplete_AppliesToNoChangePredictionsOnly()
    {
        // G1 has a value in one of the two old samples.
        var (table, resolved) = Members(("G1", new[] { 20, double.NaN, 20.5, 20.7 }));
        var setup = SetContrastSetup.Build(Freeze(Row("G1")).Parts.Single().Contrasts.Single(), FourSamples, AgeSpec);
        var rows = new Dictionary<string, DifferentialResult> { ["G1"] = Stat1("G1", "fitted") };

        var unchanged = SetMemberStateFromRow.Apply(resolved, SetPrediction.Unchanged, setup, table, rows).Single();
        var up = SetMemberStateFromRow.Apply(resolved, SetPrediction.Up, setup, table, rows).Single();

        Assert.That(unchanged.State, Is.EqualTo(SetMemberStateFromRow.Incomplete));
        Assert.That(unchanged.SamplesWithValueNumerator, Is.EqualTo(1));
        Assert.That(unchanged.SamplesWithValueDenominator, Is.EqualTo(2));
        Assert.That(unchanged.InEstimate, Is.True);
        Assert.That(up.State, Is.EqualTo(SetMemberStateFromRow.Tested));
    }

    [TestCase(1, 4, true)]
    [TestCase(2, 4, false)]
    [TestCase(2, 5, true)]
    [TestCase(3, 5, false)]
    [TestCase(1, 1, true)]
    public void ASideIsIncomplete_BelowHalfOrBelowTwo(int withValue, int size, bool incomplete) =>
        Assert.That(SetMemberStateFromRow.IncompleteSide(withValue, size), Is.EqualTo(incomplete));

    [Test]
    public void AGroupWithNoStat1Row_IsRefused()
    {
        var (table, resolved) = Members(("G1", All4));
        var setup = SetContrastSetup.Build(Freeze(Row("G1")).Parts.Single().Contrasts.Single(), FourSamples, AgeSpec);
        Assert.That(() => SetMemberStateFromRow.Apply(resolved, SetPrediction.Up, setup, table, new Dictionary<string, DifferentialResult>()),
            Throws.ArgumentException.With.Message.Contains("no STAT1 row"));
    }

    // ── every prediction's test, against R's pt and qt ───────────────────────

    private const double Qt95Df10 = 1.8124611228116756;   // qt(0.95, 10)
    private const double PUpT25 = 0.015723422118304388;   // pt(2.5, 10, lower.tail = FALSE)
    private const double PtMinus4 = 0.0012591663123683451; // pt(-4, 10)
    private const double PtMinus5 = 0.00026866680137822624; // pt(-5, 10)

    private static SetEstimate Est(double b, double se, double df = 10) =>
        new() { Log2Effect = b, StandardError = se, Df = df, DfMethod = "satterthwaite", ModelUsed = SetModel.SetMixedModel };

    [Test]
    public void Up_IsOneSidedAgainstZero_AndContradictedWhenTheIntervalIsBelowIt()
    {
        var t = SetVerdict.Test(SetPrediction.Up, Est(0.5, 0.2));
        Assert.That(t.PValue, Is.EqualTo(PUpT25).Within(1e-12));
        Assert.That(t.Ci90Low, Is.EqualTo(0.5 - Qt95Df10 * 0.2).Within(1e-12));
        Assert.That(t.Ci90High, Is.EqualTo(0.5 + Qt95Df10 * 0.2).Within(1e-12));
        Assert.That(t.Contradicted, Is.False);
        Assert.That(t.Bound, Is.EqualTo("0"));
        Assert.That(t.SmallestMarginPassed, Is.Null);
        Assert.That(SetVerdict.Test(SetPrediction.Up, Est(-0.5, 0.2)).Contradicted, Is.True);
        Assert.That(SetVerdict.Test(SetPrediction.Up, Est(-0.3, 0.2)).Contradicted, Is.False);
    }

    [Test]
    public void Down_MirrorsUp()
    {
        var t = SetVerdict.Test(SetPrediction.Down, Est(-0.5, 0.2));
        Assert.That(t.PValue, Is.EqualTo(PUpT25).Within(1e-12));
        Assert.That(SetVerdict.Test(SetPrediction.Down, Est(0.5, 0.2)).Contradicted, Is.True);
    }

    [Test]
    public void Unchanged_IsTost_TheLargerOfTwoOneSidedPs()
    {
        var t = SetVerdict.Test(SetPrediction.Unchanged, Est(0.1, 0.1));
        Assert.That(t.PValue, Is.EqualTo(PtMinus4).Within(1e-12));          // max(pt(-6), pt(-4))
        Assert.That(t.Bound, Is.EqualTo("±0.5"));
        Assert.That(t.Margin, Is.EqualTo(0.5));
        Assert.That(t.SmallestMarginPassed, Is.EqualTo(0.1 + Qt95Df10 * 0.1).Within(1e-12));
        Assert.That(t.Contradicted, Is.False);
        Assert.That(SetVerdict.Test(SetPrediction.Unchanged, Est(1.0, 0.1)).Contradicted, Is.True);
        Assert.That(SetVerdict.Test(SetPrediction.Unchanged, Est(-1.0, 0.1)).Contradicted, Is.True);
        Assert.That(SetVerdict.Test(SetPrediction.Unchanged, Est(0.6, 0.1)).Contradicted, Is.False);
        // The smallest margin is the interval's larger end in absolute value, on either side.
        Assert.That(SetVerdict.Test(SetPrediction.Unchanged, Est(-0.1, 0.1)).SmallestMarginPassed,
            Is.EqualTo(0.1 + Qt95Df10 * 0.1).Within(1e-12));
    }

    [Test]
    public void NotUpAndNotDown_AreOneSidedAgainstTheMargin()
    {
        var notUp = SetVerdict.Test(SetPrediction.NotUp, Est(0.1, 0.1));
        Assert.That(notUp.PValue, Is.EqualTo(PtMinus4).Within(1e-12));
        Assert.That(notUp.Bound, Is.EqualTo("+0.5"));
        Assert.That(notUp.SmallestMarginPassed, Is.EqualTo(0.1 + Qt95Df10 * 0.1).Within(1e-12));
        Assert.That(SetVerdict.Test(SetPrediction.NotUp, Est(1.0, 0.1)).Contradicted, Is.True);
        Assert.That(SetVerdict.Test(SetPrediction.NotUp, Est(-0.4, 0.1)).SmallestMarginPassed, Is.EqualTo(0));

        var notDown = SetVerdict.Test(SetPrediction.NotDown, Est(-0.1, 0.1));
        Assert.That(notDown.PValue, Is.EqualTo(PtMinus4).Within(1e-12));
        Assert.That(notDown.Bound, Is.EqualTo("-0.5"));
        Assert.That(notDown.SmallestMarginPassed, Is.EqualTo(0.1 + Qt95Df10 * 0.1).Within(1e-12));
        Assert.That(SetVerdict.Test(SetPrediction.NotDown, Est(-1.0, 0.1)).Contradicted, Is.True);
    }

    [Test]
    public void AttenuatedDown_TestsBothBounds_ReadingB()
    {
        // m̄ = -2, f = 0.5, so f·m̄ = -1.
        var t = SetVerdict.Test(SetPrediction.AttenuatedDown, Est(-0.5, 0.1), mrnaMean: -2);
        Assert.That(t.PValue, Is.EqualTo(PtMinus5).Within(1e-12));          // max(pt(-5), pt(-5))
        Assert.That(t.Bound, Is.EqualTo("0;-1"));
        Assert.That(t.Attenuation, Is.EqualTo(0.5));
        Assert.That(t.SmallestMarginPassed, Is.EqualTo((-0.5 - Qt95Df10 * 0.1) / -2).Within(1e-12));
        Assert.That(t.Contradicted, Is.False);
        // No decrease refutes it, and so does a fall beyond f·m̄.
        Assert.That(SetVerdict.Test(SetPrediction.AttenuatedDown, Est(0.3, 0.1), mrnaMean: -2).Contradicted, Is.True);
        Assert.That(SetVerdict.Test(SetPrediction.AttenuatedDown, Est(-1.5, 0.1), mrnaMean: -2).Contradicted, Is.True);
        Assert.That(SetVerdict.Test(SetPrediction.AttenuatedDown, Est(0.1, 0.1), mrnaMean: -2).SmallestMarginPassed, Is.Null);
    }

    [Test]
    public void AttenuatedUp_MirrorsAttenuatedDown()
    {
        var t = SetVerdict.Test(SetPrediction.AttenuatedUp, Est(0.5, 0.1), mrnaMean: 2);
        Assert.That(t.PValue, Is.EqualTo(PtMinus5).Within(1e-12));
        Assert.That(t.SmallestMarginPassed, Is.EqualTo((0.5 + Qt95Df10 * 0.1) / 2).Within(1e-12));
        Assert.That(SetVerdict.Test(SetPrediction.AttenuatedUp, Est(-0.3, 0.1), mrnaMean: 2).Contradicted, Is.True);
        Assert.That(SetVerdict.Test(SetPrediction.AttenuatedUp, Est(1.5, 0.1), mrnaMean: 2).Contradicted, Is.True);
    }

    [Test]
    public void AFrozenMarginAndFraction_ReplaceTheDefaults()
    {
        var t = SetVerdict.Test(SetPrediction.Unchanged, Est(0.1, 0.1), margin: 0.25);
        Assert.That(t.Bound, Is.EqualTo("±0.25"));
        Assert.That(SetVerdict.Test(SetPrediction.AttenuatedDown, Est(-0.5, 0.1), attenuation: 0.25, mrnaMean: -2).Bound,
            Is.EqualTo("0;-0.5"));
    }

    [Test]
    public void InfiniteDf_UsesTheNormal()
    {
        var t = SetVerdict.Test(SetPrediction.Up, Est(0.5, 0.2, double.PositiveInfinity));
        Assert.That(t.PValue, Is.EqualTo(0.0062096653257761349).Within(1e-12)); // pnorm(-2.5)
        Assert.That(t.Ci90High - 0.5, Is.EqualTo(1.6448536269514715 * 0.2).Within(1e-12));
    }

    [TestCase(SetPrediction.AttenuatedDown, 1.0)]
    [TestCase(SetPrediction.AttenuatedDown, 0.0)]
    [TestCase(SetPrediction.AttenuatedDown, double.NaN)]
    [TestCase(SetPrediction.AttenuatedDown, null)]
    [TestCase(SetPrediction.AttenuatedUp, -1.0)]
    [TestCase(SetPrediction.AttenuatedUp, 0.0)]
    [TestCase(SetPrediction.AttenuatedUp, null)]
    public void AnAttenuatedPredictionWithoutMrnaOnItsSide_IsRefused(SetPrediction prediction, double? mrnaMean) =>
        Assert.That(() => SetVerdict.Test(prediction, Est(-0.5, 0.1), mrnaMean: mrnaMean),
            Throws.ArgumentException.With.Message.Contains("on its side of 0"));

    [TestCase(double.NaN, 0.1, 10.0)]
    [TestCase(double.PositiveInfinity, 0.1, 10.0)]
    [TestCase(0.5, 0.0, 10.0)]
    [TestCase(0.5, double.NaN, 10.0)]
    [TestCase(0.5, 0.1, 0.0)]
    public void AnEstimateWithoutAFiniteEffectAPositiveSeAndPositiveDf_IsRefused(double b, double se, double df) =>
        Assert.That(() => SetVerdict.Test(SetPrediction.Up, Est(b, se, df)),
            Throws.ArgumentException.With.Message.Contains("finite effect"));

    [TestCase(0.0, null)]
    [TestCase(-0.5, null)]
    [TestCase(double.PositiveInfinity, null)]
    [TestCase(null, 0.0)]
    [TestCase(null, 1.0)]
    [TestCase(null, double.NaN)]
    public void AMarginOrFractionOutOfRange_IsRefused(double? margin, double? attenuation) =>
        Assert.That(() => SetVerdict.Test(SetPrediction.AttenuatedDown, Est(-0.5, 0.1), margin, attenuation, mrnaMean: -2),
            Throws.InstanceOf<ArgumentOutOfRangeException>());

    [Test]
    public void AnUnknownPrediction_IsRefused() =>
        Assert.That(() => SetVerdict.Test((SetPrediction)99, Est(0.5, 0.2)), Throws.InstanceOf<ArgumentOutOfRangeException>());

    // ── a part on a contrast: singles, the minimum, the guards, the verdict words ──

    private static (SetHypothesisPart Part, SetContrastOutcome Outcome) RunPart(string prediction, string kind,
        (string Group, double[] Values, string Status)[] groups, string mrna = "", double? stat1Effect = null, bool? robust = null)
    {
        var freeze = Freeze(groups.Select(g => Row(g.Group, prediction, kind, mrna: mrna)).ToArray());
        var part = freeze.Parts.Single();
        var contrast = part.Contrasts.Single();
        var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
        var (table, resolved) = Members(groups.Select(g => (g.Group, g.Values)).ToArray());
        var rows = groups.ToDictionary(g => g.Group, g => Stat1(g.Group, g.Status, stat1Effect, stat1Effect is null ? null : 0.1, stat1Effect is null ? null : 12));
        return (part, robust is null
            ? SetTest.Run(part, contrast, setup, resolved, table, rows)
            : SetTest.Run(part, contrast, setup, resolved, table, rows, robust.Value));
    }

    private static double[] Shifted(double d) => new[] { 20 + d, 21 + d, 20.5, 20.7 };

    /// <summary>Three members of two peptides each over o1, o2, y1, y2: a set the model fits in a moment.</summary>
    private static (SetContrastSetup Setup, ObservationTable Table, List<string> Groups) SmallSet()
    {
        var setup = SetContrastSetup.Build(Freeze(Row("G1"), Row("G2"), Row("G3")).Parts.Single().Contrasts.Single(), FourSamples, AgeSpec);
        var (table, _) = Members(("G1", Shifted(1.0)), ("G2", Shifted(1.3)), ("G3", Shifted(0.8)));
        return (setup, table, new List<string> { "G1", "G2", "G3" });
    }

    [TestCase(null, true)]
    [TestCase(false, false)]
    public void ASet_IsFittedWithTheRobustSettingItIsGiven(bool? robust, bool expected)
    {
        var (_, outcome) = RunPart("up", "set",
            new[] { ("G1", Shifted(1.0), "fitted"), ("G2", Shifted(1.3), "fitted"), ("G3", Shifted(0.8), "fitted") }, robust: robust);
        var estimate = outcome.Test!.Estimate;
        Assert.That(estimate.Robust, Is.EqualTo(expected));
        Assert.That(estimate.Weights.Any(v => v.Weight < 1), Is.EqualTo(expected));
    }

    [TestCase(true)]
    [TestCase(false)]
    public void ASingle_CarriesItsStat1RowsRobustSetting(bool robust)
    {
        var estimate = SetModel.FromStat1Row(Stat1("G1", "fitted", 0.4, 0.1, 12) with { Robust = robust });
        Assert.That(estimate.Robust, Is.EqualTo(robust));
        Assert.That(estimate.RobustFellBack, Is.Null, "a STAT1 row does not carry the flag until M5b");
        Assert.That(estimate.Weights, Is.Empty);
    }

    [Test]
    public void ASingle_TakesItsStat1Row()
    {
        var (_, outcome) = RunPart("up", "single", new[] { ("G1", All4, "fitted") }, stat1Effect: 0.4);
        Assert.That(outcome.Untested, Is.Null);
        Assert.That(outcome.Test!.Estimate.ModelUsed, Is.EqualTo(SetModel.Stat1Row));
        Assert.That(outcome.Test.Estimate.Log2Effect, Is.EqualTo(0.4));
        Assert.That(outcome.Test.Estimate.Df, Is.EqualTo(12));
    }

    [Test]
    public void ASingleWithoutATestedMember_IsNotTestable_AndOneWhoseRowFailed_IsNotRun()
    {
        Assert.That(RunPart("up", "single", new[] { ("G1", All4, "absent_in_denominator") }).Outcome.Verdict(null),
            Is.EqualTo((SetVerdictName.NotTestable, "no_tested_member")));
        Assert.That(RunPart("up", "single", new[] { ("G1", All4, "not_converged") }).Outcome.Verdict(null),
            Is.EqualTo((SetVerdictName.NotRun, "not_converged")));
    }

    [Test]
    public void ASetWithFewerThanThreeTestedGroups_IsNotTestable()
    {
        var (_, outcome) = RunPart("up", "set", new[] { ("G1", Shifted(1), "fitted"), ("G2", Shifted(1.2), "fitted"), ("G3", All4, "below_support:min_2_per_side") });
        Assert.That(outcome.Verdict(null), Is.EqualTo((SetVerdictName.NotTestable, "fewer_than_3_tested")));
    }

    [Test]
    public void ThreeMembersSharingOneGroupThroughTheException_CountAsOneTerm()
    {
        var part = Freeze(Row("H4c1"), Row("H4c2"), Row("H4c3"), Row("G1")).Parts.Single();
        var contrast = part.Contrasts.Single();
        var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
        var (table, _) = Members(("H4X", Shifted(1.0)), ("G1", Shifted(1.2)));
        var resolved = new[] { "H4c1", "H4c2", "H4c3" }
            .Select(m => new SetMemberAssignment { Member = m, ProteinGroup = "H4X", ByAllMembersException = true })
            .Append(new SetMemberAssignment { Member = "G1", ProteinGroup = "G1" }).ToList();
        var rows = new Dictionary<string, DifferentialResult> { ["H4X"] = Stat1("H4X", "fitted"), ["G1"] = Stat1("G1", "fitted") };

        var outcome = SetTest.Run(part, contrast, setup, resolved, table, rows);

        Assert.That(outcome.Members.Count(m => m.InEstimate), Is.EqualTo(4));
        Assert.That(outcome.Verdict(null), Is.EqualTo((SetVerdictName.NotTestable, "fewer_than_3_tested")));
    }

    [Test]
    public void ASetOfThree_IsFitted_AndUpCountsOneSideMembersWithoutFoldingThem()
    {
        var (_, outcome) = RunPart("up", "set", new[]
        {
            ("G1", Shifted(1.0), "fitted"), ("G2", Shifted(1.3), "fitted"), ("G3", Shifted(0.8), "fitted"),
            ("G4", All4, "absent_in_denominator"), ("G5", All4, "absent_in_numerator"), ("G6", All4, "absent_in_denominator"),
        });
        Assert.That(outcome.Untested, Is.Null);
        Assert.That(outcome.Test!.Estimate.Terms, Is.EqualTo(3));
        Assert.That(outcome.OneSideAgree, Is.EqualTo(2));
        Assert.That(outcome.OneSideDisagree, Is.EqualTo(1));
        Assert.That(outcome.BlockedBy, Is.Empty);
    }

    [TestCase("unchanged", "absent_in_denominator", true)]
    [TestCase("unchanged", "absent_in_numerator", true)]
    [TestCase("not_up", "absent_in_denominator", true)]
    [TestCase("not_up", "absent_in_numerator", false)]
    [TestCase("not_down", "absent_in_numerator", true)]
    [TestCase("not_down", "absent_in_denominator", false)]
    public void ANoChangeGuard_BlocksSupportedOnTheSidesTheContractNames(string prediction, string oneSide, bool blocks)
    {
        var (_, outcome) = RunPart(prediction, "set", new[]
        {
            ("G1", Shifted(0.05), "fitted"), ("G2", Shifted(-0.05), "fitted"), ("G3", Shifted(0.0), "fitted"), ("G4", All4, oneSide),
        });
        Assert.That(outcome.BlockedBy, Is.EqualTo(blocks ? new[] { "G4" } : Array.Empty<string>()));
        if (blocks)
            Assert.That(outcome.Verdict(0.001), Is.EqualTo((SetVerdictName.NotInformative, "blocked_by:G4")));
    }

    private static readonly (string, double[], string)[] IncompleteG3 =
    {
        ("G1", Shifted(0.05), "fitted"), ("G2", Shifted(-0.05), "fitted"), ("G3", new[] { 20, double.NaN, 20.5, 20.7 }, "fitted"),
    };

    [Test]
    public void AnIncompleteMember_StaysInTheEstimate_AndBlocksSupported()
    {
        var (_, outcome) = RunPart("unchanged", "set", IncompleteG3);
        Assert.That(outcome.Members.Single(m => m.Assignment.Member == "G3").State, Is.EqualTo(SetMemberStateFromRow.Incomplete));
        Assert.That(outcome.Test!.Estimate.Terms, Is.EqualTo(3));
        Assert.That(outcome.BlockedBy, Is.EqualTo(new[] { "G3" }));
    }

    [Test]
    public void ARobustSetWhoseReweightingNeverSettles_FallsBackToOneRound_AsR()
    {
        // GR-22 amended 2026-11-03 (was: not_run, not_converged). Two peptides per member and nearly no noise between them:
        // (1|sample:member) can take one peptide's values exactly, so reweighting drives the other peptide's weights toward
        // 0 and never settles; R's loop on the same 22 values cycles for its 100 rounds. The set falls back to msqrob2's one
        // round, is fitted and tested on it, and says so. R's one-round fit, with those weights in the lmer call
        // (make_settest_model_fixtures.R, "incomplete"), matches.
        var r = Tsv("settest_model_incomplete_robust_result.tsv").Single();
        var expected = Tsv("settest_model_incomplete_robust_weights.tsv").ToDictionary(c => (c[0], c[1], c[2]), c => D(c[3]));
        Assert.That(bool.Parse(r[8]), Is.False, "R never settles");

        var (_, outcome) = RunPart("unchanged", "set", IncompleteG3);

        Assert.That(outcome.Untested, Is.Null, "tested, not not_run");
        var estimate = outcome.Test!.Estimate;
        Assert.That(estimate.Robust, Is.True);
        Assert.That(estimate.RobustFellBack, Is.True);
        Assert.That(estimate.Weights.Select(v => (v.ProteinGroup, v.Peptide, v.Sample)), Is.EquivalentTo(expected.Keys));
        foreach (var v in estimate.Weights)
            Assert.That(v.Weight, Is.EqualTo(expected[(v.ProteinGroup, v.Peptide, v.Sample)]).Within(Fit),
                $"one-round weight {v.Peptide} {v.Sample}");
        Assert.That(estimate.Weights.Count(v => v.Weight < 1), Is.EqualTo(int.Parse(r[9])), "values down-weighted, as R");
        Close(estimate.Log2Effect, r[0], Fit, "one-round estimate");
        Close(estimate.StandardError, r[1], Fit, "one-round se");
        Close(estimate.Df, r[2], Df, "one-round Satterthwaite df");
        Close(estimate.ResidualVariance, r[3], Fit, "one-round sigma^2");
        Close(estimate.SampleVariance, r[4], Fit, "one-round var (1|sample)");
        Close(estimate.SampleMemberVariance, r[5], Fit, "one-round var (1|sample:member)");
        Assert.That(estimate.Observations, Is.EqualTo(int.Parse(r[12])));
        Assert.That(outcome.Verdict(null).Item1, Is.Not.EqualTo(SetVerdictName.NotRun), "a fallback set gets a verdict");
    }

    [Test]
    public void AttenuatedDown_TakesMbarOverTheMembersInTheEstimate()
    {
        var part = Freeze(Row("G1", "attenuated_down", mrna: "-2"), Row("G2", "attenuated_down", mrna: "-1"),
            Row("G3", "attenuated_down", mrna: "-3"), Row("G4", "attenuated_down", mrna: "-10")).Parts.Single();
        var contrast = part.Contrasts.Single();
        var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
        var (table, resolved) = Members(("G1", Shifted(-0.4)), ("G2", Shifted(-0.5)), ("G3", Shifted(-0.3)), ("G4", All4));
        var rows = new Dictionary<string, DifferentialResult>
        {
            ["G1"] = Stat1("G1", "fitted"), ["G2"] = Stat1("G2", "fitted"), ["G3"] = Stat1("G3", "fitted"),
            ["G4"] = Stat1("G4", "below_support:min_2_per_side"),
        };
        var outcome = SetTest.Run(part, contrast, setup, resolved, table, rows);
        Assert.That(outcome.Test!.MrnaMean, Is.EqualTo(-2));   // G4 (too_few) is left out
    }

    [Test]
    public void AttenuatedDown_WeighsASharedGroupOnceInMbar_AsInBetaBar()
    {
        // Alexander-Sol's case (#1572 review): three histone members share H4X through the all-members exception, each with
        // mRNA -6, beside G1 and G2 at -1. β̄ has three equally weighted groups, so m̄ is (-6 - 1 - 1) / 3 = -8/3, not the
        // member mean -4, and the default bound f·m̄ is -4/3, not -2.
        var part = Freeze(Row("H4c1", "attenuated_down", mrna: "-6"), Row("H4c2", "attenuated_down", mrna: "-6"),
            Row("H4c3", "attenuated_down", mrna: "-6"), Row("G1", "attenuated_down", mrna: "-1"),
            Row("G2", "attenuated_down", mrna: "-1")).Parts.Single();
        var contrast = part.Contrasts.Single();
        var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
        var (table, _) = Members(("H4X", Shifted(-0.4)), ("G1", Shifted(-0.5)), ("G2", Shifted(-0.3)));
        var resolved = new[] { "H4c1", "H4c2", "H4c3" }
            .Select(m => new SetMemberAssignment { Member = m, ProteinGroup = "H4X", ByAllMembersException = true })
            .Append(new SetMemberAssignment { Member = "G1", ProteinGroup = "G1" })
            .Append(new SetMemberAssignment { Member = "G2", ProteinGroup = "G2" }).ToList();
        var rows = new Dictionary<string, DifferentialResult>
        {
            ["H4X"] = Stat1("H4X", "fitted"), ["G1"] = Stat1("G1", "fitted"), ["G2"] = Stat1("G2", "fitted"),
        };

        var outcome = SetTest.Run(part, contrast, setup, resolved, table, rows);

        Assert.That(outcome.Untested, Is.Null);
        Assert.That(outcome.Test!.Estimate.Terms, Is.EqualTo(3));
        Assert.That(outcome.Test.MrnaMean!.Value, Is.EqualTo(-8.0 / 3).Within(1e-12));
        Assert.That(outcome.Test.Bound, Is.EqualTo("0;" + (0.5 * (-8.0 / 3)).ToString("R", CultureInfo.InvariantCulture)));
    }

    [Test]
    public void AttenuatedDownWhoseMrnaRises_IsNotRun()
    {
        var (_, outcome) = RunPart("attenuated_down", "set", new[]
        {
            ("G1", Shifted(-0.4), "fitted"), ("G2", Shifted(-0.5), "fitted"), ("G3", Shifted(-0.3), "fitted"),
        }, mrna: "1.5");
        Assert.That(outcome.Verdict(null), Is.EqualTo((SetVerdictName.NotRun, "mrna_not_usable:direction")));
    }

    [Test]
    public void TheVerdictWords_ComeInTheContractsOrder()
    {
        var (_, up) = RunPart("up", "set", new[] { ("G1", Shifted(1.0), "fitted"), ("G2", Shifted(1.3), "fitted"), ("G3", Shifted(0.8), "fitted") });
        Assert.That(up.Verdict(0.01).Verdict, Is.EqualTo(SetVerdictName.Supported));
        Assert.That(up.Verdict(0.2), Is.EqualTo((SetVerdictName.NotInformative, "bh_adjusted_p>0.05")));
        Assert.That(up.Verdict(null), Is.EqualTo((SetVerdictName.NotInformative, "not_in_family")));
        // A clear rise, far outside the interval's width on four samples, contradicts down.
        var (_, down) = RunPart("down", "set", new[] { ("G1", Shifted(3.0), "fitted"), ("G2", Shifted(3.3), "fitted"), ("G3", Shifted(2.8), "fitted") });
        Assert.That(down.Test!.Ci90Low, Is.GreaterThan(0));
        Assert.That(down.Test.Contradicted, Is.True);
        Assert.That(down.Verdict(0.01), Is.EqualTo((SetVerdictName.Contradicted, "uncorrected")));
    }

    [Test]
    public void ANotRunSetup_IsNotRunWithItsReason()
    {
        var part = Freeze(Row("G1", contexts: "diet=chow")).Parts.Single();
        var setup = SetContrastSetup.Build(part.Contrasts.Single(), FourSamples, AgeSpec);
        var (table, resolved) = Members(("G1", All4));
        var outcome = SetTest.Run(part, part.Contrasts.Single(), setup, resolved, table, new Dictionary<string, DifferentialResult>());
        Assert.That(outcome.Verdict(null), Is.EqualTo((SetVerdictName.NotRun, "not_in_design:context_factor=diet")));
    }

    [Test]
    public void MarginsThatDifferBetweenMemberRows_AreRefused()
    {
        string header = Header + "\tmargin_log2";
        var part = SetTestFreeze.Read(Encoding.UTF8.GetBytes(header + "\n" + Row("G1", "unchanged") + "\t0.3\n" +
            Row("G2", "unchanged") + "\t0.4\n" + Row("G3", "unchanged") + "\t0.3\n"), "t").Parts.Single();
        var setup = SetContrastSetup.Build(part.Contrasts.Single(), FourSamples, AgeSpec);
        var (table, resolved) = Members(("G1", All4), ("G2", All4), ("G3", All4));
        var rows = resolved.ToDictionary(r => r.ProteinGroup!, r => Stat1(r.ProteinGroup!, "fitted"));
        Assert.That(() => SetTest.Run(part, part.Contrasts.Single(), setup, resolved, table, rows),
            Throws.ArgumentException.With.Message.Contains("margin_log2"));
    }

    [Test]
    public void Down_CountsOneSideMembersTheOtherWayRound()
    {
        var (_, outcome) = RunPart("down", "set", new[]
        {
            ("G1", Shifted(-1.0), "fitted"), ("G2", Shifted(-1.3), "fitted"), ("G3", Shifted(-0.8), "fitted"),
            ("G4", All4, "absent_in_numerator"), ("G5", All4, "absent_in_denominator"), ("G6", All4, "absent_in_numerator"),
        });
        Assert.That(outcome.Untested, Is.Null);
        Assert.That(outcome.Test!.Estimate.Terms, Is.EqualTo(3));
        // Seen only at the reference level is a fall, so it agrees with down.
        Assert.That(outcome.OneSideAgree, Is.EqualTo(2));
        Assert.That(outcome.OneSideDisagree, Is.EqualTo(1));
        Assert.That(outcome.BlockedBy, Is.Empty);
    }

    [Test]
    public void AMemberWithoutAGroup_KeepsItsState_AndStaysOutOfTheEstimate()
    {
        var part = Freeze(Row("G1"), Row("G2"), Row("G3"), Row("Gx")).Parts.Single();
        var contrast = part.Contrasts.Single();
        var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
        var (table, resolved) = Members(("G1", Shifted(1.0)), ("G2", Shifted(1.3)), ("G3", Shifted(0.8)));
        var all = resolved.Append(new SetMemberAssignment { Member = "Gx", State = SetMemberState.NoGroup }).ToList();
        var rows = resolved.ToDictionary(r => r.ProteinGroup!, r => Stat1(r.ProteinGroup!, "fitted"));

        var outcome = SetTest.Run(part, contrast, setup, all, table, rows);

        var gx = outcome.Members.Single(m => m.Assignment.Member == "Gx");
        Assert.That(gx.State, Is.EqualTo(SetMemberState.NoGroup));
        Assert.That(gx.Stat1Row, Is.Null);
        Assert.That(gx.InEstimate, Is.False);
        Assert.That(outcome.Untested, Is.Null);
        Assert.That(outcome.Test!.Estimate.Terms, Is.EqualTo(3));
    }

    [Test]
    public void ASetModelThatCannotBeFormed_IsNotRunNamingWhy()
    {
        // Three one-peptide groups on four samples: (1|sample:member) has one level per value, so it cannot be told from
        // the residual (too_few_groups, as lme4 refuses it), and the contrast is not run.
        var part = Freeze(Row("G1"), Row("G2"), Row("G3")).Parts.Single();
        var contrast = part.Contrasts.Single();
        var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
        var groups = new[] { (Group: "G1", Values: Shifted(1.0)), (Group: "G2", Values: Shifted(1.3)), (Group: "G3", Values: Shifted(0.8)) };
        var table = Table(Four, groups.SelectMany(g => Four.Select((s, i) => (g.Group, g.Group + "PEPA", s, g.Values[i]))));
        var resolved = groups.Select(g => new SetMemberAssignment { Member = g.Group, ProteinGroup = g.Group }).ToList();
        var rows = groups.ToDictionary(g => g.Group, g => Stat1(g.Group, "fitted"));

        Assert.That(SetModel.FitSet(setup, groups.Select(g => g.Group).ToList(), table, out string? notRun), Is.Null);
        Assert.That(notRun, Is.EqualTo("not_converged:too_few_groups"));
        Assert.That(SetTest.Run(part, contrast, setup, resolved, table, rows).Verdict(null),
            Is.EqualTo((SetVerdictName.NotRun, "not_converged:too_few_groups")));
    }

    [Test]
    public void AMemberInTheEstimateWithoutMrna_MakesAnAttenuatedPartNotRun_AndOneOutsideItDoesNot()
    {
        SetContrastOutcome Run(string g4Status)
        {
            var part = Freeze(Row("G1", "attenuated_down", mrna: "-2"), Row("G2", "attenuated_down", mrna: "-1"),
                Row("G3", "attenuated_down", mrna: "-3"), Row("G4", "attenuated_down", mrna: "")).Parts.Single();
            var contrast = part.Contrasts.Single();
            var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
            var (table, resolved) = Members(("G1", Shifted(-0.4)), ("G2", Shifted(-0.5)), ("G3", Shifted(-0.3)), ("G4", Shifted(-0.2)));
            var rows = new Dictionary<string, DifferentialResult>
            {
                ["G1"] = Stat1("G1", "fitted"), ["G2"] = Stat1("G2", "fitted"), ["G3"] = Stat1("G3", "fitted"), ["G4"] = Stat1("G4", g4Status),
            };
            return SetTest.Run(part, contrast, setup, resolved, table, rows);
        }

        Assert.That(Run("fitted").Verdict(null), Is.EqualTo((SetVerdictName.NotRun, "mrna_not_usable:missing")));
        var withoutG4 = Run("below_support:min_2_per_side");
        Assert.That(withoutG4.Untested, Is.Null);
        Assert.That(withoutG4.Test!.MrnaMean, Is.EqualTo(-2));
    }

    [Test]
    public void AnAttenuatedSetWithNoTestedMember_IsNotTestable_WithoutNeedingMrna()
    {
        var (_, outcome) = RunPart("attenuated_down", "set", new[]
        {
            ("G1", All4, "below_support:min_2_per_side"), ("G2", All4, "absent_in_both"), ("G3", All4, "below_support:min_2_per_side"),
        });
        Assert.That(outcome.Verdict(null), Is.EqualTo((SetVerdictName.NotTestable, "fewer_than_3_tested")));
    }

    [Test]
    public void AttenuatedUpWhoseMrnaFalls_IsNotRun()
    {
        var (_, outcome) = RunPart("attenuated_up", "set", new[]
        {
            ("G1", Shifted(0.4), "fitted"), ("G2", Shifted(0.5), "fitted"), ("G3", Shifted(0.3), "fitted"),
        }, mrna: "-1.5");
        Assert.That(outcome.Verdict(null), Is.EqualTo((SetVerdictName.NotRun, "mrna_not_usable:direction")));
    }

    [Test]
    public void ASingleResolvedToTwoMembers_IsRefused()
    {
        var part = Freeze(Row("G1", kind: "single")).Parts.Single();
        var contrast = part.Contrasts.Single();
        var setup = SetContrastSetup.Build(contrast, FourSamples, AgeSpec);
        var (table, resolved) = Members(("G1", All4), ("G2", All4));
        var rows = resolved.ToDictionary(r => r.ProteinGroup!, r => Stat1(r.ProteinGroup!, "fitted", 0.4, 0.1, 12));
        Assert.That(() => SetTest.Run(part, contrast, setup, resolved, table, rows),
            Throws.ArgumentException.With.Message.Contains("tested as a single but lists 2 members"));
    }
}
