using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// What <see cref="AnalysisDesign"/> promises beyond the matrix values: column names, level order, default and explicit
/// contrasts, estimability and confounding, strata, refusals, and determinism (STAT1 M1; GR-1, GR-5, GR-14;
/// <c>DEF-DIFF-STRATUM</c>, <c>DEF-DIFF-STATUS</c>).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class AnalysisDesignTests
{
    private static AnalysisSample S(string id, params (string Factor, string Level)[] factors) =>
        new(id) { Factors = factors.ToDictionary(f => f.Factor, f => f.Level) };

    private static DesignSpec Spec(params (string Name, string? Reference)[] factors) =>
        new() { Factors = factors.Select(f => new FactorSpec(f.Name, f.Reference)).ToArray() };

    private static AnalysisDesign FromFixture(string name, DesignSpec spec) =>
        AnalysisDesign.Create(AnalysisDesignReferenceTests.ReadSamples(name).Samples, spec);

    private static readonly DesignSpec TwoFactors = Spec(("age_group", "young"), ("treatment", "vehicle"));

    [Test]
    public void ColumnsAreNamedFactorEqualsLevel()
    {
        var design = FromFixture("interaction", TwoFactors with { Interactions = new[] { ("age_group", "treatment") } });
        Assert.That(design.ColumnNames, Is.EqualTo(new[]
            { "intercept", "age_group=old", "treatment=drug", "age_group=old:treatment=drug" }));
    }

    [Test]
    public void BatchAndCovariateColumnsAreNamed()
    {
        var spec = Spec(("condition", "ctrl")) with
        {
            Covariates = new[] { new CovariateSpec("age", "years", Scale: 10, ScaledUnit: "decade", Centre: 50) }
        };
        var design = FromFixture("factor_batch_covariate", spec);
        Assert.That(design.ColumnNames, Is.EqualTo(new[] { "intercept", "condition=treated", "batch=b2", "age (per decade)" }));
    }

    [Test]
    public void ReferenceComesFirstThenOrdinalOrder()
    {
        var design = FromFixture("nonalpha_reference", Spec(("age_group", "young")));
        Assert.That(design.Levels("age_group"), Is.EqualTo(new[] { "young", "middle", "old" }));
    }

    [Test]
    public void IndividualIsCarriedButNotAFixedColumn()
    {
        var design = FromFixture("paired", Spec(("time", "pre")));
        Assert.That(design.ColumnNames, Is.EqualTo(new[] { "intercept", "time=post" }));
        Assert.That(design.Samples.Select(s => s.Individual).Distinct().Count(), Is.EqualTo(4));
    }

    [Test]
    public void DefaultContrastsAreEachLevelAgainstItsReference()
    {
        var contrasts = FromFixture("two_factors", TwoFactors).DefaultContrasts();
        Assert.That(contrasts.Select(c => c.Label), Is.EqualTo(new[]
            { "age_group=old vs age_group=young", "treatment=drug vs treatment=vehicle" }));
        Assert.That(contrasts.Select(c => c.Id), Is.EqualTo(new[] { "c1", "c2" }));
        Assert.That(contrasts[0].Weights, Is.EqualTo(new[] { 0.0, 1, 0 }));
        Assert.That(contrasts[1].Weights, Is.EqualTo(new[] { 0.0, 0, 1 }));
        Assert.That(contrasts[0].Numerator, Is.EqualTo("age_group=old"));
        Assert.That(contrasts[0].Denominator, Is.EqualTo("age_group=young"));
    }

    [Test]
    public void ThreeLevelsGiveTwoDefaultContrasts()
    {
        var contrasts = FromFixture("lost_sample", Spec(("condition", "ctrl"))).DefaultContrasts();
        Assert.That(contrasts.Select(c => c.Label), Is.EqualTo(new[]
            { "condition=x vs condition=ctrl", "condition=y vs condition=ctrl" }));
    }

    [Test]
    public void NoReferenceComparesEveryPairAndWarns()
    {
        var design = FromFixture("no_reference", Spec(("condition", null)));
        Assert.That(design.Warnings, Has.Some.Contains("No reference level").And.Some.Contains("condition"));
        var contrasts = design.DefaultContrasts();
        Assert.That(contrasts.Select(c => c.Label), Is.EqualTo(new[]
            { "condition=B vs condition=A", "condition=C vs condition=A", "condition=C vs condition=B" }));
        Assert.That(contrasts[2].Weights, Is.EqualTo(new[] { 0.0, -1, 1 }));
    }

    [Test]
    public void BetweenTwoNonReferenceLevels()
    {
        var contrast = FromFixture("lost_sample", Spec(("condition", "ctrl"))).Between("condition", "y", "x");
        Assert.That(contrast.Weights, Is.EqualTo(new[] { 0.0, -1, 1 }));
        Assert.That(contrast.Label, Is.EqualTo("condition=y vs condition=x"));
    }

    [Test]
    public void BetweenAgainstTheReferenceIsOneCoefficient()
    {
        var contrast = FromFixture("lost_sample", Spec(("condition", "ctrl"))).Between("condition", "x", "ctrl");
        Assert.That(contrast.Weights, Is.EqualTo(new[] { 0.0, 1, 0 }));
    }

    [Test]
    public void BetweenAnUnknownLevelThrowsAndNamesTheLevels()
    {
        var design = FromFixture("lost_sample", Spec(("condition", "ctrl")));
        var e = Assert.Throws<ArgumentException>(() => design.Between("condition", "z", "ctrl"));
        Assert.That(e!.Message, Does.Contain("z").And.Contain("ctrl, x, y"));
    }

    [Test]
    public void SlopeIsPerScaledUnitAndStatesItsScale()
    {
        var spec = new DesignSpec
        {
            Covariates = new[] { new CovariateSpec("age", "years", Scale: 10, ScaledUnit: "decade", Centre: 50) }
        };
        var design = FromFixture("age_slope", spec);
        var slope = design.Slope("age");
        Assert.That(slope.Weights, Is.EqualTo(new[] { 0.0, 1 }));
        Assert.That(slope.Covariate, Is.EqualTo("age"));
        Assert.That(slope.CovariateUnit, Is.EqualTo("decade"));
        Assert.That(slope.CovariateScale, Is.EqualTo("years/10, centred at 50 years"));
        Assert.That(slope.Numerator, Is.Null);
        Assert.That(design.DefaultContrasts().Single().Label, Is.EqualTo("age (per decade)"));
    }

    [Test]
    public void AnUnstatedCentreIsTheStratumMean()
    {
        var spec = Spec(("sex", "female")) with
        {
            Covariates = new[] { new CovariateSpec("age", "years", Scale: 10, ScaledUnit: "decade") }
        };
        Assert.That(FromFixture("age_default_centre", spec).CovariateCentre("age"), Is.EqualTo(14.25).Within(1e-12));
    }

    [Test]
    public void CustomContrastByColumnName()
    {
        var design = FromFixture("interaction", TwoFactors with { Interactions = new[] { ("age_group", "treatment") } });
        var c = design.Custom("i1", "drug effect in old", new Dictionary<string, double>
            { ["treatment=drug"] = 1, ["age_group=old:treatment=drug"] = 1 });
        Assert.That(c.Weights, Is.EqualTo(new[] { 0.0, 0, 1, 1 }));
        Assert.Throws<ArgumentException>(() => design.Custom("i2", "x", new Dictionary<string, double> { ["nope"] = 1 }));
    }

    [Test]
    public void ConfoundedWithBatchIsNotEstimable()
    {
        var design = FromFixture("confounded", Spec(("condition", "A")));
        Assert.That(design.Rank, Is.EqualTo(2));
        Assert.That(design.AliasedTerms, Is.EquivalentTo(new[] { "condition", "batch" }));
        Assert.That(design.Estimability(design.DefaultContrasts().Single()), Is.EqualTo("not_estimable:confounded"));
    }

    [Test]
    public void BalancedOverBatchIsEstimable()
    {
        var design = FromFixture("batch", Spec(("condition", "A")));
        Assert.That(design.AliasedTerms, Is.Empty);
        Assert.That(design.Estimability(design.DefaultContrasts().Single()), Is.EqualTo("estimable"));
    }

    [Test]
    public void AnEmptyInteractionCellIsRankDeficientButMainEffectsStayEstimable()
    {
        // No sample is old AND on drug, so the interaction column is all zero.
        var samples = new[]
        {
            S("s1", ("age_group", "young"), ("treatment", "vehicle")), S("s2", ("age_group", "young"), ("treatment", "drug")),
            S("s3", ("age_group", "old"), ("treatment", "vehicle")), S("s4", ("age_group", "young"), ("treatment", "vehicle")),
            S("s5", ("age_group", "young"), ("treatment", "drug")), S("s6", ("age_group", "old"), ("treatment", "vehicle")),
        };
        var design = AnalysisDesign.Create(samples, TwoFactors with { Interactions = new[] { ("age_group", "treatment") } });
        var interaction = design.Custom("i", "interaction", new Dictionary<string, double> { ["age_group=old:treatment=drug"] = 1 });
        Assert.That(design.Estimability(interaction), Is.EqualTo("not_estimable:rank_deficient"));
        Assert.That(design.DefaultContrasts().Select(design.Estimability), Is.All.EqualTo("estimable"));
    }

    [Test]
    public void AFactorWithOneLevelInTheStratumIsDroppedAndSaysSo()
    {
        var samples = new[]
        {
            S("s1", ("condition", "A"), ("sex", "female")), S("s2", ("condition", "B"), ("sex", "female")),
            S("s3", ("condition", "A"), ("sex", "female")), S("s4", ("condition", "B"), ("sex", "female")),
        };
        var design = AnalysisDesign.Create(samples, Spec(("condition", "A"), ("sex", "female")));
        Assert.That(design.ColumnNames, Is.EqualTo(new[] { "intercept", "condition=B" }));
        Assert.That(design.Warnings, Has.Some.Contains("sex"));
        Assert.That(design.NotRun.Select(n => n.Reason), Has.Some.EqualTo("not_in_stratum"));
    }

    [Test]
    public void AReferenceWithNoSampleInTheStratumRunsNoContrastAgainstIt()
    {
        var samples = new[]
        {
            S("s1", ("age_group", "old")), S("s2", ("age_group", "middle")),
            S("s3", ("age_group", "old")), S("s4", ("age_group", "middle")),
        };
        var design = AnalysisDesign.Create(samples, Spec(("age_group", "young")));
        Assert.That(design.DefaultContrasts().Where(c => c.Denominator == "age_group=young"), Is.Empty);
        Assert.That(design.NotRun, Has.Some.Matches<NotRunContrast>(n =>
            n.Label.Contains("age_group=young") && n.Reason == "not_in_stratum"));
        Assert.That(design.Warnings, Has.Some.Contains("young"));
    }

    [Test]
    public void SamplesAreSplitByStratum()
    {
        AnalysisSample T(string id, string tissue, string condition) =>
            S(id, ("condition", condition)) with { Strata = new Dictionary<string, string> { ["organism part"] = tissue } };
        var samples = new[]
        {
            T("s1", "liver", "A"), T("s2", "muscle", "A"), T("s3", "liver", "B"),
            T("s4", "muscle", "B"), T("s5", "liver", "A"), T("s6", "muscle", "B"),
        };
        var designs = AnalysisDesign.ByStratum(samples, Spec(("condition", "A")));
        Assert.That(designs.Select(d => d.Stratum), Is.EqualTo(new[] { "organism part=liver", "organism part=muscle" }));
        Assert.That(designs[0].Samples.Select(s => s.SampleId), Is.EqualTo(new[] { "s1", "s3", "s5" }));
    }

    [Test]
    public void StratumKeyIsSortedEscapedAndAllWhenEmpty()
    {
        Assert.That(AnalysisDesign.StratumKey(new Dictionary<string, string>
            { ["organism part"] = "liver", ["cell type"] = "hep;a=b%" }),
            Is.EqualTo("cell type=hep%3Ba%3Db%25;organism part=liver"));
        Assert.That(AnalysisDesign.StratumKey(new Dictionary<string, string>()), Is.EqualTo("all"));
    }

    [Test]
    public void ResultDoesNotDependOnInputOrder()
    {
        var samples = AnalysisDesignReferenceTests.ReadSamples("interaction").Samples;
        var spec = TwoFactors with { Interactions = new[] { ("age_group", "treatment") } };
        var a = AnalysisDesign.Create(samples, spec).Matrix;
        var b = AnalysisDesign.Create(Enumerable.Reverse(samples).ToList(), spec).Matrix;
        Assert.That(b, Is.EqualTo(a));
    }

    [Test]
    public void MatrixIsACopy()
    {
        var design = FromFixture("balanced", Spec(("condition", "A")));
        design.Matrix[0, 0] = 99;
        Assert.That(design.Matrix[0, 0], Is.EqualTo(1));
    }

    [Test]
    public void ADuplicateSampleIsRefused()
    {
        var e = Assert.Throws<ArgumentException>(() => AnalysisDesign.Create(
            new[] { S("s1", ("condition", "A")), S("s1", ("condition", "B")) }, Spec(("condition", "A"))));
        Assert.That(e!.Message, Does.Contain("s1"));
    }

    [Test]
    public void AMissingFactorValueIsRefusedAndNamesTheSample()
    {
        var e = Assert.Throws<ArgumentException>(() => AnalysisDesign.Create(
            new[] { S("s1", ("condition", "A")), S("s2") }, Spec(("condition", "A"))));
        Assert.That(e!.Message, Does.Contain("s2").And.Contain("condition"));
    }

    [Test]
    public void ANonFiniteCovariateIsRefused()
    {
        var samples = new[]
        {
            new AnalysisSample("s1") { Covariates = new Dictionary<string, double> { ["age"] = 30 } },
            new AnalysisSample("s2") { Covariates = new Dictionary<string, double> { ["age"] = double.NaN } },
        };
        var e = Assert.Throws<ArgumentException>(() => AnalysisDesign.Create(samples,
            new DesignSpec { Covariates = new[] { new CovariateSpec("age", "years") } }));
        Assert.That(e!.Message, Does.Contain("s2").And.Contain("age"));
    }

    [Test]
    public void AnInteractionWithAnUnknownFactorIsRefused()
    {
        Assert.Throws<ArgumentException>(() => AnalysisDesign.Create(
            new[] { S("s1", ("condition", "A")), S("s2", ("condition", "B")) },
            Spec(("condition", "A")) with { Interactions = new[] { ("condition", "sex") } }));
    }

    [Test]
    public void BatchOnSomeSamplesOnlyIsRefused()
    {
        var samples = new[] { S("s1", ("condition", "A")) with { Batch = "1" }, S("s2", ("condition", "B")) };
        var e = Assert.Throws<ArgumentException>(() => AnalysisDesign.Create(samples, Spec(("condition", "A"))));
        Assert.That(e!.Message, Does.Contain("s2").And.Contain("batch"));
    }
}
