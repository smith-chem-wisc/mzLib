using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;
using StatisticalModels;

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

    /// <summary>
    /// With no scaled unit named, a slope is per <c>Scale</c> stored units, and its unit says so: <c>Scale: 10</c> on
    /// years is an effect per 10 years, not per year. Written in a comma-decimal culture, so the number is invariant.
    /// </summary>
    [TestCase(10, "10 years", "years/10, centred at 50 years")]
    [TestCase(2.5, "2.5 years", "years/2.5, centred at 50 years")]
    [TestCase(1, "years", "years/1, centred at 50 years")]
    [SetCulture("de-DE")]
    public void WithoutAScaledUnitTheSlopeIsPerScaleStoredUnits(double scale, string unit, string scaleText)
    {
        var spec = new DesignSpec { Covariates = new[] { new CovariateSpec("age", "years", Scale: scale, Centre: 50) } };
        var design = FromFixture("age_slope", spec);
        var slope = design.Slope("age");
        Assert.That(design.ColumnNames.Last(), Is.EqualTo($"age (per {unit})"));
        Assert.That(slope.Label, Is.EqualTo($"age (per {unit})"));
        Assert.That(slope.CovariateUnit, Is.EqualTo(unit));
        Assert.That(slope.CovariateScale, Is.EqualTo(scaleText));
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
    public void ANearlyRedundantCovariateIsNotFullRankBecauseTheFitRefusesIt()
    {
        // b is a plus or minus 1e-8: its part orthogonal to the intercept and a is about 4e-9 of its own length, under
        // LinearModel's 1e-7, so LinearModel.Fit refuses the design; but the smallest singular value is about 2e-9 of the
        // largest, over the old 1e-10 cut-off, which called every contrast estimable.
        var samples = Enumerable.Range(0, 8).Select(i => new AnalysisSample($"s{i}")
        {
            Covariates = new Dictionary<string, double> { ["a"] = i, ["b"] = i + (i % 2 == 0 ? 1e-8 : -1e-8) },
        }).ToList();
        var design = AnalysisDesign.Create(samples, new DesignSpec
        {
            Covariates = new[] { new CovariateSpec("a", "units"), new CovariateSpec("b", "units") },
        });

        Assert.That(LinearModel.HasFullColumnRank(design.Matrix), Is.False);
        Assert.Throws<ArgumentException>(() => LinearModel.Fit(new double[1, 8], design.Matrix));
        Assert.That(design.Rank, Is.EqualTo(2));
        Assert.That(design.AliasedTerms, Is.EquivalentTo(new[] { "a", "b" }));
        Assert.That(design.Estimability(design.Slope("a")), Is.EqualTo("not_estimable:confounded"));
    }

    [Test]
    public void ACovariateOnAnIntensityScaleBesideTheInterceptIsFullRank()
    {
        // Centred values near ±1e11 beside an intercept of length sqrt(12): the smallest singular value is about 3e-11
        // of the largest, so the old cut-off called the slope not estimable, while LinearModel's test, relative to each
        // column's own length, fits it.
        var samples = Enumerable.Range(0, 12).Select(i => new AnalysisSample($"s{i:00}")
        {
            Covariates = new Dictionary<string, double> { ["summed_intensity"] = 1e11 * (1 + 0.1 * i) },
        }).ToList();
        var design = AnalysisDesign.Create(samples, new DesignSpec
        {
            Covariates = new[] { new CovariateSpec("summed_intensity", "intensity") },
        });

        Assert.That(LinearModel.HasFullColumnRank(design.Matrix), Is.True);
        Assert.That(design.Rank, Is.EqualTo(2));
        Assert.That(design.AliasedTerms, Is.Empty);
        Assert.That(design.Estimability(design.Slope("summed_intensity")), Is.EqualTo("estimable"));
    }

    /// <summary>
    /// Review round 4 (@pcruzparri): beside one real redundancy (a covariate that does not vary) the rank was sized by an
    /// SVD cut-off at 1e-10 of the largest singular value, which the intensity-scale covariate falls under, so the
    /// condition was called not estimable while <see cref="AnalysisDesign.FittedColumns"/> kept it. One rule now decides
    /// both.
    /// </summary>
    [Test]
    public void AnIntensityScaleCovariateBesideAConstantOneLeavesTheConditionEstimable()
    {
        var samples = Enumerable.Range(0, 12).Select(i => new AnalysisSample($"s{i:00}")
        {
            Factors = new Dictionary<string, string> { ["condition"] = i % 2 == 0 ? "A" : "B" },
            Covariates = new Dictionary<string, double> { ["summed_intensity"] = 1e11 * (1 + 0.1 * i), ["bmi"] = 22.5 },
        }).ToList();
        var design = AnalysisDesign.Create(samples, Spec(("condition", "A")) with
        {
            Covariates = new[] { new CovariateSpec("summed_intensity", "intensity"), new CovariateSpec("bmi", "kg/m2") },
        });

        Assert.That(design.FittedColumns.Select(j => design.ColumnNames[j]),
            Is.EqualTo(new[] { "intercept", "condition=B", "summed_intensity (per intensity)" }));
        Assert.That(design.Rank, Is.EqualTo(design.FittedColumns.Count));
        Assert.That(design.AliasedTerms, Is.Empty, "the constant covariate is redundant on its own");
        var contrast = design.Between("condition", "B", "A");
        Assert.That(design.Estimability(contrast), Is.EqualTo(AnalysisDesign.Estimable));
        Assert.That(design.Estimability(design.Slope("summed_intensity")), Is.EqualTo(AnalysisDesign.Estimable));
        Assert.That(design.Estimability(design.Slope("bmi")), Is.EqualTo(AnalysisDesign.NotEstimableRankDeficient));
        Assert.That(design.FittedWeights(contrast), Is.EqualTo(new[] { 0.0, 1.0, 0.0 }));
    }

    /// <summary>
    /// Review round 4 (@pcruzparri): the other side of the same disagreement. Two covariates equal to within 1e-8, which
    /// the fit cannot separate, were called estimable once a constant covariate was added beside them, and the slope of
    /// the left-out one got weights that select nothing.
    /// </summary>
    [Test]
    public void NearlyEqualCovariatesBesideAConstantOneStayNotEstimable()
    {
        var samples = Enumerable.Range(0, 8).Select(i => new AnalysisSample($"s{i}")
        {
            Covariates = new Dictionary<string, double>
                { ["a"] = i, ["b"] = i + (i % 2 == 0 ? 1e-8 : -1e-8), ["bmi"] = 22.5 },
        }).ToList();
        var design = AnalysisDesign.Create(samples, new DesignSpec
        {
            Covariates = new[] { new CovariateSpec("a", "units"), new CovariateSpec("b", "units"), new CovariateSpec("bmi", "kg/m2") },
        });

        Assert.That(design.FittedColumns.Select(j => design.ColumnNames[j]), Is.EqualTo(new[] { "intercept", "a (per units)" }));
        Assert.That(design.Rank, Is.EqualTo(2));
        Assert.That(design.AliasedTerms, Is.EquivalentTo(new[] { "a", "b" }));
        Assert.That(design.Estimability(design.Slope("a")), Is.EqualTo(AnalysisDesign.NotEstimableConfounded));
        Assert.That(design.Estimability(design.Slope("b")), Is.EqualTo(AnalysisDesign.NotEstimableConfounded));
        Assert.Throws<ArgumentException>(() => design.FittedWeights(design.Slope("b")));
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

    /// <summary>
    /// Ten (or six) samples, half condition A in batch b1 and half B in b2, all of one age: condition is confounded with
    /// batch, and age does not vary. <paramref name="centre"/> is the stated centre, or none.
    /// </summary>
    private static AnalysisDesign ConstantAgeBesideAConfoundedCondition(double age, int n, double? centre = null) =>
        AnalysisDesign.Create(Enumerable.Range(0, n).Select(i => new AnalysisSample($"s{i:D2}")
        {
            Factors = new Dictionary<string, string> { ["condition"] = i < n / 2 ? "A" : "B" },
            Batch = i < n / 2 ? "b1" : "b2",
            Covariates = new Dictionary<string, double> { ["age"] = age },
        }).ToList(), Spec(("condition", "A")) with { Covariates = new[] { new CovariateSpec("age", "months", Centre: centre) } });

    /// <summary>
    /// Review round 3 (@pcruzparri): with no stated centre a constant covariate is centred at its mean, and the mean of
    /// ten 18.2s is not exactly 18.2, so its column was 3.6e-15, not 0. Beside condition confounded with batch, the slope
    /// came back confounded and age was named among the aliased terms. A covariate that does not vary is rank deficient
    /// on its own, and its centre is that value.
    /// </summary>
    [TestCase(18.2, 10)]
    [TestCase(24.0, 10)]
    [TestCase(0.1, 6)]
    [TestCase(0.1, 10)]
    [TestCase(0.7, 6)]
    [TestCase(0.7, 10)]
    [TestCase(1.1, 6)]
    [TestCase(1.1, 10)]
    [TestCase(0.3, 10)]
    public void AConstantCovariateBesideAConfoundedConditionIsRankDeficient(double age, int n)
    {
        var design = ConstantAgeBesideAConfoundedCondition(age, n);

        Assert.That(design.Estimability(design.Slope("age")), Is.EqualTo(AnalysisDesign.NotEstimableRankDeficient));
        Assert.That(design.Estimability(design.Between("condition", "B", "A")),
            Is.EqualTo(AnalysisDesign.NotEstimableConfounded));
        Assert.That(design.AliasedTerms, Is.EqualTo(new[] { "condition", "batch" }));
        Assert.That(design.CovariateCentre("age"), Is.EqualTo(age), "the mean of equal values is that value");
    }

    /// <summary>
    /// The reason must not depend on which basis of a two-dimensional null space the SVD returns. A stated centre a
    /// rounding error off the values (the floating-point mean of ten 18.2s) leaves the same 3.6e-15 column whatever the
    /// default centre does; the reason and the aliased terms come from the projection onto the null space.
    /// </summary>
    [Test]
    public void TheReasonDoesNotDependOnTheNullSpaceBasis()
    {
        double offCentre = Enumerable.Repeat(18.2, 10).Average();
        Assert.That(offCentre, Is.Not.EqualTo(18.2), "the fixture needs the rounding error");

        var design = ConstantAgeBesideAConfoundedCondition(18.2, 10, offCentre);

        Assert.That(design.Estimability(design.Slope("age")), Is.EqualTo(AnalysisDesign.NotEstimableRankDeficient));
        Assert.That(design.AliasedTerms, Is.EqualTo(new[] { "condition", "batch" }));
    }

    /// <summary>A design of full column rank is fitted as it is.</summary>
    [Test]
    public void AFullRankDesignIsFittedOnEveryColumn()
    {
        var design = FromFixture("interaction", TwoFactors with { Interactions = new[] { ("age_group", "treatment") } });

        Assert.That(design.FittedColumns, Is.EqualTo(Enumerable.Range(0, design.ColumnNames.Count)));
        Assert.That(design.FittedMatrix, Is.EqualTo(design.Matrix));
    }

    /// <summary>
    /// Review round 3 (@pcruzparri): a rank-deficient design called some contrasts estimable, but <c>LinearModel.Fit</c>
    /// refuses its <see cref="AnalysisDesign.Matrix"/>. The design now leaves the later redundant columns out of
    /// <see cref="AnalysisDesign.FittedMatrix"/>, and gives each estimable contrast its weights on the kept columns. Its
    /// estimate equals the one from the design without the redundant term, built independently.
    /// </summary>
    [Test]
    public void AnEstimableContrastOfARankDeficientDesignIsFittedOnTheKeptColumns()
    {
        // No sample is old AND on drug, so the interaction column is all zero.
        var samples = new[]
        {
            S("s1", ("age_group", "young"), ("treatment", "vehicle")), S("s2", ("age_group", "young"), ("treatment", "drug")),
            S("s3", ("age_group", "old"), ("treatment", "vehicle")), S("s4", ("age_group", "young"), ("treatment", "vehicle")),
            S("s5", ("age_group", "young"), ("treatment", "drug")), S("s6", ("age_group", "old"), ("treatment", "vehicle")),
        };
        var design = AnalysisDesign.Create(samples, TwoFactors with { Interactions = new[] { ("age_group", "treatment") } });
        var additive = AnalysisDesign.Create(samples, TwoFactors);
        Assert.Throws<ArgumentException>(() => LinearModel.Fit(new double[1, 6], design.Matrix), "the full matrix is refused");

        Assert.That(design.FittedColumns.Select(j => design.ColumnNames[j]), Is.EqualTo(additive.ColumnNames));
        var y = new double[,] { { 20.1, 21.7, 19.4, 20.6, 22.3, 18.9 } };
        var fit = LinearModel.Fit(y, design.FittedMatrix);
        var reference = LinearModel.Fit(y, additive.Matrix);
        foreach (var (contrast, independent) in design.DefaultContrasts().Zip(additive.DefaultContrasts()))
        {
            Assert.That(design.Estimability(contrast), Is.EqualTo(AnalysisDesign.Estimable));
            double[] w = design.FittedWeights(contrast);
            double estimate = w.Select((x, j) => x * fit.Coefficient(0, j)).Sum();
            double expected = independent.Weights.Select((x, j) => x * reference.Coefficient(0, j)).Sum();
            Assert.That(estimate, Is.EqualTo(expected).Within(1e-10), contrast.Id);
        }
    }

    /// <summary>The same when a covariate does not vary: the condition is fitted without it.</summary>
    [Test]
    public void AConditionBesideAConstantCovariateIsFittedWithoutIt()
    {
        var samples = Enumerable.Range(0, 6).Select(i => new AnalysisSample($"s{i}")
        {
            Factors = new Dictionary<string, string> { ["condition"] = i % 2 == 0 ? "A" : "B" },
            Covariates = new Dictionary<string, double> { ["age"] = 18.2 },
        }).ToList();
        var design = AnalysisDesign.Create(samples, Spec(("condition", "A")) with { Covariates = new[] { new CovariateSpec("age", "months") } });
        var withoutAge = AnalysisDesign.Create(samples, Spec(("condition", "A")));

        Assert.That(design.FittedColumns.Select(j => design.ColumnNames[j]), Is.EqualTo(withoutAge.ColumnNames));
        var y = new double[,] { { 20.1, 21.7, 19.4, 20.6, 22.3, 18.9 } };
        var contrast = design.Between("condition", "B", "A");
        var fit = LinearModel.Fit(y, design.FittedMatrix);
        double estimate = design.FittedWeights(contrast).Select((x, j) => x * fit.Coefficient(0, j)).Sum();
        Assert.That(estimate, Is.EqualTo(LinearModel.Fit(y, withoutAge.Matrix).Coefficient(0, 1)).Within(1e-10));
    }

    /// <summary>A contrast that is not estimable has no weights on the kept columns: its estimate would depend on which were kept.</summary>
    [Test]
    public void ANonEstimableContrastHasNoFittedWeights()
    {
        var design = ConstantAgeBesideAConfoundedCondition(24, 10);

        var e = Assert.Throws<ArgumentException>(() => design.FittedWeights(design.Slope("age")));
        Assert.That(e!.Message, Does.Contain(AnalysisDesign.NotEstimableRankDeficient));
        Assert.Throws<ArgumentException>(() => design.FittedWeights(design.Between("condition", "B", "A")));
    }
}
