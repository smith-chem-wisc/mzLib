using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// <see cref="NormalizationCheck"/> (<c>DEF-DIFF-NORMCHECK</c>): its thresholds, scores and verdicts, and the runner's
/// check of each comparison: where it runs, where it is not checked and why, and that it changes no result. The empirical
/// null itself is checked against locfdr in <c>Test.StatisticalModels.EmpiricalNullTests</c>, and the whole check against
/// the third NORMCHECK attempt's 3,273 comparisons in QuantProject (results/normcheck/parity).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class NormalizationCheckTests
{
    private static double Gaussian(Random rng) =>
        Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

    /// <summary>n proteins whose effects are N(shift, se²), se in [0.1, 0.3], df 6.</summary>
    private static List<(double, double, double)> Proteins(int n, double shift, int seed)
    {
        var rng = new Random(seed);
        var list = new List<(double, double, double)>();
        for (int i = 0; i < n; i++)
        {
            double se = 0.1 + 0.2 * rng.NextDouble();
            list.Add((shift + se * Gaussian(rng), se, 6));
        }
        return list;
    }

    [Test]
    public void TheThresholdsAreTheYardsticksAndLargerGroupsTakeTheOneForFour()
    {
        Assert.That(NormalizationCheck.Threshold(2), Is.EqualTo(0.410702939399882));
        Assert.That(NormalizationCheck.Threshold(3), Is.EqualTo(0.363279703570296));
        Assert.That(NormalizationCheck.Threshold(4), Is.EqualTo(0.429775562501974));
        Assert.That(NormalizationCheck.Threshold(9), Is.EqualTo(NormalizationCheck.Threshold(4)));
        Assert.Throws<ArgumentOutOfRangeException>(() => NormalizationCheck.Threshold(1));
    }

    [Test]
    public void ADeltaOnTheThresholdDoesNotWarnAndOneJustAboveDoes()
    {
        foreach (int k in new[] { 2, 3, 4, 7 })
        {
            double t = NormalizationCheck.Threshold(k);
            Assert.That(NormalizationCheck.Warns(t, k), Is.False);
            Assert.That(NormalizationCheck.Warns(-t, k), Is.False);
            Assert.That(NormalizationCheck.Warns(Math.BitIncrement(t), k), Is.True);
            Assert.That(NormalizationCheck.Warns(-Math.BitIncrement(t), k), Is.True);
        }
    }

    [TestCase(0.5, 3)]
    [TestCase(-1.7, 6)]
    [TestCase(2.9, 12)]
    [TestCase(-4.2, 40)]
    public void TheScoreIsTheNormalQuantileOfTheTProbability(double t, double df)
    {
        // Efron's z = qnorm(pt(t, df)), here with t = effect / se and se 0.2.
        double expected = Normal.InvCDF(0, 1, StudentT.CDF(0, 1, df, t));
        Assert.That(NormalizationCheck.Score(0.2 * t, 0.2, df), Is.EqualTo(expected).Within(1e-10));
        Assert.That(NormalizationCheck.Score(-0.2 * t, 0.2, df), Is.EqualTo(-expected).Within(1e-10));
    }

    [Test]
    public void TheScoreIsZeroAtZeroUsesTheNormalAtInfiniteDfAndIsClampedAtEight()
    {
        Assert.That(NormalizationCheck.Score(0, 0.3, 5), Is.EqualTo(0));
        Assert.That(NormalizationCheck.Score(0.4, 0.2, double.PositiveInfinity), Is.EqualTo(2).Within(1e-12));
        // t = 300: with 100 df its tail is far below Phi(-8); with 4 df it is about 1e-10, so z is about 6.2, not clamped.
        Assert.That(NormalizationCheck.Score(30, 0.1, 100), Is.EqualTo(NormalizationCheck.ScoreBound));
        Assert.That(NormalizationCheck.Score(-30, 0.1, 100), Is.EqualTo(-NormalizationCheck.ScoreBound));
        Assert.That(NormalizationCheck.Score(30, 0.1, 4), Is.EqualTo(Normal.InvCDF(0, 1, StudentT.CDF(0, 1, 4, 300))).Within(1e-6));
        Assert.That(NormalizationCheck.Score(1e6, 1e-6, 1e6), Is.EqualTo(NormalizationCheck.ScoreBound));
    }

    [Test]
    public void AnUnshiftedSetHasNoWarningAndAShiftedOneWarns()
    {
        var flat = NormalizationCheck.Check("all", "msms_only", "c1", Proteins(2000, 0, 20261124), 3);
        Assert.That(flat.Status, Is.EqualTo(NormalizationCheck.NoWarning));
        Assert.That(flat.NotCheckedReason, Is.Null);
        Assert.That(Math.Abs(flat.Delta0!.Value), Is.LessThan(0.1));
        Assert.That(flat.Sigma0, Is.EqualTo(1).Within(0.1));
        Assert.That(flat.Threshold, Is.EqualTo(NormalizationCheck.Threshold(3)));
        Assert.That(flat.FittedProteins, Is.EqualTo(2000));
        Assert.That(flat.Yardstick, Is.EqualTo(NormalizationCheck.Yardstick));

        // Every protein moved by 0.15 log2, about 0.75 standard errors: delta0 near 0.75, above 0.36.
        var shifted = NormalizationCheck.Check("all", "msms_only", "c1", Proteins(2000, 0.15, 20261124), 3);
        Assert.That(shifted.Status, Is.EqualTo(NormalizationCheck.Warning));
        Assert.That(shifted.Delta0, Is.EqualTo(0.75).Within(0.1));
    }

    [Test]
    public void EHatIsDeltaTimesTheMedianStandardError()
    {
        var proteins = Proteins(501, 0.1, 7);
        var check = NormalizationCheck.Check("all", null, "c1", proteins, 4);
        var ses = proteins.Select(p => p.Item2).OrderBy(v => v).ToList();
        Assert.That(check.EHat, Is.EqualTo(check.Delta0!.Value * ses[250]).Within(1e-12));
    }

    [Test]
    public void ProteinsWithoutAFiniteEffectSeOrDfAreNotScoredAndFewerThan200AreNotChecked()
    {
        var proteins = Proteins(199, 0, 3);
        proteins.Add((double.NaN, 0.2, 6));
        proteins.Add((0.1, 0, 6));
        proteins.Add((0.1, 0.2, double.PositiveInfinity));
        var check = NormalizationCheck.Check("all", "msms_only", "c1", proteins, 3);
        Assert.That(check.Status, Is.EqualTo(NormalizationCheck.NotChecked));
        Assert.That(check.NotCheckedReason, Is.EqualTo(NormalizationCheck.FewerThanMinimumFitted));
        Assert.That(check.FittedProteins, Is.EqualTo(199));
        Assert.That(new[] { check.Delta0, check.Sigma0, check.EHat, check.Threshold }, Has.All.Null);
        Assert.That(NormalizationCheck.Check("all", "msms_only", "c1", Proteins(200, 0, 3), 3).Status, Is.Not.EqualTo(NormalizationCheck.NotChecked));
    }

    [Test]
    public void AGroupOfOneAndScoresThatCannotBeFittedAreNotChecked()
    {
        Assert.That(NormalizationCheck.Check("all", null, "c1", Proteins(300, 0, 5), 1).NotCheckedReason,
            Is.EqualTo(NormalizationCheck.GroupOfOne));
        // Every protein the same: the scores tie, the window has no width, and no null can be fitted.
        var same = Enumerable.Repeat((0.1, 0.2, 6.0), 300);
        Assert.That(NormalizationCheck.Check("all", null, "c1", same, 3).NotCheckedReason, Is.EqualTo(NormalizationCheck.FitFailed));
    }

    // ---- the runner ----

    private static readonly string[] Samples = { "a1", "a2", "a3", "b1", "b2", "b3" };

    /// <summary>
    /// A label-free 3 v 3 set: <paramref name="proteins"/> protein groups of 3 peptides, noise 0.2 log2. In group A,
    /// <paramref name="changedShare"/> of the proteins (the first ones) are raised by 1 log2, which moves the normalization.
    /// </summary>
    private static DifferentialRequest Request(int proteins, double changedShare, int seed, QuantBasis basis = QuantBasis.MsmsOnly)
    {
        var rng = new Random(seed);
        var peptides = new List<PeptideRunValues>();
        var groups = new List<ProteinGroupInfo>();
        int changed = (int)Math.Round(changedShare * proteins);
        for (int p = 0; p < proteins; p++)
        {
            string group = $"P{p:000}";
            groups.Add(new ProteinGroupInfo(group, 0, false, false));
            for (int k = 0; k < 3; k++)
            {
                double level = 18 + 6 * rng.NextDouble();
                peptides.Add(new PeptideRunValues($"{group}PEP{k}", $"{group}PEP{k}", new[] { group }, Samples.ToDictionary(s => s,
                    s => new RunValue(Math.Pow(2, level + (p < changed && s[0] == 'a' ? 1 : 0) + 0.2 * Gaussian(rng)), ObservationState.Quantified))));
            }
        }
        var runs = Samples.Select(s => new ObservationRun(s, s)).ToList();
        return new DifferentialRequest
        {
            Table = ObservationTable.Build(runs, Samples, peptides, groups, basis),
            Samples = Samples.Select(s => new AnalysisSample(s)
            {
                Factors = new Dictionary<string, string> { ["grp"] = s[..1].ToUpperInvariant() },
                Covariates = new Dictionary<string, double> { ["age"] = 20 + 5 * Array.IndexOf(Samples, s) },
            }).ToList(),
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("grp", "B") } },
            Contrasts = new[] { DifferentialContrastRequest.Between("grp", "A", "B") },
        };
    }

    private static DifferentialNormalizationCheck OnlyCheck(DifferentialAnalysisResult result) =>
        result.Metadata.NormalizationChecks.Single();

    [Test]
    public void TheRunnerWarnsWhenManyProteinsChangeOneWayAndNotOnANullSet()
    {
        var quiet = DifferentialAnalysis.Run(Request(240, 0, 20261124), maxThreads: 4);
        var check = OnlyCheck(quiet);
        Assert.That(check.Status, Is.EqualTo(NormalizationCheck.NoWarning));
        Assert.That(check.Stratum, Is.EqualTo(quiet.Rows[0].Stratum));
        Assert.That(check.QuantBasis, Is.EqualTo("msms_only"));
        Assert.That(check.Contrast, Is.EqualTo(quiet.Rows[0].ContrastId));
        Assert.That(check.GroupSize, Is.EqualTo(3));
        Assert.That(check.FittedProteins, Is.EqualTo(quiet.Rows.Count(r => DifferentialStatus.IsFitted(r.Status))));

        // 20% of the proteins up by 1 log2 on one side. Each sample's median difference from the reference then moves by
        // about 0.06, up in A and down in B, so the unchanged proteins' A - B effect is about -0.12, some -1.2 standard
        // errors: delta0 near -1.2. (At 40% the changed proteins reach the quartiles, the fitted null spans both clusters
        // and delta0 takes the wrong sign; the empirical null assumes most proteins do not change.)
        var shifted = DifferentialAnalysis.Run(Request(240, 0.2, 20261124), maxThreads: 4);
        var warned = OnlyCheck(shifted);
        Assert.That(warned.Status, Is.EqualTo(NormalizationCheck.Warning));
        Assert.That(warned.Delta0, Is.LessThan(-NormalizationCheck.Threshold(3)), "the unchanged proteins are pushed down");

        // The check changes nothing: no decision line names it, the analysis id is the settings' and inputs', and the
        // unchanged proteins keep the shift the normalization gave them.
        Assert.That(shifted.Decisions.Select(d => d.Rule), Has.None.Contain("normalization_check"));
        Assert.That(shifted.Metadata.AnalysisId, Is.EqualTo(DifferentialMetadata.ComputeAnalysisId(shifted.Metadata.Inputs, shifted.Metadata.Settings)));
        double unchanged = shifted.Rows.Where(r => string.CompareOrdinal(r.FeatureId, "P048") >= 0).Select(r => r.Log2Effect!.Value).Average();
        Assert.That(unchanged, Is.EqualTo(-0.12).Within(0.05), "the effects are reported as fitted, not corrected");

        string report = StatisticalReportWriter.Build(shifted.Metadata, shifted.Rows, shifted.Diagnostics, shifted.Decisions).Markdown;
        Assert.That(report, Does.Contain("The data were normalized."));
        Assert.That(report, Does.Contain("**Normalization check.** 1 of 1 comparison checked: 1 with a warning; 0 not checked."));
        Assert.That(report, Does.Contain("looks off. The unchanged proteins' scores centre at -"));
        Assert.That(report, Does.Contain("The normalization check warned for this comparison (section 3)."));
        string quietReport = StatisticalReportWriter.Build(quiet.Metadata, quiet.Rows, quiet.Diagnostics, quiet.Decisions).Markdown;
        Assert.That(quietReport, Does.Contain("1 of 1 comparison checked: 0 with a warning; 0 not checked."));
        Assert.That(quietReport, Does.Not.Contain("looks off").And.Not.Contain("warned for this comparison"));
    }

    [Test]
    public void TheGroupSizeIsTheSmallerSide()
    {
        // Three samples in A and two in B: b3 is left out of the analysis, so the smaller side has 2 and the larger 3.
        string[] kept = { "a1", "a2", "a3", "b1", "b2" };
        var request = Request(40, 0, 11);
        request = request with { Samples = request.Samples.Where(s => kept.Contains(s.SampleId)).ToList() };
        var check = OnlyCheck(DifferentialAnalysis.Run(request, maxThreads: 2));
        Assert.That(check.GroupSize, Is.EqualTo(2));
    }

    [Test]
    public void FewerThan200FittedProteinsAreNotChecked()
    {
        var result = DifferentialAnalysis.Run(Request(40, 0, 11), maxThreads: 2);
        var check = OnlyCheck(result);
        Assert.That(check.NotCheckedReason, Is.EqualTo(NormalizationCheck.FewerThanMinimumFitted));
        Assert.That(check.FittedProteins, Is.EqualTo(40));
        string report = StatisticalReportWriter.Build(result.Metadata, result.Rows, result.Diagnostics, result.Decisions).Markdown;
        Assert.That(report, Does.Contain("fewer than 200 fitted proteins (40)."));
    }

    [Test]
    public void SettingsItWasNotCalibratedOnAreNamed()
    {
        string Reason(DifferentialRequest request) => OnlyCheck(DifferentialAnalysis.Run(request, maxThreads: 2)).NotCheckedReason!;

        Assert.That(Reason(Request(40, 0, 11) with { Settings = new DifferentialSettings { Robust = false } }),
            Is.EqualTo("not_calibrated:robust_off"));
        Assert.That(Reason(Request(40, 0, 11) with { Settings = new DifferentialSettings { Normalization = SharedPeptideNormalization.None } }),
            Is.EqualTo("not_calibrated:none"));
        Assert.That(Reason(Request(40, 0, 11, QuantBasis.MbrKept)), Is.EqualTo("not_calibrated:mbr_kept"));
        Assert.That(Reason(Request(40, 0, 11) with { QuantStyle = DifferentialConstants.DiaStyle }), Is.EqualTo("not_calibrated:dia"));
        Assert.That(Reason(Request(40, 0, 11) with
        {
            // 40 declared groups of 3 peptides: 120 reference peptides, above the floor of 100, so the stratum is shifted.
            Settings = new DifferentialSettings
            {
                Normalization = "background_set:steady", BackgroundSetGroups = Enumerable.Range(0, 40).Select(p => $"P{p:000}").ToArray(),
            },
        }), Is.EqualTo("not_calibrated:background_set"));
        // Three declared groups give 9 reference peptides, below the floor: the stratum is not shifted, its entry says none,
        // and so do the check and the report.
        var below = DifferentialAnalysis.Run(Request(40, 0, 11) with
        {
            Settings = new DifferentialSettings { Normalization = "background_set:steady", BackgroundSetGroups = new[] { "P000", "P001", "P002" } },
        }, 2);
        Assert.That(OnlyCheck(below).NotCheckedReason, Is.EqualTo("not_calibrated:none"));
        string belowReport = StatisticalReportWriter.Build(below.Metadata, below.Rows, below.Diagnostics, below.Decisions).Markdown;
        Assert.That(belowReport, Does.Contain("The data were not normalized.").And.Not.Contain("The data were normalized."));

        var none = DifferentialAnalysis.Run(Request(40, 0, 11) with { Settings = new DifferentialSettings { Normalization = SharedPeptideNormalization.None } }, 2);
        string report = StatisticalReportWriter.Build(none.Metadata, none.Rows, none.Diagnostics, none.Decisions).Markdown;
        Assert.That(report, Does.Contain("The data were not normalized.").And.Not.Contain("The data were normalized."));
        Assert.That(report, Does.Contain("this run used `none`."));
    }

    [Test]
    public void ASlopeIsNotChecked()
    {
        var request = Request(40, 0, 11) with
        {
            Spec = new DesignSpec { Covariates = new[] { new CovariateSpec("age", "years") } },
            Contrasts = new[] { DifferentialContrastRequest.Slope("age") },
        };
        var check = OnlyCheck(DifferentialAnalysis.Run(request, maxThreads: 2));
        Assert.That(check.NotCheckedReason, Is.EqualTo(NormalizationCheck.Slope));
        Assert.That(check.GroupSize, Is.Null);
    }
}
