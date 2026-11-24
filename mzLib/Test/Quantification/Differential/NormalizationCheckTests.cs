using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// <see cref="NormalizationCheck"/>: its empirical null against R's locfdr 1.1.8 (normalization_check_null.tsv,
/// PROVENANCE_normalization_check.txt; R is never run here), its verdicts, and the runner and report around it.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class NormalizationCheckTests
{
    private static List<Dictionary<string, string>> Fixture()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential", "ReferenceData",
            "normalization_check_null.tsv");
        var lines = File.ReadAllLines(path).Where(l => l.Length > 0).ToArray();
        string[] header = lines[0].Split('\t');
        return lines.Skip(1).Select(l => header.Zip(l.Split('\t')).ToDictionary(p => p.First, p => p.Second)).ToList();
    }

    private static double D(string s) => s == "NA" ? double.NaN : double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    /// <summary>A recipe's z-values, as the R script builds them: parts "n,mean,sd" (mean + sd Φ⁻¹((i − 0.5)/n), or n
    /// copies of mean when sd is 0) separated by ';', and "|k" to repeat each value k times.</summary>
    private static double[] Build(string recipe)
    {
        string[] halves = recipe.Split('|');
        int each = halves.Length > 1 ? int.Parse(halves[1], CultureInfo.InvariantCulture) : 1;
        return halves[0].Split(';').Select(p => p.Split(',').Select(D).ToArray())
            .SelectMany(v => Enumerable.Range(1, (int)v[0]).Select(i => v[2] == 0 ? v[1] : v[1] + v[2] * Normal.InvCDF(0, 1, (i - 0.5) / v[0])))
            .SelectMany(x => Enumerable.Repeat(x, each)).ToArray();
    }

    private static IEnumerable<TestCaseData> Scenarios() =>
        Fixture().Select(r => new TestCaseData(r["scenario"]).SetName("TheNullMatchesLocfdr_" + r["scenario"]));

    [TestCaseSource(nameof(Scenarios))]
    public void TheNullMatchesLocfdr(string scenario)
    {
        var r = Fixture().Single(x => x["scenario"] == scenario);
        var fit = NormalizationCheck.FitNull(Build(r["recipe"]));

        // The exact maximum within the simplex's precision; locfdr's own answer within its stopping rule (about 1e-5,
        // more as sigma grows). Where locfdr stops with "ML estimation failed" (clamped_both_25) the maximum exists all
        // the same, and it is found.
        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.Value.Delta, Is.EqualTo(D(r["mle_delta0"])).Within(1e-6));
        Assert.That(fit.Value.Sigma, Is.EqualTo(D(r["mle_sigma0"])).Within(1e-6));
        if (r["locfdr_delta0"] == "NA") return;
        double tolerance = 2e-5 * Math.Max(1, fit.Value.Sigma);
        Assert.That(fit.Value.Delta, Is.EqualTo(D(r["locfdr_delta0"])).Within(tolerance));
        Assert.That(fit.Value.Sigma, Is.EqualTo(D(r["locfdr_sigma0"])).Within(tolerance));
    }

    /// <summary>Proteins whose scores are the recipe's z-values: effect z, two-sided p 2Φ(−|z|).</summary>
    private static List<(double, double)> Proteins(string recipe) =>
        Build(recipe).Select(z => (z, 2 * Normal.CDF(0, 1, -Math.Abs(z)))).ToList();

    [Test]
    public void TheVerdictWarnsAboveTheThresholdForTheGroupSizeAndOtherwiseSaysWhy()
    {
        Assert.That(new[] { 2, 3, 4, 9 }.Select(NormalizationCheck.Threshold),
            Is.EqualTo(new[] { 0.410702939399882, 0.363279703570296, 0.429775562501974, 0.429775562501974 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => NormalizationCheck.Threshold(1));

        Assert.That(NormalizationCheck.Centre(Proteins("2000,0.25,1")), Is.EqualTo(0.25).Within(1e-5));
        Assert.That(NormalizationCheck.Verdict(Proteins("2000,0.25,1"), 3), Is.EqualTo("no warning (centre 0.25, threshold 0.363)"));
        Assert.That(NormalizationCheck.Verdict(Proteins("2000,0.25,1"), 2), Is.EqualTo("no warning (centre 0.25, threshold 0.411)"));
        Assert.That(NormalizationCheck.Verdict(Proteins("2000,-0.4,1"), 3), Does.StartWith("warning: the unchanged proteins' scores centre at -0.4,"));
        Assert.That(NormalizationCheck.Verdict(Proteins("199,0,1"), 3), Is.EqualTo("not checked: 199 fitted proteins, fewer than 200"));
        Assert.That(NormalizationCheck.Verdict(Proteins("300,0,1"), 1), Is.EqualTo("not checked: a group of one sample"));
        Assert.That(NormalizationCheck.Verdict(Proteins("300,0.5,0"), 3), Is.EqualTo("not checked: no empirical null could be fitted"));
        // Scores are clamped at 8, as R's were. With 30% of the proteins at p = 0 the upper quartile is a clamped score, so a
        // score of -Φ⁻¹(ε) ≈ 8.1 instead would widen the window and move δ0.
        var extreme = Proteins("1400,0,1").Concat(Enumerable.Repeat((1.0, 0.0), 600)).ToList();
        Assert.That(NormalizationCheck.Centre(extreme), Is.EqualTo(NormalizationCheck.FitNull(Build("1400,0,1;600,8,0"))!.Value.Delta).Within(1e-6));
    }

    // ---- the runner and the report ----

    private static readonly string[] Samples = { "a1", "a2", "a3", "b1", "b2", "b3" };

    private static double Gaussian(Random rng) => Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

    /// <summary>A label-free set of 240 protein groups of 3 peptides, noise 0.2 log2, with the first
    /// <paramref name="changedShare"/> of them 1 log2 higher in group A.</summary>
    private static DifferentialRequest Request(double changedShare, params string[] samples)
    {
        var rng = new Random(20261124);
        var peptides = new List<PeptideRunValues>();
        var groups = new List<ProteinGroupInfo>();
        for (int p = 0; p < 240; p++)
        {
            groups.Add(new ProteinGroupInfo($"P{p:000}", 0, false, false));
            for (int k = 0; k < 3; k++)
            {
                double level = 18 + 6 * rng.NextDouble();
                peptides.Add(new PeptideRunValues($"P{p:000}PEP{k}", $"P{p:000}PEP{k}", new[] { $"P{p:000}" }, Samples.ToDictionary(s => s,
                    s => new RunValue(Math.Pow(2, level + (p < changedShare * 240 && s[0] == 'a' ? 1 : 0) + 0.2 * Gaussian(rng)), ObservationState.Quantified))));
            }
        }
        var runs = Samples.Select(s => new ObservationRun(s, s)).ToList();
        return new DifferentialRequest
        {
            Table = ObservationTable.Build(runs, Samples, peptides, groups, QuantBasis.MsmsOnly),
            Samples = (samples.Length > 0 ? samples : Samples).Select(s => new AnalysisSample(s)
                { Factors = new Dictionary<string, string> { ["grp"] = s[..1].ToUpperInvariant() } }).ToList(),
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("grp", "B") } },
            Contrasts = new[] { DifferentialContrastRequest.Between("grp", "A", "B") },
        };
    }

    private static string Report(DifferentialAnalysisResult r) =>
        StatisticalReportWriter.Build(r.Metadata, r.Rows, r.Diagnostics, r.Decisions).Markdown;

    [Test]
    public void TheRunnerChecksEachComparisonAndTheReportSaysSo()
    {
        // 3 v 2 (b3 left out), nothing changed: no warning, judged against the threshold for groups of 2.
        var quiet = DifferentialAnalysis.Run(Request(0, "a1", "a2", "a3", "b1", "b2"), maxThreads: 4);
        Assert.That(quiet.Metadata.Normalization.Single().Checks!["c1"], Does.StartWith("no warning (centre ").And.EndWith("threshold 0.411)"));
        Assert.That(Report(quiet), Does.Contain("log2 scale. The data were normalized.")
            .And.Contain("This assumes that most of the reference peptides do not change")
            .And.Contain("Normalization check: `c1` no warning")
            .And.Contain("a comparison without a warning can still be off."));

        // 20% of the proteins 1 log2 up in A: each sample's median difference moves by about 0.06, up in A and down in B,
        // so the unchanged proteins' A − B effect is about −0.12, some −1.2 standard errors.
        var shifted = DifferentialAnalysis.Run(Request(0.2), maxThreads: 4);
        string verdict = shifted.Metadata.Normalization.Single().Checks!["c1"];
        Assert.That(verdict, Does.StartWith("warning: the unchanged proteins' scores centre at -1."));
        Assert.That(Encoding.UTF8.GetString(DifferentialMetadataWriter.ToBytes(shifted.Metadata, shifted.Rows)),
            Does.Contain("\"normalization_check\": {").And.Contain("\"c1\": \"warning:"));
        Assert.That(Report(shifted), Does.Contain("Normalization check: `c1` warning: the unchanged proteins' scores centre at -1."));

        // Where it was not calibrated it does not run, and the report says so.
        var robustOff = DifferentialAnalysis.Run(Request(0.2) with { Settings = new DifferentialSettings { Robust = false } }, 4);
        Assert.That(robustOff.Metadata.Normalization.Single().Checks, Is.Null);
        Assert.That(Report(robustOff), Does.Contain("The normalization check was not run: it is calibrated only for label-free"));
        var none = DifferentialAnalysis.Run(Request(0.2) with { Settings = new DifferentialSettings { Normalization = SharedPeptideNormalization.None } }, 4);
        Assert.That(Report(none), Does.Contain("log2 scale. The data were not normalized.").And.Not.Contain("normalization check"));
    }
}
