using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The STAT1 runner end to end against R (<c>ReferenceData/make_runner_fixtures.R</c>, which says what R does and why):
/// a two-stratum label-free dataset with a protein missing a group, one missing the reference, a peptide seen in one group
/// only, a protein with one value per sample (least squares, ruling 16), one-peptide proteins, a shared peptide, and decoy,
/// entrapment and undescribed groups. Every row R fits is fitted by the runner with the same effect, means, SE, df, p and
/// adjusted p, and the runner fits nothing else; the global shifts agree; the groups with no rows have none. R is never
/// run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialAnalysisReferenceTests
{
    /// <summary>Effects and means, absolute log2.</summary>
    private const double Effect = 1e-6;
    /// <summary>SE, p and adjusted p, relative.</summary>
    private const double Inference = 1e-4;
    /// <summary>Satterthwaite df, relative (numerical derivatives, as lmerTest's).</summary>
    private const double Df = 1e-3;

    private static List<Dictionary<string, string>> Read(string file)
    {
        var lines = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential",
            "ReferenceData", file)).Where(l => l.Length > 0).ToList();
        var header = lines[0].Split('\t');
        return lines.Skip(1).Select(l => header.Zip(l.Split('\t')).ToDictionary(z => z.First, z => z.Second)).ToList();
    }

    private static double Num(string s) => s == "NaN" ? double.NaN : double.Parse(s, CultureInfo.InvariantCulture);

    private static DifferentialAnalysisResult Run()
    {
        var values = Read("runner_values.tsv");
        var samples = values.GroupBy(v => v["sample"]).OrderBy(g => g.Key, StringComparer.Ordinal)
            .Select(g => new AnalysisSample(g.Key)
            {
                Factors = new Dictionary<string, string> { ["group"] = g.First()["group"] },
                Strata = new Dictionary<string, string> { ["tissue"] = g.First()["tissue"] },
            }).ToList();
        var peptides = values.GroupBy(v => v["peptide"]).Select(g => new PeptideRunValues(g.Key, g.Key,
            g.First()["groups"].Split(';'), g.ToDictionary(v => v["sample"],
                v => new RunValue(Math.Pow(2, Num(v["y"])), ObservationState.Quantified)))).ToList();
        var groups = Enumerable.Range(1, 34).Select(i => $"P{i:00}").Concat(Enumerable.Range(1, 5).Select(i => $"S{i:00}"))
            .Select(g => new ProteinGroupInfo(g, 0, false, false))
            .Append(new ProteinGroupInfo("DECOY", 0, true, false))
            .Append(new ProteinGroupInfo("ENTRAP", 0, false, false, IsEntrapment: true));
        var runs = samples.Select(s => new ObservationRun(s.SampleId, s.SampleId)).ToList();
        return DifferentialAnalysis.Run(new DifferentialRequest
        {
            Table = ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, groups, QuantBasis.MsmsOnly),
            Samples = samples,
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("group", "A") } },
            Contrasts = new[]
            {
                DifferentialContrastRequest.Between("group", "B", "A"), DifferentialContrastRequest.Between("group", "C", "A"),
                DifferentialContrastRequest.Between("group", "C", "B"),
            },
            Settings = new DifferentialSettings { Robust = false, Normalization = SharedPeptideNormalization.None },
        }, maxThreads: 1);
    }

    [Test]
    public void EveryRowRFitsIsFittedWithTheSameNumbersAndNothingElseIs()
    {
        var result = Run();
        var expected = Read("runner_results.tsv");
        var fitted = result.Rows.Where(r => DifferentialStatus.IsFitted(r.Status))
            .ToDictionary(r => (r.Stratum, r.FeatureId, r.ContrastId));

        Assert.That(fitted.Keys, Is.EquivalentTo(expected.Select(e => (e["stratum"], e["protein"], e["contrast"]))),
            "the runner fits exactly the rows R fits");
        foreach (var e in expected)
        {
            var r = fitted[(e["stratum"], e["protein"], e["contrast"])];
            string at = $"{e["stratum"]} {e["protein"]} {e["contrast"]}";
            Assert.That(r.Status, Is.EqualTo(e["status"]), at);
            Close(r.Log2Effect!.Value, e["log2"], Effect, false, at + " log2");
            Close(r.MeanNumerator!.Value, e["mean_num"], Effect, false, at + " mean numerator");
            Close(r.MeanDenominator!.Value, e["mean_den"], Effect, false, at + " mean denominator");
            Close(r.StandardError!.Value, e["se"], Inference, true, at + " se");
            Close(r.Df!.Value, e["df"], Df, true, at + " df");
            Close(r.PValue!.Value, e["p"], Inference, true, at + " p");
            Close(r.PAdjusted!.Value, e["p_adj"], Inference, true, at + " adjusted p");
            if (e["model"] == "lm" && e["status"] == DifferentialStatus.Fitted)
                Assert.That(r.DfMethod, Is.EqualTo(DifferentialConstants.DfResidual), at + " least squares");
        }
    }

    /// <summary>
    /// The report's mean-variance diagnostic against R: one entry per stratum and model family, its proteins those R fits in
    /// that family, each with R's residual variance (<c>sigma(m)^2</c>), the family's prior (<c>squeezeVar</c>'s
    /// <c>var.prior</c>; none with fewer than 3 fitted) and the mean of the values it was fitted on.
    /// </summary>
    [Test]
    public void TheMeanVarianceDiagnosticHasRsVariancesAndPriors()
    {
        var result = Run();
        var models = new Dictionary<string, string>
        {
            ["mixed"] = DifferentialConstants.MixedModel, ["reduced"] = DifferentialConstants.MixedModelReduced,
            ["t"] = DifferentialConstants.ModeratedT,
        };
        var expected = Read("runner_results.tsv").GroupBy(e => (e["stratum"], e["family"], e["protein"])).Select(g => g.First()).ToList();
        Assert.That(result.Diagnostics.Fits.Select(f => (f.Stratum, f.ModelUsed)),
            Is.EquivalentTo(expected.Select(e => (e["stratum"], models[e["family"]])).Distinct()), "one entry per stratum and family");
        foreach (var fit in result.Diagnostics.Fits)
        {
            Assert.That(fit.QuantBasis, Is.EqualTo("msms_only"));
            var mine = expected.Where(e => e["stratum"] == fit.Stratum && models[e["family"]] == fit.ModelUsed).ToList();
            Assert.That(fit.FeatureIds, Is.EquivalentTo(mine.Select(e => e["protein"])), $"{fit.Stratum} {fit.ModelUsed}");
            Assert.That(fit.FeatureIds, Is.Ordered.Using((IComparer<string>)StringComparer.Ordinal));
            foreach (var e in mine)
            {
                int i = fit.FeatureIds.ToList().IndexOf(e["protein"]);
                string at = $"{fit.Stratum} {fit.ModelUsed} {e["protein"]}";
                Close(fit.ResidualVariance[i], e["s2"], Inference, true, at + " residual variance");
                Close(fit.AverageLog2[i], e["ave"], 1e-9, false, at + " average log2");
                if (e["s2_prior"] == "NaN") Assert.That(fit.PriorVariance[i], Is.NaN, at + " no prior");
                else Close(fit.PriorVariance[i], e["s2_prior"], Inference, true, at + " prior variance");
            }
        }
    }

    /// <summary>
    /// The report's sample-correlation values: per stratum, every sample's log2 value (normalization is off here) on every
    /// peptide unique to a group in the row set, so not the shared peptide, nor the decoy's, the entrapment group's or the
    /// undescribed group's.
    /// </summary>
    [Test]
    public void TheSampleMatrixHoldsEverySamplesValueOnTheRowSetsUniquePeptides()
    {
        var result = Run();
        var values = Read("runner_values.tsv");
        string[] excluded = { "DECOY", "ENTRAP", "UNDESC" };

        Assert.That(result.Diagnostics.Samples.Select(s => s.Stratum), Is.EquivalentTo(new[] { "tissue=brain", "tissue=liver" }));
        foreach (var matrix in result.Diagnostics.Samples)
        {
            string tissue = matrix.Stratum["tissue=".Length..];
            var mine = values.Where(v => v["tissue"] == tissue).ToList();
            var peptides = mine.Where(v => !v["groups"].Contains(';') && !excluded.Contains(v["groups"])).Select(v => v["peptide"])
                .Distinct().Order(StringComparer.Ordinal).ToList();
            Assert.That(matrix.QuantBasis, Is.EqualTo("msms_only"));
            Assert.That(matrix.SampleIds, Is.EqualTo(mine.Select(v => v["sample"]).Distinct().Order(StringComparer.Ordinal)));
            var cell = mine.ToDictionary(v => (v["sample"], v["peptide"]), v => Num(v["y"]));
            for (int s = 0; s < matrix.SampleIds.Count; s++)
            {
                Assert.That(matrix.Values[s], Has.Count.EqualTo(peptides.Count), "the peptides, in order");
                for (int k = 0; k < peptides.Count; k++)
                {
                    double v = matrix.Values[s][k];
                    if (cell.TryGetValue((matrix.SampleIds[s], peptides[k]), out double y))
                        Assert.That(v, Is.EqualTo(y).Within(1e-9), $"{matrix.SampleIds[s]} {peptides[k]}");
                    else Assert.That(v, Is.NaN, $"{matrix.SampleIds[s]} {peptides[k]} has no value");
                }
            }
        }
    }

    [Test]
    public void TheGlobalShiftsAgreeAndTheExcludedGroupsHaveNoRows()
    {
        var result = Run();

        foreach (var e in Read("runner_shifts.tsv"))
        {
            var shift = result.Metadata.Contrasts.Single(c => c.Id == e["contrast"]).GlobalShift!.Single(s => s.Stratum == e["stratum"]);
            Close(shift.Value, e["value"], 1e-12, false, $"{e["stratum"]} {e["contrast"]} shift");
            Assert.That(shift.Peptides.ToString(CultureInfo.InvariantCulture), Is.EqualTo(e["peptides"]));
        }
        Assert.That(result.Rows.Select(r => r.FeatureId).Distinct(), Has.None.AnyOf("DECOY", "ENTRAP", "UNDESC"));
        Assert.That(result.Rows.Where(r => r.FeatureId == "P31" && r.ContrastId != "c1").Select(r => r.Status),
            Is.All.EqualTo(DifferentialStatus.AbsentInNumerator));
        Assert.That(result.Rows.Where(r => r.FeatureId == "P32" && r.ContrastId != "c3").Select(r => r.Status),
            Is.All.EqualTo(DifferentialStatus.AbsentInDenominator));
        Assert.That(result.Decisions.Where(d => d.Rule == DifferentialDecisionRules.RandomTermsDropped).Select(d => d.FeatureId).Distinct(),
            Is.EquivalentTo(new[] { "P34" }));
    }

    /// <summary>
    /// <c>ReferenceData/make_runner_batch_fixtures.R</c>: one stratum, groups A and B, batches b1-b3 (batch fitted), with a
    /// protein never measured in b3 (B09), one never in the reference batch b1 (B10), one in b2 only (B11), and one measured
    /// in b3 in group B only (B12, which still has every batch).
    /// </summary>
    private static DifferentialAnalysisResult RunBatch()
    {
        var values = Read("runner_batch_values.tsv");
        var samples = values.GroupBy(v => v["sample"]).OrderBy(g => g.Key, StringComparer.Ordinal)
            .Select(g => new AnalysisSample(g.Key)
            {
                Factors = new Dictionary<string, string> { ["group"] = g.First()["group"] },
                Batch = g.First()["batch"],
            }).ToList();
        var peptides = values.GroupBy(v => v["peptide"]).Select(g => new PeptideRunValues(g.Key, g.Key,
            g.First()["groups"].Split(';'), g.ToDictionary(v => v["sample"],
                v => new RunValue(Math.Pow(2, Num(v["y"])), ObservationState.Quantified)))).ToList();
        var groups = Enumerable.Range(1, 12).Select(i => new ProteinGroupInfo($"B{i:00}", 0, false, false));
        var runs = samples.Select(s => new ObservationRun(s.SampleId, s.SampleId)).ToList();
        return DifferentialAnalysis.Run(new DifferentialRequest
        {
            Table = ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, groups, QuantBasis.MsmsOnly),
            Samples = samples,
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("group", "A") }, FitBatch = true },
            Contrasts = new[] { DifferentialContrastRequest.Between("group", "B", "A") },
            Settings = new DifferentialSettings { Robust = false, Normalization = SharedPeptideNormalization.None },
        }, maxThreads: 1);
    }

    [Test]
    public void WithBatchesEveryRowMatchesRIncludingTheMeansOfProteinsMissingABatch()
    {
        var result = RunBatch();
        var expected = Read("runner_batch_results.tsv");
        var fitted = result.Rows.Where(r => DifferentialStatus.IsFitted(r.Status)).ToDictionary(r => r.FeatureId);

        Assert.That(fitted.Keys, Is.EquivalentTo(expected.Select(e => e["protein"])), "the runner fits exactly the rows R fits");
        foreach (var e in expected)
        {
            var r = fitted[e["protein"]];
            string at = $"{e["protein"]} (batches {e["batches_seen"]})";
            Close(r.Log2Effect!.Value, e["log2"], Effect, false, at + " log2");
            Assert.That(r.MeanNumerator, Is.Not.Null, at + " has a numerator mean");
            Close(r.MeanNumerator!.Value, e["mean_num"], Effect, false, at + " mean numerator");
            Close(r.MeanDenominator!.Value, e["mean_den"], Effect, false, at + " mean denominator");
            Close(r.StandardError!.Value, e["se"], Inference, true, at + " se");
            Close(r.Df!.Value, e["df"], Df, true, at + " df");
            Close(r.PValue!.Value, e["p"], Inference, true, at + " p");
            Assert.That(r.MeanNumerator!.Value - r.MeanDenominator!.Value, Is.EqualTo(r.Log2Effect!.Value).Within(1e-9),
                at + " effect = numerator mean - denominator mean");
        }
    }

    [Test]
    public void AProteinMissingABatchHasItsMeansAveragedOverTheBatchesItWasMeasuredIn()
    {
        var result = RunBatch();

        var over = result.Decisions.Where(d => d.Rule == DifferentialDecisionRules.MeansOverSeenLevels)
            .ToDictionary(d => d.FeatureId!, d => d.Facts);
        Assert.That(over, Is.EquivalentTo(new Dictionary<string, string>
        {
            ["B09"] = "no value in batch=b3",
            ["B10"] = "no value in batch=b1",
            ["B11"] = "no value in batch=b1, batch=b3",
        }), "only the proteins never measured in a batch, naming it; B12 was measured in b3 (group B only)");
        Assert.That(result.Decisions.Where(d => d.Rule == DifferentialDecisionRules.MeansOverSeenLevels).Select(d => d.ContrastId),
            Is.All.EqualTo("c1"));
        Assert.That(result.Decisions.Select(d => d.Rule), Has.None.EqualTo(DifferentialDecisionRules.MeansNotEstimable));
        Assert.That(DifferentialDecisionRules.Ruling(DifferentialDecisionRules.MeansOverSeenLevels), Is.EqualTo("M5b-18"));
    }

    private static void Close(double actual, string expectedText, double tolerance, bool relative, string where)
    {
        double expected = Num(expectedText);
        double scale = relative ? Math.Max(Math.Abs(expected), 1e-300) : 1;
        Assert.That(Math.Abs(actual - expected) / scale, Is.LessThanOrEqualTo(tolerance), $"{where}: {actual:R} vs {expected:R}");
    }
}
