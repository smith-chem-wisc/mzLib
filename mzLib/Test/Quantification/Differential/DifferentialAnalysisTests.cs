using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The STAT1 runner end to end on a small label-free set with known effects (A = 0, B = +1, C = +2 in log2): every
/// status rule of ruling 2 in its order, the reference-absent comparison (ruling 1), the one-peptide path (ruling 6), the
/// row set (ruling 12), and rows the results writer accepts.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialAnalysisTests
{
    private static readonly string[] Ids = { "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9" };
    private static string GroupOf(string id) => (Array.IndexOf(Ids, id) / 3) switch { 0 => "A", 1 => "B", _ => "C" };
    private static double Effect(string group) => group switch { "A" => 0, "B" => 1, _ => 2 };

    private static PeptideRunValues Peptide(Random rng, string sequence, string groups, double offset, IEnumerable<string> samples) =>
        new(sequence, sequence, groups.Split(';'), samples.ToDictionary(s => s,
            s => new RunValue(Math.Pow(2, 20 + offset + Effect(GroupOf(s)) + 0.05 * Gaussian(rng)), ObservationState.Quantified)));

    private static double Gaussian(Random rng) =>
        Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

    private static DifferentialRequest Request()
    {
        var rng = new Random(20261110);
        var all = Ids;
        var notC = Ids.Where(s => GroupOf(s) != "C").ToArray();
        var notA = Ids.Where(s => GroupOf(s) != "A").ToArray();
        var peptides = new List<PeptideRunValues>();
        for (int k = 0; k < 3; k++) peptides.Add(Peptide(rng, $"P1PEP{k}", "P1", k * 0.5, all));
        for (int k = 0; k < 3; k++) peptides.Add(Peptide(rng, $"P2PEP{k}", "P2", k * 0.5, notC));
        for (int k = 0; k < 3; k++) peptides.Add(Peptide(rng, $"P3PEP{k}", "P3", k * 0.5, notA));
        peptides.Add(Peptide(rng, "P4PEP", "P4", 0, all));
        // A second one-peptide protein so the t-test's moderation has more than one feature.
        peptides.Add(Peptide(rng, "P6PEP", "P6", 0.3, all));
        peptides.Add(Peptide(rng, "SHAREDPEP", "P5;P1", 0, all));
        peptides.Add(Peptide(rng, "DECOYPEP", "DECOY", 0, all));
        peptides.Add(Peptide(rng, "WEAKPEP", "WEAK", 0, all));
        var groups = new[]
        {
            new ProteinGroupInfo("P1", 0, false, false), new ProteinGroupInfo("P2", 0, false, false),
            new ProteinGroupInfo("P3", 0, false, false), new ProteinGroupInfo("P4", 0, false, false),
            new ProteinGroupInfo("P5", 0, false, false), new ProteinGroupInfo("P6", 0, false, false),
            new ProteinGroupInfo("DECOY", 0, true, false), new ProteinGroupInfo("WEAK", 0.05, false, false),
        };
        var runs = Ids.Select(s => new ObservationRun(s, s)).ToList();
        var table = ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, groups, QuantBasis.MsmsOnly);
        return new DifferentialRequest
        {
            Table = table,
            Samples = Ids.Select(s => new AnalysisSample(s) { Factors = new Dictionary<string, string> { ["group"] = GroupOf(s) } }).ToList(),
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("group") } },
            Contrasts = new[]
            {
                DifferentialContrastRequest.Between("group", "B", "A"), DifferentialContrastRequest.Between("group", "C", "A"),
                DifferentialContrastRequest.Between("group", "C", "B"),
            },
            Settings = new DifferentialSettings { Normalization = SharedPeptideNormalization.None },
        };
    }

    private static DifferentialResult Row(DifferentialAnalysisResult r, string feature, string contrast) =>
        r.Rows.Single(x => x.FeatureId == feature && x.ContrastId == contrast);

    [Test]
    public void EveryStatusRuleAppliesInItsOrderAndFittedEffectsRecoverTheTruth()
    {
        var result = DifferentialAnalysis.Run(Request(), maxThreads: 1);

        Assert.That(result.Rows.Select(r => r.FeatureId).Distinct(), Is.EquivalentTo(new[] { "P1", "P2", "P3", "P4", "P5", "P6" }),
            "the decoy and the group above 1% q get no rows");

        // P1: every comparison fitted, near the truth.
        Assert.That(Row(result, "P1", "c1").Status, Is.EqualTo(DifferentialStatus.Fitted));
        Assert.That(Row(result, "P1", "c1").Log2Effect, Is.EqualTo(1).Within(0.15));
        Assert.That(Row(result, "P1", "c2").Log2Effect, Is.EqualTo(2).Within(0.15));
        Assert.That(Row(result, "P1", "c3").Log2Effect, Is.EqualTo(1).Within(0.15));

        // P2 lacks C: B vs A fitted; the comparisons needing C are absent on the side that lacks it.
        Assert.That(Row(result, "P2", "c1").Status, Is.EqualTo(DifferentialStatus.Fitted));
        Assert.That(Row(result, "P2", "c1").Log2Effect, Is.EqualTo(1).Within(0.15));
        Assert.That(Row(result, "P2", "c2").Status, Is.EqualTo(DifferentialStatus.AbsentInNumerator));
        Assert.That(Row(result, "P2", "c3").Status, Is.EqualTo(DifferentialStatus.AbsentInNumerator));

        // P3 lacks the reference A: only C vs B, and it is C minus B (1), not the surviving column read as C vs A (2).
        Assert.That(Row(result, "P3", "c1").Status, Is.EqualTo(DifferentialStatus.AbsentInDenominator));
        Assert.That(Row(result, "P3", "c2").Status, Is.EqualTo(DifferentialStatus.AbsentInDenominator));
        Assert.That(Row(result, "P3", "c3").Status, Is.EqualTo(DifferentialStatus.Fitted));
        Assert.That(Row(result, "P3", "c3").Log2Effect, Is.EqualTo(1).Within(0.15));

        // P4: one usable peptide, the moderated t-test.
        Assert.That(Row(result, "P4", "c2").Status, Is.EqualTo(DifferentialStatus.FittedSinglePeptide));
        Assert.That(Row(result, "P4", "c2").ModelUsed, Is.EqualTo(DifferentialConstants.ModeratedT));
        Assert.That(Row(result, "P4", "c2").Log2Effect, Is.EqualTo(2).Within(0.2));

        // P5: only a shared peptide.
        Assert.That(Row(result, "P5", "c1").Status, Is.EqualTo(DifferentialStatus.BelowSupport(DifferentialConstants.NoUniquePeptides)));
        Assert.That(Row(result, "P5", "c1").Log2Effect, Is.Null);
    }

    /// <summary>
    /// Model means (ruling 3): each side's estimated mean averages the protein's peptides with equal weight, so P1 (peptides
    /// at +0, +0.5 and +1 over 20; groups at 0, +1 and +2) has means near 20.5, 21.5 and 22.5, P3 with its reference absent
    /// still has C near 22.5 and B near 21.5, and the one-peptide P4 near 20 and 22. On every fitted row the numerator mean
    /// minus the denominator mean is the effect; a row with no fit has no means.
    /// </summary>
    [Test]
    public void TheMeansAverageThePeptidesAndTheirDifferenceIsTheEffect()
    {
        var result = DifferentialAnalysis.Run(Request(), maxThreads: 1);

        var fitted = result.Rows.Where(r => DifferentialStatus.IsFitted(r.Status)).ToList();
        Assert.That(fitted, Has.Count.GreaterThan(5));
        foreach (var r in fitted)
            Assert.That(r.MeanNumerator!.Value - r.MeanDenominator!.Value, Is.EqualTo(r.Log2Effect!.Value).Within(1e-9),
                $"{r.FeatureId} {r.ContrastId}");
        Assert.That(result.Rows.Where(r => !DifferentialStatus.IsFitted(r.Status)),
            Is.All.Matches<DifferentialResult>(r => r.MeanNumerator is null && r.MeanDenominator is null));

        Assert.That(Row(result, "P1", "c1").MeanNumerator, Is.EqualTo(21.5).Within(0.1));
        Assert.That(Row(result, "P1", "c1").MeanDenominator, Is.EqualTo(20.5).Within(0.1));
        Assert.That(Row(result, "P1", "c2").MeanNumerator, Is.EqualTo(22.5).Within(0.1));
        Assert.That(Row(result, "P3", "c3").MeanNumerator, Is.EqualTo(22.5).Within(0.1));
        Assert.That(Row(result, "P3", "c3").MeanDenominator, Is.EqualTo(21.5).Within(0.1));
        Assert.That(Row(result, "P4", "c2").MeanNumerator, Is.EqualTo(22).Within(0.1));
        Assert.That(Row(result, "P4", "c2").MeanDenominator, Is.EqualTo(20).Within(0.1));
    }

    /// <summary>
    /// A protein never seen in batch b2 still compares its groups (batch b2 is dropped from its design, ruling 1), but a mean
    /// at the stratum's batch proportions cannot come from its data: both means are left empty, the effect is reported,
    /// and the decisions record says so for each of its contrasts. A protein seen in both batches has its means.
    /// </summary>
    [Test]
    public void AProteinMissingFromABatchKeepsItsEffectAndItsEmptyMeansAreRecorded()
    {
        var rng = new Random(7);
        string BatchOf(string id) => Array.IndexOf(Ids, id) % 3 == 0 ? "b2" : "b1";
        var inB1 = Ids.Where(s => BatchOf(s) == "b1").ToArray();
        var peptides = new List<PeptideRunValues>();
        foreach (var protein in new[] { "FULL1", "FULL2" })
            for (int k = 0; k < 3; k++) peptides.Add(Peptide(rng, $"{protein}PEP{k}", protein, k * 0.5, Ids));
        for (int k = 0; k < 3; k++) peptides.Add(Peptide(rng, $"ONLYB1PEP{k}", "ONLYB1", k * 0.5, inB1));
        var groups = new[] { "FULL1", "FULL2", "ONLYB1" }.Select(g => new ProteinGroupInfo(g, 0, false, false)).ToArray();
        var runs = Ids.Select(s => new ObservationRun(s, s)).ToList();
        var request = Request() with
        {
            Table = ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, groups, QuantBasis.MsmsOnly),
            Samples = Ids.Select(s => new AnalysisSample(s)
            {
                Factors = new Dictionary<string, string> { ["group"] = GroupOf(s) },
                Batch = BatchOf(s),
            }).ToList(),
        };

        var result = DifferentialAnalysis.Run(request, maxThreads: 1);

        var only = result.Rows.Where(r => r.FeatureId == "ONLYB1").ToList();
        Assert.That(only.Select(r => r.Status), Is.All.EqualTo(DifferentialStatus.Fitted));
        Assert.That(only.Select(r => r.Log2Effect), Is.All.Not.Null);
        Assert.That(only, Is.All.Matches<DifferentialResult>(r => r.MeanNumerator is null && r.MeanDenominator is null));
        var recorded = result.Decisions.Where(d => d.Rule == DifferentialDecisionRules.MeansNotEstimable).ToList();
        Assert.That(recorded.Select(d => (d.FeatureId, d.ContrastId)),
            Is.EquivalentTo(new[] { ("ONLYB1", "c1"), ("ONLYB1", "c2"), ("ONLYB1", "c3") }));
        Assert.That(recorded[0].Facts, Does.Contain("batch=b2"));
        Assert.That(result.Rows.Where(r => r.FeatureId == "FULL1"), Is.All.Matches<DifferentialResult>(r => r.MeanNumerator is not null));
    }

    [Test]
    public void TheRowsPassTheWritersValidationAndOneFamilyCoversBothModels()
    {
        var result = DifferentialAnalysis.Run(Request(), maxThreads: 1);

        Assert.DoesNotThrow(() => DifferentialResultWriter.Write(TextWriter.Null, result.Rows, DifferentialHeaderStyle.Machine));
        var c2 = result.Rows.Where(r => r.ContrastId == "c2" && DifferentialStatus.IsFitted(r.Status)).ToList();
        Assert.That(c2.Select(r => r.ModelUsed).Distinct(), Is.EquivalentTo(new[] { DifferentialConstants.MixedModel, DifferentialConstants.ModeratedT }));
        Assert.That(c2.Select(r => r.FamilySize), Is.All.EqualTo(c2.Count));
    }

    [Test]
    public void TheDecisionsRecordNamesEachRuleThatChangedSomething()
    {
        var result = DifferentialAnalysis.Run(Request(), maxThreads: 1);

        Assert.That(result.Decisions.Select(d => d.Rule), Is.All.Matches<string>(DifferentialDecisionRules.IsValid));
        var rowSet = result.Decisions.Single(d => d.Rule == DifferentialDecisionRules.RowSet);
        Assert.That(rowSet.Facts, Does.Contain("decoy 1").And.Contain("above_1pct_q 1"));
        Assert.That(result.Decisions.Any(d => d.Rule == DifferentialDecisionRules.ReducedDesign && d.FeatureId == "P3"), Is.True);
        Assert.That(result.Decisions.Any(d => d.Rule == DifferentialDecisionRules.AbsentSide && d.FeatureId == "P2" && d.ContrastId == "c2"), Is.True);
        Assert.That(result.Decisions.Any(d => d.Rule == DifferentialDecisionRules.SinglePeptide && d.FeatureId == "P4"), Is.True);
        Assert.That(result.Decisions.Any(d => d.Rule == DifferentialDecisionRules.NoUniquePeptides && d.FeatureId == "P5"), Is.True);
    }

    [Test]
    public void TheOutputDoesNotDependOnTheThreadCount()
    {
        string Written(int threads)
        {
            var text = new StringWriter();
            DifferentialResultWriter.Write(text, DifferentialAnalysis.Run(Request(), threads).Rows, DifferentialHeaderStyle.Machine);
            return text.ToString();
        }

        // Compared as written: records hold lists, which compare by reference.
        Assert.That(Written(4), Is.EqualTo(Written(1)));
    }
}
