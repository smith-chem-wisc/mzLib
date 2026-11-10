using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
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

    internal static DifferentialRequest Request()
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

        // covariates_fitted: the group is the comparison's own factor, so a full peptide model shows its random sample term
        // only; the t-test row and a row with no fit show nothing. design_sha256: one stratum, one hash, on every row.
        Assert.That(Row(result, "P1", "c1").CovariatesFitted, Is.EqualTo("sample(random)"));
        Assert.That(Row(result, "P4", "c2").CovariatesFitted, Is.Null);
        Assert.That(Row(result, "P5", "c1").CovariatesFitted, Is.Null);
        Assert.That(result.Rows.Select(r => r.DesignSha256).Distinct().Single(), Does.Match("^[0-9a-f]{64}$"));
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
        Assert.That(result.Rows.Where(r => r.FeatureId == "FULL1").Select(r => r.CovariatesFitted), Is.All.EqualTo("batch;sample(random)"));
        Assert.That(only.Select(r => r.CovariatesFitted), Is.All.EqualTo("sample(random)"), "batch left out of its own fit");
    }

    /// <summary>
    /// robust_fell_back (ruling 8, DATA-DEFINITIONS v3.21): a robust mixed-model row says whether its reweighting settled
    /// (here every one did: false); a moderated_t row is never reweighted and an evidence-only row has no fit, so both leave
    /// it empty; with robust weighting off no row carries it.
    /// </summary>
    [Test]
    public void RobustFellBackIsSetOnRobustMixedModelRowsOnly()
    {
        var result = DifferentialAnalysis.Run(Request(), maxThreads: 1);

        var mixed = result.Rows.Where(r => r.Status == DifferentialStatus.Fitted).ToList();
        Assert.That(mixed, Is.Not.Empty);
        Assert.That(mixed.Select(r => r.RobustFellBack), Is.All.EqualTo(false));
        Assert.That(result.Rows.Where(r => r.Status != DifferentialStatus.Fitted).Select(r => r.RobustFellBack), Is.All.Null);
        Assert.That(result.Decisions.Any(d => d.Rule == DifferentialDecisionRules.RobustFellBack), Is.False);

        var plain = DifferentialAnalysis.Run(Request() with { Settings = Request().Settings with { Robust = false } }, maxThreads: 1);
        Assert.That(plain.Rows.Select(r => r.RobustFellBack), Is.All.Null);
        Assert.That(plain.Rows.Where(r => r.Status == DifferentialStatus.Fitted).Select(r => r.Robust), Is.All.EqualTo(false));
    }

    /// <summary>
    /// Ruling 16: a protein with two peptides but one value per sample (X in s1, s2, s4, s7; Y in the rest) cannot have its
    /// sample term estimated, so it is fitted with no random term (least squares, residual df) instead of being refused, and
    /// the decisions record says so. It is moderated in a family of its own (ruling 17): alone there, it has no prior.
    /// </summary>
    [Test]
    public void AProteinWithOneValuePerSampleIsFittedWithoutItsSampleTerm()
    {
        var request = Request();
        var rng = new Random(5);
        var peptides = new List<PeptideRunValues>
        {
            Peptide(rng, "SPLITX", "SPLIT", 0, new[] { "s1", "s2", "s4", "s7" }),
            Peptide(rng, "SPLITY", "SPLIT", 0.5, new[] { "s3", "s5", "s6", "s8", "s9" }),
        };
        var runs = Ids.Select(s => new ObservationRun(s, s)).ToList();
        var table = ObservationTable.Build(runs, runs.Select(r => r.FileName),
            request.Table!.Peptides.Select((pep, p) => new PeptideRunValues(pep.FullSequence, pep.BaseSequence, pep.ProteinGroups,
                Ids.Select((id, s) => (id, v: request.Table.Log2Intensity(p, s))).Where(t => !double.IsNaN(t.v))
                    .ToDictionary(t => t.id, t => new RunValue(Math.Pow(2, t.v), ObservationState.Quantified)))).Concat(peptides),
            request.Table.ProteinGroups.Values.Append(new ProteinGroupInfo("SPLIT", 0, false, false)), QuantBasis.MsmsOnly);

        var result = DifferentialAnalysis.Run(request with { Table = table }, maxThreads: 1);

        var split = result.Rows.Where(r => r.FeatureId == "SPLIT").ToList();
        Assert.That(split.Select(r => r.Status), Is.All.EqualTo(DifferentialStatus.Fitted));
        Assert.That(split.Select(r => r.ModelUsed), Is.All.EqualTo(DifferentialConstants.MixedModelReduced));
        Assert.That(Row(result, "P1", "c1").ModelUsed, Is.EqualTo(DifferentialConstants.MixedModel));
        var reducedModel = result.Metadata.Models.Single(m => m.ModelUsed == DifferentialConstants.MixedModelReduced);
        Assert.That(reducedModel.Prior!.Df, Is.EqualTo(0), "one protein in the family: no prior (ruling 5)");
        Assert.That(reducedModel.Formula, Does.Contain("only the random terms each protein's values can estimate"));
        Assert.That(result.Metadata.Models.Single(m => m.ModelUsed == DifferentialConstants.MixedModel).Prior!.Features, Is.EqualTo(3),
            "P1, P2 and P3, not SPLIT");
        Assert.That(split.Select(r => r.DfMethod), Is.All.EqualTo(DifferentialConstants.DfResidual));
        Assert.That(Row(result, "SPLIT", "c2").Log2Effect, Is.EqualTo(2).Within(0.2));
        Assert.That(Row(result, "P1", "c1").DfMethod, Is.EqualTo(DifferentialConstants.DfSatterthwaite));
        var dropped = result.Decisions.Single(d => d.Rule == DifferentialDecisionRules.RandomTermsDropped);
        Assert.That(dropped.FeatureId, Is.EqualTo("SPLIT"));
        Assert.That(dropped.Facts, Does.Contain("(1|sample)").And.Contain("9 value(s)"));
        Assert.That(dropped.Effect, Does.Contain("no random term"));
        Assert.That(split.Select(r => r.CovariatesFitted), Is.All.Null, "no term besides the group, no random term");
    }

    /// <summary>
    /// Ruling 14: the fewest samples per side is a setting (default 2), named in the status, written into the metadata and
    /// so into the analysis id. With 4, every comparison of these 3-sample groups is below support. Fewer than 2 is refused.
    /// </summary>
    [Test]
    public void TheMinimumSamplesPerSideIsASetting()
    {
        var request = Request();
        var standard = DifferentialAnalysis.Run(request, maxThreads: 1);
        var four = DifferentialAnalysis.Run(request with { Settings = request.Settings with { MinimumSamplesPerSide = 4 } }, maxThreads: 1);

        Assert.That(Row(standard, "P1", "c1").Status, Is.EqualTo(DifferentialStatus.Fitted));
        Assert.That(Row(four, "P1", "c1").Status, Is.EqualTo(DifferentialStatus.BelowSupport("min_4_per_side")));
        Assert.That(four.Rows.Where(r => r.FeatureId == "P1").Select(r => r.Status), Is.All.EqualTo("below_support:min_4_per_side"));
        Assert.That(four.Decisions.Any(d => d.Rule == DifferentialDecisionRules.MinPerSide && d.FeatureId == "P1"), Is.True);
        Assert.That(four.Metadata.AnalysisId, Is.Not.EqualTo(standard.Metadata.AnalysisId));
        Assert.That(Encoding.UTF8.GetString(DifferentialMetadataWriter.ToBytes(four.Metadata, four.Rows)),
            Does.Contain("\"minimum_samples_per_side\": 4").And.Contain("\"normalization_floor\": 100"));
        Assert.Throws<ArgumentException>(() => DifferentialAnalysis.Run(
            request with { Settings = request.Settings with { MinimumSamplesPerSide = 1 } }, maxThreads: 1));
        Assert.Throws<ArgumentException>(() => DifferentialAnalysis.Run(
            request with { Settings = request.Settings with { NormalizationFloor = 0 } }, maxThreads: 1));
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
