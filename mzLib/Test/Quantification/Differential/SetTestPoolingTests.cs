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
/// SETTEST1 S3 (QuantProject <c>DEF-DIFF-SET-POOL</c>, <c>-FAMILY</c>, <c>-POSTHOC</c>, <c>-MATURATION</c>): which dataset
/// contrasts build a part's primary result; the pooled row (DerSimonian–Laird τ², Hartung–Knapp standard error on k − 1
/// df); the aging and maturation families; Benjamini–Bogomolov on the pooled datasets of a supported part; and the
/// "supported in at least 2" p. One scenario against R (<c>ReferenceData/make_settest_pooling_fixtures.R</c>, metafor and
/// base R, never run here), then one hand-made case per rule.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SetTestPoolingTests
{
    private const double Tight = 1e-8;

    private const string Header =
        "hypothesis_id\tkind\tprediction\tmember\tspecies\tdataset\tfactor\ttest_level\treference_level\tcontexts\tmrna_log2" +
        "\tindependence\tpost_hoc\tprimary\tnote";

    // ── fixtures ─────────────────────────────────────────────────────────────

    private static (string[] Header, List<string[]> Rows) Tsv(string file)
    {
        var lines = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential",
            "ReferenceData", file)).Where(l => l.Length > 0).ToArray();
        return (lines[0].Split('\t'), lines.Skip(1).Select(l => l.Split('\t')).ToList());
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static void Close(double actual, double expected, string what)
    {
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(Tight * Math.Max(Math.Abs(expected), 1e-300)),
            $"{what}: {actual:R} vs {expected:R}");
    }

    private static SetTestFreeze Freeze(params string[] rows) =>
        SetTestFreeze.Read(Encoding.UTF8.GetBytes(string.Join('\n', new[] { Header }.Concat(rows)) + "\n"), "t");

    private static string Row(string part, string prediction, string dataset, string independence = "independent",
        string postHoc = "no", string primary = "", string note = "", string test = "old", string member = "",
        string mrna = "") =>
        string.Join('\t', part, "single", prediction, member == "" ? part + "_m" : member, "mouse", dataset, "age group", test,
            "young", "", mrna, independence, postHoc, primary, note);

    private static SetHypothesisPart Part(SetTestFreeze freeze, string id) => freeze.Parts.Single(p => p.HypothesisId == id);

    private static SetDatasetContrast Contrast(SetHypothesisPart part, string dataset, int index = 0) =>
        part.Contrasts.Where(c => c.Dataset == dataset).ElementAt(index);

    private static SetEstimate Estimate(double b, double se, double df) => new()
    {
        Log2Effect = b, StandardError = se, Df = df, DfMethod = "satterthwaite", ModelUsed = "set_mixed_model",
    };

    private static SetContrastOutcome Tested(SetHypothesisPart part, SetDatasetContrast contrast, double b, double se,
        double df = 10, double? margin = null, double? mrna = null, params string[] blockedBy) => new()
    {
        Part = part,
        Contrast = contrast,
        Setup = new SetContrastSetup(),
        Test = SetVerdict.Test(part.Prediction, Estimate(b, se, df), margin, null, mrna),
        BlockedBy = blockedBy,
    };

    private static SetContrastOutcome Untested(SetHypothesisPart part, SetDatasetContrast contrast, string verdict, string reason) => new()
    {
        Part = part, Contrast = contrast, Setup = new SetContrastSetup(), Untested = verdict, UntestedReason = reason,
    };

    private static SetDatasetResult Row(SetPoolingResult result, string part, string dataset, int index = 0) =>
        result.Parts.Single(p => p.Part.HypothesisId == part).Datasets.Where(d => d.Outcome.Contrast.Dataset == dataset).ElementAt(index);

    private static SetPooledResult Pooled(SetPoolingResult result, string part) =>
        result.Parts.Single(p => p.Part.HypothesisId == part).Pooled!;

    // ── one scenario against R ───────────────────────────────────────────────

    /// <summary>
    /// Seven parts: HA up on 3 datasets, HB down on 2, HD unchanged on 2 and HE attenuated_down on 2 pool; HC not_up has one
    /// dataset; HG up keeps one of its two (the other reproduces); HF's one contrast is a maturation contrast. R gives each
    /// pooled row, both families' adjusted p's, R = 3 of m = 6, the level 0.025, and HA's and HD's datasets under
    /// Benjamini–Bogomolov (HA's third, at 0.035, fails the level but would pass at 0.05). HE's m̄ is the pooling-weight
    /// average of its datasets' (ruled by the user 2026-11-14).
    /// </summary>
    [Test]
    public void TheScenario_MatchesR()
    {
        var (dh, datasets) = Tsv("settest_pooling_datasets.tsv");
        string Cell(string[] r, string c) => r[Array.IndexOf(dh, c)];
        var freeze = Freeze(datasets.Select(r => Row(Cell(r, "part"), Cell(r, "prediction"), Cell(r, "dataset"),
            independence: Cell(r, "independence"),
            note: Cell(r, "kind") == "maturation" ? "MATURATION CONTRAST: P30 against 12 months" : "",
            mrna: Cell(r, "mrna"))).ToArray());
        var outcomes = datasets.Select(r =>
        {
            var part = Part(freeze, Cell(r, "part"));
            var o = Tested(part, Contrast(part, Cell(r, "dataset")), D(Cell(r, "b")), D(Cell(r, "se")), D(Cell(r, "df")),
                mrna: Cell(r, "mrna") == "" ? null : D(Cell(r, "mrna")));
            Close(o.Test!.PValue, D(Cell(r, "p")), $"{Cell(r, "part")} {Cell(r, "dataset")} p (S2's rules, against R)");
            return o;
        }).ToList();

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes);

        Assert.That(result.Freeze, Is.EqualTo("t"));
        Assert.That(result.TestsSha256, Is.EqualTo(freeze.Sha256));
        Assert.That(result.Basis, Is.EqualTo(QuantBasis.MsmsOnly));

        var (ph, pooledRows) = Tsv("settest_pooling_pooled.tsv");
        Assert.That(result.Parts.Count(p => p.Pooled is not null), Is.EqualTo(pooledRows.Count));
        foreach (var r in pooledRows)
        {
            string part = r[Array.IndexOf(ph, "part")];
            double Out(string c) => D(r[Array.IndexOf(ph, c)]);
            var pooled = Pooled(result, part);
            Assert.That(pooled.Datasets, Has.Count.EqualTo((int)Out("k")), part);
            Assert.That(pooled.Estimate!.DfMethod, Is.EqualTo(SetPoolingName.PooledDfMethod));
            Assert.That(pooled.Estimate.ModelUsed, Is.EqualTo(SetPoolingName.PooledModel));
            Assert.That(pooled.Estimate.Df, Is.EqualTo(Out("df")), part);
            Close(pooled.Estimate.Log2Effect, Out("estimate"), $"{part} estimate");
            Close(pooled.Estimate.StandardError, Out("se"), $"{part} se");
            Assert.That(pooled.Tau2, Is.EqualTo(Out("tau2")).Within(1e-12), $"{part} tau2");
            Close(pooled.Test!.PValue, Out("p"), $"{part} p");
            Close(pooled.Test.Ci90Low, Out("ci_low"), $"{part} ci low");
            Close(pooled.Test.Ci90High, Out("ci_high"), $"{part} ci high");
            Assert.That(pooled.Test.Contradicted, Is.EqualTo(Out("contradicted") == 1), $"{part} contradicted");
            Assert.That(pooled.DatasetsAgree, Is.EqualTo((int)Out("n_agree")), $"{part} n_datasets_agree");
            Close(pooled.PSupportedInTwo, Out("p_in_2"), $"{part} p supported in 2");
            if (r[Array.IndexOf(ph, "mrna_mean")] != "") Close(pooled.Test.MrnaMean!.Value, Out("mrna_mean"), $"{part} m̄");
        }

        var (fh, family) = Tsv("settest_pooling_family.tsv");
        foreach (var r in family)
        {
            string part = r[Array.IndexOf(fh, "part")], dataset = r[Array.IndexOf(fh, "dataset")];
            string fam = r[Array.IndexOf(fh, "family")];
            bool supported = r[Array.IndexOf(fh, "supported")] == "1";
            if (dataset == "")
            {
                var pooled = Pooled(result, part);
                Assert.That(pooled.Family, Is.EqualTo(fam), part);
                Assert.That(pooled.FamilySize, Is.EqualTo(6), part);
                Close(pooled.PAdjusted!.Value, D(r[Array.IndexOf(fh, "p_adjusted")]), $"{part} pooled adjusted p");
                Assert.That(pooled.Verdict == SetVerdictName.Supported, Is.EqualTo(supported), $"{part} pooled verdict");
            }
            else
            {
                var row = Row(result, part, dataset);
                Assert.That(row.IsPrimary && row.Family == fam, $"{part} {dataset} is a primary of {fam}");
                Assert.That(row.ContrastKind, Is.EqualTo(fam));
                Close(row.PAdjusted!.Value, D(r[Array.IndexOf(fh, "p_adjusted")]), $"{part} {dataset} adjusted p");
                Assert.That(row.Verdict == SetVerdictName.Supported, Is.EqualTo(supported), $"{part} {dataset} verdict");
            }
        }

        var (sh, summary) = Tsv("settest_pooling_summary.tsv");
        Assert.That(result.AgingFamilySize, Is.EqualTo((int)D(summary[0][Array.IndexOf(sh, "aging_family_size")])));
        Assert.That(result.AgingSupported, Is.EqualTo((int)D(summary[0][Array.IndexOf(sh, "aging_supported")])));
        Close(result.BogomolovLevel, D(summary[0][Array.IndexOf(sh, "bb_level")]), "Benjamini–Bogomolov level");
        Assert.That(result.MaturationFamilySize, Is.EqualTo((int)D(summary[0][Array.IndexOf(sh, "maturation_family_size")])));

        var (bh, bogomolov) = Tsv("settest_pooling_bogomolov.tsv");
        Assert.That(bogomolov, Has.Count.EqualTo(5));
        foreach (var r in bogomolov)
        {
            string part = r[Array.IndexOf(bh, "part")], dataset = r[Array.IndexOf(bh, "dataset")];
            var row = Row(result, part, dataset);
            double level = D(r[Array.IndexOf(bh, "level")]), adjusted = D(r[Array.IndexOf(bh, "adjusted")]);
            Assert.That(row.Pooled && !row.InFamily, $"{part} {dataset} is pooled, outside the family");
            Close(row.PAdjusted!.Value, adjusted, $"{part} {dataset} Benjamini–Bogomolov adjusted p");
            Close(row.BogomolovLevel!.Value, level, $"{part} {dataset} level");
            Assert.That((row.Verdict, row.VerdictDetail), Is.EqualTo(adjusted <= level
                ? (SetVerdictName.Supported, SetPoolingName.BogomolovPassed)
                : (SetVerdictName.NotInformative, SetPoolingName.BogomolovFailed)), $"{part} {dataset}");
        }

        // Pooled datasets of a part that is not supported carry no adjusted p.
        foreach (var dataset in new[] { "D1", "D4" })
        {
            var row = Row(result, "HB", dataset);
            Assert.That(row.Pooled && row.PAdjusted is null && row.BogomolovLevel is null);
            Assert.That((row.Verdict, row.VerdictDetail), Is.EqualTo((SetVerdictName.NotInformative, "not_in_family")));
        }
        Assert.That(Row(result, "HG", "D2").LeftOut, Is.EqualTo(SetPoolingName.Reproduction));
        Assert.That(Row(result, "HG", "D2").VerdictDetail, Is.EqualTo(SetPoolingName.Reproduction));
    }

    // ── which contrasts build a primary result ───────────────────────────────

    [Test]
    public void EachContrastLeftOut_SaysWhy_AndBuildsNothing()
    {
        var freeze = Freeze(
            Row("H1", "up", "D1"),
            Row("H1", "up", "D2", independence: "probably reproduces [102]"),
            Row("H1", "up", "D3", independence: "weak (not a test)"),
            Row("H1", "up", "D4", postHoc: "yes"),
            Row("H1", "up", "D5"), Row("H1", "up", "D5", test: "middle"));
        var part = Part(freeze, "H1");
        var outcomes = part.Contrasts.Select(c => Tested(part, c, 1.0, 0.2)).ToList();

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes);

        var rows = result.Parts.Single().Datasets;
        Assert.That(rows.Select(r => r.LeftOut), Is.EqualTo(new[]
        {
            null, SetPoolingName.Reproduction, SetPoolingName.NotPooled, SetPoolingName.PostHocPrefix + "H1_m",
            SetPoolingName.SeveralContrastsPerDataset, SetPoolingName.SeveralContrastsPerDataset,
        }));
        Assert.That(rows.Skip(1).Select(r => r.VerdictDetail), Is.EqualTo(rows.Skip(1).Select(r => r.LeftOut)));
        Assert.That(rows.Skip(1).All(r => !r.Pooled && !r.IsPrimary && !r.InFamily && r.PAdjusted is null));
        Assert.That(result.Parts.Single().Pooled, Is.Null, "one contrast is left, so it is the primary itself");
        Assert.That(rows[0].IsPrimary && rows[0].Family == SetPoolingName.Aging && rows[0].FamilySize == 1);
        Assert.That(result.AgingFamilySize, Is.EqualTo(1));
    }

    /// <summary>A post-hoc contrast names only the members seen before the freeze (v3.13).</summary>
    [Test]
    public void APostHocContrast_NamesTheMembersSeenBeforeTheFreeze()
    {
        var freeze = Freeze(Row("H1", "up", "D1", member: "A"), Row("H1", "up", "D1", member: "B", postHoc: "yes"),
            Row("H1", "up", "D1", member: "C", postHoc: "yes"));
        var part = Part(freeze, "H1");

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, new[] { Tested(part, part.Contrasts.Single(), 2.0, 0.1) });

        var row = result.Parts.Single().Datasets.Single();
        Assert.That(row.LeftOut, Is.EqualTo(SetPoolingName.PostHocPrefix + "B;C"));
        Assert.That((row.Verdict, row.InFamily), Is.EqualTo((SetVerdictName.NotInformative, false)));
        Assert.That(result.AgingFamilySize, Is.EqualTo(0));
    }

    [Test]
    public void TheContrastMarkedPrimary_IsTheOneItsDatasetPools()
    {
        var freeze = Freeze(
            Row("H1", "up", "D1", primary: "no"), Row("H1", "up", "D1", primary: "yes", test: "middle"),
            Row("H1", "up", "D2"),
            Row("H1", "up", "D3", primary: "no"));
        var part = Part(freeze, "H1");
        var outcomes = part.Contrasts.Select((c, i) => Tested(part, c, 1.0 + 0.1 * i, 0.2)).ToList();

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes);

        Assert.That(Row(result, "H1", "D1", 0).LeftOut, Is.EqualTo(SetPoolingName.NotPrimary));
        Assert.That(Row(result, "H1", "D1", 1).Pooled);
        Assert.That(Row(result, "H1", "D2").Pooled);
        Assert.That(Row(result, "H1", "D3").LeftOut, Is.EqualTo(SetPoolingName.NotPrimary), "a lone contrast marked no");
        Assert.That(Pooled(result, "H1").Datasets.Select(d => d.Contrast),
            Is.EqualTo(new[] { Contrast(part, "D1", 1), Contrast(part, "D2") }));
    }

    [Test]
    public void TwoContrastsOfOneDatasetBothMarkedPrimary_AreBothLeftOut()
    {
        var freeze = Freeze(Row("H1", "up", "D1", primary: "yes"), Row("H1", "up", "D1", primary: "yes", test: "middle"));
        var part = Part(freeze, "H1");

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, part.Contrasts.Select(c => Tested(part, c, 1.0, 0.2)).ToList());

        Assert.That(result.Parts.Single().Datasets.Select(d => d.LeftOut),
            Is.EqualTo(new[] { SetPoolingName.SeveralContrastsPerDataset, SetPoolingName.SeveralContrastsPerDataset }));
        Assert.That(result.AgingFamilySize, Is.EqualTo(0));
    }

    [Test]
    public void AMaturationContrast_NeverPools_AndIsItsOwnPrimaryInItsOwnFamily()
    {
        var freeze = Freeze(Row("H1", "down", "D1"), Row("H1", "down", "D2"),
            Row("H1", "down", "D2", test: "P30", note: "MATURATION CONTRAST: juvenile"),
            Row("H1", "down", "D3", note: "MATURATION CONTRAST: juvenile"));
        var part = Part(freeze, "H1");
        var outcomes = part.Contrasts.Select((c, i) => Tested(part, c, -1.5 - 0.1 * i, 0.2)).ToList();

        var result = SetPooling.Run(freeze, QuantBasis.MbrKept, outcomes);

        Assert.That(Pooled(result, "H1").Datasets.Select(d => d.Contrast.Dataset), Is.EqualTo(new[] { "D1", "D2" }));
        var maturation = result.Parts.Single().Datasets.Where(d => d.ContrastKind == SetPoolingName.Maturation).ToList();
        Assert.That(maturation, Has.Count.EqualTo(2), "D2's maturation contrast is not 'several per dataset' with its aging one");
        Assert.That(maturation.All(d => d.IsPrimary && !d.Pooled && d.Family == SetPoolingName.Maturation && d.FamilySize == 2));
        Assert.That(result.MaturationFamilySize, Is.EqualTo(2));
        Assert.That(result.AgingFamilySize, Is.EqualTo(1), "only the pooled aging row");
        Assert.That(maturation.All(d => d.Verdict == SetVerdictName.Supported));
        Assert.That(result.AgingSupported, Is.EqualTo(1), "a maturation contrast never counts toward R");
    }

    [Test]
    public void AnUntestedContrast_LeavesThePool_SoTheOneLeftIsThePrimary()
    {
        var freeze = Freeze(Row("H1", "up", "D1"), Row("H1", "up", "D2"));
        var part = Part(freeze, "H1");
        var outcomes = new[]
        {
            Tested(part, Contrast(part, "D1"), 1.0, 0.2),
            Untested(part, Contrast(part, "D2"), SetVerdictName.NotTestable, "fewer_than_3_tested"),
        };

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes);

        Assert.That(result.Parts.Single().Pooled, Is.Null);
        Assert.That(Row(result, "H1", "D1").IsPrimary && Row(result, "H1", "D1").InFamily);
        Assert.That((Row(result, "H1", "D2").Verdict, Row(result, "H1", "D2").VerdictDetail),
            Is.EqualTo((SetVerdictName.NotTestable, "fewer_than_3_tested")));
        Assert.That(Row(result, "H1", "D2").LeftOut, Is.Null, "nothing in the freeze leaves it out; it was just not tested");
    }

    [Test]
    public void APartWithNothingTested_HasNoPrimary_AndNothingInTheFamily()
    {
        var freeze = Freeze(Row("H1", "up", "D1"), Row("H2", "up", "D1"));
        var h1 = Part(freeze, "H1");
        var h2 = Part(freeze, "H2");

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, new[]
        {
            Untested(h1, h1.Contrasts.Single(), SetVerdictName.NotRun, "spans_strata"),
            Tested(h2, h2.Contrasts.Single(), 1.0, 0.2),
        });

        Assert.That(result.AgingFamilySize, Is.EqualTo(1));
        Assert.That(Row(result, "H1", "D1").InFamily, Is.False);
        Assert.That(Row(result, "H2", "D1").PAdjusted, Is.EqualTo(Row(result, "H2", "D1").Outcome.Test!.PValue).Within(1e-15),
            "a family of one adjusts nothing");
    }

    // ── the pooled row ───────────────────────────────────────────────────────

    /// <summary>
    /// Ruled by the user 2026-11-14: a guard on any pooled dataset blocks the pooled row's <c>supported</c>, naming
    /// <c>dataset:member</c>. The row stays in the family (it has a p), but does not count toward R.
    /// </summary>
    [Test]
    public void AGuardOnAPooledDataset_BlocksThePooledRow()
    {
        var freeze = Freeze(Row("H1", "unchanged", "D1"), Row("H1", "unchanged", "D2"));
        var part = Part(freeze, "H1");
        var outcomes = new[]
        {
            Tested(part, Contrast(part, "D1"), 0.01, 0.05),
            Tested(part, Contrast(part, "D2"), -0.01, 0.05, blockedBy: new[] { "Col1a1", "Col3a1" }),
        };

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes);

        var pooled = Pooled(result, "H1");
        Assert.That(pooled.PAdjusted, Is.LessThan(0.05), "the interval alone would pass");
        Assert.That((pooled.Verdict, pooled.VerdictDetail),
            Is.EqualTo((SetVerdictName.NotInformative, "blocked_by:D2:Col1a1;D2:Col3a1")));
        Assert.That(pooled.Family, Is.EqualTo(SetPoolingName.Aging));
        Assert.That(result.AgingSupported, Is.EqualTo(0));
        Assert.That(result.Parts.Single().Datasets.All(d => d.BogomolovLevel is null));
    }

    [Test]
    public void APooledIntervalOnTheRefutedSide_IsContradicted_Uncorrected()
    {
        var freeze = Freeze(Row("H1", "up", "D1"), Row("H1", "up", "D2"), Row("H1", "up", "D3"));
        var part = Part(freeze, "H1");
        var outcomes = new[] { -1.2, -1.0, -1.1 }.Select((b, i) => Tested(part, part.Contrasts[i], b, 0.1)).ToList();

        var pooled = Pooled(SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes), "H1");

        Assert.That((pooled.Verdict, pooled.VerdictDetail), Is.EqualTo((SetVerdictName.Contradicted, "uncorrected")));
        Assert.That(pooled.DatasetsAgree, Is.EqualTo(0));
    }

    [Test]
    public void EqualEstimates_LeaveThePooledRowNotRun_AndOutOfTheFamily()
    {
        var freeze = Freeze(Row("H1", "up", "D1"), Row("H1", "up", "D2"));
        var part = Part(freeze, "H1");

        var result = SetPooling.Run(freeze, QuantBasis.MsmsOnly,
            part.Contrasts.Select((c, i) => Tested(part, c, 0.8, 0.1 + 0.1 * i)).ToList());

        var pooled = Pooled(result, "H1");
        Assert.That((pooled.Verdict, pooled.VerdictDetail), Is.EqualTo((SetVerdictName.NotRun, SetPoolingName.PooledStandardErrorZero)));
        Assert.That(pooled.Estimate is null && pooled.Test is null && pooled.Family is null && pooled.PAdjusted is null);
        Assert.That(result.AgingFamilySize, Is.EqualTo(0));
    }

    [Test]
    public void PooledDatasetsTestingDifferentMargins_AreRefused()
    {
        var freeze = Freeze(Row("H1", "unchanged", "D1"), Row("H1", "unchanged", "D2"));
        var part = Part(freeze, "H1");
        var outcomes = new[]
        {
            Tested(part, Contrast(part, "D1"), 0.0, 0.1, margin: 0.5),
            Tested(part, Contrast(part, "D2"), 0.0, 0.2, margin: 0.3),
        };

        Assert.That(() => SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes),
            Throws.ArgumentException.With.Message.Contains("margins are frozen per hypothesis"));
    }

    [Test]
    public void SeenOneSideCounts_AreSummedOverThePooledDatasets()
    {
        var freeze = Freeze(Row("H1", "up", "D1"), Row("H1", "up", "D2"));
        var part = Part(freeze, "H1");
        var outcomes = new[]
        {
            Tested(part, Contrast(part, "D1"), 1.0, 0.2) with { OneSideAgree = 2, OneSideDisagree = 1 },
            Tested(part, Contrast(part, "D2"), 1.4, 0.3) with { OneSideAgree = 1 },
        };

        var pooled = Pooled(SetPooling.Run(freeze, QuantBasis.MsmsOnly, outcomes), "H1");

        Assert.That((pooled.OneSideAgree, pooled.OneSideDisagree), Is.EqualTo((3, 1)));
    }

    // ── refusals ─────────────────────────────────────────────────────────────

    [Test]
    public void AnOutcomeFromAnotherFreeze_OrSuppliedTwice_IsRefused()
    {
        var freeze = Freeze(Row("H1", "up", "D1"));
        var other = Freeze(Row("H1", "up", "D1"));
        var part = Part(freeze, "H1");
        var stranger = Part(other, "H1");

        Assert.That(() => SetPooling.Run(freeze, QuantBasis.MsmsOnly, new[] { Tested(stranger, stranger.Contrasts.Single(), 1, 0.2) }),
            Throws.ArgumentException.With.Message.Contains("is not a part of freeze t"));
        Assert.That(() => SetPooling.Run(freeze, QuantBasis.MsmsOnly, new[] { Tested(part, stranger.Contrasts.Single(), 1, 0.2) }),
            Throws.ArgumentException.With.Message.Contains("is not a contrast of H1 up"));
        var once = Tested(part, part.Contrasts.Single(), 1, 0.2);
        Assert.That(() => SetPooling.Run(freeze, QuantBasis.MsmsOnly, new[] { once, once }),
            Throws.ArgumentException.With.Message.Contains("is supplied twice"));
        Assert.That(() => SetPooling.Run(freeze, (QuantBasis)9, new[] { once }), Throws.TypeOf<ArgumentOutOfRangeException>());
    }
}
