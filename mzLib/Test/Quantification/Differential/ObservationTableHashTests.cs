using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The observation tables' canonical sha256 (STAT1 M5b ruling 10), which the runner records as its first input so the
/// analysis id follows the data: the same table gives the same hash whatever the input order or culture, and anything the
/// analysis reads (a value, a state, a peptide's groups, a group's q-value or flags, a warning) changes it.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ObservationTableHashTests
{
    private static readonly string[] Ids = { "s1", "s2", "s3", "s4" };

    private static ObservationTable Table(QuantBasis basis = QuantBasis.MsmsOnly, bool reverse = false,
        double s1Value = 1000.5, ObservationState s1State = ObservationState.Quantified, double weakQ = 0.05,
        bool p2Decoy = false, string p1Groups = "P1")
    {
        var peptides = new List<PeptideRunValues>
        {
            new("PEPA", "PEPA", p1Groups.Split(';'), new Dictionary<string, RunValue>
            {
                ["s1"] = new(s1Value, s1State), ["s2"] = new(2000.25, ObservationState.Quantified),
                ["s3"] = new(1500.125, ObservationState.Quantified),
            }),
            new("PEPB", "PEPB", new[] { "P2" }, Ids.ToDictionary(s => s, s => new RunValue(3000.75, ObservationState.Quantified))),
            new("PEPC", "PEPC", new[] { "WEAK" }, Ids.ToDictionary(s => s, s => new RunValue(123.4, ObservationState.Quantified))),
        };
        if (reverse) peptides.Reverse();
        var groups = new[]
        {
            new ProteinGroupInfo("P1", 0, false, false, Genes: "GENE1"), new ProteinGroupInfo("P2", 0, p2Decoy, false),
            new ProteinGroupInfo("WEAK", weakQ, false, false),
        };
        var runs = Ids.Select(s => new ObservationRun(s, s)).ToList();
        return ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, groups, basis);
    }

    [Test]
    public void TheSameTableGivesTheSameHashWhateverTheInputOrderAndCulture()
    {
        string hash = Table().Sha256();
        Assert.That(hash, Does.Match("^[0-9a-f]{64}$"));
        Assert.That(Table(reverse: true).Sha256(), Is.EqualTo(hash), "peptides given in another order");

        var culture = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = new CultureInfo("de-DE");
            Assert.That(Table().Sha256(), Is.EqualTo(hash), "a culture with a decimal comma");
        }
        finally
        {
            CultureInfo.CurrentCulture = culture;
        }

        string text = Table().ToCanonicalText();
        Assert.That(text, Does.StartWith("observation_table\tlabel_free\tmsms_only\nsamples\ts1\ts2\ts3\ts4\n"));
        Assert.That(text, Does.Contain("group\tP1\t0\t0\t0\t0\tGENE1\n"));
        Assert.That(text, Does.Contain("group\tP2\t0\t0\t0\t0\t\\N\n"), "no genes cell is not an empty one");
        Assert.That(text, Does.Contain("\tNaN:NotDetected\n"), "PEPA has no value in s4");
        Assert.That(text, Does.Not.Contain(","));
    }

    [Test]
    public void AnythingTheAnalysisReadsChangesTheHash()
    {
        string hash = Table().Sha256();
        var changed = new Dictionary<string, string>
        {
            ["a value"] = Table(s1Value: 1000.75).Sha256(),
            ["a state with the same value"] = Table(QuantBasis.MbrKept, s1State: ObservationState.MbrTransferred).Sha256(),
            ["the basis"] = Table(QuantBasis.MbrKept).Sha256(),
            ["a group's q-value"] = Table(weakQ: 0.02).Sha256(),
            ["a group's decoy flag"] = Table(p2Decoy: true).Sha256(),
            ["a peptide's groups"] = Table(p1Groups: "P1;P2").Sha256(),
        };
        foreach (var (what, other) in changed)
            Assert.That(other, Is.Not.EqualTo(hash), what);
        Assert.That(changed.Values.Distinct().Count(), Is.EqualTo(changed.Count), "every change gives its own hash");
        Assert.That(Table(QuantBasis.MbrKept, s1State: ObservationState.MbrTransferred).Sha256(),
            Is.Not.EqualTo(Table(QuantBasis.MbrKept).Sha256()), "under one basis, only the state differs");
    }

    [Test]
    public void ATextFieldCannotForgeAnotherFieldOrLine()
    {
        Assert.That(CanonicalText.Escape("a\tb\nc\\d\re"), Is.EqualTo(@"a\tb\nc\\d\re"));
        Assert.That(CanonicalText.Escape(null), Is.EqualTo(@"\N"));
        Assert.That(CanonicalText.Escape(@"\N"), Is.EqualTo(@"\\N"), "a text that reads like null is not null");
        Assert.That(CanonicalText.Number(double.NaN), Is.EqualTo("NaN"));
        Assert.That(CanonicalText.Number(0.1 + 0.2), Is.EqualTo("0.30000000000000004"), "round-trip precision");
    }

    /// <summary>
    /// The runner records the table it analyses as its first input, with no file name and the table's hash, before the
    /// caller's own inputs; a different table gives a different analysis id, and a caller's own <c>observation_table</c>
    /// input is refused.
    /// </summary>
    [Test]
    public void TheRunnerRecordsTheTableFirstAndItsAnalysisIdFollowsTheData()
    {
        static DifferentialRequest Request(ObservationTable table) => new()
        {
            Table = table,
            Samples = Ids.Select(s => new AnalysisSample(s)
            {
                Factors = new Dictionary<string, string> { ["group"] = s is "s1" or "s2" ? "A" : "B" },
            }).ToList(),
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("group", "A") } },
        };
        var caller = new DifferentialInput("quant_peptides", "AllQuantifiedPeptides.tsv", new string('a', 64));
        var result = DifferentialAnalysis.Run(Request(Table()) with { Inputs = new[] { caller } }, maxThreads: 1);

        Assert.That(result.Metadata.Inputs, Is.EqualTo(new[]
        {
            new DifferentialInput(DifferentialConstants.ObservationTableRole, "", Table().Sha256()), caller,
        }));
        Assert.That(result.Rows, Is.Not.Empty);
        Assert.That(result.Rows.Select(r => r.AnalysisId).Distinct(), Is.EqualTo(new[] { result.Metadata.AnalysisId }));

        string same = DifferentialAnalysis.Run(Request(Table()), maxThreads: 1).Metadata.AnalysisId;
        Assert.That(DifferentialAnalysis.Run(Request(Table(reverse: true)), maxThreads: 1).Metadata.AnalysisId, Is.EqualTo(same));
        Assert.That(DifferentialAnalysis.Run(Request(Table(s1Value: 1000.75)), maxThreads: 1).Metadata.AnalysisId,
            Is.Not.EqualTo(same), "one value different, the same settings and no caller inputs");

        var own = new DifferentialInput(DifferentialConstants.ObservationTableRole, "peptides.tsv", new string('b', 64));
        Assert.That(() => DifferentialAnalysis.Run(Request(Table()) with { Inputs = new[] { own } }),
            Throws.ArgumentException.With.Message.Contains("records the observation table itself"));
    }
}
