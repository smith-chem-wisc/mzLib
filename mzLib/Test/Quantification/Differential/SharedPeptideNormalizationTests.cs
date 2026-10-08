using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// GR-2 / GR-12 / GR-25 (<c>DEF-DIFF-NORM</c>): each sample's shift is its median log2 difference from the peptide's
/// mean over the samples, measured on the peptides with a value in every sample (fallback: in at least half), never on
/// decoy, contaminant or entrapment peptides.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SharedPeptideNormalizationTests
{
    private static readonly string[] SampleIds = { "S1", "S2", "S3", "S4" };
    private static readonly double[] Offsets = { 0.5, -0.25, 1.0, 0.0 };

    /// <summary>log2 intensity = 20 + p/100 + the sample's offset, for each sample index given.</summary>
    private static PeptideRunValues Shifted(string sequence, string group, int p, IEnumerable<int> samples,
        Func<int, double>? extra = null, ObservationState state = ObservationState.Quantified) =>
        new(sequence, sequence, new[] { group }, samples.ToDictionary(s => SampleIds[s],
            s => new RunValue(Math.Pow(2, 20 + p / 100.0 + Offsets[s] + (extra?.Invoke(s) ?? 0)), state)));

    private static IEnumerable<PeptideRunValues> EverySample(int count, string prefix = "PEP", string group = "P1") =>
        Enumerable.Range(0, count).Select(p => Shifted($"{prefix}{p:D4}", group, p, Enumerable.Range(0, 4)));

    private static ObservationTable Table(IEnumerable<PeptideRunValues> peptides, QuantBasis basis = QuantBasis.MsmsOnly,
        IEnumerable<ProteinGroupInfo>? groups = null)
    {
        var runs = SampleIds.Select(s => new ObservationRun(s, s)).ToList();
        return ObservationTable.Build(runs, SampleIds, peptides, groups ?? Array.Empty<ProteinGroupInfo>(), basis);
    }

    private static void AssertShifts(NormalizationResult result, double[] expected)
    {
        for (int s = 0; s < SampleIds.Length; s++)
            Assert.That(result.PerSampleShift[SampleIds[s]], Is.EqualTo(expected[s]).Within(1e-9), SampleIds[s]);
    }

    /// <summary>The known truth: each sample's offset minus the mean offset (0.3125).</summary>
    private static readonly double[] TrueShifts = Offsets.Select(o => o - Offsets.Average()).ToArray();

    [Test]
    public void RecoversKnownShiftsOnPeptidesInEverySample()
    {
        var result = SharedPeptideNormalization.Measure(Table(EverySample(120)));

        Assert.That(result.Setting, Is.EqualTo(SharedPeptideNormalization.SharedPeptideMedian));
        Assert.That(result.ReferenceSetSize, Is.EqualTo(120));
        Assert.That(result.MinimumReferenceSetSize, Is.EqualTo(100));
        AssertShifts(result, TrueShifts);
        Assert.That(result.Warnings, Is.Empty);
    }

    [Test]
    public void PeptidesMissingFromASampleDoNotMoveTheShifts()
    {
        // 200 peptides missing from S1 and 3 log2 units high in S2: on every value they would drag S2 up.
        var partial = Enumerable.Range(0, 200).Select(p =>
            Shifted($"GAP{p:D4}", "P1", p, new[] { 1, 2, 3 }, s => s == 1 ? 3 : 0));

        var result = SharedPeptideNormalization.Measure(Table(EverySample(120).Concat(partial)));

        Assert.That(result.ReferenceSetSize, Is.EqualTo(120));
        AssertShifts(result, TrueShifts);
    }

    [TestCase("contaminant")]
    [TestCase("decoy")]
    [TestCase("entrapment")]
    public void DecoyContaminantAndEntrapmentPeptidesStayOutOfTheReference(string kind)
    {
        var group = new ProteinGroupInfo("X1", 0.001, IsDecoy: kind == "decoy", IsContaminant: kind == "contaminant",
            IsEntrapment: kind == "entrapment");
        var excluded = Enumerable.Range(0, 300).Select(p =>
            Shifted($"XXX{p:D4}", "X1", p, Enumerable.Range(0, 4), s => s == 2 ? 4 : 0));
        var table = Table(EverySample(120).Concat(excluded), groups: new[] { group, new ProteinGroupInfo("P1", 0, false, false) });

        var result = SharedPeptideNormalization.Measure(table);

        Assert.That(result.ReferenceSetSize, Is.EqualTo(120));
        AssertShifts(result, TrueShifts);
        Assert.That(SharedPeptideNormalization.IsReferenceEligible(table, table.Peptides.First(p => p.FullSequence.StartsWith("XXX"))),
            Is.False);
    }

    [Test]
    public void APeptideSharedWithAContaminantGroupStaysOut()
    {
        var groups = new[] { new ProteinGroupInfo("P1", 0, false, false), new ProteinGroupInfo("C1", 0, false, true) };
        var shared = new PeptideRunValues("SHARED", "SHARED", new[] { "P1", "C1" },
            new Dictionary<string, RunValue> { ["S1"] = new(5, ObservationState.Quantified) });
        var table = Table(new[] { shared }, groups: groups);

        Assert.That(SharedPeptideNormalization.IsReferenceEligible(table, table.Peptides.Single()), Is.False);
    }

    [Test]
    public void FallsBackToTheHalfSetBelowTheFloor()
    {
        // 99 peptides in every sample, 40 in two of the four: below GR-25's 100, so the half set (139) is used.
        var halves = Enumerable.Range(0, 40).Select(p => Shifted($"HALF{p:D4}", "P1", p, new[] { 0, 2 }));
        var thirds = Enumerable.Range(0, 30).Select(p => Shifted($"ONE{p:D4}", "P1", p, new[] { 3 }));

        var below = SharedPeptideNormalization.Measure(Table(EverySample(99).Concat(halves).Concat(thirds)));
        Assert.That(below.Setting, Is.EqualTo(SharedPeptideNormalization.SharedPeptideMedianHalf));
        Assert.That(below.ReferenceSetSize, Is.EqualTo(139), "a peptide in one sample of four is not in the half set");

        var at = SharedPeptideNormalization.Measure(Table(EverySample(100).Concat(halves)));
        Assert.That(at.Setting, Is.EqualTo(SharedPeptideNormalization.SharedPeptideMedian));
        Assert.That(at.ReferenceSetSize, Is.EqualTo(100));
    }

    [Test]
    public void TheHalfSetUsesEachPeptidesMeanOverTheSamplesThatHaveIt()
    {
        // Two peptides in S1 and S2 only; with no every-sample peptide the half set is used.
        // Differences S1 - S2 are 1 and 3, so S1's ratios are 0.5 and 1.5 (median 1) and S2's -0.5 and -1.5 (median -1).
        var peptides = new[]
        {
            new PeptideRunValues("A", "A", new[] { "P1" }, new Dictionary<string, RunValue>
                { ["S1"] = new(Math.Pow(2, 11), ObservationState.Quantified), ["S2"] = new(Math.Pow(2, 10), ObservationState.Quantified) }),
            new PeptideRunValues("B", "B", new[] { "P1" }, new Dictionary<string, RunValue>
                { ["S1"] = new(Math.Pow(2, 13), ObservationState.Quantified), ["S2"] = new(Math.Pow(2, 10), ObservationState.Quantified) }),
        };

        var result = SharedPeptideNormalization.Measure(Table(peptides));

        Assert.That(result.Setting, Is.EqualTo(SharedPeptideNormalization.SharedPeptideMedianHalf));
        Assert.That(result.PerSampleShift["S1"], Is.EqualTo(1).Within(1e-12), "the mean of the middle two");
        Assert.That(result.PerSampleShift["S2"], Is.EqualTo(-1).Within(1e-12));
        Assert.That(result.PerSampleShift["S3"], Is.NaN);
        Assert.That(result.PerSampleShift["S4"], Is.NaN);
        Assert.That(result.Warnings, Has.Count.EqualTo(1));
        Assert.That(result.Warnings[0], Does.Contain("S3").And.Contain("S4"));
    }

    [Test]
    public void ApplySubtractsEachSamplesShiftAndLeavesMissingAlone()
    {
        var peptides = EverySample(120).Append(Shifted("ZZZ", "P1", 0, new[] { 0 })).ToList();
        var table = Table(peptides);
        var result = SharedPeptideNormalization.Measure(table);

        var normalized = SharedPeptideNormalization.Apply(table, result);

        for (int p = 0; p < table.Peptides.Count; p++)
        for (int s = 0; s < table.Samples.Count; s++)
        {
            double expected = table.Log2Intensity(p, s) - result.PerSampleShift[table.Samples[s]];
            Assert.That(normalized.Log2Intensity(p, s), double.IsNaN(expected) ? Is.NaN : Is.EqualTo(expected));
            Assert.That(normalized.State(p, s), Is.EqualTo(table.State(p, s)));
        }

        // After the shift, the four samples agree on every reference peptide.
        int first = 0;
        var values = Enumerable.Range(0, 4).Select(s => normalized.Log2Intensity(first, s)).ToList();
        Assert.That(values.Max() - values.Min(), Is.LessThan(1e-9));
    }

    [Test]
    public void ASampleWithoutAShiftIsLeftAsItIs()
    {
        var table = Table(EverySample(120));
        var shifts = SampleIds.ToDictionary(s => s, s => s == "S2" ? double.NaN : 1.0);
        var result = new NormalizationResult(SharedPeptideNormalization.SharedPeptideMedian, 120, 100, shifts, Array.Empty<string>());

        var normalized = SharedPeptideNormalization.Apply(table, result);

        Assert.That(normalized.Log2Intensity(0, 1), Is.EqualTo(table.Log2Intensity(0, 1)));
        Assert.That(normalized.Log2Intensity(0, 0), Is.EqualTo(table.Log2Intensity(0, 0) - 1));
    }

    [Test]
    public void TheReferenceSetFollowsTheBasis()
    {
        // 50 peptides are transferred into S4: they are in every sample only when transfers are kept.
        var transferred = Enumerable.Range(0, 50).Select(p => new PeptideRunValues($"MBR{p:D4}", $"MBR{p:D4}", new[] { "P1" },
            Enumerable.Range(0, 4).ToDictionary(s => SampleIds[s], s => new RunValue(Math.Pow(2, 20 + Offsets[s]),
                s == 3 ? ObservationState.MbrTransferred : ObservationState.Quantified))));
        var peptides = EverySample(100).Concat(transferred).ToList();

        Assert.That(SharedPeptideNormalization.Measure(Table(peptides, QuantBasis.MsmsOnly)).ReferenceSetSize, Is.EqualTo(100));
        Assert.That(SharedPeptideNormalization.Measure(Table(peptides, QuantBasis.MbrKept)).ReferenceSetSize, Is.EqualTo(150));
    }

    [Test]
    public void AStratumIsMeasuredOnItsOwnSamples()
    {
        var table = Table(EverySample(120));
        var stratum = table.ForSamples(new[] { "S1", "S2" });

        var result = SharedPeptideNormalization.Measure(stratum);

        Assert.That(result.PerSampleShift.Keys, Is.EquivalentTo(new[] { "S1", "S2" }));
        Assert.That(result.PerSampleShift["S1"], Is.EqualTo((Offsets[0] - Offsets[1]) / 2).Within(1e-9));
    }
}
