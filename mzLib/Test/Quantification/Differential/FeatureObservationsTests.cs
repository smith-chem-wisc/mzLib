using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The runner's per-protein values from either table (STAT1 M5b): usable levels (unique to the group, valued in at least 2
/// samples, ruling 6), the random-effect labels (ruling 13), and the row set (ruling 12).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class FeatureObservationsTests
{
    private static PeptideRunValues Peptide(string sequence, string groups, params (string File, double Intensity, ObservationState State)[] runs) =>
        new(sequence, sequence, groups.Split(';'), runs.ToDictionary(r => r.File, r => new RunValue(r.Intensity, r.State)));

    private static ObservationTable LabelFree(params PeptideRunValues[] peptides)
    {
        var runs = new[] { "s1", "s2", "s3", "s4" }.Select(s => new ObservationRun(s, s)).ToList();
        return ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, Array.Empty<ProteinGroupInfo>(), QuantBasis.MbrKept);
    }

    private static Dictionary<string, AnalysisSample> Samples(string? individualOfS1 = null, params string[] ids) =>
        ids.ToDictionary(id => id, id => new AnalysisSample(id) { Individual = id == "s1" ? individualOfS1 : null }, StringComparer.Ordinal);

    [Test]
    public void LabelFreeKeepsUniquePeptidesValuedInTwoSamplesAndFlagsTransfers()
    {
        const ObservationState q = ObservationState.Quantified, mbr = ObservationState.MbrTransferred;
        var table = LabelFree(
            Peptide("AAA", "G1", ("s1", 100, q), ("s2", 200, mbr), ("s3", 300, q)),
            Peptide("CCC", "G1", ("s1", 100, q)), // unique but in one sample only: not usable
            Peptide("DDD", "G1;G2", ("s1", 100, q), ("s2", 100, q)), // shared: never a level of either group
            Peptide("EEE", "G2", ("s2", 50, q), ("s4", 60, q)));

        var features = FeatureObservationBuilder.FromTable(table, Samples(null, "s1", "s2", "s3", "s4"), minimumSamples: 2);

        var g1 = features["G1"];
        Assert.That(g1.UniqueLevels, Is.EqualTo(2));
        Assert.That(g1.NuisanceLevels, Is.EqualTo(new[] { "AAA" }));
        Assert.That(g1.Values.Select(v => (v.SampleId, v.IsMbr)), Is.EqualTo(new[] { ("s1", false), ("s2", true), ("s3", false) }));
        Assert.That(g1.Values.Select(v => v.RandomSample), Is.EqualTo(new[] { "s1", "s2", "s3" }));
        Assert.That(g1.Values.Select(v => v.RandomIndividual), Is.All.Null);
        Assert.That(g1.SamplesWithValue(new[] { "s1", "s4" }), Is.EqualTo(1));
        Assert.That(features["G2"].NuisanceLevels, Is.EqualTo(new[] { "EEE" }));
    }

    [Test]
    public void LabelFreeRestrictsToTheStratumsSamplesAndCarriesTheIndividual()
    {
        const ObservationState q = ObservationState.Quantified;
        var table = LabelFree(Peptide("AAA", "G1", ("s1", 100, q), ("s2", 200, q), ("s3", 300, q)));

        var features = FeatureObservationBuilder.FromTable(table, Samples("mouse7", "s1", "s2"), minimumSamples: 2);

        Assert.That(features["G1"].SampleIds, Is.EqualTo(new[] { "s1", "s2" }), "s3 is in another stratum");
        Assert.That(features["G1"].Values[0].RandomIndividual, Is.EqualTo("mouse7"));
        Assert.That(FeatureObservationBuilder.FromTable(table, Samples(null, "s1"), 2)["G1"].Values, Is.Empty,
            "valued in one stratum sample only: no usable level");
    }

    private static IsobaricQuantSampleInfo Channel(string file, int plex, string channel, string condition, string? name,
        int techrep = 0, bool reference = false) =>
        new($@"C:\data\{file}.raw", condition, 0, techrep, 0, plex, channel, 126.127726, reference) { SampleName = name };

    [Test]
    public void IsobaricLevelsAreIonsInTheirPlexRunWithTheChannelAsIndividualOrTheAnimal()
    {
        var design = new SampleExperimentalDesign();
        foreach (int run in new[] { 1, 2 })
        {
            var file = $"r{run}";
            design.Add($@"C:\data\{file}.raw", new ISampleInfo[]
            {
                Channel(file, 1, "126", "A", "a", run), Channel(file, 1, "127N", "B", "b", run),
                Channel(file, 1, "127C", "pool", "pool", run, reference: true),
            });
        }
        var psms = new[] { "r1", "r2" }.Select(f =>
            new IsobaricPsm(f, 1, "PEPTIDE", "PEPTIDE", 2, new[] { "G1" }, new double[] { 100, 200, 300 }));
        var table = IsobaricObservationTable.Build(design, psms, Array.Empty<ProteinGroupInfo>());
        var samples = new Dictionary<string, AnalysisSample> { ["a"] = new("a"), ["b"] = new("b") };

        var g1 = FeatureObservationBuilder.FromIsobaric(table, samples, minimumSamples: 2)["G1"];

        Assert.That(g1.NuisanceLevels, Is.EqualTo(new[] { "PEPTIDE/2@plex1.run1", "PEPTIDE/2@plex1.run2" }));
        Assert.That(g1.Values.Select(v => (v.SampleId, v.RandomSample, v.RandomIndividual)), Is.EqualTo(new[]
        {
            ("a", "a@plex1.run1", "a"), ("b", "b@plex1.run1", "b"), ("a", "a@plex1.run2", "a"), ("b", "b@plex1.run2", "b"),
        }));

        var animals = new Dictionary<string, AnalysisSample> { ["a"] = new("a") { Individual = "m1" }, ["b"] = new("b") { Individual = "m2" } };
        var withAnimals = FeatureObservationBuilder.FromIsobaric(table, animals, 2)["G1"];
        Assert.That(withAnimals.Values.Select(v => (v.RandomSample, v.RandomIndividual)).Distinct(),
            Is.EqualTo(new[] { ("a", "m1"), ("b", "m2") }), "with real individuals: individual, then the channel's sample");
    }

    [TestCase("G1", null)]
    [TestCase("UNKNOWN", "undescribed")]
    [TestCase("DECOY", "decoy")]
    [TestCase("ENTRAP", "entrapment")]
    [TestCase("CONTAM", "contaminant")]
    [TestCase("WEAK", "above_1pct_q")]
    public void TheRowSetKeepsTargetsAtOnePercentAndNamesEveryOtherReason(string group, string? reason)
    {
        var groups = new[]
        {
            new ProteinGroupInfo("G1", 0.01, false, false), new ProteinGroupInfo("DECOY", 0.0, true, false),
            new ProteinGroupInfo("ENTRAP", 0.0, false, false, IsEntrapment: true), new ProteinGroupInfo("CONTAM", 0.0, false, true),
            new ProteinGroupInfo("WEAK", 0.02, false, false),
        }.ToDictionary(g => g.Name, StringComparer.Ordinal);

        Assert.That(FeatureObservationBuilder.RowSetExclusion(group, groups), Is.EqualTo(reason));
    }
}
