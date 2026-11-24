using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// Precursor rows in the label-free observation table (QuantProject dia 073, DEF-DIFF-ROWSET for <c>dia</c>): a row may
/// carry a precursor charge, so one peptide at two charges is two rows, each its own level in the model; a table carries
/// a charge on every row or on none; a table without charges reads, orders and hashes exactly as before.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ObservationTablePrecursorTests
{
    private static readonly string[] SampleIds = { "s1", "s2", "s3", "s4" };

    private static PeptideRunValues Precursor(string sequence, int? charge, string groups = "G1", double level = 100,
        string? baseSequence = null) =>
        new(sequence, baseSequence ?? sequence, groups.Split(';'),
            SampleIds.ToDictionary(s => s, s => new RunValue(level * (1 + Array.IndexOf(SampleIds, s)), ObservationState.Quantified)),
            charge);

    private static ObservationTable Table(params PeptideRunValues[] rows)
    {
        var runs = SampleIds.Select(s => new ObservationRun(s, s)).ToList();
        return ObservationTable.Build(runs, SampleIds, rows, new[] { new ProteinGroupInfo("G1", 0, false, false) }, QuantBasis.MsmsOnly);
    }

    private static Dictionary<string, AnalysisSample> Samples() =>
        SampleIds.ToDictionary(id => id, id => new AnalysisSample(id), StringComparer.Ordinal);

    [Test]
    public void OnePeptideAtTwoChargesIsTwoRowsOrderedBySequenceThenCharge()
    {
        var table = Table(Precursor("BBB", 2, level: 30), Precursor("AAA", 3, level: 20), Precursor("AAA", 2, level: 10));

        Assert.That(table.Peptides.Select(p => p.FullSequence), Is.EqualTo(new[] { "AAA", "AAA", "BBB" }));
        Assert.That(Enumerable.Range(0, 3).Select(table.Charge), Is.EqualTo(new int?[] { 2, 3, 2 }));
        Assert.That(table.Log2Intensity(0, 0), Is.EqualTo(Math.Log2(10)).Within(1e-12), "AAA/2 keeps its own values");
        Assert.That(table.Log2Intensity(1, 0), Is.EqualTo(Math.Log2(20)).Within(1e-12), "AAA/3 keeps its own values");
    }

    [Test]
    public void APrecursorListedTwiceIsRefused()
    {
        var ex = Assert.Throws<ArgumentException>(() => Table(Precursor("AAA", 2), Precursor("AAA", 2, level: 50)));

        Assert.That(ex!.Message, Does.Contain("'AAA'").And.Contain("charge 2").And.Contain("twice"));
    }

    [Test]
    public void ARowWithAChargeAndARowWithoutOneAreRefusedTogether()
    {
        var ex = Assert.Throws<ArgumentException>(() => Table(Precursor("AAA", 2), Precursor("BBB", null)));

        Assert.That(ex!.Message, Does.Contain("'AAA'").And.Contain("'BBB'"));
    }

    [Test]
    public void AChargeBelowOneIsRefused()
    {
        var ex = Assert.Throws<ArgumentException>(() => Table(Precursor("AAA", 0)));

        Assert.That(ex!.Message, Does.Contain("'AAA'").And.Contain("at least 1"));
    }

    [Test]
    public void ThePrecursorsOfOnePeptideMustAgreeOnItsBaseSequenceAndGroups()
    {
        Assert.That(Assert.Throws<ArgumentException>(() => Table(Precursor("AAA", 2, groups: "G1"), Precursor("AAA", 3, groups: "G2")))!
            .Message, Does.Contain("'AAA'").And.Contain("'G1'").And.Contain("'G2'"));
        Assert.That(Assert.Throws<ArgumentException>(() => Table(Precursor("AAA", 2, baseSequence: "AAA"), Precursor("AAA", 3, baseSequence: "AAB")))!
            .Message, Does.Contain("'AAB'"));
    }

    [Test]
    public void APrecursorRowIsWrittenWithItsChargeAndTheChargeChangesTheHash()
    {
        string text = Table(Precursor("AAA", 2)).ToCanonicalText();

        Assert.That(text, Does.Contain("\nprecursor\tAAA\tAAA\tG1\t2\t"));
        Assert.That(text, Does.Not.Contain("\npeptide\t"));
        Assert.That(Table(Precursor("AAA", 3)).Sha256(), Is.Not.EqualTo(Table(Precursor("AAA", 2)).Sha256()));
        Assert.That(Table(Precursor("AAA", null)).ToCanonicalText(), Does.Contain("\npeptide\tAAA\tAAA\tG1\t"),
            "a table without charges is written as before");
        Assert.That(Table(Precursor("AAA", null)).Charge(0), Is.Null);
    }

    [Test]
    public void TheChargeSurvivesASubsetOfSamplesANewScaleAndTheLongForm()
    {
        var table = Table(Precursor("AAA", 3), Precursor("AAA", 2));

        var subset = table.ForSamples(new[] { "s1", "s3" });
        Assert.That(new[] { subset.Charge(0), subset.Charge(1) }, Is.EqualTo(new int?[] { 2, 3 }));
        var shifted = SharedPeptideNormalization.Apply(table, new NormalizationResult(SharedPeptideNormalization.SharedPeptideMedian, 2, 1,
            SampleIds.ToDictionary(s => s, _ => 0.5), Array.Empty<string>()));
        Assert.That(new[] { shifted.Charge(0), shifted.Charge(1) }, Is.EqualTo(new int?[] { 2, 3 }));
        Assert.That(table.Observations().Where(o => o.SampleId == "s1").Select(o => o.Charge), Is.EqualTo(new int?[] { 2, 3 }));
    }

    [Test]
    public void EachPrecursorIsALevelOfItsOwnAndCountsOnceAsAPeptide()
    {
        var table = Table(Precursor("AAA", 2), Precursor("AAA", 3), Precursor("CCC", 2));

        var g1 = FeatureObservationBuilder.FromTable(table, Samples(), minimumSamples: 2)["G1"];

        Assert.That(g1.NuisanceLevels, Is.EqualTo(new[] { "AAA/2", "AAA/3", "CCC/2" }));
        Assert.That(g1.UniqueLevels, Is.EqualTo(3));
        Assert.That(g1.PeptideCount, Is.EqualTo(2), "two peptides, three precursors");
        Assert.That(g1.Values.Select(v => v.Peptide).Distinct(), Is.EqualTo(new[] { "AAA", "CCC" }));
    }
}
