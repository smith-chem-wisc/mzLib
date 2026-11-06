using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The isobaric observation table (REQ-QP-15 PLAN 6f rule 1): one value per peptide ion per plex run from its match with
/// the largest reporter sum, ions with 1 or 2 channels measured dropped, 0 read as missing, and the design's samples.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class IsobaricObservationTableTests
{
    /// <summary>A four-channel plex: two sample channels, a pooled reference, and a channel the plex does not use.</summary>
    private static ISampleInfo[] Plex(string file, int plex, int techrep = 0, int fraction = 0) => new ISampleInfo[]
    {
        Channel(file, plex, "126", "A", $"a{plex}", techrep, fraction),
        Channel(file, plex, "127N", "B", $"b{plex}", techrep, fraction),
        Channel(file, plex, "127C", "pool", "pool", techrep, fraction, reference: true),
        Channel(file, plex, "128N", "", null, techrep, fraction),
    };

    private static IsobaricQuantSampleInfo Channel(string file, int plex, string channel, string condition, string? name,
        int techrep = 0, int fraction = 0, bool reference = false) =>
        new($@"C:\data\{file}.raw", condition, 0, techrep, fraction, plex, channel, 126.127726, reference) { SampleName = name };

    private static SampleExperimentalDesign Design(params ISampleInfo[][] files)
    {
        var design = new SampleExperimentalDesign();
        foreach (var samples in files)
            design.Add(samples[0].FullFilePathWithExtension, samples);
        return design;
    }

    private static IsobaricPsm Psm(string file, int scan, string sequence, int charge, params double[] intensities) =>
        new(file, scan, sequence, sequence, charge, new[] { "P1" }, intensities);

    private static IsobaricObservationTable Build(SampleExperimentalDesign design, params IsobaricPsm[] psms) =>
        IsobaricObservationTable.Build(design, psms, Array.Empty<ProteinGroupInfo>());

    private static double Cell(IsobaricObservationTable table, int row, string sample) =>
        table.Log2Intensity(row, table.Samples.ToList().IndexOf(sample));

    [Test]
    public void TheMatchWithTheLargestReporterSumGivesTheIonsValuesOverAllFractions()
    {
        var design = Design(Plex("f1", 1, fraction: 1), Plex("f2", 1, fraction: 2));
        var table = Build(design,
            Psm("f1", 10, "PEPTIDE", 2, 100, 200, 300, 0),
            Psm("f2", 20, "PEPTIDE", 2, 400, 100, 300, 0), // sum 800 > 600: this one, though channel b is lower
            Psm("f1", 11, "PEPTIDE", 2, 50, 60, 70, 0));

        Assert.That(table.Rows, Has.Count.EqualTo(1));
        Assert.That(table.Rows[0].FileName, Is.EqualTo("f2"));
        Assert.That(table.Rows[0].ScanNumber, Is.EqualTo(20));
        Assert.That(Cell(table, 0, "a1"), Is.EqualTo(Math.Log2(400)));
        Assert.That(Cell(table, 0, "b1"), Is.EqualTo(Math.Log2(100)));
    }

    [Test]
    public void ChargesAreSeparateIonsAndTechnicalReplicatesAreSeparatePlexRuns()
    {
        var design = Design(Plex("t1", 1, techrep: 1), Plex("t2", 1, techrep: 2), Plex("p2", 2));
        var table = Build(design,
            Psm("t2", 1, "PEPTIDE", 2, 1, 2, 3, 4),
            Psm("t1", 1, "PEPTIDE", 3, 1, 2, 3, 4),
            Psm("t1", 2, "PEPTIDE", 2, 1, 2, 3, 4),
            Psm("p2", 1, "AAA", 2, 1, 2, 3, 4));

        Assert.That(table.Rows.Select(r => (r.Peptide.FullSequence, r.Charge, r.Plex, r.TechnicalReplicate)), Is.EqualTo(new[]
        {
            ("AAA", 2, 2, 0), ("PEPTIDE", 2, 1, 1), ("PEPTIDE", 2, 1, 2), ("PEPTIDE", 3, 1, 1),
        }));
    }

    [Test]
    public void AReporterWrittenAsZeroIsMissingNotZero()
    {
        var table = Build(Design(Plex("f", 1)), Psm("f", 1, "PEPTIDE", 2, 0, 8, 4, 2));

        int a = table.Samples.ToList().IndexOf("a1"), b = table.Samples.ToList().IndexOf("b1");
        Assert.That(table.Log2Intensity(0, a), Is.NaN);
        Assert.That(table.State(0, a), Is.EqualTo(ObservationState.IdentifiedNotQuantified));
        Assert.That(table.Log2Intensity(0, b), Is.EqualTo(3));
        Assert.That(table.State(0, b), Is.EqualTo(ObservationState.Quantified));
    }

    [Test]
    public void ReferenceAndUnusedChannelsHoldNoSampleButCountTowardsTheFewChannelRule()
    {
        var table = Build(Design(Plex("f", 1)),
            Psm("f", 1, "THREE", 2, 0, 8, 4, 2), // three measured, one in a sample channel: kept
            Psm("f", 2, "TWO", 2, 16, 8, 0, 0), // two measured: dropped
            Psm("f", 3, "ONE", 2, 0, 0, 0, 5), // one measured: dropped
            Psm("f", 4, "NONE", 2, 0, 0, 0, 0)); // none measured: dropped, not an all-missing row

        Assert.That(table.Samples, Is.EqualTo(new[] { "a1", "b1" }));
        Assert.That(table.Rows.Select(r => r.Peptide.FullSequence), Is.EqualTo(new[] { "THREE" }));
        Assert.That(table.IonsWithFewChannels, Is.EqualTo(3));
        Assert.That(table.Observations().Select(o => o.SampleId), Is.EqualTo(new[] { "a1", "b1" }));
    }

    [Test]
    public void TheFewChannelRuleJudgesTheChosenMatchNotTheBestChannelsOfAll()
    {
        // The larger-sum match has two channels measured, so the ion is dropped although another match had four.
        var table = Build(Design(Plex("f", 1)),
            Psm("f", 1, "PEPTIDE", 2, 1, 1, 1, 1),
            Psm("f", 2, "PEPTIDE", 2, 100, 100, 0, 0));

        Assert.That(table.Rows, Is.Empty);
        Assert.That(table.IonsWithFewChannels, Is.EqualTo(1));
    }

    [Test]
    public void ARowHasValuesOnlyForTheSamplesOfItsPlex()
    {
        var table = Build(Design(Plex("p1", 1), Plex("p2", 2)), Psm("p2", 1, "PEPTIDE", 2, 2, 4, 8, 16));

        Assert.That(table.Samples, Is.EqualTo(new[] { "a1", "a2", "b1", "b2" }));
        Assert.That(table.SamplesOfPlex(2).Select(s => table.Samples[s]), Is.EqualTo(new[] { "a2", "b2" }));
        Assert.That(table.SamplesOfPlex(3), Is.Empty);
        Assert.That(table.PlexOf(0), Is.EqualTo(1));
        Assert.That(Cell(table, 0, "a1"), Is.NaN);
        Assert.That(table.State(0, 0), Is.EqualTo(ObservationState.NotDetected));
        Assert.That(Cell(table, 0, "a2"), Is.EqualTo(1));
        Assert.That(Cell(table, 0, "b2"), Is.EqualTo(2));
        Assert.That(table.Observations().Select(o => (o.SampleId, o.Log2Intensity)), Is.EqualTo(new[] { ("a2", 1.0), ("b2", 2.0) }));
    }

    [Test]
    public void ATieGoesToTheEarlierFileThenTheLowerScanWhateverTheInputOrder()
    {
        var design = Design(Plex("f1", 1, fraction: 1), Plex("f2", 1, fraction: 2));
        var psms = new[]
        {
            Psm("f2", 1, "PEPTIDE", 2, 1, 2, 3, 4),
            Psm("f1", 9, "PEPTIDE", 2, 4, 3, 2, 1),
            Psm("f1", 5, "PEPTIDE", 2, 2, 2, 3, 3),
        };

        foreach (var order in new[] { psms, psms.Reverse().ToArray(), new[] { psms[2], psms[0], psms[1] } })
        {
            var table = Build(design, order);
            Assert.That((table.Rows[0].FileName, table.Rows[0].ScanNumber), Is.EqualTo(("f1", 5)));
            Assert.That(Cell(table, 0, "a1"), Is.EqualTo(1));
        }
    }

    [Test]
    public void ACalibratedFileNameMatchesItsDesignFileButNotAlongsideTheUncalibratedName()
    {
        var design = Design(Plex("f", 1));
        Assert.That(Build(design, Psm("f-calib", 1, "PEPTIDE", 2, 1, 2, 3, 4)).Rows[0].FileName, Is.EqualTo("f-calib"));

        var ex = Assert.Throws<ArgumentException>(() =>
            Build(design, Psm("f-calib", 1, "PEPTIDE", 2, 1, 2, 3, 4), Psm("f", 2, "OTHER", 2, 1, 2, 3, 4)));
        Assert.That(ex!.Message, Does.Contain("both 'f-calib' and 'f'"));
    }

    [Test]
    public void RowsAreOrderedBySequenceOrdinallyThenCharge()
    {
        var table = Build(Design(Plex("f", 1)),
            Psm("f", 1, "b", 2, 1, 2, 3, 4), Psm("f", 2, "B", 3, 1, 2, 3, 4), Psm("f", 3, "B", 2, 1, 2, 3, 4));

        Assert.That(table.Rows.Select(r => (r.Peptide.FullSequence, r.Charge)), Is.EqualTo(new[] { ("B", 2), ("B", 3), ("b", 2) }));
    }

    [Test]
    public void OnePeptideIsOneObjectWithItsGroupsSortedAndUndescribedGroupsAreWarned()
    {
        var psms = new[]
        {
            new IsobaricPsm("f", 1, "PEPTIDE", "PEPTIDE", 2, new[] { "P2", "P1", "P2" }, new double[] { 1, 2, 3, 4 }),
            new IsobaricPsm("f", 2, "PEPTIDE", "PEPTIDE", 3, new[] { "P1", "P2" }, new double[] { 1, 2, 3, 4 }),
        };
        var table = IsobaricObservationTable.Build(Design(Plex("f", 1)), psms,
            new[] { new ProteinGroupInfo("P1", 0.001, false, false) });

        Assert.That(table.Rows[0].Peptide, Is.SameAs(table.Rows[1].Peptide));
        Assert.That(table.Rows[0].Peptide.ProteinGroups, Is.EqualTo(new[] { "P1", "P2" }));
        Assert.That(table.Warnings.Single(), Does.Contain("1 protein group(s)").And.Contain("P2"));
    }

    private static readonly object[] Refused =
    {
        new object[] { Psm("g", 1, "PEPTIDE", 2, 1, 2, 3, 4), "does not have" },
        new object[] { Psm("f", 1, "PEPTIDE", 2, 1, 2, 3), "3 reporter intensities" },
        new object[] { Psm("f", 1, "PEPTIDE", 2, 1, -2, 3, 4), "-2" },
        new object[] { Psm("f", 1, "PEPTIDE", 2, 1, double.NaN, 3, 4), "NaN" },
        new object[] { Psm("f", 1, "PEPTIDE", 2, 1, double.PositiveInfinity, 3, 4), "finite" },
        new object[] { Psm("f", 1, "PEPTIDE", 0, 1, 2, 3, 4), "charge 0" },
        new object[] { Psm("f", 1, "", 2, 1, 2, 3, 4), "no full sequence" },
        new object[] { new IsobaricPsm("f", 1, "PEPTIDE", "PEPTIDE", 2, new[] { "P1;P2" }, new double[] { 1, 2, 3, 4 }), "';'" },
    };

    [TestCaseSource(nameof(Refused))]
    public void AMatchThatCannotBeReadIsRefusedByName(IsobaricPsm psm, string message)
    {
        var ex = Assert.Throws<ArgumentException>(() => Build(Design(Plex("f", 1)), psm));
        Assert.That(ex!.Message, Does.Contain(message));
    }

    [Test]
    public void APeptideWhoseMatchesDisagreeOnItsGroupsOrBaseSequenceIsRefused()
    {
        var design = Design(Plex("f", 1));
        var first = new IsobaricPsm("f", 1, "PEP[x]TIDE", "PEPTIDE", 2, new[] { "P1" }, new double[] { 1, 2, 3, 4 });

        Assert.Throws<ArgumentException>(() => Build(design, first, first with { ScanNumber = 2, ProteinGroups = new[] { "P2" } }));
        Assert.Throws<ArgumentException>(() => Build(design, first, first with { ScanNumber = 2, BaseSequence = "PEPTIDF" }));
    }

    [Test]
    public void ADesignTheTableCannotReadIsRefused()
    {
        var mixed = new[] { Channel("f", 1, "126", "A", "a"), Channel("f", 2, "127N", "B", "b") };
        Assert.That(Assert.Throws<ArgumentException>(() => Build(Design(mixed)))!.Message,
            Does.Contain("not all in one plex"));

        Assert.That(Assert.Throws<ArgumentException>(() => Build(Design(Plex("f", 1), Plex("g", 1))))!.Message,
            Does.Contain("are both plex 1"));

        var labelFree = Design(new ISampleInfo[] { new SpectraFileInfo(@"C:\data\lf.raw", "A", 0, 0, 0) });
        Assert.That(Assert.Throws<ArgumentException>(() => Build(labelFree))!.Message, Does.Contain("not an isobaric channel"));
    }
}
