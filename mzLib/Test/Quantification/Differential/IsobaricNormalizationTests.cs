using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The global channel-median normalization (REQ-QP-15 PLAN 6f rule 3, <c>DEF-DIFF-NORM</c>): each channel's median in each
/// plex run moved to the median of those medians, never measured on decoy, contaminant or entrapment rows.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class IsobaricNormalizationTests
{
    /// <summary>A plex run of three sample channels (s{plex}a, s{plex}b, s{plex}c) and a pooled reference.</summary>
    private static ISampleInfo[] Plex(string file, int plex, int techrep = 0) => new[] { "126", "127N", "127C", "128N" }
        .Select((c, k) => (ISampleInfo)new IsobaricQuantSampleInfo($@"C:\data\{file}.raw", k == 0 ? "A" : k < 3 ? "B" : "pool", 0,
            techrep, 0, plex, c, 126.127726, k == 3) { SampleName = k < 3 ? $"s{plex}{(char)('a' + k)}" : "pool" })
        .ToArray();

    private static SampleExperimentalDesign Design(params ISampleInfo[][] files)
    {
        var design = new SampleExperimentalDesign();
        foreach (var samples in files)
            design.Add(samples[0].FullFilePathWithExtension, samples);
        return design;
    }

    /// <summary>Ion i in a file: log2 intensity 10 + i + the channel's offset in that file; the reference 1024.</summary>
    private static IEnumerable<IsobaricPsm> Ions(string file, int count, double[] offsets, string group = "P1", string prefix = "PEP") =>
        Enumerable.Range(0, count).Select(i => new IsobaricPsm(file, i, $"{prefix}{i:D3}", $"{prefix}{i:D3}", 2, new[] { group },
            offsets.Select(o => Math.Pow(2, 10 + i + o)).Append(1024).ToArray()));

    /// <summary>The ion number a row's sequence ends in.</summary>
    private static int Ion(IsobaricIonRow row) => int.Parse(row.Peptide.FullSequence[^3..]);

    private static double Shift(IsobaricNormalizationResult result, string sample, int techrep = 0) =>
        result.Shifts.Single(s => s.SampleId == sample && s.TechnicalReplicate == techrep).Shift;

    [Test]
    public void EachChannelIsMovedToTheMedianOfEveryPlexRunsChannelMedians()
    {
        var design = Design(Plex("p1", 1), Plex("p2", 2));
        double[] o1 = { 0.5, -0.25, 1.0 }, o2 = { 0.0, 2.0, -1.0 };
        var table = IsobaricObservationTable.Build(design, Ions("p1", 5, o1).Concat(Ions("p2", 5, o2)),
            Array.Empty<ProteinGroupInfo>());

        var result = IsobaricNormalization.Measure(table);

        // Every channel's median is 12 + its offset; the median of the six offsets {-1, -0.25, 0, 0.5, 1, 2} is 0.25.
        Assert.That(result.Setting, Is.EqualTo(IsobaricNormalization.GlobalChannelMedian));
        Assert.That(result.Warnings, Is.Empty);
        Assert.That(result.Shifts.Select(s => (s.Plex, s.SampleId, s.Values)),
            Is.EqualTo(new[] { (1, "s1a", 5), (1, "s1b", 5), (1, "s1c", 5), (2, "s2a", 5), (2, "s2b", 5), (2, "s2c", 5) }));
        var expected = o1.Concat(o2).Select(o => o - 0.25).ToArray();
        Assert.That(result.Shifts.Select(s => s.Shift), Is.EqualTo(expected).Within(1e-12));

        var normalized = IsobaricNormalization.Apply(table, result);
        for (int r = 0; r < normalized.Rows.Count; r++)
            foreach (int s in normalized.SamplesOfPlex(normalized.Rows[r].Plex))
                Assert.That(normalized.Log2Intensity(r, s), Is.EqualTo(10.25 + Ion(normalized.Rows[r])).Within(1e-12));
    }

    [Test]
    public void EachTechnicalReplicateOfAPlexIsShiftedOnItsOwn()
    {
        var design = Design(Plex("t1", 1, techrep: 1), Plex("t2", 1, techrep: 2));
        var table = IsobaricObservationTable.Build(design,
            Ions("t1", 3, new[] { 0.0, 1.0, 2.0 }).Concat(Ions("t2", 3, new[] { 3.0, 4.0, 5.0 })), Array.Empty<ProteinGroupInfo>());

        var result = IsobaricNormalization.Measure(table);

        // Offsets {0..5}: median 2.5.
        Assert.That(Shift(result, "s1a", 1), Is.EqualTo(-2.5).Within(1e-12));
        Assert.That(Shift(result, "s1a", 2), Is.EqualTo(0.5).Within(1e-12));
        var normalized = IsobaricNormalization.Apply(table, result);
        Assert.That(Enumerable.Range(0, normalized.Rows.Count).SelectMany(r =>
                normalized.SamplesOfPlex(1).Select(s => normalized.Log2Intensity(r, s) - Ion(normalized.Rows[r]))),
            Has.All.EqualTo(12.5).Within(1e-12));
    }

    [Test]
    public void DecoyContaminantAndEntrapmentRowsDoNotEnterTheMedians()
    {
        var design = Design(Plex("p1", 1));
        var psms = Ions("p1", 3, new[] { 0.0, 0.0, 0.0 })
            .Concat(Ions("p1", 3, new[] { 50.0, 0.0, 0.0 }, "DECOY_P9", "DEC"))
            .Concat(Ions("p1", 3, new[] { 50.0, 0.0, 0.0 }, "CON_P8", "CON"))
            .Concat(Ions("p1", 3, new[] { 50.0, 0.0, 0.0 }, "ENT_P7", "ENT"));
        var groups = new[]
        {
            new ProteinGroupInfo("DECOY_P9", 1, IsDecoy: true, IsContaminant: false),
            new ProteinGroupInfo("CON_P8", 0, IsDecoy: false, IsContaminant: true),
            new ProteinGroupInfo("ENT_P7", 0, IsDecoy: false, IsContaminant: false, IsEntrapment: true),
        };
        var table = IsobaricObservationTable.Build(design, psms, groups);

        var result = IsobaricNormalization.Measure(table);

        Assert.That(result.Shifts.Select(s => s.Shift), Has.All.EqualTo(0).Within(1e-12));
        Assert.That(result.Shifts.Select(s => s.Values), Has.All.EqualTo(3));
    }

    [Test]
    public void TheMedianOfAnEvenCountIsTheMeanOfTheMiddleTwo()
    {
        var design = Design(Plex("p1", 1));
        var psms = new[]
        {
            new IsobaricPsm("p1", 1, "A", "A", 2, new[] { "P1" }, new[] { Math.Pow(2, 10), 1024, 1024, 1024 }),
            new IsobaricPsm("p1", 2, "B", "B", 2, new[] { "P1" }, new[] { Math.Pow(2, 13), 1024, 1024, 1024 }),
        };
        var result = IsobaricNormalization.Measure(IsobaricObservationTable.Build(design, psms, Array.Empty<ProteinGroupInfo>()));

        // Channel medians 11.5, 10, 10; target 10.
        Assert.That(result.Shifts.Select(s => s.Shift), Is.EqualTo(new[] { 1.5, 0, 0 }).Within(1e-12));
    }

    [Test]
    public void AChannelWithNoValueInAPlexRunHasNoShiftIsWarnedAndIsLeftAsItIs()
    {
        var design = Design(Plex("p1", 1));
        var psms = Enumerable.Range(0, 3).Select(i =>
            new IsobaricPsm("p1", i, $"P{i}", $"P{i}", 2, new[] { "P1" }, new[] { Math.Pow(2, 10 + i), Math.Pow(2, 11 + i), 0, 1024 }))
            .Append(new IsobaricPsm("p1", 9, "Z", "Z", 2, new[] { "DECOY" }, new[] { 2.0, 4.0, 8.0, 16.0 }));
        var table = IsobaricObservationTable.Build(design, psms, new[] { new ProteinGroupInfo("DECOY", 1, true, false) });

        var result = IsobaricNormalization.Measure(table);

        Assert.That(Shift(result, "s1c"), Is.NaN);
        Assert.That(Shift(result, "s1a"), Is.EqualTo(-0.5).Within(1e-12));
        Assert.That(result.Warnings.Single(), Does.Contain("'s1c' in plex 1, technical replicate 0"));
        var normalized = IsobaricNormalization.Apply(table, result);
        int decoyRow = normalized.Rows.ToList().FindIndex(r => r.Peptide.FullSequence == "Z");
        Assert.That(normalized.Log2Intensity(decoyRow, 2), Is.EqualTo(3), "an unshifted channel keeps its value");
        Assert.That(normalized.Log2Intensity(0, 2), Is.NaN, "a missing value stays missing");
        Assert.That(normalized.State(0, 2), Is.EqualTo(ObservationState.IdentifiedNotQuantified));
    }

    [Test]
    public void NothingToMeasureIsSettingNoneAndAnEmptyTableIsNotShifted()
    {
        var design = Design(Plex("p1", 1));
        var decoyOnly = IsobaricObservationTable.Build(design,
            new[] { new IsobaricPsm("p1", 1, "Z", "Z", 2, new[] { "DECOY" }, new[] { 2.0, 4.0, 8.0, 16.0 }) },
            new[] { new ProteinGroupInfo("DECOY", 1, true, false) });
        var result = IsobaricNormalization.Measure(decoyOnly);
        Assert.That(result.Setting, Is.EqualTo(IsobaricNormalization.None));
        Assert.That(result.Warnings.Single(), Does.Contain("nothing is shifted"));

        var empty = IsobaricNormalization.Measure(IsobaricObservationTable.Build(design, Array.Empty<IsobaricPsm>(),
            Array.Empty<ProteinGroupInfo>()));
        Assert.That((empty.Setting, empty.Shifts.Count, empty.Warnings.Count), Is.EqualTo((IsobaricNormalization.None, 0, 0)));
    }

    [Test]
    public void ApplyRefusesAResultWithoutAShiftForAChannelOfAPlexRun()
    {
        var table = IsobaricObservationTable.Build(Design(Plex("p1", 1)), Ions("p1", 3, new[] { 0.0, 1.0, 2.0 }),
            Array.Empty<ProteinGroupInfo>());
        var result = IsobaricNormalization.Measure(table);
        var partial = result with { Shifts = result.Shifts.Where(s => s.SampleId != "s1b").ToList() };

        var ex = Assert.Throws<ArgumentException>(() => IsobaricNormalization.Apply(table, partial));
        Assert.That(ex!.Message, Does.Contain("'s1b' in plex 1"));
    }
}
