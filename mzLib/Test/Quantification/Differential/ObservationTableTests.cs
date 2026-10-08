using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The LFQ observation table of STAT1 M5: fractions summed before the log, technical replicates averaged on the log
/// scale, match-between-runs values kept apart by basis, samples without values kept, states never collapsed, and
/// MetaMorpheus's protein groups with the cell rule for unique peptides (GR-23).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ObservationTableTests
{
    private static RunValue Msms(double intensity) => new(intensity, ObservationState.Quantified);
    private static RunValue Mbr(double intensity) => new(intensity, ObservationState.MbrTransferred);
    private static RunValue NoValue(ObservationState state) => new(0, state);

    private static PeptideRunValues Peptide(string sequence, string groups, params (string Run, RunValue Value)[] values) =>
        new(sequence, sequence, groups.Split(';', StringSplitOptions.RemoveEmptyEntries),
            values.ToDictionary(v => v.Run, v => v.Value));

    private static ObservationTable Build(IReadOnlyList<ObservationRun> runs, QuantBasis basis,
        params PeptideRunValues[] peptides) =>
        ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, Array.Empty<ProteinGroupInfo>(), basis);

    private static double Cell(ObservationTable table, string sequence, string sample) =>
        table.Log2Intensity(Index(table, sequence), table.Samples.ToList().IndexOf(sample));

    private static ObservationState StateOf(ObservationTable table, string sequence, string sample) =>
        table.State(Index(table, sequence), table.Samples.ToList().IndexOf(sample));

    private static int Index(ObservationTable table, string sequence) =>
        table.Peptides.Select(p => p.FullSequence).ToList().IndexOf(sequence);

    [Test]
    public void FractionsOfOneReplicateAreSummedBeforeTheLog()
    {
        var runs = new[] { new ObservationRun("f1", "S1", Fraction: 1), new ObservationRun("f2", "S1", Fraction: 2) };
        var table = Build(runs, QuantBasis.MsmsOnly, Peptide("PEPA", "P1", ("f1", Msms(100)), ("f2", Msms(300))));

        Assert.That(table.Samples, Is.EqualTo(new[] { "S1" }));
        Assert.That(Cell(table, "PEPA", "S1"), Is.EqualTo(Math.Log2(400)));
        Assert.That(StateOf(table, "PEPA", "S1"), Is.EqualTo(ObservationState.Quantified));
    }

    [Test]
    public void TechnicalReplicatesAreAveragedOnTheLog2Scale()
    {
        var runs = new[]
        {
            new ObservationRun("t1", "S1", TechnicalReplicate: 1), new ObservationRun("t2", "S1", TechnicalReplicate: 2),
        };
        var table = Build(runs, QuantBasis.MsmsOnly, Peptide("PEPA", "P1", ("t1", Msms(100)), ("t2", Msms(400))));

        // log2 200, the log of the geometric mean; a linear mean would give log2 250.
        Assert.That(Cell(table, "PEPA", "S1"), Is.EqualTo((Math.Log2(100) + Math.Log2(400)) / 2));
        Assert.That(Cell(table, "PEPA", "S1"), Is.Not.EqualTo(Math.Log2(250)).Within(1e-6));
    }

    [Test]
    public void FractionsAreSummedWithinEachReplicateThenReplicatesAveraged()
    {
        var runs = new[]
        {
            new ObservationRun("t1f1", "S1", 1, 1), new ObservationRun("t1f2", "S1", 2, 1),
            new ObservationRun("t2f1", "S1", 1, 2), new ObservationRun("t2f2", "S1", 2, 2),
        };
        var table = Build(runs, QuantBasis.MsmsOnly, Peptide("PEPA", "P1",
            ("t1f1", Msms(100)), ("t1f2", Msms(300)), ("t2f1", Msms(50)), ("t2f2", Msms(150))));

        Assert.That(Cell(table, "PEPA", "S1"), Is.EqualTo((Math.Log2(400) + Math.Log2(200)) / 2));
    }

    [Test]
    public void AReplicateWithoutAValueIsLeftOutOfTheMean()
    {
        var runs = new[]
        {
            new ObservationRun("t1", "S1", TechnicalReplicate: 1), new ObservationRun("t2", "S1", TechnicalReplicate: 2),
        };
        var table = Build(runs, QuantBasis.MsmsOnly,
            Peptide("PEPA", "P1", ("t1", Msms(100)), ("t2", NoValue(ObservationState.NotDetected))));

        Assert.That(Cell(table, "PEPA", "S1"), Is.EqualTo(Math.Log2(100)));
    }

    [Test]
    public void MsmsOnlyExcludesEveryMbrValue()
    {
        var runs = new[]
        {
            new ObservationRun("a", "S1"),
            new ObservationRun("b1", "S2", Fraction: 1), new ObservationRun("b2", "S2", Fraction: 2),
            new ObservationRun("c", "S3"),
        };
        var peptides = new[]
        {
            Peptide("PEPA", "P1", ("a", Mbr(500)), ("b1", Msms(100)), ("b2", Mbr(300)), ("c", Msms(70))),
            Peptide("PEPB", "P1", ("a", Msms(20)), ("b1", Mbr(40)), ("c", Mbr(60))),
        };

        var msmsOnly = Build(runs, QuantBasis.MsmsOnly, peptides);
        Assert.That(Cell(msmsOnly, "PEPA", "S1"), Is.NaN);
        Assert.That(StateOf(msmsOnly, "PEPA", "S1"), Is.EqualTo(ObservationState.MbrTransferred));
        Assert.That(Cell(msmsOnly, "PEPA", "S2"), Is.EqualTo(Math.Log2(100)), "only the MS/MS fraction is summed");
        Assert.That(StateOf(msmsOnly, "PEPA", "S2"), Is.EqualTo(ObservationState.Quantified));
        Assert.That(msmsOnly.Observations().Where(o => !double.IsNaN(o.Log2Intensity)).Select(o => o.State),
            Is.All.EqualTo(ObservationState.Quantified), "no value under msms_only carries a transfer");

        var mbrKept = Build(runs, QuantBasis.MbrKept, peptides);
        Assert.That(Cell(mbrKept, "PEPA", "S1"), Is.EqualTo(Math.Log2(500)));
        Assert.That(StateOf(mbrKept, "PEPA", "S1"), Is.EqualTo(ObservationState.MbrTransferred));
        Assert.That(Cell(mbrKept, "PEPA", "S2"), Is.EqualTo(Math.Log2(400)));
        Assert.That(StateOf(mbrKept, "PEPA", "S2"), Is.EqualTo(ObservationState.MbrTransferred),
            "a value with any transferred part is a transfer");
        Assert.That(mbrKept.Basis, Is.EqualTo(QuantBasis.MbrKept));
    }

    [Test]
    public void ASampleWithNoValueIsKeptAndReported()
    {
        var runs = new[] { new ObservationRun("a", "S1"), new ObservationRun("b", "S2"), new ObservationRun("c", "S3") };
        var peptides = new[]
        {
            Peptide("PEPA", "P1", ("a", Msms(100)), ("c", Mbr(80))),
            Peptide("PEPB", "P1", ("a", Msms(200)), ("b", NoValue(ObservationState.IdentifiedNotQuantified))),
        };

        var msmsOnly = Build(runs, QuantBasis.MsmsOnly, peptides);
        Assert.That(msmsOnly.Samples, Is.EqualTo(new[] { "S1", "S2", "S3" }), "no sample is dropped");
        Assert.That(msmsOnly.SamplesWithoutValues, Is.EqualTo(new[] { "S2", "S3" }));
        Assert.That(Cell(msmsOnly, "PEPA", "S2"), Is.NaN);
        Assert.That(StateOf(msmsOnly, "PEPA", "S2"), Is.EqualTo(ObservationState.NotDetected));
        Assert.That(StateOf(msmsOnly, "PEPB", "S2"), Is.EqualTo(ObservationState.IdentifiedNotQuantified));

        var mbrKept = Build(runs, QuantBasis.MbrKept, peptides);
        Assert.That(mbrKept.SamplesWithoutValues, Is.EqualTo(new[] { "S2" }));
    }

    [Test]
    public void ACellWithoutAValueKeepsItsMostInformativeState()
    {
        var runs = new[]
        {
            new ObservationRun("a1", "S1", Fraction: 1), new ObservationRun("a2", "S1", Fraction: 2),
            new ObservationRun("b1", "S2", Fraction: 1), new ObservationRun("b2", "S2", Fraction: 2),
            new ObservationRun("c1", "S3", Fraction: 1), new ObservationRun("c2", "S3", Fraction: 2),
        };
        var table = Build(runs, QuantBasis.MsmsOnly, Peptide("PEPA", "P1",
            ("a1", NoValue(ObservationState.NotDetected)), ("a2", NoValue(ObservationState.IdentifiedNotQuantified)),
            ("b1", NoValue(ObservationState.AmbiguousPeak)), ("b2", NoValue(ObservationState.IdentifiedNotQuantified)),
            ("c1", Mbr(10)), ("c2", NoValue(ObservationState.AmbiguousPeak))));

        Assert.That(StateOf(table, "PEPA", "S1"), Is.EqualTo(ObservationState.IdentifiedNotQuantified));
        Assert.That(StateOf(table, "PEPA", "S2"), Is.EqualTo(ObservationState.AmbiguousPeak));
        Assert.That(StateOf(table, "PEPA", "S3"), Is.EqualTo(ObservationState.MbrTransferred));
        Assert.That(table.Observations().Select(o => o.Log2Intensity), Is.All.NaN, "a state without a value is never 0");
    }

    [Test]
    public void ARunTheProducerOmitsIsNotDetected()
    {
        var runs = new[] { new ObservationRun("a", "S1"), new ObservationRun("b", "S2") };
        var table = Build(runs, QuantBasis.MsmsOnly, Peptide("PEPA", "P1", ("a", Msms(5))));

        Assert.That(StateOf(table, "PEPA", "S2"), Is.EqualTo(ObservationState.NotDetected));
        Assert.That(Cell(table, "PEPA", "S2"), Is.NaN);
    }

    [TestCase("P1", "P1")]
    [TestCase("P1|P2", "P1|P2")]
    [TestCase("P1;P2", null)]
    [TestCase("UNDEFINED", null)]
    [TestCase("P1;UNDEFINED", null)]
    [TestCase("", null)]
    public void APeptideIsUniqueWhenItsCellNamesOneGroup(string cell, string? unique)
    {
        var runs = new[] { new ObservationRun("a", "S1") };
        var table = Build(runs, QuantBasis.MsmsOnly, Peptide("PEPA", cell, ("a", Msms(5))));

        Assert.That(table.Peptides.Single().UniqueProteinGroup, Is.EqualTo(unique));
    }

    [Test]
    public void GroupNamesAreDistinctAndOrdinalSorted()
    {
        var runs = new[] { new ObservationRun("a", "S1") };
        var table = Build(runs, QuantBasis.MsmsOnly, Peptide("PEPA", "b;B;a;b", ("a", Msms(5))));

        Assert.That(table.Peptides.Single().ProteinGroups, Is.EqualTo(new[] { "B", "a", "b" }));
    }

    [Test]
    public void ACalibratedColumnMatchesItsDesignFile()
    {
        var runs = new[] { new ObservationRun("A", "S1"), new ObservationRun("B", "S2") };
        var peptides = new[] { Peptide("PEPA", "P1", ("A-calib", Msms(7)), ("B", Msms(9))) };

        var table = ObservationTable.Build(runs, new[] { "A-calib", "B" }, peptides, Array.Empty<ProteinGroupInfo>(),
            QuantBasis.MsmsOnly);

        Assert.That(Cell(table, "PEPA", "S1"), Is.EqualTo(Math.Log2(7)));
        Assert.That(Cell(table, "PEPA", "S2"), Is.EqualTo(Math.Log2(9)));
    }

    [Test]
    public void AProducerColumnTheDesignDoesNotHaveIsRefused()
    {
        var runs = new[] { new ObservationRun("A", "S1") };
        var ex = Assert.Throws<ArgumentException>(() => ObservationTable.Build(runs, new[] { "A", "Z" },
            Array.Empty<PeptideRunValues>(), Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("'Z'"));
    }

    [Test]
    public void ADesignFileTheProducerDidNotWriteIsRefused()
    {
        var runs = new[] { new ObservationRun("A", "S1"), new ObservationRun("B", "S2") };
        var ex = Assert.Throws<ArgumentException>(() => ObservationTable.Build(runs, new[] { "A" },
            Array.Empty<PeptideRunValues>(), Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("'B'"));
    }

    [Test]
    public void TwoColumnsClaimingOneFileAreRefused()
    {
        var runs = new[] { new ObservationRun("A", "S1") };
        var ex = Assert.Throws<ArgumentException>(() => ObservationTable.Build(runs, new[] { "A", "A-calib" },
            Array.Empty<PeptideRunValues>(), Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("'A'"));
    }

    [Test]
    public void TwoRunsWithOneFileNameAreRefused()
    {
        var runs = new[] { new ObservationRun("A", "S1"), new ObservationRun("A", "S2") };
        var ex = Assert.Throws<ArgumentException>(() => Build(runs, QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("'A'"));
    }

    [Test]
    public void TwoRunsForOneFractionOfOneReplicateAreRefused()
    {
        var runs = new[] { new ObservationRun("A", "S1", 1, 1), new ObservationRun("B", "S1", 1, 1) };
        var ex = Assert.Throws<ArgumentException>(() => Build(runs, QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("S1"));
    }

    [Test]
    public void APeptideListedTwiceIsRefused()
    {
        var runs = new[] { new ObservationRun("a", "S1") };
        var ex = Assert.Throws<ArgumentException>(() => Build(runs, QuantBasis.MsmsOnly,
            Peptide("PEPA", "P1", ("a", Msms(5))), Peptide("PEPA", "P1", ("a", Msms(6)))));
        Assert.That(ex!.Message, Does.Contain("PEPA"));
    }

    [TestCase(0.0)]
    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void AValueThatIsNotFiniteAndPositiveIsRefused(double intensity)
    {
        var runs = new[] { new ObservationRun("a", "S1") };
        var ex = Assert.Throws<ArgumentException>(() =>
            Build(runs, QuantBasis.MbrKept, Peptide("PEPA", "P1", ("a", Mbr(intensity)))));
        Assert.That(ex!.Message, Does.Contain("PEPA").And.Contain("'a'"));
    }

    [Test]
    public void AGroupNameContainingASemicolonIsRefused()
    {
        var runs = new[] { new ObservationRun("a", "S1") };
        var peptide = new PeptideRunValues("PEPA", "PEPA", new[] { "P1;P2" },
            new Dictionary<string, RunValue> { ["a"] = Msms(5) });
        var ex = Assert.Throws<ArgumentException>(() => Build(runs, QuantBasis.MsmsOnly, peptide));
        Assert.That(ex!.Message, Does.Contain("P1;P2"));
    }

    [Test]
    public void TheTableDoesNotDependOnInputOrder()
    {
        var runs = new[]
        {
            new ObservationRun("a1", "S1", 1), new ObservationRun("a2", "S1", 2), new ObservationRun("b", "S2"),
        };
        var peptides = new[]
        {
            Peptide("PEPB", "P2;P1", ("a1", Msms(0.1)), ("a2", Msms(0.2)), ("b", Mbr(3))),
            Peptide("PEPA", "P1", ("a1", Msms(1e10)), ("a2", Msms(1)), ("b", Msms(7))),
        };

        var forward = Build(runs, QuantBasis.MbrKept, peptides);
        var reversed = ObservationTable.Build(runs.Reverse(), runs.Reverse().Select(r => r.FileName), peptides.Reverse(),
            Array.Empty<ProteinGroupInfo>(), QuantBasis.MbrKept);

        Assert.That(reversed.Peptides.Select(p => p.FullSequence), Is.EqualTo(new[] { "PEPA", "PEPB" }));
        Assert.That(reversed.Runs, Is.EqualTo(forward.Runs));
        Assert.That(reversed.Observations().Select(o => (o.Peptide.FullSequence, o.SampleId, o.Log2Intensity, o.State)),
            Is.EqualTo(forward.Observations().Select(o => (o.Peptide.FullSequence, o.SampleId, o.Log2Intensity, o.State))));
    }

    [Test]
    public void ForSamplesKeepsEveryPeptideAndOnlyTheNamedSamples()
    {
        var runs = new[] { new ObservationRun("a", "S1"), new ObservationRun("b", "S2"), new ObservationRun("c", "S3") };
        var table = Build(runs, QuantBasis.MsmsOnly,
            Peptide("PEPA", "P1", ("a", Msms(2))), Peptide("PEPB", "P1", ("b", Msms(4)), ("c", Msms(8))));

        var subset = table.ForSamples(new[] { "S3", "S1" });

        Assert.That(subset.Samples, Is.EqualTo(new[] { "S1", "S3" }));
        Assert.That(subset.Peptides.Select(p => p.FullSequence), Is.EqualTo(new[] { "PEPA", "PEPB" }));
        Assert.That(Cell(subset, "PEPB", "S3"), Is.EqualTo(3));
        Assert.That(subset.Runs.Select(r => r.FileName), Is.EqualTo(new[] { "a", "c" }));
        Assert.That(subset.SamplesWithoutValues, Is.Empty);
        Assert.Throws<ArgumentException>(() => table.ForSamples(new[] { "S9" }));
    }

    [Test]
    public void GroupsTheProteinTableDoesNotDescribeAreWarned()
    {
        var runs = new[] { new ObservationRun("a", "S1") };
        var groups = new[] { new ProteinGroupInfo("P1", 0.001, false, false) };
        var table = ObservationTable.Build(runs, new[] { "a" }, new[]
        {
            Peptide("PEPA", "P1;P9", ("a", Msms(5))), Peptide("PEPB", "UNDEFINED", ("a", Msms(6))),
        }, groups, QuantBasis.MsmsOnly);

        Assert.That(table.ProteinGroups.Keys, Is.EqualTo(new[] { "P1" }));
        Assert.That(table.Warnings, Has.Count.EqualTo(1));
        Assert.That(table.Warnings[0], Does.Contain("P9").And.Not.Contain("UNDEFINED"));
    }

    [Test]
    public void AProteinGroupDescribedTwiceDifferentlyIsRefused()
    {
        var runs = new[] { new ObservationRun("a", "S1") };
        var same = new[] { new ProteinGroupInfo("P1", 0.001, false, false), new ProteinGroupInfo("P1", 0.001, false, false) };
        Assert.DoesNotThrow(() => ObservationTable.Build(runs, new[] { "a" }, Array.Empty<PeptideRunValues>(), same,
            QuantBasis.MsmsOnly));

        var different = new[] { new ProteinGroupInfo("P1", 0.001, false, false), new ProteinGroupInfo("P1", 0.5, false, false) };
        var ex = Assert.Throws<ArgumentException>(() => ObservationTable.Build(runs, new[] { "a" },
            Array.Empty<PeptideRunValues>(), different, QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("P1"));
    }

    [Test]
    public void AGroupsAccessionsAreItsNameSplitAndOrdinalSorted()
    {
        Assert.That(new ProteinGroupInfo("Q2|P1|a3", 0, false, false).Accessions, Is.EqualTo(new[] { "P1", "Q2", "a3" }));
    }

    [Test]
    public void RunsFromALabelFreeDesignUseTheQuantifiedTablesSampleLabel()
    {
        var runs = ObservationRun.FromSpectraFiles(new[]
        {
            new SpectraFileInfo(@"C:\data\run7.mzML", "ctrl", biorep: 0, techrep: 1, fraction: 2),
        });

        Assert.That(runs.Single(), Is.EqualTo(new ObservationRun("run7", "ctrl_1", Fraction: 2, TechnicalReplicate: 1)));
    }

    [TestCase(QuantBasis.MsmsOnly, "msms_only")]
    [TestCase(QuantBasis.MbrKept, "mbr_kept")]
    public void BasesHaveTheirDefinitionNames(QuantBasis basis, string name) =>
        Assert.That(basis.ToMachineName(), Is.EqualTo(name));
}
