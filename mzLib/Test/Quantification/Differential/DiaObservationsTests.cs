using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// DIA precursors into the observation table (QuantProject dia 067, 072, 076, 078): one row per precursor, never summed
/// over charges; a value is <c>Quantified</c>, below the floor and not measured are <c>IdentifiedNotQuantified</c>, and a
/// precursor that fails either q-value in a run, or has no record there, is <c>NotDetected</c>; no cell is ever 0; decoy
/// and entrapment precursors never enter, whatever groups they name.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaObservationsTests
{
    private static readonly string[] Files = { "run1", "run2", "run3" };

    private static ObservationRun[] OneRunPerSample() => Files.Select(f => new ObservationRun(f, "s_" + f)).ToArray();

    private static readonly ProteinGroupInfo[] Groups =
    {
        new("G1", 0.001, false, false),
        new("G2", 0.002, false, false),
    };

    private static DiaPrecursorRun Run(string file, string sequence = "PEPTIDEK", int charge = 2, double quantity = 1000,
        PrecursorQuantityKind kind = PrecursorQuantityKind.Value, string groups = "G1", double runQ = 0.001,
        double globalQ = 0.001, bool isDecoy = false, bool isEntrapment = false, string? baseSequence = null,
        bool isContaminant = false) =>
        new(file, sequence, baseSequence ?? sequence, charge, groups.Split(';'),
            kind == PrecursorQuantityKind.Value ? quantity : double.NaN, kind, runQ, globalQ, isDecoy, isEntrapment, isContaminant);

    private static ObservationTable Build(IEnumerable<DiaPrecursorRun> precursors, ObservationRun[]? runs = null,
        double qValueThreshold = DiaObservations.DefaultQValueThreshold) =>
        DiaObservations.Build(runs ?? OneRunPerSample(), precursors, Groups, qValueThreshold);

    private static int Row(ObservationTable table, string sequence, int charge) =>
        Enumerable.Range(0, table.Peptides.Count).Single(p => table.Peptides[p].FullSequence == sequence && table.Charge(p) == charge);

    private static int Sample(ObservationTable table, string file) => table.Samples.ToList().IndexOf("s_" + file);

    [Test]
    public void EachKindMapsToItsStateAndNoCellIsZero()
    {
        var table = Build(new[]
        {
            Run("run1", quantity: 1000),
            Run("run2", kind: PrecursorQuantityKind.BelowFloor),
            Run("run3", kind: PrecursorQuantityKind.NotLooked),
            Run("run1", sequence: "OTHERK", quantity: 50),
        });

        int p = Row(table, "PEPTIDEK", 2);
        Assert.That(table.State(p, Sample(table, "run1")), Is.EqualTo(ObservationState.Quantified));
        Assert.That(table.Log2Intensity(p, Sample(table, "run1")), Is.EqualTo(Math.Log2(1000)).Within(1e-12));
        Assert.That(table.State(p, Sample(table, "run2")), Is.EqualTo(ObservationState.IdentifiedNotQuantified), "below the floor");
        Assert.That(table.State(p, Sample(table, "run3")), Is.EqualTo(ObservationState.IdentifiedNotQuantified), "not measured");
        int other = Row(table, "OTHERK", 2);
        Assert.That(table.State(other, Sample(table, "run2")), Is.EqualTo(ObservationState.NotDetected), "no record in the run");

        var noValue = table.Observations().Where(o => o.State != ObservationState.Quantified).ToList();
        Assert.That(noValue, Has.Count.EqualTo(4), "the fixture has four cells without a value");
        Assert.That(noValue.Select(o => o.Log2Intensity), Is.All.NaN);
        Assert.That(table.Basis, Is.EqualTo(QuantBasis.MsmsOnly), "DIA has no match-between-runs values (078)");
    }

    /// <summary>A kind added to <see cref="PrecursorQuantityKind"/> fails here until the adapter maps it.</summary>
    [Test]
    public void EveryPrecursorQuantityKindIsMappedAndAnUndefinedOneIsRefused()
    {
        foreach (PrecursorQuantityKind kind in Enum.GetValues<PrecursorQuantityKind>())
        {
            var table = Build(new[] { Run("run1", kind: kind) });
            Assert.That(table.State(0, Sample(table, "run1")), Is.Not.EqualTo(ObservationState.NotDetected), kind.ToString());
        }

        var ex = Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1") with { Kind = (PrecursorQuantityKind)99 } }));
        Assert.That(ex!.Message, Does.Contain("'PEPTIDEK'").And.Contain("99"));
    }

    [Test]
    public void EachPrecursorIsItsOwnRowNeverSummedOverCharges()
    {
        var table = Build(new[]
        {
            Run("run1", charge: 2, quantity: 100), Run("run1", charge: 3, quantity: 300),
            Run("run2", charge: 3, quantity: 30),
        });

        Assert.That(table.Peptides, Has.Count.EqualTo(2));
        Assert.That(Enumerable.Range(0, 2).Select(table.Charge), Is.EqualTo(new int?[] { 2, 3 }));
        Assert.That(table.Log2Intensity(Row(table, "PEPTIDEK", 2), Sample(table, "run1")), Is.EqualTo(Math.Log2(100)).Within(1e-12));
        Assert.That(table.Log2Intensity(Row(table, "PEPTIDEK", 3), Sample(table, "run1")), Is.EqualTo(Math.Log2(300)).Within(1e-12));
        Assert.That(table.State(Row(table, "PEPTIDEK", 2), Sample(table, "run2")), Is.EqualTo(ObservationState.NotDetected));
    }

    [Test]
    public void ARunThatFailsEitherQValueIsNotDetectedThere()
    {
        var precursors = new[]
        {
            Run("run1", runQ: 0.01, globalQ: 0.01),
            Run("run2", runQ: 0.02),
            Run("run3", globalQ: 0.02),
        };

        var table = Build(precursors);
        Assert.That(table.State(0, Sample(table, "run1")), Is.EqualTo(ObservationState.Quantified), "at the threshold passes");
        Assert.That(table.State(0, Sample(table, "run2")), Is.EqualTo(ObservationState.NotDetected), "run q above 1%");
        Assert.That(table.State(0, Sample(table, "run3")), Is.EqualTo(ObservationState.NotDetected), "global q above 1%");

        var looser = Build(precursors, qValueThreshold: 0.05);
        Assert.That(Files.Select(f => looser.State(0, Sample(looser, f))), Is.All.EqualTo(ObservationState.Quantified));

        Assert.That(Build(new[] { Run("run1", runQ: 0.5) }).Peptides, Is.Empty, "a precursor passing in no run has no row");
    }

    /// <summary>
    /// A decoy, entrapment or contaminant precursor that names a target group must not reach it: groups are excluded by
    /// their own flags, but a precursor's group name is the producer's, and a permuted decoy can carry its owner's name.
    /// </summary>
    [Test]
    public void DecoyEntrapmentAndContaminantPrecursorsNeverEnterEvenUnderATargetGroup()
    {
        var decoy = Run("run1", sequence: "KEDITPEP", isDecoy: true);
        var entrapment = Run("run2", sequence: "PEPTLDEK", isEntrapment: true);
        var contaminant = Run("run3", sequence: "KERATINK", isContaminant: true);
        Assert.That(new[] { decoy, entrapment, contaminant }.SelectMany(p => p.ProteinGroups), Is.All.EqualTo("G1"),
            "each names the target group");
        Assert.That(Groups.Single(g => g.Name == "G1"), Has.Property("IsDecoy").False.And.Property("IsContaminant").False);

        var table = Build(new[] { Run("run1"), Run("run2"), Run("run3"), decoy, entrapment, contaminant });

        Assert.That(table.Peptides.Select(p => p.FullSequence), Is.EqualTo(new[] { "PEPTIDEK" }));
    }

    [Test]
    public void APrecursorTwiceInOneRunIsRefused()
    {
        var ex = Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1"), Run("run1", quantity: 5) }));

        Assert.That(ex!.Message, Does.Contain("'PEPTIDEK'").And.Contain("charge 2").And.Contain("'run1'").And.Contain("twice"));
    }

    [Test]
    public void APrecursorsRunsMustAgreeOnItsBaseSequenceAndGroups()
    {
        Assert.That(Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1", groups: "G1"), Run("run2", groups: "G2") }))!
            .Message, Does.Contain("'PEPTIDEK'").And.Contain("'G1'").And.Contain("'G2'").And.Contain("'run2'"));
        Assert.That(Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1"), Run("run2", baseSequence: "PEPTIDEKK") }))!
            .Message, Does.Contain("'PEPTIDEK'").And.Contain("'PEPTIDEKK'"));
        Assert.That(Build(new[] { Run("run1", groups: "G1;G2"), Run("run2", groups: "G2;G1") }).Peptides, Has.Count.EqualTo(1),
            "the same groups in another order agree");
    }

    [TestCase(PrecursorQuantityKind.Value, 0.0)]
    [TestCase(PrecursorQuantityKind.Value, -1.0)]
    [TestCase(PrecursorQuantityKind.Value, double.NaN)]
    [TestCase(PrecursorQuantityKind.Value, double.PositiveInfinity)]
    [TestCase(PrecursorQuantityKind.BelowFloor, 5.0)]
    [TestCase(PrecursorQuantityKind.NotLooked, 0.0)]
    public void AValueMustBeFiniteAndAboveZeroAndAKindWithoutOneCarriesNoNumber(PrecursorQuantityKind kind, double quantity)
    {
        var ex = Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1") with { Kind = kind, Quantity = quantity } }));

        Assert.That(ex!.Message, Does.Contain("'PEPTIDEK'").And.Contain("'run1'").And.Contain(kind.ToString()));
    }

    [TestCase(double.NaN, 0.001)]
    [TestCase(0.001, double.NaN)]
    [TestCase(-0.1, 0.001)]
    [TestCase(0.001, 1.5)]
    public void AQValueThatIsNotBetweenZeroAndOneIsRefused(double runQ, double globalQ)
    {
        var ex = Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1", runQ: runQ, globalQ: globalQ) }));

        Assert.That(ex!.Message, Does.Contain("'PEPTIDEK'").And.Contain("'run1'").And.Contain("q-value"));
    }

    [TestCase(0.0)]
    [TestCase(1.5)]
    [TestCase(double.NaN)]
    public void AThresholdOutsideZeroToOneIsRefused(double threshold)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Build(new[] { Run("run1") }, qValueThreshold: threshold));
        Assert.That(Build(new[] { Run("run1", runQ: 0.9, globalQ: 1) }, qValueThreshold: 1).Peptides, Has.Count.EqualTo(1),
            "a threshold of 1 keeps every target");
    }

    [Test]
    public void TechnicalReplicatesAreAveragedOnTheLogScaleAndFractionsAreRefused()
    {
        var replicates = new[] { new ObservationRun("run1", "s1", 0, 1), new ObservationRun("run2", "s1", 0, 2), new ObservationRun("run3", "s3") };
        var table = Build(new[] { Run("run1", quantity: 100), Run("run2", quantity: 400) }, replicates);
        Assert.That(table.Log2Intensity(0, table.Samples.ToList().IndexOf("s1")), Is.EqualTo((Math.Log2(100) + Math.Log2(400)) / 2).Within(1e-12));

        var fractions = new[] { new ObservationRun("run1", "s1", 1), new ObservationRun("run2", "s1", 2), new ObservationRun("run3", "s3") };
        var ex = Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1") }, fractions));
        Assert.That(ex!.Message, Does.Contain("'run1'").And.Contain("'run2'").And.Contain("fraction"));
    }

    [Test]
    public void ADesignRunWithNothingReportedIsASampleWithoutValues()
    {
        var table = Build(new[] { Run("run1"), Run("run2", runQ: 0.2) });

        Assert.That(table.SamplesWithoutValues, Is.EqualTo(new[] { "s_run2", "s_run3" }));
    }

    /// <summary>
    /// A run that adds nothing to a sample whose other runs do would drop out of the replicate average with no record
    /// anywhere (QuantProject 086), so it is refused, naming it and its sample; a sample with nothing in any run is kept.
    /// </summary>
    [Test]
    public void ARunLostInsideASampleIsRefusedButAWholeSampleWithNothingIsKept()
    {
        var replicates = new[] { new ObservationRun("run1", "s1", 0, 1), new ObservationRun("run2", "s1", 0, 2), new ObservationRun("run3", "s3") };

        foreach (var precursors in new[]
                 {
                     new[] { Run("run1"), Run("run3") },
                     new[] { Run("run1"), Run("run2", runQ: 0.2), Run("run3") },
                     new[] { Run("run1"), Run("run2", isDecoy: true), Run("run3") },
                 })
        {
            var ex = Assert.Throws<ArgumentException>(() => Build(precursors, replicates));
            Assert.That(ex!.Message, Does.Contain("'run2'").And.Contain("'s1'").And.Contain("'run1'"));
        }

        var nothingInSampleOne = Build(new[] { Run("run3") }, replicates);
        Assert.That(nothingInSampleOne.SamplesWithoutValues, Is.EqualTo(new[] { "s1" }));
    }

    [Test]
    public void APrecursorInARunTheDesignDoesNotHaveIsRefused()
    {
        var ex = Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1.raw") }));

        Assert.That(ex!.Message, Does.Contain("'run1.raw'"));
    }

    [Test]
    public void ANullPrecursorOrRunIsRefused()
    {
        Assert.That(Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1"), null! }))!.Message, Does.Contain("precursor is null"));
        Assert.That(Assert.Throws<ArgumentException>(() => Build(new[] { Run("run1") }, OneRunPerSample().Append(null!).ToArray()))!
            .Message, Does.Contain("run is null"));
    }

    /// <summary>MetaMorpheus quantifies calibrated files, so a record's file may carry <c>-calib</c>; it is that design run.</summary>
    [Test]
    public void ACalibratedFileIsItsDesignRun()
    {
        var table = Build(new[] { Run("run1-calib", quantity: 100), Run("run2-calib", quantity: 200) });

        Assert.That(table.Log2Intensity(0, Sample(table, "run1")), Is.EqualTo(Math.Log2(100)).Within(1e-12));
        Assert.That(table.Log2Intensity(0, Sample(table, "run2")), Is.EqualTo(Math.Log2(200)).Within(1e-12));
        Assert.That(table.SamplesWithoutValues, Is.EqualTo(new[] { "s_run3" }));

        var replicates = new[] { new ObservationRun("run1", "s1", 0, 1), new ObservationRun("run2", "s1", 0, 2), new ObservationRun("run3", "s3") };
        var averaged = Build(new[] { Run("run1-calib", quantity: 100), Run("run2-calib", quantity: 400) }, replicates);
        Assert.That(averaged.Log2Intensity(0, averaged.Samples.ToList().IndexOf("s1")),
            Is.EqualTo((Math.Log2(100) + Math.Log2(400)) / 2).Within(1e-12), "both replicates of s1 report, under their calibrated names");
    }

    [Test]
    public void APrecursorWithoutGroupsHasNone()
    {
        var table = Build(new[] { Run("run1") with { ProteinGroups = null! } });

        Assert.That(table.Peptides.Single().ProteinGroups, Is.Empty);
        Assert.That(table.Peptides.Single().UniqueProteinGroup, Is.Null);
    }

    [Test]
    public void TheTableDoesNotDependOnInputOrder()
    {
        var precursors = new List<DiaPrecursorRun>();
        foreach (string file in Files)
            for (int k = 0; k < 4; k++)
                precursors.Add(Run(file, sequence: $"PEP{k}K", charge: 2 + k % 2, quantity: 100 * (k + 1) + file.Length,
                    groups: k < 2 ? "G1" : "G2"));
        precursors.Add(Run("run2", sequence: "PEP0K", charge: 3, kind: PrecursorQuantityKind.BelowFloor));

        string forward = Build(precursors).Sha256();
        precursors.Reverse();
        Assert.That(Build(precursors).Sha256(), Is.EqualTo(forward));
    }

    /// <summary>
    /// The adapter's table runs as a <c>dia</c> request (078): built under <c>msms_only</c>, so the runner does not refuse it,
    /// and one peptide at two charges is two levels and one peptide.
    /// </summary>
    [Test]
    public void TheRunnerFitsTheAdaptersTableAsADiaRequest()
    {
        string[] samples = { "a1", "a2", "a3", "b1", "b2", "b3" };
        var runs = samples.Select(s => new ObservationRun("f_" + s, s)).ToArray();
        // Ten unchanged groups (U) and ten that are one log2 unit higher in B (C), each with peptide A at charges 2 and 3
        // and peptide B at charge 2.
        string[] groups = Enumerable.Range(0, 20).Select(g => g < 10 ? $"U{g:D2}" : $"C{g:D2}").ToArray();
        var rng = new Random(20261124);
        var precursors = new List<DiaPrecursorRun>();
        foreach (string group in groups)
            foreach (var (sequence, charge) in new[] { ($"{group}AK", 2), ($"{group}AK", 3), ($"{group}BK", 2) })
                foreach (string s in samples)
                {
                    double effect = group[0] == 'C' && s[0] == 'b' ? 1 : 0;
                    double log2 = 20 + charge + effect + 0.1 * (rng.NextDouble() - 0.5);
                    precursors.Add(Run("f_" + s, sequence: sequence, charge: charge, quantity: Math.Pow(2, log2), groups: group));
                }

        var result = DifferentialAnalysis.Run(new DifferentialRequest
        {
            Table = DiaObservations.Build(runs, precursors, groups.Select(g => new ProteinGroupInfo(g, 0.001, false, false))),
            Samples = samples.Select(s => new AnalysisSample(s) { Factors = new Dictionary<string, string> { ["mix"] = s[..1].ToUpperInvariant() } }).ToList(),
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("mix", "A") } },
            Contrasts = new[] { DifferentialContrastRequest.Between("mix", "B", "A") },
            QuantStyle = "dia",
            Settings = new DifferentialSettings { Normalization = SharedPeptideNormalization.None },
        }, maxThreads: 1);

        Assert.That(result.Rows, Has.Count.EqualTo(20));
        Assert.That(result.Rows.Select(r => r.Status), Is.All.EqualTo(DifferentialStatus.Fitted));
        Assert.That(result.Rows.Select(r => r.NPeptides), Is.All.EqualTo(2), "two peptides, three precursors");
        Assert.That(result.Rows.Select(r => r.NObservations), Is.All.EqualTo(18));
        Assert.That(result.Rows.Where(r => r.FeatureId[0] == 'C').Select(r => r.Log2Effect!.Value), Is.All.EqualTo(1).Within(0.1));
        Assert.That(result.Rows.Where(r => r.FeatureId[0] == 'U').Select(r => r.Log2Effect!.Value), Is.All.EqualTo(0).Within(0.1));
    }
}
