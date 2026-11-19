using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The global shift (GR-27; STAT1 M5b ruling 9): the median over peptides (isobaric: ions in their plex runs) of the
/// numerator side's mean log2 minus the denominator side's, on values before normalization; decoy, contaminant and
/// entrapment peptides never count (9a); a lost run, or a channel with no value in a plex run, is skipped (9b, 9c).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class GlobalShiftTests
{
    private static readonly string[] Numerator = { "n1", "n2" };
    private static readonly string[] Denominator = { "d1", "d2", "d3" };

    private static PeptideRunValues Peptide(string sequence, string group, params (string Sample, double Log2)[] values) =>
        new(sequence, sequence, new[] { group }, values.ToDictionary(v => v.Sample,
            v => new RunValue(Math.Pow(2, v.Log2), ObservationState.Quantified)));

    /// <summary>
    /// Three peptides count: 2.0, 1.2 and 1.0, so the shift is 1.2 over 3. One is missing a numerator sample, one belongs
    /// to a decoy, one to a contaminant (each would move the median). d3, a run with no value at all, is skipped: else no
    /// peptide is valued in every denominator sample and there is no shift.
    /// </summary>
    private static ObservationTable LabelFreeTable()
    {
        var peptides = new List<PeptideRunValues>
        {
            Peptide("PA", "P1", ("n1", 22), ("n2", 22.2), ("d1", 20), ("d2", 20.2)),
            Peptide("PB", "P1", ("n1", 21), ("n2", 21.4), ("d1", 20), ("d2", 20)),
            Peptide("PC", "P2", ("n1", 25), ("n2", 25), ("d1", 24), ("d2", 24)),
            Peptide("PMISSING", "P2", ("n1", 30), ("d1", 20), ("d2", 20)),
            Peptide("PDECOY", "DECOY", ("n1", 30), ("n2", 30), ("d1", 20), ("d2", 20)),
            Peptide("PCONT", "CONT", ("n1", 30), ("n2", 30), ("d1", 20), ("d2", 20)),
            Peptide("PTRAP", "TRAP", ("n1", 30), ("n2", 30), ("d1", 20), ("d2", 20)),
        };
        var groups = new[]
        {
            new ProteinGroupInfo("P1", 0, false, false), new ProteinGroupInfo("P2", 0, false, false),
            new ProteinGroupInfo("DECOY", 0, true, false), new ProteinGroupInfo("CONT", 0, false, true),
            new ProteinGroupInfo("TRAP", 0, false, false, IsEntrapment: true),
        };
        var runs = Numerator.Concat(Denominator).Select(s => new ObservationRun(s, s)).ToList();
        return ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, groups, QuantBasis.MsmsOnly);
    }

    [Test]
    public void LabelFreeIsTheMedianOverPeptidesValuedOnEveryKeptSample()
    {
        var shift = GlobalShift.LabelFree(LabelFreeTable(), Numerator, Denominator, "all", "msms_only");

        Assert.That(shift.Value, Is.EqualTo(1.2).Within(1e-12));
        Assert.That(shift.Peptides, Is.EqualTo(3));
        Assert.That((shift.Stratum, shift.QuantBasis), Is.EqualTo(("all", "msms_only")));
    }

    [Test]
    public void WithNoPeptideOnBothSidesTheShiftIsNaNOverNone()
    {
        var shift = GlobalShift.LabelFree(LabelFreeTable(), new[] { "d3" }, Denominator, "all", "msms_only");

        Assert.That(shift.Value, Is.NaN);
        Assert.That(shift.Peptides, Is.EqualTo(0));
    }

    /// <summary>
    /// Two plexes of five channels (A, A, B, B, pool). Plex 1: ions at +1 and +2; an ion missing an A channel and a decoy
    /// ion (+10) do not count. Plex 2: channel 128N has no value at all there, so its B side is 127C alone, and its ion
    /// (+3) counts. The shift is the median of 1, 2, 3 = 2 over 3 ions.
    /// </summary>
    [Test]
    public void IsobaricIsTheMedianOverIonsWithinTheirPlexRuns()
    {
        string[] labels = { "126", "127N", "127C", "128N", "128C" };
        var design = new SampleExperimentalDesign();
        foreach (int plex in new[] { 1, 2 })
            design.Add($@"C:\data\plex{plex}.raw", labels.Select((label, i) => (ISampleInfo)new IsobaricQuantSampleInfo(
                $@"C:\data\plex{plex}.raw", i == 4 ? "pool" : i < 2 ? "A" : "B", 0, 0, 0, plex, label, 126 + i, isReferenceChannel: i == 4)
            { SampleName = i == 4 ? $"pool{plex}" : $"{(i < 2 ? "a" : "b")}{plex}_{i % 2}" }).ToArray());
        double P(double log2) => log2 == 0 ? 0 : Math.Pow(2, log2);
        IsobaricPsm Ion(int plex, int scan, string sequence, string group, params double[] log2) =>
            new($"plex{plex}", scan, sequence, sequence, 2, new[] { group }, log2.Select(P).ToArray());
        var psms = new[]
        {
            Ion(1, 1, "IONA", "P1", 20, 20, 21, 21, 20),
            Ion(1, 2, "IONB", "P1", 18, 18.2, 20.1, 20.1, 20),
            Ion(1, 3, "IONMISSING", "P1", 20, 0, 25, 25, 20),
            Ion(1, 4, "IONDECOY", "DECOY", 20, 20, 30, 30, 20),
            Ion(2, 5, "IONC", "P2", 19, 19, 22, 0, 20),
            Ion(2, 6, "IOND", "P2", 19, 19.4, 0, 0, 20),
        };
        var groups = new[]
        {
            new ProteinGroupInfo("P1", 0, false, false), new ProteinGroupInfo("P2", 0, false, false),
            new ProteinGroupInfo("DECOY", 0, true, false),
        };
        var table = IsobaricObservationTable.Build(design, psms, groups);
        var samples = AnalysisSample.FromExperimentalDesign(design);
        string[] Side(string level) => samples.Where(s => s.Factors[AnalysisSample.ConditionFactor] == level).Select(s => s.SampleId).ToArray();

        var shift = GlobalShift.Isobaric(table, Side("B"), Side("A"), "all");

        Assert.That(shift.Value, Is.EqualTo(2).Within(1e-12));
        Assert.That(shift.Peptides, Is.EqualTo(3));
        Assert.That(shift.QuantBasis, Is.Null);
    }

    /// <summary>
    /// End to end: the runner writes each two-sided contrast's shift into the metadata, per stratum, on raw values. In the
    /// label-free fixture every peptide valued in all A and B samples sits 1 higher in B (P1's, P2's, P4's, P6's, the shared
    /// one's and the weak group's: 10; the decoy's does not count, P3's lack A).
    /// </summary>
    [Test]
    public void TheRunnerWritesEachContrastsShiftIntoTheMetadata()
    {
        var result = DifferentialAnalysis.Run(DifferentialAnalysisTests.Request(), maxThreads: 1);

        var c1 = result.Metadata.Contrasts.Single(c => c.Id == "c1").GlobalShift!.Single();
        Assert.That(c1.Value, Is.EqualTo(1).Within(0.1));
        Assert.That(c1.Peptides, Is.EqualTo(10));
        Assert.That((c1.Stratum, c1.QuantBasis), Is.EqualTo((AnalysisDesign.AllSamples, "msms_only")));
    }

    /// <summary>
    /// The shift is taken BEFORE normalization (GR-27). B is loaded 1 log2 higher than A on every one of 120 peptides (40
    /// proteins), so the shared-peptide normalization moves every B sample down by about 1 and every fitted effect comes
    /// out near 0; the shift still reads about 1, which is what it is for: showing what the normalization took out.
    /// </summary>
    [Test]
    public void TheShiftIsMeasuredBeforeNormalization()
    {
        var request = LoadedRequest();

        var result = DifferentialAnalysis.Run(request, maxThreads: 1);

        Assert.That(result.Metadata.Normalization.Single().Setting, Is.EqualTo(SharedPeptideNormalization.SharedPeptideMedian));
        Assert.That(result.Rows.Where(r => r.Status == DifferentialStatus.Fitted).Select(r => r.Log2Effect!.Value),
            Is.All.EqualTo(0).Within(0.05), "normalization took B's loading out");
        var shift = result.Metadata.Contrasts.Single().GlobalShift!.Single();
        Assert.That(shift.Value, Is.EqualTo(1).Within(0.05));
        Assert.That(shift.Peptides, Is.EqualTo(120));
    }

    /// <summary>
    /// The normalization floor is a setting (GR-25, ruling 14; default 100): with 200, the 120 shared peptides are too few,
    /// no sample is shifted, and B's loading stays in every effect.
    /// </summary>
    [Test]
    public void TheNormalizationFloorIsASetting()
    {
        var request = LoadedRequest();

        var result = DifferentialAnalysis.Run(request with { Settings = request.Settings with { NormalizationFloor = 200 } }, maxThreads: 1);

        Assert.That(result.Metadata.Normalization.Single().Setting, Is.EqualTo(SharedPeptideNormalization.None));
        Assert.That(result.Rows.Where(r => r.Status == DifferentialStatus.Fitted).Select(r => r.Log2Effect!.Value),
            Is.All.EqualTo(1).Within(0.05));
    }

    /// <summary>B loaded 1 log2 higher than A on every one of 120 peptides of 40 proteins, three samples a side.</summary>
    private static DifferentialRequest LoadedRequest()
    {
        var rng = new Random(11);
        string[] ids = { "a1", "a2", "a3", "b1", "b2", "b3" };
        var peptides = new List<PeptideRunValues>();
        for (int g = 0; g < 40; g++)
            for (int k = 0; k < 3; k++)
            {
                double offset = 0.1 * g + 0.4 * k;
                peptides.Add(new PeptideRunValues($"G{g}PEP{k}", $"G{g}PEP{k}", new[] { $"G{g}" }, ids.ToDictionary(s => s,
                    s => new RunValue(Math.Pow(2, 20 + offset + (s[0] == 'b' ? 1 : 0) + 0.05 * (rng.NextDouble() - 0.5)), ObservationState.Quantified))));
            }
        var runs = ids.Select(s => new ObservationRun(s, s)).ToList();
        var request = new DifferentialRequest
        {
            Table = ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides,
                Enumerable.Range(0, 40).Select(g => new ProteinGroupInfo($"G{g}", 0, false, false)), QuantBasis.MsmsOnly),
            Samples = ids.Select(s => new AnalysisSample(s) { Factors = new Dictionary<string, string> { ["group"] = s[0] == 'a' ? "A" : "B" } }).ToList(),
            Spec = new DesignSpec { Factors = new[] { new FactorSpec("group", "A") } },
            Contrasts = new[] { DifferentialContrastRequest.Between("group", "B", "A") },
        };
        return request;
    }
}
