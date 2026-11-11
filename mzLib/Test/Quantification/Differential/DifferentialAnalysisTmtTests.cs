using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The STAT1 runner on an isobaric table (model A, REQ-QP-15 PLAN 6h): two plexes of seven channels (three A, three B, a
/// pooled reference), several ions per protein, B = A + 1 in log2 for three proteins against 100 unchanged ones, each plex
/// run at its own level. The ions in their plex runs absorb the plex levels, the reference channel is no sample, and the
/// effect is recovered after normalization, which needs the unchanged background to align on.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialAnalysisTmtTests
{
    private static readonly string[] Labels = { "126", "127N", "127C", "128N", "128C", "129N", "129C" };

    private static (IsobaricObservationTable Table, IReadOnlyList<AnalysisSample> Samples) Data(int seed)
    {
        var rng = new Random(seed);
        var design = new SampleExperimentalDesign();
        foreach (int plex in new[] { 1, 2 })
        {
            string file = $"plex{plex}";
            var channels = Labels.Select((label, i) => (ISampleInfo)new IsobaricQuantSampleInfo($@"C:\data\{file}.raw",
                i == 6 ? "pool" : i < 3 ? "A" : "B", 0, 0, 0, plex, label, 126 + i, isReferenceChannel: i == 6)
            { SampleName = i == 6 ? $"pool{plex}" : $"{(i < 3 ? "a" : "b")}{plex}_{i % 3}" }).ToArray();
            design.Add($@"C:\data\{file}.raw", channels);
        }
        var psms = new List<IsobaricPsm>();
        int scan = 0;
        var proteins = new[] { "P1", "P2", "P3" }.Concat(Enumerable.Range(1, 100).Select(i => $"BG{i:000}")).ToArray();
        foreach (string protein in proteins)
            for (int ion = 0; ion < 4; ion++)
                foreach (int plex in new[] { 1, 2 })
                {
                    double ionLevel = 10 + 6 * Level(protein, ion) + plex * 0.4; // each ion in each plex run at its own level
                    var values = Labels.Select((_, i) => i == 6 ? Math.Pow(2, ionLevel + 0.5)
                        : Math.Pow(2, ionLevel + (i >= 3 && protein.StartsWith('P') ? 1 : 0) + 0.05 * Gaussian(rng))).ToArray();
                    psms.Add(new IsobaricPsm($"plex{plex}", ++scan, $"{protein}ION{ion}", $"{protein}ION{ion}", 2,
                        new[] { protein }, values));
                }
        var groups = proteins.Select(g => new ProteinGroupInfo(g, 0, false, false));
        var table = IsobaricObservationTable.Build(design, psms, groups);
        return (table, AnalysisSample.FromExperimentalDesign(design));
    }

    /// <summary>A fixed level in [0, 1) per protein and ion, spread continuously as real intensities are.</summary>
    private static double Level(string protein, int ion) => new Random(protein.Aggregate(17, (h, c) => unchecked(h * 31 + c)) * 7 + ion).NextDouble(); // a fixed seed, never string.GetHashCode

    private static double Gaussian(Random rng) =>
        Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

    private static DifferentialAnalysisResult Run(int threads = 1)
    {
        var (table, samples) = Data(20261110);
        return DifferentialAnalysis.Run(new DifferentialRequest
        {
            IsobaricTable = table,
            Samples = samples,
            Spec = new DesignSpec { Factors = new[] { new FactorSpec(AnalysisSample.ConditionFactor, "A") } },
            Contrasts = new[] { DifferentialContrastRequest.Between(AnalysisSample.ConditionFactor, "B", "A") },
            ReporterAcquisition = "ms3_sps",
        }, threads);
    }

    /// <summary>
    /// The isobaric table's canonical sha256 (ruling 10): the same data rebuilt give the same hash, one changed cell another,
    /// and the runner records it as its first input, so the analysis id follows the data.
    /// </summary>
    [Test]
    public void TheIsobaricTableIsTheFirstInputAndItsHashFollowsTheData()
    {
        var (table, _) = Data(20261110);
        string hash = table.Sha256();
        Assert.That(Data(20261110).Table.Sha256(), Is.EqualTo(hash), "the same data, rebuilt");
        var changed = table.WithLog2((r, s) => r == 0 && s == table.SamplesOfPlex(table.Rows[0].Plex)[0]
            ? table.Log2Intensity(r, s) + 0.001 : table.Log2Intensity(r, s));
        Assert.That(changed.Sha256(), Is.Not.EqualTo(hash), "one cell different");

        string text = table.ToCanonicalText();
        Assert.That(text, Does.StartWith("observation_table\tisobaric\nsamples\t"));
        Assert.That(text, Does.Contain("\nion\tBG001ION0\tBG001ION0\tBG001\t2\t1\t0\tplex1\t"));
        Assert.That(text, Does.Contain("\nions_with_few_channels\t0\n"));

        var input = Run().Metadata.Inputs.First();
        Assert.That(input, Is.EqualTo(new DifferentialInput(DifferentialConstants.ObservationTableRole, "", hash)));
    }

    /// <summary>
    /// The global shift (ruling 9c) on the raw table: every ion of every plex run is valued in all six sample channels, so
    /// all 824 count (103 proteins x 4 ions x 2 plexes), and with 100 of 103 proteins unchanged the median is near 0.
    /// </summary>
    [Test]
    public void TheGlobalShiftIsTakenPerIonWithinItsPlexRun()
    {
        var shift = Run().Metadata.Contrasts.Single().GlobalShift!.Single();

        Assert.That(shift.Peptides, Is.EqualTo(824));
        Assert.That(shift.Value, Is.EqualTo(0).Within(0.05));
        Assert.That(shift.QuantBasis, Is.Null);
    }

    [Test]
    public void ModelARecoversTheEffectAcrossPlexesWithTheReferenceLeftOut()
    {
        var result = Run();

        Assert.That(result.Rows, Has.Count.EqualTo(103));
        Assert.That(result.Rows.Where(r => r.FeatureId.StartsWith("BG")).Select(r => r.Log2Effect!.Value),
            Is.All.EqualTo(0).Within(0.1), "the unchanged background");
        foreach (var row in result.Rows.Where(r => r.FeatureId.StartsWith('P')))
        {
            Assert.That(row.Status, Is.EqualTo(DifferentialStatus.Fitted), row.FeatureId);
            Assert.That(row.Log2Effect, Is.EqualTo(1).Within(0.1), row.FeatureId);
            Assert.That(row.QuantStyle, Is.EqualTo("tmt"));
            Assert.That(row.ReporterAcquisition, Is.EqualTo("ms3_sps"));
            Assert.That(row.QuantBasis, Is.Null);
            Assert.That((row.NSamplesNumerator, row.NSamplesDenominator), Is.EqualTo((6, 6)), "three per group per plex; no reference");
            Assert.That(row.NPeptides, Is.EqualTo(8), "four ions in each of two plex runs");
        }
        var model = result.Metadata.Models.Single();
        Assert.That(model.Formula, Does.Contain("ion_in_run").And.Contain("(1|sample/run)"));
        Assert.That(model.Formula, Does.Not.Contain("batch"), "plex is never a batch column for model A");
        Assert.That(model.Prior!.Variance, Is.Not.Null.And.GreaterThan(0), "the prior variance, from whichever part has a fit");
        Assert.That(result.Metadata.Normalization.Single().Setting, Is.EqualTo(IsobaricNormalization.GlobalChannelMedian));
    }

    [Test]
    public void TheIsobaricRowsPassTheWritersValidationAndDoNotDependOnTheThreadCount()
    {
        string Written(int threads)
        {
            var text = new StringWriter();
            DifferentialResultWriter.Write(text, Run(threads).Rows, DifferentialHeaderStyle.Machine);
            return text.ToString();
        }

        Assert.That(Written(4), Is.EqualTo(Written(1)));
    }
}
