using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The runner on a DIA-style table and the normalization settings (QuantProject DATA-DEFINITIONS v3.29; dia 072-074): a
/// <c>dia</c> row has no quant basis and no MBR count, counts peptides rather than precursors, and fits each precursor as
/// its own level; a setting the runner does not know is refused; and <c>background_set:&lt;name&gt;</c> shifts each sample
/// on the protein groups the caller declares.
/// </summary>
/// <remarks>
/// The data: six samples (A and B, three each), each with its own loading offset. Twenty "human" groups are unchanged and
/// forty "yeast" groups are one log2 unit higher in B. Each group has four peptides, every one at charge 2 and the even
/// ones also at charge 3, so six precursors. ONE has a single peptide at charges 2 and 3. Yeast outnumbers human, so the
/// shared-peptide median aligns the samples on yeast and moves human, and only the declared human set gets both right.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialDiaTests
{
    private static readonly string[] Ids = { "a1", "a2", "a3", "b1", "b2", "b3" };
    private static readonly double[] Loading = { 0.3, -0.2, 0.1, 0.4, -0.1, 0.0 };
    private static readonly string[] Human = Enumerable.Range(0, 20).Select(g => $"H{g:D2}").ToArray();
    private static readonly string[] Yeast = Enumerable.Range(0, 40).Select(g => $"Y{g:D2}").ToArray();

    private static double Gaussian(Random rng) =>
        Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

    private static PeptideRunValues Precursor(Random rng, string sequence, string group, int charge, double level, double effectInB) =>
        new(sequence, sequence, new[] { group }, Ids.ToDictionary(s => s, s =>
        {
            int i = Array.IndexOf(Ids, s);
            double log2 = 20 + level + Loading[i] + (i >= 3 ? effectInB : 0) + 0.05 * Gaussian(rng);
            return new RunValue(Math.Pow(2, log2), ObservationState.Quantified);
        }), charge);

    private static ObservationTable Table(QuantBasis basis = QuantBasis.MsmsOnly)
    {
        var rng = new Random(20261116);
        var rows = new List<PeptideRunValues>();
        foreach (var (groups, effect) in new[] { (Human, 0.0), (Yeast, 1.0) })
            foreach (string g in groups)
                for (int k = 0; k < 4; k++)
                    foreach (int z in k % 2 == 0 ? new[] { 2, 3 } : new[] { 2 })
                        rows.Add(Precursor(rng, $"{g}PEP{k}", g, z, 0.7 * k + 0.2 * z, effect));
        rows.Add(Precursor(rng, "ONEPEP", "ONE", 2, 0.5, 0));
        rows.Add(Precursor(rng, "ONEPEP", "ONE", 3, 0.9, 0));
        var groupInfo = Human.Concat(Yeast).Append("ONE").Select(g => new ProteinGroupInfo(g, 0, false, false));
        var runs = Ids.Select(s => new ObservationRun(s, s)).ToList();
        return ObservationTable.Build(runs, Ids, rows, groupInfo, basis);
    }

    private static DifferentialRequest Request(DifferentialSettings settings, QuantBasis basis = QuantBasis.MsmsOnly) => new()
    {
        Table = Table(basis),
        Samples = Ids.Select(s => new AnalysisSample(s) { Factors = new Dictionary<string, string> { ["mix"] = s[..1].ToUpperInvariant() } }).ToList(),
        Spec = new DesignSpec { Factors = new[] { new FactorSpec("mix", "A") } },
        Contrasts = new[] { DifferentialContrastRequest.Between("mix", "B", "A") },
        QuantStyle = "dia",
        Settings = settings,
    };

    private static DifferentialSettings Background(params string[] extra) => new()
    {
        Normalization = "background_set:human",
        BackgroundSetGroups = Human.Concat(extra).ToArray(),
    };

    private static double MedianEffect(DifferentialAnalysisResult result, IEnumerable<string> groups)
    {
        var effects = result.Rows.Where(r => groups.Contains(r.FeatureId)).Select(r => r.Log2Effect!.Value).Order().ToArray();
        return (effects[(effects.Length - 1) / 2] + effects[effects.Length / 2]) / 2;
    }

    [Test]
    public void ADiaRowHasNoQuantBasisAndNoMbrCountAndCountsPeptidesNotPrecursors()
    {
        var result = DifferentialAnalysis.Run(Request(new DifferentialSettings { Normalization = SharedPeptideNormalization.None }), maxThreads: 1);

        var y00 = result.Rows.Single(r => r.FeatureId == "Y00");
        Assert.That(y00.Status, Is.EqualTo(DifferentialStatus.Fitted));
        Assert.That(y00.QuantStyle, Is.EqualTo("dia"));
        Assert.That(y00.QuantBasis, Is.Null, "DEF-DIFF-COLUMNS: empty for a style with no MBR");
        Assert.That(y00.NMbrValues, Is.Null);
        Assert.That(y00.ReporterAcquisition, Is.Null);
        Assert.That(y00.NPeptides, Is.EqualTo(4), "four peptides, six precursors");
        Assert.That(y00.NObservations, Is.EqualTo(36), "six precursors in six samples");
        Assert.That(result.Decisions.Select(d => d.QuantBasis), Is.All.Null);
        // The column registry names quant_basis as a column; no entry carries it as a key or a basis value.
        Assert.That(Encoding.UTF8.GetString(DifferentialMetadataWriter.ToBytes(result.Metadata, result.Rows)),
            Does.Not.Contain("\"quant_basis\": ").And.Not.Contain("msms_only").And.Not.Contain("mbr_kept"));

        // The model follows the levels, as for TMT ions: one peptide at two charges is two levels, so it is fitted with the
        // mixed model, not as a single peptide, while n_peptides is 1.
        var one = result.Rows.Single(r => r.FeatureId == "ONE");
        Assert.That(one.Status, Is.EqualTo(DifferentialStatus.Fitted));
        Assert.That(one.NPeptides, Is.EqualTo(1));
        Assert.That(one.ModelUsed, Is.Not.EqualTo(DifferentialConstants.ModeratedT));
    }

    [Test]
    public void ADiaTableBuiltWithMatchBetweenRunsValuesIsRefused()
    {
        var ex = Assert.Throws<ArgumentException>(() => DifferentialAnalysis.Run(
            Request(new DifferentialSettings { Normalization = SharedPeptideNormalization.None }, QuantBasis.MbrKept), maxThreads: 1));

        Assert.That(ex!.Message, Does.Contain("'dia'").And.Contain("mbr_kept"));
    }

    [TestCase("shared_peptide_mediann")]
    [TestCase("shared_peptide_median_half")]
    [TestCase("global_channel_median")]
    [TestCase("background_set:")]
    [TestCase("background_set:two words")]
    public void ANormalizationSettingTheRunnerDoesNotKnowIsRefused(string setting)
    {
        var request = Request(new DifferentialSettings { Normalization = setting });

        var ex = Assert.Throws<ArgumentException>(() => DifferentialAnalysis.Run(request, maxThreads: 1));

        Assert.That(ex!.Message, Does.Contain($"'{setting}'"));
    }

    [Test]
    public void ABackgroundSetNeedsItsGroupsAndGroupsNeedABackgroundSet()
    {
        foreach (var groups in new[] { null, Array.Empty<string>(), new[] { "H00", " " } })
            Assert.That(Assert.Throws<ArgumentException>(() => DifferentialAnalysis.Run(Request(new DifferentialSettings
                { Normalization = "background_set:human", BackgroundSetGroups = groups }), maxThreads: 1))!.Message,
                Does.Contain("'background_set:human'"));
        Assert.That(Assert.Throws<ArgumentException>(() => DifferentialAnalysis.Run(
            Request(new DifferentialSettings { BackgroundSetGroups = Human }), maxThreads: 1))!.Message, Does.Contain("BackgroundSetGroups"));
    }

    /// <summary>
    /// Yeast outnumbers human, so the shared-peptide median aligns the samples on yeast: human then reads about -1 and yeast
    /// about 0. Declaring the human groups as the background set recovers both: human 0, yeast +1.
    /// </summary>
    [Test]
    public void TheBackgroundSetRecoversBothSpeciesWhereTheSharedPeptideMedianDoesNot()
    {
        var shared = DifferentialAnalysis.Run(Request(new DifferentialSettings()), maxThreads: 1);
        var background = DifferentialAnalysis.Run(Request(Background("HGONE")), maxThreads: 1);

        Assert.That(MedianEffect(shared, Human), Is.LessThan(-0.5), "aligned on the changing species");
        Assert.That(MedianEffect(background, Human), Is.EqualTo(0).Within(0.05));
        Assert.That(MedianEffect(background, Yeast), Is.EqualTo(1).Within(0.05));
        Assert.That(background.Rows.Select(r => r.Normalization), Is.All.EqualTo("background_set:human"));

        var normalization = background.Metadata.Normalization.Single();
        Assert.That(normalization.Setting, Is.EqualTo("background_set:human"));
        Assert.That(normalization.QuantBasis, Is.Null);
        Assert.That(normalization.ReferenceSetSize, Is.EqualTo(120), "the human precursors; ONE is not declared");
        Assert.That(normalization.DeclaredGroups, Is.EqualTo(21));
        Assert.That(normalization.DeclaredGroupsAbsent, Is.EqualTo(new[] { "HGONE" }));
    }

    [Test]
    public void TheDeclaredSetIsRecordedAndEntersTheAnalysisId()
    {
        var result = DifferentialAnalysis.Run(Request(Background("HGONE")), maxThreads: 1);

        string json = Encoding.UTF8.GetString(DifferentialMetadataWriter.ToBytes(result.Metadata, result.Rows));
        Assert.That(json, Does.Contain("\"normalization\": \"background_set:human\"")
            .And.Contain("\"background_set_sha256\": \"").And.Contain("\"declared_groups\": 21")
            .And.Contain("\"declared_groups_absent\": [").And.Contain("\"HGONE\""));
        var reversed = Background("HGONE") with { BackgroundSetGroups = Background("HGONE").BackgroundSetGroups!.Reverse().ToArray() };
        Assert.That(DifferentialMetadata.ComputeAnalysisId(result.Metadata.Inputs, reversed), Is.EqualTo(result.Metadata.AnalysisId),
            "the same set in another order is the same analysis");
        Assert.That(DifferentialMetadata.ComputeAnalysisId(result.Metadata.Inputs, Background()), Is.Not.EqualTo(result.Metadata.AnalysisId));
        string report = StatisticalReportWriter.Build(result.Metadata, result.Rows, result.Diagnostics, result.Decisions).Markdown;
        Assert.That(report, Does.Contain("declared background set `human` (120 peptides; 21 protein groups declared, 1 of them not " +
            "in the table, listed in the metadata file)."));
    }

    [Test]
    public void SettingsWithoutABackgroundSetWriteNoBackgroundKeys()
    {
        var result = DifferentialAnalysis.Run(Request(new DifferentialSettings()), maxThreads: 1);

        string json = Encoding.UTF8.GetString(DifferentialMetadataWriter.ToBytes(result.Metadata, result.Rows));
        Assert.That(json, Does.Not.Contain("background_set").And.Not.Contain("declared_groups"));
    }
}
