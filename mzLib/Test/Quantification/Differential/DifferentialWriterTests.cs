using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The differential results table and its metadata file, against QuantProject's <c>DEF-DIFF-*</c> v1 contract
/// (DATA-DEFINITIONS v3.7, and v3.8's global shift per contrast and stratum): both header names for every column, in
/// order; missing as an empty cell; round-trip doubles; the status vocabulary and its "evidence only" rule; the fixed row
/// order (features never interleave); family sizes; the metadata's key order, analysis id and determinism; the global
/// shift's place; and byte-for-byte golden files (STAT1 milestone M4).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialWriterTests
{
    /// <summary>The 56 columns of DEF-DIFF-COLUMNS v1 (DATA-DEFINITIONS v3.7), machine and readable names, in order.</summary>
    private static readonly (string Machine, string Human)[] Contract =
    {
        ("definition_id", "Definition ID"), ("analysis_id", "Analysis ID"), ("quant_style", "Quantification Style"),
        ("reporter_acquisition", "Reporter Ion Acquisition"), ("stratum", "Stratum"), ("grain", "Grain"),
        ("feature_id", "Feature"), ("feature_accessions", "Accessions"), ("genes", "Genes"), ("quantity", "Quantity"),
        ("effect_type", "Effect Type"), ("quant_basis", "Quantification Basis"),
        ("contrast_id", "Contrast ID"), ("contrast_label", "Contrast"), ("numerator", "Numerator"),
        ("denominator", "Denominator"), ("covariate", "Covariate"), ("covariate_unit", "Covariate Unit"),
        ("covariate_scale", "Covariate Scale"),
        ("log2_effect", "Log2 Effect"), ("natural_effect", "Effect (Natural Scale)"), ("natural_unit", "Natural Unit"),
        ("delta_percentage_points", "Change (Percentage Points)"), ("mean_numerator", "Numerator Mean"),
        ("mean_denominator", "Denominator Mean"),
        ("se", "Standard Error"), ("ci_low", "CI Lower"), ("ci_high", "CI Upper"), ("ci_level", "CI Level"),
        ("statistic", "Test Statistic"), ("df", "Degrees of Freedom"), ("df_method", "DF Method"), ("p_value", "P-Value"),
        ("p_adjusted", "Adjusted P-Value"), ("adjustment_method", "Adjustment Method"), ("family_size", "Family Size"),
        ("pep", "Posterior Error Probability"), ("bayesian_fdr", "Bayesian False Discovery Rate"),
        ("bayes_factor", "Bayes Factor"), ("null_width", "Null Hypothesis Width"),
        ("n_samples_numerator", "Numerator Samples With Value"), ("n_samples_denominator", "Denominator Samples With Value"),
        ("n_samples_with_value", "Samples With Value"), ("n_samples_total", "Samples In Contrast"), ("n_peptides", "Peptides"),
        ("n_observations", "Peptide-Sample Values"), ("n_mbr_values", "MBR Values"),
        ("status", "Status"), ("status_detail", "Status Detail"), ("method", "Method"), ("model_used", "Model Used"),
        ("robust", "Robust Weighting"), ("normalization", "Normalization"), ("covariates_fitted", "Covariates Fitted"),
        ("design_sha256", "Design SHA-256"), ("method_version", "Method Version"),
    };

    private const string Analysis = "0123456789abcdef";

    private static DifferentialResult Fitted(string feature, string contrast = "c1", string? basis = "msms_only",
        string stratum = "all", double log2 = 0.5, int familySize = 1) => new()
    {
        DefinitionId = "QuantProject:DEF-DIFF-ABUNDANCE v1", AnalysisId = Analysis, QuantStyle = "lfq", Stratum = stratum,
        Grain = "protein_group", FeatureId = feature, FeatureAccessions = new[] { feature }, Quantity = "abundance",
        EffectType = "abundance_log2_ratio", QuantBasis = basis, ContrastId = contrast, ContrastLabel = "age=old vs age=young",
        Numerator = "age=old", Denominator = "age=young", Log2Effect = log2, NaturalEffect = Math.Pow(2, log2),
        NaturalUnit = "ratio", MeanNumerator = 20.5, MeanDenominator = 20.0, StandardError = 0.2, CiLow = log2 - 0.4,
        CiHigh = log2 + 0.4, CiLevel = 0.95, Statistic = log2 / 0.2, Df = 12.5, DfMethod = "satterthwaite_moderated",
        PValue = 0.03, PAdjusted = 0.06, AdjustmentMethod = "benjamini_hochberg", FamilySize = familySize,
        NSamplesNumerator = 6, NSamplesDenominator = 6, NSamplesWithValue = 12, NSamplesTotal = 12, NPeptides = 4,
        NObservations = 44, NMbrValues = 0, Status = DifferentialStatus.Fitted, Method = "moderated",
        ModelUsed = "peptide_mixed_model", Robust = true, Normalization = "shared_peptide_median",
        CovariatesFitted = "sex;individual(random);sample(random)", DesignSha256 = new string('a', 64), MethodVersion = "diff-1.0",
    };

    private static DifferentialResult EvidenceOnly(string feature, string status, string contrast = "c1") => Fitted(feature, contrast) with
    {
        Log2Effect = null, NaturalEffect = null, MeanNumerator = null, MeanDenominator = null, StandardError = null, CiLow = null,
        CiHigh = null, CiLevel = null, Statistic = null, Df = null, DfMethod = null, PValue = null, PAdjusted = null,
        AdjustmentMethod = null, FamilySize = null, Status = status, StatusDetail = "counts 0 vs 6", NSamplesNumerator = 0,
    };

    private static string[] Lines(string text) => text.Split('\n');

    private static string WriteToString(IEnumerable<DifferentialResult> rows, DifferentialHeaderStyle style)
    {
        using var w = new StringWriter { NewLine = "\n" };
        DifferentialResultWriter.Write(w, rows, style);
        return w.ToString();
    }

    private static Dictionary<string, string> Row(string text, int index)
    {
        var lines = Lines(text);
        var header = lines[0].Split('\t');
        var cells = lines[1 + index].Split('\t');
        return header.Zip(cells).ToDictionary(z => z.First, z => z.Second);
    }

    [Test]
    public void EveryColumnHasBothNamesInContractOrder()
    {
        var all = DifferentialColumns.All;
        Assert.That(all.Select(c => (c.MachineName, c.HumanName)), Is.EqualTo(Contract));
        Assert.That(all.Select(c => c.MachineName).Distinct().Count(), Is.EqualTo(56));
        Assert.That(all.Select(c => c.HumanName).Distinct().Count(), Is.EqualTo(56));
        Assert.That(all.Select(c => c.Group).Distinct(),
            Is.EqualTo(new[] { "Identity", "Contrast", "Effect", "Confidence", "Evidence", "Status" }));
    }

    [TestCase(DifferentialHeaderStyle.Machine)]
    [TestCase(DifferentialHeaderStyle.Human)]
    public void TheHeaderUsesTheChosenStyleAndEveryRowIsAsWide(DifferentialHeaderStyle style)
    {
        var text = WriteToString(new[] { Fitted("P1"), EvidenceOnly("P2", DifferentialStatus.AbsentInDenominator) }, style);
        var lines = Lines(text);
        var expected = Contract.Select(c => style == DifferentialHeaderStyle.Machine ? c.Machine : c.Human);
        Assert.That(lines[0].Split('\t'), Is.EqualTo(expected));
        Assert.That(lines.Length, Is.EqualTo(4), "header, two rows, and the empty string after the final newline");
        Assert.That(lines[3], Is.Empty, "the file ends with a newline");
        Assert.That(lines.Take(3).All(l => l.Split('\t').Length == 56));
        Assert.That(text, Does.Not.Contain("\r"));
    }

    [Test]
    public void MissingIsAnEmptyCellNeverZeroOrNaN()
    {
        var row = Fitted("P1") with { DeltaPercentagePoints = null, Pep = null, Statistic = double.NaN, BayesFactor = double.PositiveInfinity };
        var cells = Row(WriteToString(new[] { row }, DifferentialHeaderStyle.Machine), 0);
        Assert.That(cells["delta_percentage_points"], Is.Empty);
        Assert.That(cells["pep"], Is.Empty);
        Assert.That(cells["statistic"], Is.Empty, "NaN is written empty");
        Assert.That(cells["bayes_factor"], Is.Empty, "Infinity is written empty");
        Assert.That(cells["reporter_acquisition"], Is.Empty);
        Assert.That(cells.Values, Has.None.EqualTo("NaN"));
        Assert.That(cells["n_mbr_values"], Is.EqualTo("0"), "a real 0 is written");
    }

    [Test]
    public void DoublesRoundTrip()
    {
        double[] values = { 1.0 / 3, 1e-300, 123456.78901234567, -0.1, 2.5e17, 0.1 + 0.2 };
        foreach (double v in values)
        {
            var cells = Row(WriteToString(new[] { Fitted("P1", log2: v) }, DifferentialHeaderStyle.Machine), 0);
            double back = double.Parse(cells["log2_effect"], CultureInfo.InvariantCulture);
            Assert.That(BitConverter.DoubleToInt64Bits(back), Is.EqualTo(BitConverter.DoubleToInt64Bits(v)), cells["log2_effect"]);
        }
    }

    [Test]
    public void BooleansAndListsHaveOneSpelling()
    {
        var row = Fitted("P1") with
        {
            FeatureAccessions = new[] { "Q9", "A1", "M5" }, Genes = new string?[] { "GENEQ", null, "GENEM" }, Robust = false,
        };
        var cells = Row(WriteToString(new[] { row }, DifferentialHeaderStyle.Machine), 0);
        Assert.That(cells["feature_accessions"], Is.EqualTo("A1;M5;Q9"), "ordinal-sorted");
        Assert.That(cells["genes"], Is.EqualTo(";GENEM;GENEQ"), "in the sorted accessions' order, empty where none");
        Assert.That(cells["robust"], Is.EqualTo("false"));
    }

    [Test]
    public void EveryStatusWritesAndEvidenceOnlyRowsHaveNoNumbers()
    {
        string[] statuses =
        {
            DifferentialStatus.AbsentInNumerator, DifferentialStatus.AbsentInDenominator, DifferentialStatus.AbsentInBoth,
            DifferentialStatus.BelowSupport("min_2_per_side"), DifferentialStatus.NotEstimable("confounded"),
            DifferentialStatus.NotConverged,
        };
        var rows = new List<DifferentialResult> { Fitted("A", familySize: 2), Fitted("B", familySize: 2) with { Status = DifferentialStatus.FittedSinglePeptide, ModelUsed = "moderated_t" } };
        rows.AddRange(statuses.Select((s, i) => EvidenceOnly($"Z{i}", s)));
        var text = WriteToString(rows, DifferentialHeaderStyle.Machine);
        string[] effectAndConfidence = Contract.Skip(19).Take(21).Select(c => c.Machine).Where(n => n != "natural_unit").ToArray();
        for (int i = 2; i < rows.Count; i++)
        {
            var cells = Row(text, i);
            Assert.That(statuses, Does.Contain(cells["status"]));
            foreach (var column in effectAndConfidence) Assert.That(cells[column], Is.Empty, $"{cells["status"]}: {column}");
            Assert.That(cells["n_samples_denominator"], Is.EqualTo("6"), "evidence is still written");
        }
        Assert.That(DifferentialStatus.IsFitted(DifferentialStatus.FittedSinglePeptide));
        Assert.That(DifferentialStatus.IsFitted(DifferentialStatus.BelowSupport("x")), Is.False);
    }

    [Test]
    public void ContractViolationsAreRefused()
    {
        Assert.Throws<ArgumentException>(() => WriteToString(new[] { Fitted("P1") with { Status = "maybe" } }, DifferentialHeaderStyle.Machine), "unknown status");
        Assert.Throws<ArgumentException>(() => WriteToString(new[] { EvidenceOnly("P1", DifferentialStatus.AbsentInBoth) with { PValue = 1 } }, DifferentialHeaderStyle.Machine), "a number on an evidence-only row");
        Assert.Throws<ArgumentException>(() => WriteToString(new[] { Fitted("P1") with { StandardError = null } }, DifferentialHeaderStyle.Machine), "a fitted row without SE");
        Assert.Throws<ArgumentException>(() => WriteToString(new[] { Fitted("P1"), Fitted("P1") }, DifferentialHeaderStyle.Machine), "duplicate row key");
        Assert.Throws<ArgumentException>(() => WriteToString(new[] { Fitted("P1", familySize: 3) }, DifferentialHeaderStyle.Machine), "family size disagrees with the rows");
        Assert.Throws<ArgumentException>(() => DifferentialStatus.BelowSupport(""));
    }

    [Test]
    public void TabsAndNewlinesInTextBecomeSpacesAndAreCounted()
    {
        using var w = new StringWriter { NewLine = "\n" };
        int cleaned = DifferentialResultWriter.Write(w, new[] { Fitted("P1") with { StatusDetail = "a\tb\nc" } }, DifferentialHeaderStyle.Machine);
        Assert.That(cleaned, Is.EqualTo(1));
        Assert.That(Row(w.ToString(), 0)["status_detail"], Is.EqualTo("a b c"));
    }

    [Test]
    public void RowsAreInContractOrderAndFeaturesNeverInterleave()
    {
        var rows = new List<DifferentialResult>();
        foreach (var stratum in new[] { "organism part=muscle", "organism part=liver" })
            foreach (var basis in new[] { "mbr_kept", "msms_only" })
                foreach (var contrast in new[] { "c3", "c1", "c2" })
                    foreach (var feature in new[] { "P2", "P1" })
                        rows.Add(Fitted(feature, contrast, basis, stratum, familySize: 2));
        var rng = new Random(20261008);
        var shuffled = rows.OrderBy(_ => rng.Next()).ToList();
        var ordered = DifferentialResultWriter.Order(shuffled);
        var keys = ordered.Select(r => (r.Stratum, r.QuantBasis, r.ContrastId, r.FeatureId)).ToList();
        var expected = keys.OrderBy(k => k.Stratum, StringComparer.Ordinal).ThenBy(k => k.QuantBasis, StringComparer.Ordinal)
            .ThenBy(k => k.ContrastId, StringComparer.Ordinal).ThenBy(k => k.FeatureId, StringComparer.Ordinal).ToList();
        Assert.That(keys, Is.EqualTo(expected));
        Assert.That(keys.First(), Is.EqualTo(("organism part=liver", "mbr_kept", "c1", "P1")));
        // Within one (stratum, basis, contrast) block both features appear once and adjacent.
        foreach (var block in keys.Chunk(2)) Assert.That(block.Select(k => k.FeatureId), Is.EqualTo(new[] { "P1", "P2" }));
        Assert.That(WriteToString(shuffled, DifferentialHeaderStyle.Machine), Is.EqualTo(WriteToString(rows, DifferentialHeaderStyle.Machine)),
            "input order does not change the file");
    }

    private static DifferentialMetadata Metadata() => new()
    {
        Software = new[] { new DifferentialSoftware("mzLib", "9.9.0"), new DifferentialSoftware("datarepo engine", "0.1.0") },
        Inputs = new[]
        {
            new DifferentialInput("observation_table", "observations.tsv", new string('1', 64)),
            new DifferentialInput("design", "design.tsv", new string('2', 64)),
        },
        Settings = new DifferentialSettings { Seed = 42 },
        Strata = new[]
        {
            new DifferentialStratumInfo("all", new Dictionary<string, string>(), "default_list",
                new Dictionary<string, int> { ["age=young"] = 6, ["age=old"] = 6 }, Array.Empty<DifferentialNotRun>())
            {
                IndividualsPerLevel = new Dictionary<string, int> { ["age=young"] = 6, ["age=old"] = 5 },
                Batches = new[] { "b2", "b1" },
                Fractions = 1,
                TechnicalReplicates = 1,
                ReferenceLevels = new Dictionary<string, string> { ["age"] = "young" },
                SamplesWithoutValues = new[] { "s12" },
            },
            new DifferentialStratumInfo("organism part=liver", new Dictionary<string, string> { ["organism part"] = "liver" },
                "curated_marking", new Dictionary<string, int> { ["age=young"] = 3, ["age=old"] = 3 }, Array.Empty<DifferentialNotRun>())
            {
                IndividualsPerLevel = new Dictionary<string, int> { ["age=young"] = 3, ["age=old"] = 3 },
                Batches = new[] { "b1" },
                Fractions = 12,
                TechnicalReplicates = 2,
                SamplesWithoutValues = Array.Empty<string>(),
            },
        },
        Contrasts = new[]
        {
            new DifferentialContrastInfo("c1", "age=old vs age=young", "age=old", "age=young", null, null, null,
                new Dictionary<string, double> { ["age=old"] = 1 },
                new[]
                {
                    new DifferentialGlobalShift("organism part=liver", "mbr_kept", -0.11, 512),
                    new DifferentialGlobalShift("all", "msms_only", 0.03, 588),
                    new DifferentialGlobalShift("all", "mbr_kept", 0.02, 640),
                    new DifferentialGlobalShift("organism part=liver", "msms_only", -0.09, 455),
                }),
            new DifferentialContrastInfo("c2", "age (per decade)", null, null, "age", "decade", "years/10, centred at 50 years",
                new Dictionary<string, double> { ["age"] = 1 }),
        },
        Models = new[]
        {
            new DifferentialModelInfo("all", "mbr_kept", "moderated", "peptide_mixed_model", "protein_group",
                "y ~ peptide + age + (1|individual/sample)", true, "satterthwaite_moderated", "huber 1.345, MAD about 0, to convergence",
                new DifferentialPriorInfo(DifferentialPriorInfo.MomentsLegacy, 4.2, false, 0.08, 1812)),
            new DifferentialModelInfo("all", "mbr_kept", "moderated", "moderated_t", "protein_group", "y ~ age", false,
                "moderated_residual", null, new DifferentialPriorInfo(DifferentialPriorInfo.MomentsLegacy, 3.1, true, null, 240)),
        },
        Normalization = new[]
        {
            new DifferentialNormalizationInfo("all", "msms_only", "shared_peptide_median", 790, 100,
                new Dictionary<string, double> { ["s01"] = 0.11, ["s02"] = -0.04 }, Array.Empty<string>()),
            new DifferentialNormalizationInfo("all", "mbr_kept", "shared_peptide_median", 812, 100,
                new Dictionary<string, double> { ["s01"] = 0.1, ["s02"] = -0.05 }, new[] { "sample s12 has no value" }),
        },
        DesignWarnings = new[] { "biological replicate 3 of condition old is absent" },
    };

    [Test]
    public void AnalysisIdHashesTheInputsAndSettingsOnly()
    {
        var m = Metadata();
        Assert.That(m.AnalysisId, Has.Length.EqualTo(16).And.Match("^[0-9a-f]{16}$"));
        Assert.That(m.AnalysisId, Is.EqualTo(DifferentialMetadata.ComputeAnalysisId(m.Inputs, m.Settings)));
        Assert.That((m with { DesignWarnings = Array.Empty<string>() }).AnalysisId, Is.EqualTo(m.AnalysisId), "only inputs and settings count");
        Assert.That((m with { Settings = m.Settings with { CiLevel = 0.9 } }).AnalysisId, Is.Not.EqualTo(m.AnalysisId));
        Assert.That((m with { Inputs = m.Inputs.Reverse().ToArray() }).AnalysisId, Is.Not.EqualTo(m.AnalysisId), "inputs keep their order");
    }

    [Test]
    public void MetadataHasItsKeysInOrderAndCountsTheRows()
    {
        var m = Metadata();
        var rows = new[] { Fitted("P1") with { AnalysisId = m.AnalysisId }, EvidenceOnly("P2", DifferentialStatus.AbsentInBoth) with { AnalysisId = m.AnalysisId } };
        var bytes = DifferentialMetadataWriter.ToBytes(m, rows);
        Assert.That(bytes.Take(3), Is.Not.EqualTo(new byte[] { 0xEF, 0xBB, 0xBF }), "no byte-order mark");
        var text = Encoding.UTF8.GetString(bytes);
        Assert.That(text, Does.Not.Contain("\r"));
        Assert.That(text, Does.EndWith("\n"));
        using var doc = JsonDocument.Parse(bytes);
        Assert.That(doc.RootElement.EnumerateObject().Select(p => p.Name), Is.EqualTo(new[]
        {
            "definition_version", "analysis_id", "software", "inputs", "settings", "columns", "strata", "contrasts", "models",
            "normalization", "families", "design_warnings", "row_count", "cleaned_text_cells",
        }));
        Assert.That(doc.RootElement.GetProperty("cleaned_text_cells").GetInt32(), Is.EqualTo(0));
        var root = doc.RootElement;
        Assert.That(root.GetProperty("definition_version").GetString(), Is.EqualTo("DEF-DIFF v1"));
        Assert.That(root.GetProperty("analysis_id").GetString(), Is.EqualTo(m.AnalysisId));
        var columns = root.GetProperty("columns").EnumerateArray().ToList();
        Assert.That(columns, Has.Count.EqualTo(56));
        Assert.That(columns[0].GetProperty("machine").GetString(), Is.EqualTo("definition_id"));
        Assert.That(columns[0].GetProperty("human").GetString(), Is.EqualTo("Definition ID"));
        Assert.That(root.GetProperty("row_count").GetProperty("total").GetInt32(), Is.EqualTo(2));
        Assert.That(root.GetProperty("row_count").GetProperty("by_status").GetProperty("absent_in_both").GetInt32(), Is.EqualTo(1));
        var family = root.GetProperty("families").EnumerateArray().Single();
        Assert.That(family.GetProperty("size").GetInt32(), Is.EqualTo(1));
        Assert.That(text, Does.Not.Match(@"\d{4}-\d{2}-\d{2}T"), "no timestamps");
        Assert.That(DifferentialMetadataWriter.ToBytes(m, rows), Is.EqualTo(bytes), "the same inputs give the same bytes");
    }

    private static JsonElement Root(DifferentialMetadata m, out JsonDocument doc)
    {
        doc = JsonDocument.Parse(DifferentialMetadataWriter.ToBytes(m, new[] { Fitted("P1") with { AnalysisId = m.AnalysisId } }));
        return doc.RootElement;
    }

    private static void AssertRefused(DifferentialMetadata m, string messagePart) =>
        Assert.That(() => DifferentialMetadataWriter.ToBytes(m, new[] { Fitted("P1") with { AnalysisId = m.AnalysisId } }),
            Throws.ArgumentException.With.Message.Contains(messagePart));

    [Test]
    public void TheGlobalShiftIsPerContrastStratumAndBasis()
    {
        var root = Root(Metadata(), out var doc);
        using (doc)
        {
            var contrasts = root.GetProperty("contrasts").EnumerateArray().ToList();
            var shifts = contrasts[0].GetProperty("global_shift").EnumerateArray().ToList();
            Assert.That(shifts.Select(s => (s.GetProperty("stratum").GetString(), s.GetProperty("quant_basis").GetString())), Is.EqualTo(new[]
            {
                ("all", "mbr_kept"), ("all", "msms_only"), ("organism part=liver", "mbr_kept"), ("organism part=liver", "msms_only"),
            }), "one entry per stratum and basis, in ordinal order");
            Assert.That(shifts[0].EnumerateObject().Select(p => p.Name), Is.EqualTo(new[] { "stratum", "quant_basis", "value", "peptides" }));
            Assert.That(shifts[1].GetProperty("value").GetDouble(), Is.EqualTo(0.03), "each basis has its own value");
            Assert.That(shifts[1].GetProperty("peptides").GetInt32(), Is.EqualTo(588));
            Assert.That(contrasts[1].TryGetProperty("global_shift", out _), Is.False, "a slope has no two sides, so no global shift");
            Assert.That(root.GetProperty("normalization").EnumerateArray().Any(n => n.TryGetProperty("global_shift", out _)), Is.False,
                "the global shift belongs to a contrast, not to a stratum's normalization");
        }
    }

    [Test]
    public void AGlobalShiftWithoutABasisLeavesTheKeyOutAndWithNoPeptidesIsNull()
    {
        var m = Metadata();
        var tmt = m.Contrasts[0] with { GlobalShift = new[] { new DifferentialGlobalShift("all", null, double.NaN, 0) } };
        var root = Root(m with { Contrasts = new[] { tmt, m.Contrasts[1] } }, out var doc);
        using (doc)
        {
            var shift = root.GetProperty("contrasts")[0].GetProperty("global_shift").EnumerateArray().Single();
            Assert.That(shift.TryGetProperty("quant_basis", out _), Is.False, "a style with no basis leaves the key out");
            Assert.That(shift.GetProperty("value").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(shift.GetProperty("peptides").GetInt32(), Is.EqualTo(0));
        }
    }

    [Test]
    public void AGlobalShiftMustNameAListedStratumATwoSidedContrastAndAppearOnce()
    {
        var m = Metadata();
        var unlisted = m.Contrasts[0] with { GlobalShift = new[] { new DifferentialGlobalShift("organism part=brain", "mbr_kept", 0.1, 10) } };
        AssertRefused(m with { Contrasts = new[] { unlisted, m.Contrasts[1] } }, "'organism part=brain'");
        var slope = m.Contrasts[1] with { GlobalShift = new[] { new DifferentialGlobalShift("all", "mbr_kept", 0.1, 10) } };
        AssertRefused(m with { Contrasts = new[] { m.Contrasts[0], slope } }, "'c2'");
        var twice = m.Contrasts[0] with
        {
            GlobalShift = new[] { new DifferentialGlobalShift("all", "mbr_kept", 0.1, 10), new DifferentialGlobalShift("all", "mbr_kept", 0.2, 11) },
        };
        AssertRefused(m with { Contrasts = new[] { twice, m.Contrasts[1] } }, "twice");
    }

    [Test]
    public void StrataCarryTheDesignFactsTheReportSummarizes()
    {
        var root = Root(Metadata(), out var doc);
        using (doc)
        {
            var strata = root.GetProperty("strata").EnumerateArray().ToList();
            Assert.That(strata[0].EnumerateObject().Select(p => p.Name), Is.EqualTo(new[]
            {
                "stratum", "factors", "source", "samples_per_level", "individuals_per_level", "batches", "fractions",
                "technical_replicates", "reference_levels", "samples_without_values", "not_run",
            }));
            Assert.That(strata[0].GetProperty("individuals_per_level").GetProperty("age=old").GetInt32(), Is.EqualTo(5));
            Assert.That(strata[0].GetProperty("batches").EnumerateArray().Select(b => b.GetString()), Is.EqualTo(new[] { "b1", "b2" }),
                "ordinal order");
            Assert.That(strata[0].GetProperty("reference_levels").GetProperty("age").GetString(), Is.EqualTo("young"));
            Assert.That(strata[0].GetProperty("samples_without_values").EnumerateArray().Single().GetString(), Is.EqualTo("s12"));
            Assert.That(strata[1].GetProperty("fractions").GetInt32(), Is.EqualTo(12));
            Assert.That(strata[1].GetProperty("technical_replicates").GetInt32(), Is.EqualTo(2));
            Assert.That(strata[1].TryGetProperty("reference_levels", out _), Is.False, "no curated reference: the key is left out");
        }
        var bare = Metadata() with
        {
            Strata = new[]
            {
                new DifferentialStratumInfo("all", new Dictionary<string, string>(), "default_list", new Dictionary<string, int>(),
                    Array.Empty<DifferentialNotRun>()),
                Metadata().Strata[1],
            },
        };
        var bareRoot = Root(bare, out var bareDoc);
        using (bareDoc)
            Assert.That(bareRoot.GetProperty("strata")[0].EnumerateObject().Select(p => p.Name),
                Is.EqualTo(new[] { "stratum", "factors", "source", "samples_per_level", "not_run" }), "facts not given are left out");
    }

    [Test]
    public void ModelsArePerStratumBasisMethodAndModelWithTheirPrior()
    {
        var root = Root(Metadata(), out var doc);
        using (doc)
        {
            var models = root.GetProperty("models").EnumerateArray().ToList();
            Assert.That(models.Select(x => x.GetProperty("model_used").GetString()), Is.EqualTo(new[] { "moderated_t", "peptide_mixed_model" }),
                "ordinal order on (stratum, basis, method, model used, grain)");
            Assert.That(models[1].EnumerateObject().Select(p => p.Name), Is.EqualTo(new[]
            {
                "stratum", "quant_basis", "method", "model_used", "grain", "formula", "reml", "df_method", "robust_rule", "prior",
            }));
            var constant = models[1].GetProperty("prior");
            Assert.That(constant.EnumerateObject().Select(p => p.Name), Is.EqualTo(new[] { "estimator", "df", "trended", "variance", "features" }));
            Assert.That(constant.GetProperty("estimator").GetString(), Is.EqualTo("moments_legacy"));
            Assert.That(constant.GetProperty("variance").GetDouble(), Is.EqualTo(0.08));
            Assert.That(constant.GetProperty("features").GetInt32(), Is.EqualTo(1812));
            var trended = models[0].GetProperty("prior");
            Assert.That(trended.GetProperty("trended").GetBoolean(), Is.True);
            Assert.That(trended.TryGetProperty("variance", out _), Is.False, "a trended prior has a variance per feature, not one number");
            Assert.That(models[0].TryGetProperty("robust_rule", out _), Is.False);
        }
    }

    [Test]
    public void AModelMustBeListedOnceWithAConsistentPrior()
    {
        var m = Metadata();
        var model = m.Models[0];
        AssertRefused(m with { Models = new[] { model, model } }, "twice");
        AssertRefused(m with { Models = new[] { model with { Stratum = "organism part=brain" } } }, "'organism part=brain'");
        AssertRefused(m with { Models = new[] { model with { Prior = model.Prior! with { Trended = true } } } }, "trended");
        AssertRefused(m with { Models = new[] { model with { Prior = model.Prior! with { Variance = null } } } }, "variance");
        AssertRefused(m with { Models = new[] { model with { Prior = model.Prior! with { Estimator = "lowess" } } } }, "'lowess'");
    }

    [Test]
    public void NormalizationIsPerStratumAndBasis()
    {
        var root = Root(Metadata(), out var doc);
        using (doc)
        {
            var entries = root.GetProperty("normalization").EnumerateArray().ToList();
            Assert.That(entries.Select(e => e.GetProperty("quant_basis").GetString()), Is.EqualTo(new[] { "mbr_kept", "msms_only" }),
                "ordinal order on (stratum, basis)");
            Assert.That(entries[0].EnumerateObject().Select(p => p.Name), Is.EqualTo(new[]
            {
                "stratum", "quant_basis", "setting", "reference_set_size", "minimum_reference_set_size", "per_sample_shift", "warnings",
            }));
            Assert.That(entries[0].GetProperty("reference_set_size").GetInt32(), Is.EqualTo(812));
            Assert.That(entries[1].GetProperty("reference_set_size").GetInt32(), Is.EqualTo(790), "each basis has its own reference set");
            Assert.That(entries[0].GetProperty("minimum_reference_set_size").GetInt32(), Is.EqualTo(100));
            Assert.That(entries[0].GetProperty("warnings").EnumerateArray().Single().GetString(), Is.EqualTo("sample s12 has no value"));
        }
        var m = Metadata();
        AssertRefused(m with { Normalization = new[] { m.Normalization[0], m.Normalization[0] } }, "twice");
        AssertRefused(m with { Normalization = new[] { m.Normalization[0] with { Stratum = "organism part=brain" } } }, "'organism part=brain'");
    }

    [Test]
    public void TheMetadataCountsTheCellsTheResultsTableCleaned()
    {
        var m = Metadata();
        var rows = new[]
        {
            Fitted("P1", familySize: 3) with { AnalysisId = m.AnalysisId, StatusDetail = "a\tb" },
            Fitted("P2", familySize: 3) with { AnalysisId = m.AnalysisId, ContrastLabel = "old\nvs young" },
            Fitted("P3", familySize: 3) with { AnalysisId = m.AnalysisId },
        };
        int written;
        using (var w = new StringWriter { NewLine = "\n" })
            written = DifferentialResultWriter.Write(w, rows, DifferentialHeaderStyle.Machine);
        using var doc = JsonDocument.Parse(DifferentialMetadataWriter.ToBytes(m, rows));
        Assert.That(written, Is.EqualTo(2));
        Assert.That(doc.RootElement.GetProperty("cleaned_text_cells").GetInt32(), Is.EqualTo(written));
    }

    [Test]
    public void RowsMustBelongToTheMetadatasAnalysis()
    {
        var m = Metadata();
        Assert.Throws<ArgumentException>(() => DifferentialMetadataWriter.ToBytes(m, new[] { Fitted("P1") }));
    }

    [Test]
    public void FilesAreWrittenWithoutABomAndMatchTheGoldenFiles()
    {
        var m = Metadata();
        var rows = GoldenRows(m.AnalysisId);
        string dir = Path.Combine(Path.GetTempPath(), "mzlib-diff-writer-" + Guid.NewGuid().ToString("N"));
        try
        {
            string path = DifferentialResultWriter.Write(dir, rows, m);
            Assert.That(Path.GetFileName(path), Is.EqualTo(DifferentialResultWriter.ResultsFileName));
            byte[] tsv = File.ReadAllBytes(path);
            byte[] json = File.ReadAllBytes(Path.Combine(dir, DifferentialResultWriter.MetadataFileName));
            Assert.That(tsv.Take(3), Is.Not.EqualTo(new byte[] { 0xEF, 0xBB, 0xBF }));
            Assert.That(tsv, Is.EqualTo(Golden("DifferentialResults.golden.tsv")), "TSV byte for byte");
            Assert.That(json, Is.EqualTo(Golden("DifferentialResults.golden.metadata.json")), "metadata byte for byte");
            using var human = new StringWriter { NewLine = "\n" };
            DifferentialResultWriter.Write(human, rows, DifferentialHeaderStyle.Human);
            Assert.That(Encoding.UTF8.GetBytes(human.ToString()), Is.EqualTo(Golden("DifferentialResults.golden.human.tsv")),
                "readable-header TSV byte for byte");
        }
        finally
        {
            if (Directory.Exists(dir)) Directory.Delete(dir, true);
        }
    }

    /// <summary>One of every status, two strata, both bases, and a slope, written once and frozen as the golden files.</summary>
    internal static List<DifferentialResult> GoldenRows(string analysis)
    {
        var rows = new List<DifferentialResult>();
        foreach (var stratum in new[] { "all", "organism part=liver" })
            foreach (var basis in new[] { "mbr_kept", "msms_only" })
            {
                rows.Add(Fitted("PG1", "c1", basis, stratum, log2: 1.0 / 3, familySize: 2) with { AnalysisId = analysis });
                rows.Add(Fitted("PG2", "c1", basis, stratum, log2: -0.75, familySize: 2) with
                {
                    AnalysisId = analysis, Status = DifferentialStatus.FittedSinglePeptide, ModelUsed = "moderated_t",
                    DfMethod = "moderated_residual", NPeptides = 1,
                });
                rows.Add(EvidenceOnly("PG3", DifferentialStatus.AbsentInNumerator) with { AnalysisId = analysis, QuantBasis = basis, Stratum = stratum });
                rows.Add(EvidenceOnly("PG4", DifferentialStatus.BelowSupport("min_2_per_side")) with { AnalysisId = analysis, QuantBasis = basis, Stratum = stratum });
                rows.Add(EvidenceOnly("PG5", DifferentialStatus.NotEstimable("confounded")) with { AnalysisId = analysis, QuantBasis = basis, Stratum = stratum });
                rows.Add(EvidenceOnly("PG6", DifferentialStatus.NotConverged) with { AnalysisId = analysis, QuantBasis = basis, Stratum = stratum });
                rows.Add(EvidenceOnly("PG7", DifferentialStatus.AbsentInDenominator) with { AnalysisId = analysis, QuantBasis = basis, Stratum = stratum });
                rows.Add(EvidenceOnly("PG8", DifferentialStatus.AbsentInBoth) with { AnalysisId = analysis, QuantBasis = basis, Stratum = stratum });
                rows.Add(Fitted("PG1", "c2", basis, stratum, log2: 0.1, familySize: 1) with
                {
                    AnalysisId = analysis, EffectType = "abundance_log2_slope", ContrastLabel = "age (per decade)", Numerator = null,
                    Denominator = null, Covariate = "age", CovariateUnit = "decade", CovariateScale = "years/10, centred at 50 years",
                    NaturalUnit = "ratio_per_decade", MeanNumerator = null, MeanDenominator = null, NSamplesNumerator = null,
                    NSamplesDenominator = null,
                });
            }
        return rows;
    }

    private static byte[] Golden(string file) =>
        File.ReadAllBytes(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential", "GoldenFiles", file));
}
