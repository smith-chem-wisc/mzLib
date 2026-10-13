using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The generated statistics report (STAT1 milestone M6, QuantProject design/STATS-FRAMEWORK.md section 6): its five sections
/// in order; a methods paragraph, design summary, results summary and journal checklist filled from the run, with a gap
/// shown as N; the diagnostic figures' numbers (left-closed p-value bins, pairwise-complete Pearson correlation, limma's
/// plotSA scale); the same bytes on every machine and culture; golden files; and the refusals.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class StatisticalReportTests
{
    private static readonly string[] Bases = { "mbr_kept", "msms_only" };

    internal static DifferentialMetadata Metadata() => new()
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
        },
        Contrasts = new[]
        {
            new DifferentialContrastInfo("c1", "age=old vs age=young", "age=old", "age=young", null, null, null,
                new Dictionary<string, double> { ["age=old"] = 1 },
                new[]
                {
                    new DifferentialGlobalShift("all", "msms_only", 0.03, 588),
                    new DifferentialGlobalShift("all", "mbr_kept", 0.02, 640),
                }),
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
            new DifferentialNormalizationInfo("all", "mbr_kept", "shared_peptide_median", 812, 100,
                new Dictionary<string, double> { ["s01"] = 0.1, ["s02"] = -0.05 }, Array.Empty<string>()),
            new DifferentialNormalizationInfo("all", "msms_only", "shared_peptide_median", 790, 100,
                new Dictionary<string, double> { ["s01"] = 0.11, ["s02"] = -0.04 }, Array.Empty<string>()),
        },
        DesignWarnings = new[] { "biological replicate 3 of condition old is absent" },
    };

    private static DifferentialResult Fitted(string analysis, string basis, int i)
    {
        double log2 = (i - 9.5) / 6.0;
        double p = Math.Pow((i + 0.5) / 20.0, 2);
        return new DifferentialResult
        {
            DefinitionId = "QuantProject:DEF-DIFF-ABUNDANCE v1", AnalysisId = analysis, QuantStyle = "lfq", Stratum = "all",
            Grain = "protein_group", FeatureId = $"PG{i:00}", FeatureAccessions = new[] { $"P{i:00}" }, Quantity = "abundance",
            EffectType = "abundance_log2_ratio", QuantBasis = basis, ContrastId = "c1", ContrastLabel = "age=old vs age=young",
            Numerator = "age=old", Denominator = "age=young", Log2Effect = log2, NaturalEffect = Math.Pow(2, log2),
            NaturalUnit = "ratio", MeanNumerator = 20 + log2 / 2, MeanDenominator = 20 - log2 / 2, StandardError = 0.2,
            CiLow = log2 - 0.4, CiHigh = log2 + 0.4, CiLevel = 0.95, Statistic = log2 / 0.2, Df = 12.5,
            DfMethod = "satterthwaite_moderated", PValue = p, PAdjusted = Math.Min(1, p * 20 / (i + 1)),
            AdjustmentMethod = "benjamini_hochberg", FamilySize = 20, NSamplesNumerator = 6, NSamplesDenominator = 6,
            NSamplesWithValue = 12, NSamplesTotal = 12, NPeptides = 4, NObservations = 44, NMbrValues = basis == "mbr_kept" ? 3 : 0,
            Status = DifferentialStatus.Fitted, Method = "moderated", ModelUsed = "peptide_mixed_model", Robust = true,
            Normalization = "shared_peptide_median", CovariatesFitted = "individual(random);sample(random)",
            DesignSha256 = new string('a', 64), MethodVersion = "diff-1.0",
        };
    }

    private static DifferentialResult NotFitted(string analysis, string basis, string feature, string status) => Fitted(analysis, basis, 0) with
    {
        FeatureId = feature, FeatureAccessions = new[] { feature }, Log2Effect = null, NaturalEffect = null, MeanNumerator = null,
        MeanDenominator = null, StandardError = null, CiLow = null, CiHigh = null, CiLevel = null, Statistic = null, Df = null,
        DfMethod = null, PValue = null, PAdjusted = null, AdjustmentMethod = null, FamilySize = null, NaturalUnit = null, Status = status,
        StatusDetail = "counts 0 vs 6", NSamplesNumerator = 0,
    };

    internal static List<DifferentialResult> Rows(string analysis)
    {
        var rows = new List<DifferentialResult>();
        foreach (var basis in Bases)
        {
            rows.AddRange(Enumerable.Range(0, 20).Select(i => Fitted(analysis, basis, i)));
            rows.Add(NotFitted(analysis, basis, "PGX1", DifferentialStatus.AbsentInNumerator));
            rows.Add(NotFitted(analysis, basis, "PGX2", DifferentialStatus.BelowSupport("min_2_per_side")));
        }
        return rows;
    }

    internal static DifferentialDiagnostics Diagnostics()
    {
        var ids = Enumerable.Range(0, 20).Select(i => $"PG{i:00}").ToArray();
        return new DifferentialDiagnostics
        {
            Fits = new[]
            {
                new DiagnosticFit("all", "mbr_kept", "moderated", "peptide_mixed_model", "protein_group", ids,
                    ids.Select((_, i) => 18 + 0.3 * i).ToArray(), ids.Select((_, i) => 0.05 + 0.01 * (i % 5)).ToArray(),
                    ids.Select(_ => 0.08).ToArray()),
                new DiagnosticFit("all", "mbr_kept", "moderated", "moderated_t", "protein_group", ids,
                    ids.Select((_, i) => 18 + 0.3 * i).ToArray(), ids.Select((_, i) => 0.12 - 0.004 * i).ToArray(),
                    ids.Select((_, i) => 0.1 - 0.002 * i).ToArray()),
            },
            Samples = new[]
            {
                new SampleLog2("all", "mbr_kept", new[] { "s01", "s02", "s03", "s04" }, new IReadOnlyList<double>[]
                {
                    new[] { 20.0, 21.0, 22.0, 23.0, 24.0, 25.0 },
                    new[] { 20.1, 21.2, 21.9, 23.1, double.NaN, 25.2 },
                    new[] { 25.0, 24.0, 23.0, 22.0, 21.0, 20.0 },
                    new[] { double.NaN, double.NaN, double.NaN, double.NaN, 22.0, 23.0 },
                }),
            },
        };
    }

    private static StatisticalReport Build(DifferentialMetadata? metadata = null, DifferentialDiagnostics? diagnostics = null,
        List<DifferentialResult>? rows = null)
    {
        var m = metadata ?? Metadata();
        return StatisticalReportWriter.Build(m, rows ?? Rows(m.AnalysisId), diagnostics ?? Diagnostics());
    }

    private static string Section(string markdown, string heading)
    {
        int start = markdown.IndexOf(heading + "\n", StringComparison.Ordinal);
        Assert.That(start, Is.GreaterThanOrEqualTo(0), heading);
        int end = markdown.IndexOf("\n## ", start + heading.Length, StringComparison.Ordinal);
        return end < 0 ? markdown[start..] : markdown[start..end];
    }

    private static List<(string Item, string Answer)> Checklist(string markdown) =>
        Section(markdown, "## 5. Reporting checklist").Split('\n')
            .Where(l => l.StartsWith("| ", StringComparison.Ordinal) && !l.StartsWith("| item ", StringComparison.Ordinal))
            .Select(l => l.Split('|').Select(c => c.Trim()).ToArray())
            .Select(c => (c[1], c[4])).ToList();

    private static string Figure(StatisticalReport report, string name) =>
        Encoding.UTF8.GetString(report.Figures.Single(f => f.Name == name).Bytes);

    private static double Pearson(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average();
        double sxy = x.Zip(y).Sum(p => (p.First - mx) * (p.Second - my));
        return sxy / Math.Sqrt(x.Sum(v => (v - mx) * (v - mx)) * y.Sum(v => (v - my) * (v - my)));
    }

    [Test]
    public void TheReportHasItsFiveSectionsInOrder()
    {
        var headings = Build().Markdown.Split('\n').Where(l => l.StartsWith("## ", StringComparison.Ordinal)).ToList();
        Assert.That(headings, Is.EqualTo(new[] { "## 1. Methods", "## 2. Design", "## 3. Results", "## 4. Diagnostics", "## 5. Reporting checklist" }));
    }

    [Test]
    public void TheMethodsSayWhatTheJournalsAskFor()
    {
        var methods = Section(Build().Markdown, "## 1. Methods");
        foreach (var phrase in new[]
                 {
                     "mzLib 9.9.0", "datarepo engine 0.1.0", "diff-1.0", "log2 scale", "812 peptides with a value in every sample",
                     "not imputed", "two-sided", "95% confidence interval", "Benjamini–Hochberg", "2 families",
                     "`y ~ peptide + age + (1|individual/sample)`", "fitted by REML", "`satterthwaite_moderated`",
                     "4.2 degrees of freedom and variance 0.08", "fitted on 1,812 features", "trended", "No significance threshold",
                     "Random seed: 42",
                 })
            Assert.That(methods, Does.Contain(phrase), phrase);
    }

    /// <summary>
    /// An infinite prior df is a real answer, not a missing one (mzLib #1453 review round 3): the Methods say "infinite", never
    /// "not recorded".
    /// </summary>
    [Test]
    public void AnInfinitePriorDfIsReportedAsInfinite()
    {
        var m = Metadata();
        m = m with { Models = m.Models.Select(x => x.Prior is { Trended: false } p ? x with { Prior = p with { Df = double.PositiveInfinity } } : x).ToArray() };

        var methods = Section(Build(m).Markdown, "## 1. Methods");

        Assert.That(methods, Does.Contain("with infinite degrees of freedom and variance 0.08"));
        Assert.That(methods, Does.Not.Contain("not recorded degrees of freedom"));
    }

    /// <summary>The half-set fallback says what was done: the set, why, and that the shifts were settled (mzLib #1455 review round 3).</summary>
    [Test]
    public void TheMethodsDescribeTheHalfSetFallback()
    {
        var m = Metadata();
        m = m with
        {
            Normalization = new[]
            {
                m.Normalization[0] with { Setting = "shared_peptide_median_half", ReferenceSetSize = 64 },
                m.Normalization[1],
            },
        };

        var methods = Section(Build(m).Markdown, "## 1. Methods");

        Assert.That(methods, Does.Contain("over the 64 peptides with a value in at least half the samples that have values, " +
            "because fewer than 100 had a value in every one; the differences were measured again on the shifted values until " +
            "the shifts stopped changing, and the shifts were centred to average 0."));
        Assert.That(methods, Does.Contain("over the 790 peptides with a value in every sample that has values"));
    }

    [Test]
    public void TheDesignGivesNPerGroupReplicatesAndEveryComparison()
    {
        var design = Section(Build().Markdown, "## 2. Design");
        foreach (var phrase in new[]
                 {
                     "| age=old | 6 | 5 |", "| age=young | 6 | 6 |", "Batches: b1, b2.", "Fractions: 1. Technical replicates: 1.",
                     "Reference level: age = young.", "Samples without values: s12.", "biological replicate 3 of condition old is absent",
                     "| c1 | age=old vs age=young | age=old | age=young |",
                 })
            Assert.That(design, Does.Contain(phrase), phrase);
    }

    [Test]
    public void TheResultsCountEachFamilysFeaturesByStatus()
    {
        var m = Metadata();
        var results = Section(Build().Markdown, "## 3. Results");
        foreach (var basis in Bases)
        {
            int below = Rows(m.AnalysisId).Count(r => r.QuantBasis == basis && r.PAdjusted < 0.05);
            Assert.That(results, Does.Contain($"| all | {basis} | protein_group | abundance | c1 | moderated | 22 | 20 | 2 | {below} |"), basis);
            Assert.That(results, Does.Contain($"| all | {basis} | c1 | moderated | absent_in_numerator | 1 |"), basis);
            Assert.That(results, Does.Contain($"| all | {basis} | c1 | moderated | below_support:min_2_per_side | 1 |"), basis);
        }
        Assert.That(results, Does.Contain("descriptive only"));
    }

    [Test]
    public void TheChecklistIsFilledFromTheRunAndShowsAGapAsN()
    {
        var checklist = Checklist(Build().Markdown);
        Assert.That(checklist, Has.Count.EqualTo(15));
        Assert.That(checklist.Where(c => c.Answer == "N").Select(c => c.Item), Is.EqualTo(new[] { "Power analysis" }));
        Assert.That(checklist.All(c => c.Answer is "Y" or "N"));

        var m = Metadata();
        string Answer(DifferentialMetadata changed, string item) => Checklist(Build(changed).Markdown).Single(c => c.Item == item).Answer;
        Assert.That(Answer(m with { Strata = new[] { m.Strata[0] with { IndividualsPerLevel = null } } }, "Biological and technical replicates"),
            Is.EqualTo("N"));
        Assert.That(Answer(m with { Models = m.Models.Select(x => x with { Prior = null }).ToArray() }, "Priors"), Is.EqualTo("N"));
        Assert.That(Answer(m with { Software = Array.Empty<DifferentialSoftware>() }, "Software and versions"), Is.EqualTo("N"));
        Assert.That(Answer(m with { Models = Array.Empty<DifferentialModelInfo>() }, "Statistical tests and sidedness"), Is.EqualTo("N"));
    }

    [Test]
    public void PValueBinsAreLeftClosedAndOneFallsInTheLastBin()
    {
        var bins = DiagnosticFigures.PValueBins(new[] { 0.0, 0.0499, 0.05, 0.15, 0.999, 1.0 });
        Assert.That(bins, Has.Length.EqualTo(20));
        Assert.That(bins[0], Is.EqualTo(2), "0 and 0.0499");
        Assert.That(bins[1], Is.EqualTo(1), "0.05 opens the second bin");
        Assert.That(bins[3], Is.EqualTo(1), "0.15 opens the fourth bin");
        Assert.That(bins[19], Is.EqualTo(2), "0.999 and 1");
        Assert.That(bins.Sum(), Is.EqualTo(6));
    }

    [Test]
    public void SampleCorrelationUsesOnlyFeaturesBothSamplesHave()
    {
        var pairs = DiagnosticFigures.Correlations(Diagnostics().Samples[0]);
        Assert.That(pairs.Select(p => (p.A, p.B)), Is.EqualTo(new[]
        {
            ("s01", "s02"), ("s01", "s03"), ("s01", "s04"), ("s02", "s03"), ("s02", "s04"), ("s03", "s04"),
        }));
        var s01s02 = pairs.Single(p => p.A == "s01" && p.B == "s02");
        Assert.That(s01s02.N, Is.EqualTo(5), "the feature s02 lacks is left out");
        Assert.That(s01s02.R, Is.EqualTo(Pearson(new[] { 20.0, 21, 22, 23, 25 }, new[] { 20.1, 21.2, 21.9, 23.1, 25.2 })).Within(1e-12));
        Assert.That(pairs.Single(p => p.A == "s01" && p.B == "s03").R, Is.EqualTo(-1).Within(1e-12));
        var s01s04 = pairs.Single(p => p.A == "s01" && p.B == "s04");
        Assert.That(s01s04.N, Is.EqualTo(2));
        Assert.That(double.IsNaN(s01s04.R), "fewer than 3 shared features leaves r empty");
    }

    [Test]
    public void MeanVarianceFollowsLimmasPlotSA()
    {
        var fit = Diagnostics().Fits[1];
        var points = DiagnosticFigures.MeanVariancePoints(fit);
        Assert.That(points, Has.Count.EqualTo(20));
        Assert.That(points[3].X, Is.EqualTo(fit.AverageLog2[3]));
        Assert.That(points[3].Y, Is.EqualTo(Math.Pow(fit.ResidualVariance[3], 0.25)).Within(1e-15), "sqrt(sigma) = s2^(1/4)");
        Assert.That(points[3].PriorY, Is.EqualTo(Math.Pow(fit.PriorVariance[3], 0.25)).Within(1e-15));
    }

    [Test]
    public void AVolcanoPointWithPZeroIsDrawnAtTheTopAndItsCellIsEmpty()
    {
        var m = Metadata();
        var rows = Rows(m.AnalysisId);
        rows[0] = rows[0] with { PValue = 0 };
        var report = Build(rows: rows);
        var table = Figure(report, "figure_02_volcano.tsv").Split('\n');
        Assert.That(table[0], Is.EqualTo("feature_id\tlog2_effect\tp_value\tneg_log10_p"));
        Assert.That(table.Single(l => l.StartsWith("PG00\t", StringComparison.Ordinal)), Does.EndWith("\t0\t"));
        Assert.That(Figure(report, "figure_02_volcano.svg"), Does.Contain("<polygon"), "a triangle at the top edge");
        Assert.That(Section(report.Markdown, "## 4. Diagnostics"), Does.Contain("1 feature with p = 0"));
    }

    [Test]
    public void TheReportAndFiguresAreTheSameOnEveryMachine()
    {
        var invariant = Build();
        var saved = CultureInfo.CurrentCulture;
        StatisticalReport german;
        try
        {
            CultureInfo.CurrentCulture = new CultureInfo("de-DE");
            german = Build();
        }
        finally
        {
            CultureInfo.CurrentCulture = saved;
        }
        Assert.That(german.Markdown, Is.EqualTo(invariant.Markdown));
        Assert.That(german.Figures.Select(f => f.Name), Is.EqualTo(invariant.Figures.Select(f => f.Name)));
        foreach (var (g, i) in german.Figures.Zip(invariant.Figures))
            Assert.That(g.Bytes, Is.EqualTo(i.Bytes), g.Name);
        Assert.That(invariant.Figures.Select(f => Encoding.UTF8.GetString(f.Bytes)), Has.None.Contain("\r"));
        Assert.That(invariant.Markdown, Does.Not.Contain("\r"));
        Assert.That(invariant.Markdown, Does.Not.Match(@"\d{4}-\d{2}-\d{2}T"), "no timestamps");
    }

    [Test]
    public void TheFiguresAreNumberedInReportOrder()
    {
        Assert.That(Build().Figures.Select(f => f.Name), Is.EqualTo(new[]
        {
            "figure_01_pvalue_histogram.svg", "figure_01_pvalue_histogram.tsv", "figure_02_volcano.svg", "figure_02_volcano.tsv",
            "figure_03_pvalue_histogram.svg", "figure_03_pvalue_histogram.tsv", "figure_04_volcano.svg", "figure_04_volcano.tsv",
            "figure_05_mean_variance.svg", "figure_05_mean_variance.tsv", "figure_06_mean_variance.svg", "figure_06_mean_variance.tsv",
            "figure_07_sample_correlation.svg", "figure_07_sample_correlation.tsv", "figure_08_global_shift.svg", "figure_08_global_shift.tsv",
        }));
    }

    [Test]
    public void InconsistentInputsAreRefused()
    {
        var m = Metadata();
        var d = Diagnostics();
        void Refused(DifferentialDiagnostics diagnostics, string part) =>
            Assert.That(() => StatisticalReportWriter.Build(m, Rows(m.AnalysisId), diagnostics), Throws.ArgumentException.With.Message.Contains(part));
        Refused(d with { Fits = new[] { d.Fits[0] with { ResidualVariance = d.Fits[0].ResidualVariance.Take(3).ToArray() } } }, "'peptide_mixed_model'");
        Refused(d with { Fits = new[] { d.Fits[0] with { Stratum = "organism part=brain" } } }, "'organism part=brain'");
        Refused(d with { Fits = new[] { d.Fits[0], d.Fits[0] } }, "twice");
        Refused(d with { Samples = new[] { d.Samples[0] with { SampleIds = new[] { "s01" } } } }, "samples");
        Refused(d with { Samples = new[] { d.Samples[0], d.Samples[0] } }, "twice");
        Assert.That(() => StatisticalReportWriter.Build(m, Rows("ffffffffffffffff"), d), Throws.ArgumentException, "rows of another analysis");
    }

    [Test]
    public void FilesAreWrittenWithoutABomAndMatchTheGoldenFiles()
    {
        var m = Metadata();
        string dir = Path.Combine(Path.GetTempPath(), "mzlib-stat-report-" + Guid.NewGuid().ToString("N"));
        try
        {
            string path = StatisticalReportWriter.Write(dir, m, Rows(m.AnalysisId), Diagnostics());
            Assert.That(Path.GetFileName(path), Is.EqualTo(StatisticalReportWriter.ReportFileName));
            byte[] report = File.ReadAllBytes(path);
            Assert.That(report.Take(3), Is.Not.EqualTo(new byte[] { 0xEF, 0xBB, 0xBF }));
            Assert.That(report, Is.EqualTo(Golden("StatisticalReport.golden.md")), "report byte for byte");
            string figures = Path.Combine(dir, StatisticalReportWriter.FiguresDirectoryName);
            Assert.That(Directory.GetFiles(figures), Has.Length.EqualTo(16));
            Assert.That(File.ReadAllBytes(Path.Combine(figures, "figure_01_pvalue_histogram.svg")),
                Is.EqualTo(Golden("StatisticalReport.figure_01_pvalue_histogram.golden.svg")), "histogram SVG byte for byte");
            Assert.That(File.ReadAllBytes(Path.Combine(figures, "figure_07_sample_correlation.tsv")),
                Is.EqualTo(Golden("StatisticalReport.figure_07_sample_correlation.golden.tsv")), "correlation table byte for byte");
        }
        finally
        {
            if (Directory.Exists(dir)) Directory.Delete(dir, true);
        }
    }

    private static byte[] Golden(string file) =>
        File.ReadAllBytes(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential", "GoldenFiles", file));
}
