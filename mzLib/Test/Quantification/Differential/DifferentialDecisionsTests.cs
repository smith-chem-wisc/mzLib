using System;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The decisions file (user requirement R-DECISIONS, <c>DEF-DIFF-DECISIONS</c>): one line per place a rule changed what
/// happened, keyed to the results table, every rule from the fixed list with its ruling, the settings once at analysis
/// scope; the same bytes for the same analysis (any thread count, any culture), a golden file, and the refusals.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DifferentialDecisionsTests
{
    private static DifferentialAnalysisResult Run(int threads = 1) => DifferentialAnalysis.Run(DifferentialAnalysisTests.Request(), threads);

    private static byte[] Bytes(DifferentialAnalysisResult r) => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows, r.Decisions);

    [Test]
    public void TheFileMatchesItsGoldenAndIsTheSameForAnyThreadCountAndCulture()
    {
        byte[] one = Bytes(Run(1));
        Assert.That(Bytes(Run(4)), Is.EqualTo(one), "thread count");
        var culture = Thread.CurrentThread.CurrentCulture;
        try
        {
            Thread.CurrentThread.CurrentCulture = new CultureInfo("de-DE");
            Assert.That(Bytes(Run(1)), Is.EqualTo(one), "culture");
        }
        finally
        {
            Thread.CurrentThread.CurrentCulture = culture;
        }
        Assert.That(one.Take(3), Is.Not.EqualTo(new byte[] { 0xEF, 0xBB, 0xBF }));
        Assert.That(one, Is.EqualTo(File.ReadAllBytes(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification",
            "Differential", "GoldenFiles", "DifferentialDecisions.golden.tsv"))), "byte for byte");
    }

    [Test]
    public void TheSettingsComeFirstOnceAndEveryLineCarriesItsRuling()
    {
        var r = Run();
        var lines = Encoding.UTF8.GetString(Bytes(r)).Split('\n').Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
        Assert.That(lines[0], Is.EqualTo(new[]
            { "analysis_id", "scope", "ruling", "rule", "stratum", "quant_basis", "feature_id", "contrast_id", "facts", "effect" }));
        Assert.That(lines[1][1..4], Is.EqualTo(new[] { "analysis", "M5b-14", "settings" }));
        Assert.That(lines[1][8], Does.Contain("at least 2 samples with a value per side").And.Contain("robust weighting on"));
        Assert.That(lines.Skip(1).Count(l => l[1] == "analysis"), Is.EqualTo(1));
        Assert.That(lines.Skip(1).Select(l => l[0]), Is.All.EqualTo(r.Metadata.AnalysisId));
        Assert.That(lines.Skip(1).Select(l => l[2]), Is.All.Matches<string>(s => s.StartsWith("M5b-")));
        var p3 = lines.Single(l => l[3] == DifferentialDecisionRules.ReducedDesign && l[6] == "P3");
        Assert.That(p3[1], Is.EqualTo(DifferentialDecisionWriter.ScopeFeature));
        var p2 = lines.Single(l => l[3] == DifferentialDecisionRules.AbsentSide && l[6] == "P2" && l[7] == "c2");
        Assert.That(p2[1], Is.EqualTo(DifferentialDecisionWriter.ScopeFeatureContrast));
        Assert.That(lines.Single(l => l[3] == DifferentialDecisionRules.RowSet)[1], Is.EqualTo(DifferentialDecisionWriter.ScopeStratum));
    }

    [Test]
    public void TheCountsGiveEachRulesFeaturesAndRowsInRuleOrder()
    {
        var counts = DifferentialDecisionWriter.Counts(Run().Decisions);

        Assert.That(counts.Select(c => DifferentialDecisionRules.Index(c.Rule)), Is.Ordered);
        Assert.That(counts.First().Rule, Is.EqualTo(DifferentialDecisionRules.Settings));
        var absent = counts.Single(c => c.Rule == DifferentialDecisionRules.AbsentSide);
        Assert.That(absent.Features, Is.EqualTo(2), "P2 lacks C, P3 lacks A");
        Assert.That(absent.Rows, Is.EqualTo(4));
    }

    [Test]
    public void ALineThatMatchesNothingOrBreaksTheRulesIsRefused()
    {
        var r = Run();
        string stratum = r.Metadata.Strata[0].Stratum;
        DifferentialDecision Line(string rule, string s, string? feature = null, string? contrast = null, string? basis = "msms_only") =>
            new(rule, s, basis, feature, contrast, "facts", "effect");

        Assert.DoesNotThrow(() => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows, new[] { Line(DifferentialDecisionRules.FitOutcome, stratum, "P1", "c1") }));
        Assert.Throws<ArgumentException>(() => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows, new[] { Line("made_up", stratum) }), "unknown rule");
        Assert.Throws<ArgumentException>(() => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows,
            new[] { Line(DifferentialDecisionRules.Settings, DifferentialDecisionWriter.AnalysisStratum, basis: "msms_only") }), "an analysis line with a key");
        Assert.Throws<ArgumentException>(() => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows, new[] { Line(DifferentialDecisionRules.RowSet, "nowhere") }), "unknown stratum");
        Assert.Throws<ArgumentException>(() => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows, new[] { Line(DifferentialDecisionRules.SinglePeptide, stratum, "NOSUCH") }), "no row for the feature");
        Assert.Throws<ArgumentException>(() => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows, new[] { Line(DifferentialDecisionRules.FitOutcome, stratum, "P1", "c9") }), "no row for the comparison");
        Assert.Throws<ArgumentException>(() => DifferentialDecisionWriter.ToBytes(r.Metadata, r.Rows,
            new[] { new DifferentialDecision(DifferentialDecisionRules.RowSet, stratum, "msms_only", null, null, " ", "effect") }), "no facts");
    }

    [Test]
    public void WritingTheResultPutsDownTheThreeTablesAndTheReport()
    {
        var r = Run();
        string dir = Path.Combine(Path.GetTempPath(), "mzlib-decisions-" + Guid.NewGuid().ToString("N"));
        try
        {
            string path = r.Write(dir);
            Assert.That(Path.GetFileName(path), Is.EqualTo(DifferentialResultWriter.ResultsFileName));
            Assert.That(File.ReadAllBytes(Path.Combine(dir, DifferentialDecisionWriter.FileName)), Is.EqualTo(Bytes(r)));
            Assert.That(File.Exists(Path.Combine(dir, DifferentialResultWriter.MetadataFileName)), Is.True);
            Assert.That(Directory.GetFiles(dir).Select(Path.GetFileName), Is.EquivalentTo(new[]
            {
                DifferentialResultWriter.ResultsFileName, DifferentialResultWriter.MetadataFileName, DifferentialDecisionWriter.FileName,
                StatisticalReportWriter.ReportFileName,
            }), "no temporary file left");
            string report = File.ReadAllText(Path.Combine(dir, StatisticalReportWriter.ReportFileName));
            Assert.That(report, Does.Contain("## 6. Decisions applied"));
            Assert.That(report, Does.Contain("Mean-variance: stratum"));
            Assert.That(report, Does.Contain("Sample correlation: stratum"));
            var figures = Directory.GetFiles(Path.Combine(dir, StatisticalReportWriter.FiguresDirectoryName)).Select(Path.GetFileName).ToList();
            Assert.That(figures.Where(f => f!.Contains("mean_variance")), Is.Not.Empty);
            Assert.That(figures.Where(f => f!.Contains("sample_correlation")), Is.Not.Empty);
        }
        finally
        {
            if (Directory.Exists(dir)) Directory.Delete(dir, true);
        }
    }
}
