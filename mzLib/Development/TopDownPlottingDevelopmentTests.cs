using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Development;

namespace Development;

[TestFixture]
[ExcludeFromCodeCoverage]
public class TopDownPlottingDevelopmentTests
{
    [Test]
    public void StaticLoadAndPlotWindow_WritesHtml()
    {
        var dataDir = GetDataDir();
        var rawPath = Path.Combine(dataDir, "02-18-20_jurkat_td_rep1_fract7.raw");
        var psmPath = Path.Combine(dataDir, "MM_1p1p4_GPTMD_Search", "AllProteoforms.psmtsv");

        Assert.That(File.Exists(rawPath), Is.True);
        Assert.That(File.Exists(psmPath), Is.True);

        var topHit = SpectrumMatchTsvReader.ReadTsv(psmPath, out _).OrderByDescending(p => p.Score).First();

        var rawFile = MsDataFileReader.GetDataFile(rawPath);
        rawFile.LoadAllStaticData();

        var scans = rawFile.GetAllScansList()
            .Where(s => s.MsnOrder == 1 && Math.Abs(s.RetentionTime - topHit.RetentionTime) <= 1.0)
            .ToArray();

        Assert.That(scans, Is.Not.Empty);

        var html = PlotlyScanWindowWriter.BuildHtml(scans, rawPath, topHit.RetentionTime, 1.0);
        Assert.That(html, Does.Contain("plotly-2.35.2.min.js"));
        Assert.That(html, Does.Contain("scattergl"));

        var output = Path.Combine(Path.GetTempPath(), $"topdown_plot_{Guid.NewGuid():N}.html");
        PlotlyScanWindowWriter.WriteHtml(scans, rawPath, topHit.RetentionTime, 1.0, output);
        Assert.That(File.Exists(output), Is.True);
    }

    private static string GetDataDir()
    {
        var dir = new DirectoryInfo(TestContext.CurrentContext.TestDirectory);
        while (dir != null)
        {
            var candidate = Path.Combine(dir.FullName, "Data");
            if (Directory.Exists(candidate))
                return candidate;
            dir = dir.Parent;
        }

        throw new DirectoryNotFoundException($"Could not find Data directory starting from: {TestContext.CurrentContext.TestDirectory}");
    }
}
