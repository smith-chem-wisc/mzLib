using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Development;
using NUnit.Framework;
using Readers;
using Test.TopDownEngine.Shared;

namespace Test.TopDownEngine.Harness;

[TestFixture]
[ExcludeFromCodeCoverage]
public class MetaMorpheusPlottingTests
{
    [Test]
    public void StaticRawWindow_CanBePassedToPlotWriter()
    {
        var psmPath = TopDownTestPaths.GetRequiredFile("MM_1p1p4_GPTMD_Search", "AllProteoforms.psmtsv");
        var rawPath = TopDownTestPaths.GetRequiredFile("02-18-20_jurkat_td_rep1_fract7.raw", "02-18-20_jurkat_td_rep2_fract7.raw");

        var topHit = SpectrumMatchTsvReader.ReadTsv(psmPath, out _).OrderByDescending(p => p.Score).First();

        var rawFile = MsDataFileReader.GetDataFile(rawPath);
        rawFile.LoadAllStaticData();

        var scans = rawFile.GetAllScansList()
            .Where(scan => scan.MsnOrder == 1 && Math.Abs(scan.RetentionTime - topHit.RetentionTime) <= 1.0)
            .ToArray();

        Assert.That(scans, Is.Not.Empty);

        var html = PlotlyScanWindowWriter.BuildHtml(scans, rawPath, topHit.RetentionTime, 1.0);
        Assert.That(html, Does.Contain("plotly-2.35.2.min.js"));
        Assert.That(html, Does.Contain(Path.GetFileName(rawPath)));
        Assert.That(html, Does.Contain("scattergl"));

        var output = Path.Combine(Path.GetTempPath(), $"topdown_plot_{Guid.NewGuid():N}.html");
        PlotlyScanWindowWriter.WriteHtml(scans, rawPath, topHit.RetentionTime, 1.0, output);
        Assert.That(File.Exists(output), Is.True);
    }
}
