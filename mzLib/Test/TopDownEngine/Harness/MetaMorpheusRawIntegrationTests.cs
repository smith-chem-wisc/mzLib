using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Test.TopDownEngine.Shared;
using MzLibUtil;

namespace Test.TopDownEngine.Harness;

[TestFixture]
[ExcludeFromCodeCoverage]
public class MetaMorpheusRawIntegrationTests
{
    [Test]
    public void HighestScoringProteoform_CanDriveRawWindowIndexing()
    {
        var psmPath = TopDownTestPaths.FindExistingFile(
            "Frac7_GPTMD_Search/Task2-TopDownSearch/AllProteoforms.psmtsv",
            "Frac7_GPTMD_Search\\Task2-TopDownSearch\\AllProteoforms.psmtsv",
            "MM_1p1p4_GPTMD_Search/AllProteoforms.psmtsv",
            "MM_1p1p4_GPTMD_Search\\AllProteoforms.psmtsv",
            "AllProteoforms.psmtsv");
        var rawPath = TopDownTestPaths.FindExistingFile("02-18-20_jurkat_td_rep1_fract7.raw", "02-18-20_jurkat_td_rep2_fract7.raw");
        if (psmPath is null || rawPath is null)
        {
            Assert.Ignore("MetaMorpheus/RAW integration fixtures are not present in local test data.");
        }

        var proteoforms = SpectrumMatchTsvReader.ReadTsv(psmPath!, out var warnings);
        Assert.That(proteoforms, Is.Not.Empty);
        Assert.That(warnings, Is.Not.Null);

        var topHit = proteoforms.OrderByDescending(p => p.Score).First();
        Assert.That(topHit.RetentionTime, Is.GreaterThanOrEqualTo(0));
        Assert.That(topHit.MonoisotopicMass, Is.GreaterThan(0));

        var rawFile = MsDataFileReader.GetDataFile(rawPath!);
        rawFile.LoadAllStaticData();

        var windowStart = Math.Max(0, topHit.RetentionTime - 1.0);
        var windowEnd = topHit.RetentionTime + 1.0;

        var scans = rawFile.GetAllScansList()
            .Where(scan => scan.MsnOrder == 1 && scan.RetentionTime >= windowStart && scan.RetentionTime <= windowEnd)
            .ToArray();

        Assert.That(scans, Is.Not.Empty);

        var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(scans);
        Assert.That(indexingEngine, Is.Not.Null);

        var apexScan = scans.OrderByDescending(s => s.TotalIonCurrent).First();
        var precursorCharge = Math.Max(1, topHit.PrecursorCharge);
        var expectedMz = topHit.MonoisotopicMass.ToMz(precursorCharge);

        var xic = indexingEngine!.GetXicByScanIndex(expectedMz, apexScan.OneBasedScanNumber - 1, new PpmTolerance(20), 2, 1.0);

        Assert.That(xic, Is.Not.Null);
    }
}
