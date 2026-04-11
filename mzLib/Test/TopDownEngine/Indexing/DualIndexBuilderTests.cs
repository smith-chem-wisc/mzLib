using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Test.TopDownEngine.Shared;
using TopDownEngine.Indexing;

namespace Test.TopDownEngine.Indexing;

[TestFixture]
[ExcludeFromCodeCoverage]
public class DualIndexBuilderTests
{
    [Test]
    public void Build_UsesSingleFineIndexForFineAndThickMzRtQueries()
    {
        string mzmlPath = TopDownTestPaths.GetRequiredFile("tiny.pwiz.1.1.mzML");
        MsDataFile file = MsDataFileReader.GetDataFile(mzmlPath);
        file.LoadAllStaticData();

        MsDataScan[] ms1Scans = file.GetMS1Scans()
            .Where(scan => scan != null && scan.MsnOrder == 1)
            .OrderBy(scan => scan.OneBasedScanNumber)
            .ToArray();

        Assert.That(ms1Scans.Length, Is.GreaterThan(1));

        MsDataScan targetScan = ms1Scans[0];
        double expectedMz = targetScan.MassSpectrum.XArray[0];
        int targetScanIndex = 0;
        double queryRetentionTime = targetScan.RetentionTime + 1e-6;

        DualIndexBuilder builder = new();
        (PeakIndexingEngine fine, ThickIndexView thick) = builder.Build(file);
        Assert.That(thick.FineIndex, Is.SameAs(fine));

        IIndexedPeak finePeak = fine.GetIndexedPeak(expectedMz, targetScanIndex, new AbsoluteTolerance(0.05));
        IIndexedPeak thickPeak = thick.GetIndexedPeak(expectedMz, targetScanIndex);

        Assert.That(finePeak, Is.Not.Null);
        Assert.That(thickPeak, Is.Not.Null);
        Assert.That(finePeak!.M, Is.EqualTo(expectedMz).Within(1e-12));
        Assert.That(thickPeak!.M, Is.EqualTo(expectedMz).Within(1e-12));
        Assert.That(finePeak.ZeroBasedScanIndex, Is.EqualTo(targetScanIndex));
        Assert.That(thickPeak.ZeroBasedScanIndex, Is.EqualTo(targetScanIndex));

        var fineXic = fine.GetXic(expectedMz, queryRetentionTime, new AbsoluteTolerance(0.05), missedScansAllowed: 0, maxPeakHalfWidth: 0.0);
        var thickXic = thick.GetXic(expectedMz, queryRetentionTime, missedScansAllowed: 0, maxPeakHalfWidth: 0.0);

        Assert.That(fineXic.Count, Is.EqualTo(1));
        Assert.That(thickXic.Count, Is.EqualTo(1));
        Assert.That(fineXic[0].M, Is.EqualTo(expectedMz).Within(1e-12));
        Assert.That(thickXic[0].M, Is.EqualTo(expectedMz).Within(1e-12));
        Assert.That(fineXic[0].ZeroBasedScanIndex, Is.EqualTo(targetScanIndex));
        Assert.That(thickXic[0].ZeroBasedScanIndex, Is.EqualTo(targetScanIndex));
    }
}
