using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownEngine.Indexing;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.TopDownEngine.Indexing;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class ThickIndexViewTests
{
    [Test]
    public void ThickIndexView_CoarsensFineBinsAcrossPointOneBoundary()
    {
        PeakIndexingEngine fineIndex = BuildFineIndexWithBoundaryCrossingPeaks();
        ThickIndexView thickIndex = new(fineIndex);

        const double targetMz = 500.05;
        List<List<IndexedMassSpectralPeak>> bins = thickIndex.GetBinsInRange(targetMz);
        var flattenedMasses = bins.SelectMany(bin => bin).Select(peak => Math.Round(peak.M, 3)).ToHashSet();

        Assert.That(flattenedMasses.Contains(500.049), Is.True);
        Assert.That(flattenedMasses.Contains(500.051), Is.True);
    }

    [Test]
    public void ThickIndexView_MatchesFineIndexPeakAndXicAtFixedPointOneWindow()
    {
        PeakIndexingEngine fineIndex = BuildFineIndexWithBoundaryCrossingPeaks();
        ThickIndexView thickIndex = new(fineIndex);

        const double targetMz = 500.05;
        const int startScanIndex = 2;

        IIndexedPeak finePeak = fineIndex.GetIndexedPeak(targetMz, zeroBasedScanIndex: 2, new AbsoluteTolerance(0.05));
        IIndexedPeak thickPeak = thickIndex.GetIndexedPeak(targetMz, zeroBasedScanIndex: 2);

        Assert.That(finePeak, Is.Not.Null);
        Assert.That(thickPeak, Is.Not.Null);
        Assert.That(thickPeak!.M, Is.EqualTo(finePeak!.M).Within(1e-12));
        Assert.That(thickPeak.ZeroBasedScanIndex, Is.EqualTo(finePeak.ZeroBasedScanIndex));

        List<IIndexedPeak> fineXic = fineIndex.GetXicByScanIndex(targetMz, startScanIndex, new AbsoluteTolerance(0.05), missedScansAllowed: 1);
        List<IIndexedPeak> thickXic = thickIndex.GetXicByScanIndex(targetMz, startScanIndex, missedScansAllowed: 1);

        Assert.That(thickXic.Count, Is.EqualTo(fineXic.Count));
        for (int i = 0; i < fineXic.Count; i++)
        {
            Assert.That(thickXic[i].ZeroBasedScanIndex, Is.EqualTo(fineXic[i].ZeroBasedScanIndex));
            Assert.That(thickXic[i].M, Is.EqualTo(fineXic[i].M).Within(1e-12));
        }
    }

    private static PeakIndexingEngine BuildFineIndexWithBoundaryCrossingPeaks()
    {
        MsDataScan[] scans =
        {
            CreateMs1Scan(1, 1.0, new[] { 500.049 }, new[] { 1000.0 }),
            CreateMs1Scan(2, 1.1, new[] { 500.051 }, new[] { 1100.0 }),
            CreateMs1Scan(3, 1.2, new[] { 500.049, 500.051 }, new[] { 1200.0, 1300.0 }),
            CreateMs1Scan(4, 1.3, new[] { 500.051 }, new[] { 1250.0 }),
            CreateMs1Scan(5, 1.4, new[] { 500.049 }, new[] { 1150.0 })
        };

        return PeakIndexingEngine.InitializeIndexingEngine(scans)!;
    }

    private static MsDataScan CreateMs1Scan(int oneBasedScanNumber, double retentionTime, double[] mzArray, double[] intensityArray)
    {
        return new MsDataScan(new MzSpectrum(mzArray, intensityArray, false), oneBasedScanNumber, 1, true,
            Polarity.Positive, retentionTime, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap,
            intensityArray.Sum(), 1.0, null, $"scan={oneBasedScanNumber}");
    }
}
