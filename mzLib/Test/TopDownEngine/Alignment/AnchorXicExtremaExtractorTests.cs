using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using TopDownEngine.Alignment;
using TopDownEngine.Indexing;

namespace Test.TopDownEngine.Alignment;

[TestFixture]
[ExcludeFromCodeCoverage]
public class AnchorXicExtremaExtractorTests
{
    [Test]
    public void ExtractExtrema_GaussianAnchorYieldsSingleMaximumNearExpectedRetentionTime()
    {
        const double anchorMz = 500.05;
        const double expectedRt = 3.0;

        ThickIndexView thickIndex = BuildGaussianAnchorIndex(anchorMz, expectedRt);
        AnchorSelector selector = new();
        AnchorBin[] anchors = selector.SelectTopK(new[] { thickIndex }, topK: 1, minFilesForAnchor: 1, intensityThreshold: 10);

        AnchorXicExtremaExtractor extractor = new();
        IReadOnlyList<AnchorFileExtrema> extremaByFile = extractor.ExtractExtrema(anchors, new[] { thickIndex }, missedScansAllowed: 2, splineRtInterval: 0.01);

        Assert.That(extremaByFile, Has.Count.EqualTo(1));

        AnchorFileExtrema fileResult = extremaByFile[0];
        var maxima = fileResult.Extrema.Where(ext => ext.Type == global::FlashLFQ.IsoTracker.ExtremumType.Maximum).ToList();

        Assert.That(maxima.Count, Is.GreaterThanOrEqualTo(1));
        var closestMaximum = maxima.OrderBy(ext => Math.Abs(ext.RetentionTime - expectedRt)).First();
        Assert.That(closestMaximum.RetentionTime, Is.EqualTo(expectedRt).Within(0.06));
    }

    private static ThickIndexView BuildGaussianAnchorIndex(double anchorMz, double centerRt)
    {
        const int scanCount = 41;
        const double rtStart = 1.0;
        const double rtStep = 0.1;
        const double sigma = 0.25;
        const double amplitude = 5000.0;

        List<MsDataScan> scans = new(scanCount);
        for (int i = 0; i < scanCount; i++)
        {
            double rt = rtStart + i * rtStep;
            double gaussian = amplitude * Math.Exp(-Math.Pow(rt - centerRt, 2) / (2 * sigma * sigma));

            double[] mzArray = { anchorMz, 700.20 };
            double[] intensityArray = { gaussian, 50.0 };
            scans.Add(CreateMs1Scan(i + 1, rt, mzArray, intensityArray));
        }

        PeakIndexingEngine fineIndex = PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
        return new ThickIndexView(fineIndex);
    }

    private static MsDataScan CreateMs1Scan(int oneBasedScanNumber, double retentionTime, double[] mzArray, double[] intensityArray)
    {
        return new MsDataScan(new MzSpectrum(mzArray, intensityArray, false), oneBasedScanNumber, 1, true,
            Polarity.Positive, retentionTime, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap,
            intensityArray.Sum(), 1.0, null, $"scan={oneBasedScanNumber}");
    }
}
