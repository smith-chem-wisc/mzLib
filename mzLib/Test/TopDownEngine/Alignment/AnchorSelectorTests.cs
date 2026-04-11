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
public class AnchorSelectorTests
{
    [Test]
    public void SelectTopK_PrefersSharedHighSignalBinsAcrossFiles()
    {
        ThickIndexView file1 = BuildIndex(
            (500.05, 2100.0),
            (700.15, 1900.0),
            (450.01, 1300.0),
            (610.33, 900.0));

        ThickIndexView file2 = BuildIndex(
            (500.05, 2200.0),
            (700.15, 1850.0),
            (620.44, 1250.0),
            (605.22, 800.0));

        ThickIndexView file3 = BuildIndex(
            (730.27, 2300.0),
            (810.88, 1800.0),
            (650.12, 700.0),
            (512.01, 650.0));

        AnchorSelector selector = new();
        AnchorBin[] anchors = selector.SelectTopK(
            new[] { file1, file2, file3 },
            topK: 2,
            minFilesForAnchor: 2,
            intensityThreshold: 1000);

        Assert.That(anchors, Has.Length.EqualTo(2));

        HashSet<int> selectedBinIndices = anchors.Select(anchor => anchor.BinIndex).ToHashSet();
        HashSet<int> expectedBinIndices = new()
        {
            GetThickBinIndex(500.05),
            GetThickBinIndex(700.15)
        };

        Assert.That(selectedBinIndices.SetEquals(expectedBinIndices), Is.True);
        Assert.That(anchors.All(anchor => anchor.FileCount == 2), Is.True);
    }

    private static ThickIndexView BuildIndex(params (double mz, double intensity)[] peaks)
    {
        MsDataScan[] scans =
        {
            CreateMs1Scan(1, 1.0, peaks.Select(p => p.mz).ToArray(), peaks.Select(p => p.intensity).ToArray()),
            CreateMs1Scan(2, 1.1, peaks.Select(p => p.mz).ToArray(), peaks.Select(p => p.intensity * 0.95).ToArray())
        };

        PeakIndexingEngine fineIndex = PeakIndexingEngine.InitializeIndexingEngine(scans)!;
        return new ThickIndexView(fineIndex);
    }

    private static int GetThickBinIndex(double mz)
    {
        int fineBinIndex = (int)Math.Round(mz * 100.0, 0);
        return fineBinIndex / 10;
    }

    private static MsDataScan CreateMs1Scan(int oneBasedScanNumber, double retentionTime, double[] mzArray, double[] intensityArray)
    {
        return new MsDataScan(new MzSpectrum(mzArray, intensityArray, false), oneBasedScanNumber, 1, true,
            Polarity.Positive, retentionTime, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap,
            intensityArray.Sum(), 1.0, null, $"scan={oneBasedScanNumber}");
    }
}
