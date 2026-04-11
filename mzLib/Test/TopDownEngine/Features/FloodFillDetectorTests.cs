using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using TopDownEngine.Features;

namespace Test.TopDownEngine.Features;

[TestFixture]
[ExcludeFromCodeCoverage]
public class FloodFillDetectorTests
{
    [TestCase(MzGrowthStrategy.SeededBinWalk)]
    [TestCase(MzGrowthStrategy.ConcentricRing)]
    public void Detect_RecoversThreeKnownFeatures_WithBoundsNearTruth(MzGrowthStrategy strategy)
    {
        PeakIndexingEngine index = BuildFeatureRichIndex();
        FloodFillDetector detector = new();

        FloodFillOptions options = new()
        {
            PeakFindingTolerance = new AbsoluteTolerance(0.005),
            MaxMissedScansAllowed = 1,
            MaxRtRange = 0.9,
            NumPeakThreshold = 5,
            MaxMzGrowthSteps = 4,
            MinAdjacentXicPeakCount = 4,
            GrowthStrategy = strategy,
            NoiseFloorIntensity = 40
        };

        IReadOnlyList<FeatureBox> boxes = detector.Detect(index, options, sourceFile: "synthetic");

        Assert.That(boxes.Count, Is.GreaterThanOrEqualTo(3));
        AssertContainsFeature(boxes, expectedMz: 500.00, expectedRt: 2.0);
        AssertContainsFeature(boxes, expectedMz: 600.00, expectedRt: 4.0);
        AssertContainsFeature(boxes, expectedMz: 700.00, expectedRt: 6.0);
    }

    [TestCase(MzGrowthStrategy.SeededBinWalk)]
    [TestCase(MzGrowthStrategy.ConcentricRing)]
    public void Detect_NoiseOnlyInput_ReturnsNoBoxes(MzGrowthStrategy strategy)
    {
        PeakIndexingEngine index = BuildNoiseOnlyIndex();
        FloodFillDetector detector = new();

        FloodFillOptions options = new()
        {
            PeakFindingTolerance = new AbsoluteTolerance(0.005),
            MaxMissedScansAllowed = 1,
            MaxRtRange = 0.6,
            NumPeakThreshold = 5,
            MaxMzGrowthSteps = 4,
            MinAdjacentXicPeakCount = 4,
            GrowthStrategy = strategy,
            NoiseFloorIntensity = 120
        };

        IReadOnlyList<FeatureBox> boxes = detector.Detect(index, options, sourceFile: "noise-only");
        Assert.That(boxes, Is.Empty);
    }

    private static void AssertContainsFeature(IReadOnlyList<FeatureBox> boxes, double expectedMz, double expectedRt)
    {
        FeatureBox match = boxes.FirstOrDefault(box =>
            box.MzRange.Minimum <= expectedMz + 0.02 &&
            box.MzRange.Maximum >= expectedMz - 0.02 &&
            box.RtRange.Minimum <= expectedRt &&
            box.RtRange.Maximum >= expectedRt);

        Assert.That(match, Is.Not.Null, $"Expected feature near m/z {expectedMz:F2}, RT {expectedRt:F2} was not recovered.");
    }

    private static PeakIndexingEngine BuildFeatureRichIndex()
    {
        List<MsDataScan> scans = new();
        int oneBasedScanNumber = 1;

        for (int i = 0; i < 91; i++)
        {
            double rt = 1.0 + (i * 0.1);
            List<(double mz, double intensity)> peaks = new()
            {
                (500.00, Gaussian(rt, 2.0, 0.18, 2000)),
                (500.01, Gaussian(rt, 2.0, 0.18, 1200)),
                (499.99, Gaussian(rt, 2.0, 0.18, 900)),

                (600.00, Gaussian(rt, 4.0, 0.2, 2200)),
                (600.01, Gaussian(rt, 4.0, 0.2, 1300)),
                (599.99, Gaussian(rt, 4.0, 0.2, 1000)),

                (700.00, Gaussian(rt, 6.0, 0.22, 1800)),
                (700.01, Gaussian(rt, 6.0, 0.22, 1100)),
                (699.99, Gaussian(rt, 6.0, 0.22, 800)),

                (450.0 + (i % 7) * 0.03, 15 + (i % 11)),
                (800.0 + (i % 5) * 0.04, 8 + (i % 9))
            };

            double[] mz = peaks.Select(p => p.mz).ToArray();
            double[] intensity = peaks.Select(p => p.intensity).ToArray();
            scans.Add(CreateMs1Scan(oneBasedScanNumber++, rt, mz, intensity));
        }

        return PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
    }

    private static PeakIndexingEngine BuildNoiseOnlyIndex()
    {
        List<MsDataScan> scans = new();
        int oneBasedScanNumber = 1;

        for (int i = 0; i < 60; i++)
        {
            double rt = 1.0 + (i * 0.1);
            double[] mz =
            {
                450.0 + (i % 5) * 0.02,
                530.0 + (i % 7) * 0.03,
                710.0 + (i % 11) * 0.01
            };

            double[] intensity =
            {
                10 + (i % 13),
                12 + (i % 9),
                7 + (i % 11)
            };

            scans.Add(CreateMs1Scan(oneBasedScanNumber++, rt, mz, intensity));
        }

        return PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
    }

    private static double Gaussian(double x, double mu, double sigma, double amplitude)
    {
        double z = (x - mu) / sigma;
        return amplitude * Math.Exp(-(z * z) / 2.0);
    }

    private static MsDataScan CreateMs1Scan(int oneBasedScanNumber, double retentionTime, double[] mzArray, double[] intensityArray)
    {
        return new MsDataScan(
            new MzSpectrum(mzArray, intensityArray, false),
            oneBasedScanNumber,
            1,
            true,
            Polarity.Positive,
            retentionTime,
            new MzRange(400, 1600),
            "f",
            MZAnalyzerType.Orbitrap,
            intensityArray.Sum(),
            1.0,
            null,
            $"scan={oneBasedScanNumber}");
    }
}
