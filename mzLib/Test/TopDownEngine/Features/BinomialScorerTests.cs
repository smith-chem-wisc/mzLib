using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using TopDownEngine.Features;

namespace Test.TopDownEngine.Features;

[TestFixture]
public class BinomialScorerTests
{
    [Test]
    public void CountPeaksInBox_UsesInclusiveLowerAndExclusiveUpperBounds_WithScanFiltering()
    {
        PeakIndexingEngine index = BuildIndex();

        int count = BinomialScorer.CountPeaksInBox(
            index,
            rtRange: new DoubleRange(1.5, 2.5),
            mzRange: new DoubleRange(500.00, 500.10));

        Assert.That(count, Is.EqualTo(2));
    }

    [Test]
    public void ScoreMatch_ClearlyMatchedBox_YieldsLowPValue()
    {
        PeakIndexingEngine index = BuildScoringIndex();
        FeatureBox donorBox = new(
            MzRange: new MzRange(499.98, 500.02),
            RtRange: new DoubleRange(1.92, 2.08),
            SeedIntensity: 1000,
            TotalIntensity: 9000,
            SourceFile: "donor",
            PeakCount: 9);

        BinomialMatchScore score = BinomialScorer.ScoreMatch(
            donorBox,
            index,
            rtWarp: rt => rt + 0.2,
            nullShifts: new[] { 0.35, 0.55 });

        Assert.That(score.KTrue, Is.GreaterThanOrEqualTo(8));
        Assert.That(score.PNull, Is.LessThan(0.25));
        Assert.That(score.PValue, Is.LessThan(0.01));
    }

    [Test]
    public void ScoreMatch_NoiseOnlyBox_YieldsHighPValue()
    {
        PeakIndexingEngine index = BuildScoringIndex();
        FeatureBox donorBox = new(
            MzRange: new MzRange(700.0, 700.04),
            RtRange: new DoubleRange(1.92, 2.08),
            SeedIntensity: 50,
            TotalIntensity: 200,
            SourceFile: "donor-noise",
            PeakCount: 9);

        BinomialMatchScore score = BinomialScorer.ScoreMatch(
            donorBox,
            index,
            rtWarp: rt => rt + 0.2,
            nullShifts: new[] { 0.35, 0.55 });

        Assert.That(score.PValue, Is.GreaterThan(0.1));
    }

    private static PeakIndexingEngine BuildIndex()
    {
        MsDataScan[] scans =
        {
            CreateMs1Scan(1, 1.0, new[] { 500.00, 500.05, 600.00 }, new[] { 1000.0, 1000.0, 100.0 }),
            CreateMs1Scan(2, 1.5, new[] { 500.00, 500.10, 600.00 }, new[] { 1000.0, 1000.0, 100.0 }),
            CreateMs1Scan(3, 2.0, new[] { 500.11, 500.05, 600.00 }, new[] { 1000.0, 1000.0, 100.0 }),
            CreateMs1Scan(4, 2.5, new[] { 500.05, 500.00, 600.00 }, new[] { 1000.0, 1000.0, 100.0 })
        };

        return PeakIndexingEngine.InitializeIndexingEngine(scans)!;
    }

    private static PeakIndexingEngine BuildScoringIndex()
    {
        List<MsDataScan> scans = new();
        int oneBasedScanNumber = 1;

        for (int i = 0; i < 70; i++)
        {
            double rt = 1.0 + (i * 0.05);
            List<(double mz, double intensity)> peaks = new()
            {
                (430.0 + ((i % 7) * 0.02), 15 + (i % 5)),
                (610.0 + ((i % 5) * 0.02), 12 + (i % 7)),
                (820.0 + ((i % 9) * 0.03), 10 + (i % 3))
            };

            if (rt >= 2.12 && rt < 2.28)
            {
                peaks.Add((499.99, 1500));
                peaks.Add((500.00, 2200));
                peaks.Add((500.01, 1400));
            }

            scans.Add(CreateMs1Scan(
                oneBasedScanNumber++,
                rt,
                peaks.Select(p => p.mz).ToArray(),
                peaks.Select(p => p.intensity).ToArray()));
        }

        return PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
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
