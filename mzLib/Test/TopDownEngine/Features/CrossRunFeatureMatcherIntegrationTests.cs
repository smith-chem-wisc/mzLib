using FlashLFQ;
using FlashLFQ.IsoTracker;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using TopDownEngine.Alignment;
using TopDownEngine.Features;

namespace Test.TopDownEngine.Features;

[TestFixture]
[ExcludeFromCodeCoverage]
public class CrossRunFeatureMatcherIntegrationTests
{
    [Test]
    public void MatchDonorBoxes_KnownFeatureAcrossThreeRuns_AllRunsMatched()
    {
        SpectraFileInfo runA = BuildFileInfo("runA");
        SpectraFileInfo runB = BuildFileInfo("runB");
        SpectraFileInfo runC = BuildFileInfo("runC");

        Dictionary<SpectraFileInfo, PeakIndexingEngine> indexes = new()
        {
            [runA] = BuildIndex(featureRtOffset: 0.00),
            [runB] = BuildIndex(featureRtOffset: 0.20),
            [runC] = BuildIndex(featureRtOffset: -0.15)
        };

        Dictionary<SpectraFileInfo, Func<double, double>> warpsToReference = BuildWarpsToReference(runA, runB, runC);
        FeatureBox donorBox = new(
            MzRange: new MzRange(499.98, 500.02),
            RtRange: new DoubleRange(1.90, 2.10),
            SeedIntensity: 2200,
            TotalIntensity: 10800,
            SourceFile: "runA",
            PeakCount: 12);

        CrossRunFeatureMatcher matcher = new();
        IReadOnlyList<FeatureGroup> groups = matcher.MatchDonorBoxes(
            donorFile: runA,
            donorBoxes: new[] { donorBox },
            indexesByFile: indexes,
            rtWarpsToReference: warpsToReference,
            nullShifts: new[] { 0.35, 0.55 },
            pValueThreshold: 0.05);

        Assert.That(groups, Has.Count.EqualTo(1));
        FeatureGroup group = groups[0];

        Assert.That(group.Matches.ContainsKey(runA), Is.True);
        Assert.That(group.Matches.ContainsKey(runB), Is.True);
        Assert.That(group.Matches.ContainsKey(runC), Is.True);
        Assert.That(group.TotalSupport, Is.EqualTo(3));
    }

    [Test]
    public void MatchDonorBoxes_DonorOnlyFeature_NoiseInOtherRuns_NotMatched()
    {
        SpectraFileInfo runA = BuildFileInfo("runA");
        SpectraFileInfo runB = BuildFileInfo("runB");
        SpectraFileInfo runC = BuildFileInfo("runC");

        Dictionary<SpectraFileInfo, PeakIndexingEngine> indexes = new()
        {
            [runA] = BuildIndexWithDonorOnlyFeature(featureRtOffset: 0.00),
            [runB] = BuildNoiseOnlyIndex(featureRtOffset: 0.20),
            [runC] = BuildNoiseOnlyIndex(featureRtOffset: -0.15)
        };

        Dictionary<SpectraFileInfo, Func<double, double>> warpsToReference = BuildWarpsToReference(runA, runB, runC);
        FeatureBox donorBox = new(
            MzRange: new MzRange(619.98, 620.02),
            RtRange: new DoubleRange(4.90, 5.10),
            SeedIntensity: 2000,
            TotalIntensity: 10000,
            SourceFile: "runA",
            PeakCount: 12);

        CrossRunFeatureMatcher matcher = new();
        IReadOnlyList<FeatureGroup> groups = matcher.MatchDonorBoxes(
            donorFile: runA,
            donorBoxes: new[] { donorBox },
            indexesByFile: indexes,
            rtWarpsToReference: warpsToReference,
            nullShifts: new[] { 0.35, 0.55 },
            pValueThreshold: 0.05);

        Assert.That(groups, Has.Count.EqualTo(1));
        FeatureGroup group = groups[0];

        Assert.That(group.Matches.ContainsKey(runA), Is.True);
        Assert.That(group.Matches.ContainsKey(runB), Is.False);
        Assert.That(group.Matches.ContainsKey(runC), Is.False);
        Assert.That(group.TotalSupport, Is.EqualTo(1));
    }

    private static Dictionary<SpectraFileInfo, Func<double, double>> BuildWarpsToReference(
        SpectraFileInfo runA,
        SpectraFileInfo runB,
        SpectraFileInfo runC)
    {
        List<SpectraFileInfo> files = new() { runA, runB, runC };
        List<AnchorFileExtrema> extrema = new();

        for (int i = 0; i < 10; i++)
        {
            double referenceRt = 1.0 + (i * 0.5);
            AnchorBin anchor = new(i, 450 + i * 3, 2, 1000);

            extrema.Add(new AnchorFileExtrema(anchor, 0, new[] { new Extremum(1000, referenceRt, ExtremumType.Maximum) }));
            extrema.Add(new AnchorFileExtrema(anchor, 1, new[] { new Extremum(1000, referenceRt + 0.20, ExtremumType.Maximum) }));
            extrema.Add(new AnchorFileExtrema(anchor, 2, new[] { new Extremum(1000, referenceRt - 0.15, ExtremumType.Maximum) }));
        }

        IdentificationFreeRtAligner aligner = new();
        return aligner.BuildRtWarps(files, extrema, referenceFileIndex: 0, extremaMatchTolerance: 0.5);
    }

    private static PeakIndexingEngine BuildIndex(double featureRtOffset)
    {
        List<MsDataScan> scans = new();
        int oneBasedScan = 1;

        for (int i = 0; i < 80; i++)
        {
            double rt = 1.0 + (i * 0.1);
            List<(double mz, double intensity)> peaks = new()
            {
                (430.0 + (i % 5) * 0.02, 12 + (i % 7)),
                (710.0 + (i % 9) * 0.03, 10 + (i % 5))
            };

            if (Math.Abs(rt - (2.0 + featureRtOffset)) <= 0.10)
            {
                peaks.Add((499.99, 1500));
                peaks.Add((500.00, 2200));
                peaks.Add((500.01, 1400));
            }

            scans.Add(CreateMs1Scan(oneBasedScan++, rt, peaks));
        }

        return PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
    }

    private static PeakIndexingEngine BuildIndexWithDonorOnlyFeature(double featureRtOffset)
    {
        List<MsDataScan> scans = new();
        int oneBasedScan = 1;

        for (int i = 0; i < 90; i++)
        {
            double rt = 1.0 + (i * 0.1);
            List<(double mz, double intensity)> peaks = new()
            {
                (430.0 + (i % 5) * 0.02, 12 + (i % 7)),
                (710.0 + (i % 9) * 0.03, 10 + (i % 5))
            };

            if (Math.Abs(rt - (5.0 + featureRtOffset)) <= 0.10)
            {
                peaks.Add((619.99, 1300));
                peaks.Add((620.00, 2000));
                peaks.Add((620.01, 1200));
            }

            scans.Add(CreateMs1Scan(oneBasedScan++, rt, peaks));
        }

        return PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
    }

    private static PeakIndexingEngine BuildNoiseOnlyIndex(double featureRtOffset)
    {
        List<MsDataScan> scans = new();
        int oneBasedScan = 1;

        for (int i = 0; i < 90; i++)
        {
            double rt = 1.0 + (i * 0.1);
            List<(double mz, double intensity)> peaks = new()
            {
                (430.0 + (i % 5) * 0.02, 12 + (i % 7)),
                (620.03 + (i % 4) * 0.02, 14 + (i % 6)),
                (710.0 + (i % 9) * 0.03, 10 + (i % 5))
            };

            if (Math.Abs(rt - (5.0 + featureRtOffset)) <= 0.10)
            {
                peaks.Add((619.90, 40));
                peaks.Add((620.10, 35));
            }

            scans.Add(CreateMs1Scan(oneBasedScan++, rt, peaks));
        }

        return PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
    }

    private static SpectraFileInfo BuildFileInfo(string stem)
    {
        return new SpectraFileInfo($"{stem}.mzML", "Cond", 1, 1, 1);
    }

    private static MsDataScan CreateMs1Scan(int oneBasedScanNumber, double retentionTime, List<(double mz, double intensity)> peaks)
    {
        return new MsDataScan(
            new MzSpectrum(peaks.Select(p => p.mz).ToArray(), peaks.Select(p => p.intensity).ToArray(), false),
            oneBasedScanNumber,
            1,
            true,
            Polarity.Positive,
            retentionTime,
            new MzRange(400, 1600),
            "f",
            MZAnalyzerType.Orbitrap,
            peaks.Sum(p => p.intensity),
            1.0,
            null,
            $"scan={oneBasedScanNumber}");
    }
}
