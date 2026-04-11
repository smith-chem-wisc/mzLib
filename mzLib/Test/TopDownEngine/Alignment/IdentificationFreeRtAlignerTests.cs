using FlashLFQ.IsoTracker;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using TopDownEngine.Alignment;

namespace Test.TopDownEngine.Alignment;

[TestFixture]
[ExcludeFromCodeCoverage]
public class IdentificationFreeRtAlignerTests
{
    [Test]
    public void BuildRtWarps_LinearDriftRecoversKnownShiftWithinOneScan()
    {
        const double slope = 1.15;
        const double intercept = 0.3;
        const double scanSpacing = 0.1;

        SpectraFileInfo referenceFile = BuildFileInfo("reference");
        SpectraFileInfo shiftedFile = BuildFileInfo("shifted");
        List<SpectraFileInfo> files = new() { referenceFile, shiftedFile };

        List<AnchorFileExtrema> extrema = new();
        for (int i = 0; i < 20; i++)
        {
            double referenceRt = 1.0 + i * 0.4;
            double shiftedRt = slope * referenceRt + intercept;
            AnchorBin anchor = new(i, 500 + i, 2, 1000);

            extrema.Add(new AnchorFileExtrema(anchor, 0, new[] { new Extremum(1000, referenceRt, ExtremumType.Maximum) }));
            extrema.Add(new AnchorFileExtrema(anchor, 1, new[] { new Extremum(1000, shiftedRt, ExtremumType.Maximum) }));
        }

        IdentificationFreeRtAligner aligner = new();
        Dictionary<SpectraFileInfo, Func<double, double>> warps = aligner.BuildRtWarps(files, extrema, extremaMatchTolerance: 4.0);

        Func<double, double> shiftedWarp = warps[shiftedFile];
        for (int i = 0; i < 20; i++)
        {
            double referenceRt = 1.0 + i * 0.4;
            double shiftedRt = slope * referenceRt + intercept;
            Assert.That(shiftedWarp(shiftedRt), Is.EqualTo(referenceRt).Within(scanSpacing));
        }
    }

    [Test]
    public void BuildRtWarps_NonLinearDriftReducesRmsErrorBelowThreshold()
    {
        SpectraFileInfo referenceFile = BuildFileInfo("reference");
        SpectraFileInfo distortedFile = BuildFileInfo("distorted");
        List<SpectraFileInfo> files = new() { referenceFile, distortedFile };

        List<AnchorFileExtrema> extrema = new();
        for (int i = 0; i < 28; i++)
        {
            double referenceRt = 1.0 + i * 0.3;
            double distortedRt = referenceRt + 0.12 * referenceRt + 0.22 * Math.Sin(0.7 * referenceRt);
            AnchorBin anchor = new(i, 650 + i, 2, 1000);

            extrema.Add(new AnchorFileExtrema(anchor, 0, new[] { new Extremum(1000, referenceRt, ExtremumType.Maximum) }));
            extrema.Add(new AnchorFileExtrema(anchor, 1, new[] { new Extremum(1000, distortedRt, ExtremumType.Maximum) }));
        }

        IdentificationFreeRtAligner aligner = new();
        Dictionary<SpectraFileInfo, Func<double, double>> warps = aligner.BuildRtWarps(files, extrema, extremaMatchTolerance: 4.0);

        Func<double, double> distortedWarp = warps[distortedFile];
        double[] referenceRts = Enumerable.Range(0, 28).Select(i => 1.0 + i * 0.3).ToArray();
        double[] distortedRts = referenceRts.Select(rt => rt + 0.12 * rt + 0.22 * Math.Sin(0.7 * rt)).ToArray();

        double baselineRms = CalculateRms(distortedRts.Zip(referenceRts, (distorted, reference) => distorted - reference));
        double alignedRms = CalculateRms(distortedRts.Zip(referenceRts, (distorted, reference) => distortedWarp(distorted) - reference));

        Assert.That(alignedRms, Is.LessThan(0.08));
        Assert.That(alignedRms, Is.LessThan(baselineRms * 0.35));
    }

    private static SpectraFileInfo BuildFileInfo(string stem)
    {
        return new SpectraFileInfo($"{stem}.mzML", "Cond", 1, 1, 1);
    }

    private static double CalculateRms(IEnumerable<double> errors)
    {
        double[] errorArray = errors.ToArray();
        return Math.Sqrt(errorArray.Select(e => e * e).Average());
    }
}
