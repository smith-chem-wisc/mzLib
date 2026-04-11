using FlashLFQ.IsoTracker;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using TopDownEngine.Indexing;

namespace TopDownEngine.Alignment;

public sealed record AnchorFileExtrema(AnchorBin Anchor, int FileIndex, IReadOnlyList<Extremum> Extrema);

public class AnchorXicExtremaExtractor
{
    public IReadOnlyList<AnchorFileExtrema> ExtractExtrema(
        IReadOnlyList<AnchorBin> anchors,
        IReadOnlyList<ThickIndexView> thickIndexes,
        int missedScansAllowed,
        double splineRtInterval = 0.02)
    {
        if (anchors == null)
        {
            throw new ArgumentNullException(nameof(anchors));
        }

        if (thickIndexes == null)
        {
            throw new ArgumentNullException(nameof(thickIndexes));
        }

        if (missedScansAllowed < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(missedScansAllowed), "missedScansAllowed must be non-negative.");
        }

        if (splineRtInterval <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(splineRtInterval), "splineRtInterval must be positive.");
        }

        List<AnchorFileExtrema> output = new(anchors.Count * thickIndexes.Count);

        for (int anchorIndex = 0; anchorIndex < anchors.Count; anchorIndex++)
        {
            AnchorBin anchor = anchors[anchorIndex];
            for (int fileIndex = 0; fileIndex < thickIndexes.Count; fileIndex++)
            {
                ThickIndexView thickIndex = thickIndexes[fileIndex];
                if (thickIndex == null)
                {
                    output.Add(new AnchorFileExtrema(anchor, fileIndex, Array.Empty<Extremum>()));
                    continue;
                }

                List<IIndexedPeak> xic = ExtractAnchorXic(thickIndex, anchor.MzCenter, missedScansAllowed);
                IReadOnlyList<Extremum> extrema = FindExtremaFromSmoothedXic(xic, splineRtInterval);
                output.Add(new AnchorFileExtrema(anchor, fileIndex, extrema));
            }
        }

        return output;
    }

    private static List<IIndexedPeak> ExtractAnchorXic(ThickIndexView thickIndex, double anchorMz, int missedScansAllowed)
    {
        var bins = thickIndex.GetBinsInRange(anchorMz);
        IIndexedPeak apexPeak = bins
            .SelectMany(bin => bin)
            .OrderByDescending(peak => peak.Intensity)
            .FirstOrDefault();

        if (apexPeak == null)
        {
            return new List<IIndexedPeak>();
        }

        return thickIndex.GetXic(anchorMz, apexPeak.RetentionTime, missedScansAllowed);
    }

    private static IReadOnlyList<Extremum> FindExtremaFromSmoothedXic(List<IIndexedPeak> xic, double splineRtInterval)
    {
        if (xic == null || xic.Count < 5)
        {
            return Array.Empty<Extremum>();
        }

        XicCubicSpline spline = new(splineRtInterval);
        float[] rtArray = xic.Select(peak => peak.RetentionTime).ToArray();
        float[] intensityArray = xic.Select(peak => peak.Intensity).ToArray();
        (double Rt, double Intensity)[] smoothed = spline.GetXicSplineData(rtArray, intensityArray, rtArray.First(), rtArray.Last());

        List<Extremum> extrema = new();
        for (int i = 1; i < smoothed.Length - 1; i++)
        {
            double previousSlope = smoothed[i].Intensity - smoothed[i - 1].Intensity;
            double nextSlope = smoothed[i + 1].Intensity - smoothed[i].Intensity;

            if (previousSlope > 0 && nextSlope <= 0)
            {
                extrema.Add(new Extremum(smoothed[i].Intensity, smoothed[i].Rt, ExtremumType.Maximum));
            }
            else if (previousSlope < 0 && nextSlope >= 0)
            {
                extrema.Add(new Extremum(smoothed[i].Intensity, smoothed[i].Rt, ExtremumType.Minimum));
            }
        }

        return extrema;
    }
}
