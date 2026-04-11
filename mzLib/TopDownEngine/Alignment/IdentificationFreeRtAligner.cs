using FlashLFQ.IsoTracker;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace TopDownEngine.Alignment;

public class IdentificationFreeRtAligner
{
    public Dictionary<SpectraFileInfo, Func<double, double>> BuildRtWarps(
        IReadOnlyList<SpectraFileInfo> files,
        IReadOnlyList<AnchorFileExtrema> extremaByFile,
        int referenceFileIndex = 0,
        double extremaMatchTolerance = 0.5,
        int minMatchedExtrema = 3,
        bool useMonotoneSpline = false)
    {
        if (files == null)
        {
            throw new ArgumentNullException(nameof(files));
        }

        if (extremaByFile == null)
        {
            throw new ArgumentNullException(nameof(extremaByFile));
        }

        if (files.Count == 0)
        {
            throw new ArgumentException("At least one file is required.", nameof(files));
        }

        if (referenceFileIndex < 0 || referenceFileIndex >= files.Count)
        {
            throw new ArgumentOutOfRangeException(nameof(referenceFileIndex), "referenceFileIndex must point to a file in the input list.");
        }

        if (extremaMatchTolerance <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(extremaMatchTolerance), "extremaMatchTolerance must be positive.");
        }

        if (minMatchedExtrema < 2)
        {
            throw new ArgumentOutOfRangeException(nameof(minMatchedExtrema), "minMatchedExtrema must be at least 2.");
        }

        Dictionary<SpectraFileInfo, Func<double, double>> warps = new(files.Count);
        SpectraFileInfo referenceFile = files[referenceFileIndex];

        for (int fileIndex = 0; fileIndex < files.Count; fileIndex++)
        {
            SpectraFileInfo file = files[fileIndex];
            if (fileIndex == referenceFileIndex)
            {
                warps[file] = static rt => rt;
                continue;
            }

            List<(double SourceRt, double TargetRt)> matchedPairs = MatchAnchorExtrema(
                extremaByFile,
                referenceFileIndex,
                fileIndex,
                extremaMatchTolerance);

            warps[file] = BuildMonotoneWarp(matchedPairs, minMatchedExtrema, useMonotoneSpline);
        }

        warps[referenceFile] = static rt => rt;
        return warps;
    }

    private static List<(double SourceRt, double TargetRt)> MatchAnchorExtrema(
        IReadOnlyList<AnchorFileExtrema> extremaByFile,
        int referenceFileIndex,
        int targetFileIndex,
        double matchTolerance)
    {
        Dictionary<AnchorBin, (List<Extremum> Reference, List<Extremum> Target)> byAnchor = new();

        foreach (AnchorFileExtrema item in extremaByFile)
        {
            if (!byAnchor.TryGetValue(item.Anchor, out (List<Extremum> Reference, List<Extremum> Target) entry))
            {
                entry = (new List<Extremum>(), new List<Extremum>());
                byAnchor[item.Anchor] = entry;
            }

            if (item.FileIndex == referenceFileIndex)
            {
                entry.Reference.AddRange(item.Extrema);
            }
            else if (item.FileIndex == targetFileIndex)
            {
                entry.Target.AddRange(item.Extrema);
            }
        }

        List<(double SourceRt, double TargetRt)> matches = new();
        foreach ((List<Extremum> referenceExtrema, List<Extremum> targetExtrema) in byAnchor.Values)
        {
            AddMatchesForType(referenceExtrema, targetExtrema, ExtremumType.Maximum, matchTolerance, matches);
            AddMatchesForType(referenceExtrema, targetExtrema, ExtremumType.Minimum, matchTolerance, matches);
        }

        return matches;
    }

    private static void AddMatchesForType(
        IReadOnlyList<Extremum> referenceExtrema,
        IReadOnlyList<Extremum> targetExtrema,
        ExtremumType type,
        double matchTolerance,
        List<(double SourceRt, double TargetRt)> output)
    {
        List<Extremum> reference = referenceExtrema
            .Where(p => p.Type == type)
            .OrderBy(p => p.RetentionTime)
            .ToList();

        List<Extremum> target = targetExtrema
            .Where(p => p.Type == type)
            .OrderBy(p => p.RetentionTime)
            .ToList();

        if (reference.Count == 0 || target.Count == 0)
        {
            return;
        }

        int referenceIndex = 0;
        for (int targetIndex = 0; targetIndex < target.Count; targetIndex++)
        {
            double targetRt = target[targetIndex].RetentionTime;

            while (referenceIndex + 1 < reference.Count
                   && Math.Abs(reference[referenceIndex + 1].RetentionTime - targetRt)
                   <= Math.Abs(reference[referenceIndex].RetentionTime - targetRt))
            {
                referenceIndex++;
            }

            double referenceRt = reference[referenceIndex].RetentionTime;
            if (Math.Abs(referenceRt - targetRt) <= matchTolerance)
            {
                output.Add((targetRt, referenceRt));
            }
        }
    }

    private static Func<double, double> BuildMonotoneWarp(
        List<(double SourceRt, double TargetRt)> matchedPairs,
        int minMatchedExtrema,
        bool useMonotoneSpline)
    {
        if (matchedPairs.Count < minMatchedExtrema)
        {
            return static rt => rt;
        }

        var knotPairs = matchedPairs
            .OrderBy(p => p.SourceRt)
            .GroupBy(p => p.SourceRt)
            .Select(g => (SourceRt: g.Key, TargetRt: g.Average(p => p.TargetRt)))
            .OrderBy(p => p.SourceRt)
            .ToList();

        if (knotPairs.Count < minMatchedExtrema)
        {
            return static rt => rt;
        }

        double[] source = knotPairs.Select(p => p.SourceRt).ToArray();
        double[] target = EnforceMonotoneIncreasing(knotPairs.Select(p => p.TargetRt).ToArray());

        if (useMonotoneSpline)
        {
            return BuildMonotoneSpline(source, target);
        }

        return BuildPiecewiseLinear(source, target);
    }

    private static double[] EnforceMonotoneIncreasing(double[] values)
    {
        List<(double Mean, int Count)> blocks = new();
        for (int i = 0; i < values.Length; i++)
        {
            blocks.Add((values[i], 1));
            while (blocks.Count >= 2)
            {
                var right = blocks[^1];
                var left = blocks[^2];
                if (left.Mean <= right.Mean)
                {
                    break;
                }

                int mergedCount = left.Count + right.Count;
                double mergedMean = (left.Mean * left.Count + right.Mean * right.Count) / mergedCount;
                blocks.RemoveAt(blocks.Count - 1);
                blocks[^1] = (mergedMean, mergedCount);
            }
        }

        double[] monotone = new double[values.Length];
        int index = 0;
        foreach (var block in blocks)
        {
            for (int i = 0; i < block.Count; i++)
            {
                monotone[index++] = block.Mean;
            }
        }

        return monotone;
    }

    private static Func<double, double> BuildPiecewiseLinear(double[] x, double[] y)
    {
        int n = x.Length;
        return rt =>
        {
            if (rt <= x[0])
            {
                return Extrapolate(rt, x[0], x[1], y[0], y[1]);
            }

            if (rt >= x[n - 1])
            {
                return Extrapolate(rt, x[n - 2], x[n - 1], y[n - 2], y[n - 1]);
            }

            int right = Array.BinarySearch(x, rt);
            if (right >= 0)
            {
                return y[right];
            }

            right = ~right;
            int left = right - 1;
            return Interpolate(rt, x[left], x[right], y[left], y[right]);
        };
    }

    private static Func<double, double> BuildMonotoneSpline(double[] x, double[] y)
    {
        int n = x.Length;
        double[] h = new double[n - 1];
        double[] delta = new double[n - 1];

        for (int i = 0; i < n - 1; i++)
        {
            h[i] = x[i + 1] - x[i];
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        double[] m = new double[n];
        m[0] = delta[0];
        m[n - 1] = delta[n - 2];

        for (int i = 1; i < n - 1; i++)
        {
            if (delta[i - 1] == 0 || delta[i] == 0 || Math.Sign(delta[i - 1]) != Math.Sign(delta[i]))
            {
                m[i] = 0;
                continue;
            }

            double w1 = 2 * h[i] + h[i - 1];
            double w2 = h[i] + 2 * h[i - 1];
            m[i] = (w1 + w2) / ((w1 / delta[i - 1]) + (w2 / delta[i]));
        }

        return rt =>
        {
            if (rt <= x[0])
            {
                return Extrapolate(rt, x[0], x[1], y[0], y[1]);
            }

            if (rt >= x[n - 1])
            {
                return Extrapolate(rt, x[n - 2], x[n - 1], y[n - 2], y[n - 1]);
            }

            int right = Array.BinarySearch(x, rt);
            if (right >= 0)
            {
                return y[right];
            }

            right = ~right;
            int left = right - 1;
            double interval = x[right] - x[left];
            double t = (rt - x[left]) / interval;
            double t2 = t * t;
            double t3 = t2 * t;

            double h00 = 2 * t3 - 3 * t2 + 1;
            double h10 = t3 - 2 * t2 + t;
            double h01 = -2 * t3 + 3 * t2;
            double h11 = t3 - t2;

            return h00 * y[left]
                + h10 * interval * m[left]
                + h01 * y[right]
                + h11 * interval * m[right];
        };
    }

    private static double Interpolate(double x, double x0, double x1, double y0, double y1)
    {
        return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    }

    private static double Extrapolate(double x, double x0, double x1, double y0, double y1)
    {
        return Interpolate(x, x0, x1, y0, y1);
    }
}
