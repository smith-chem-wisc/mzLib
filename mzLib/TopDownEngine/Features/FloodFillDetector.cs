using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace TopDownEngine.Features;

public enum MzGrowthStrategy
{
    SeededBinWalk,
    ConcentricRing
}

public sealed class FloodFillOptions
{
    public Tolerance PeakFindingTolerance { get; init; } = new AbsoluteTolerance(0.005);
    public int MaxMissedScansAllowed { get; init; } = 1;
    public double MaxRtRange { get; init; } = 0.6;
    public int NumPeakThreshold { get; init; } = 4;
    public double? CutPeakDiscriminationFactor { get; init; } = 0.6;
    public MzGrowthStrategy GrowthStrategy { get; init; } = MzGrowthStrategy.SeededBinWalk;
    public int MaxMzGrowthSteps { get; init; } = 12;
    public double MinRtOverlapFraction { get; init; } = 0.25;
    public int MinAdjacentXicPeakCount { get; init; } = 3;
    public double MaxFrontierIntensityIncreaseFraction { get; init; } = 0.1;
    public double? NoiseFloorIntensity { get; init; }
    public double NoiseFloorPercentile { get; init; } = 0.1;
}

public sealed class FloodFillDetector
{
    private const double FineBinWidth = 0.01;

    private static readonly FieldInfo IndexedPeaksField =
        typeof(IndexingEngine<IndexedMassSpectralPeak>).GetField("IndexedPeaks", BindingFlags.Instance | BindingFlags.NonPublic)
        ?? throw new MzLibException("Unable to access indexed peaks from PeakIndexingEngine.");

    public IReadOnlyList<FeatureBox> Detect(
        PeakIndexingEngine fineIndex,
        FloodFillOptions options,
        string sourceFile = "")
    {
        if (fineIndex == null)
        {
            throw new ArgumentNullException(nameof(fineIndex));
        }

        if (options == null)
        {
            throw new ArgumentNullException(nameof(options));
        }

        double noiseFloor = options.NoiseFloorIntensity
            ?? ComputeNoiseFloorFromPercentile(fineIndex, options.NoiseFloorPercentile);

        List<ExtractedIonChromatogram> xics = fineIndex.GetAllXics(
            options.PeakFindingTolerance,
            options.MaxMissedScansAllowed,
            options.MaxRtRange,
            options.NumPeakThreshold,
            options.CutPeakDiscriminationFactor);

        List<FeatureBox> boxes = new();
        foreach (ExtractedIonChromatogram xic in xics.Where(x => x?.ApexPeak != null)
            .OrderByDescending(x => x.ApexPeak.Intensity))
        {
            if (xic.ApexPeak.Intensity < noiseFloor)
            {
                continue;
            }

            HashSet<IIndexedPeak> grownPeaks = options.GrowthStrategy == MzGrowthStrategy.SeededBinWalk
                ? GrowBySeededBinWalk(fineIndex, xic, options, noiseFloor)
                : GrowByConcentricRing(fineIndex, xic, options, noiseFloor);

            if (grownPeaks.Count < options.NumPeakThreshold)
            {
                continue;
            }

            double minMz = grownPeaks.Min(p => p.M);
            double maxMz = grownPeaks.Max(p => p.M);
            double minRt = grownPeaks.Min(p => p.RetentionTime);
            double maxRt = grownPeaks.Max(p => p.RetentionTime);

            boxes.Add(new FeatureBox(
                new MzRange(minMz, maxMz),
                new DoubleRange(minRt, maxRt),
                xic.ApexPeak.Intensity,
                grownPeaks.Sum(p => p.Intensity),
                sourceFile,
                grownPeaks.Count));
        }

        return boxes;
    }

    private static HashSet<IIndexedPeak> GrowBySeededBinWalk(
        PeakIndexingEngine fineIndex,
        ExtractedIonChromatogram seedXic,
        FloodFillOptions options,
        double noiseFloor)
    {
        HashSet<IIndexedPeak> grownPeaks = new(seedXic.Peaks);
        DoubleRange seedRtRange = new(seedXic.StartRT, seedXic.EndRT);
        double apexRt = seedXic.ApexRT;

        foreach (int direction in new[] { -1, 1 })
        {
            double frontierIntensity = seedXic.ApexPeak.Intensity;

            for (int step = 1; step <= options.MaxMzGrowthSteps; step++)
            {
                double mz = seedXic.ApexPeak.M + (direction * step * FineBinWidth);
                List<IIndexedPeak> adjacentXic = fineIndex.GetXic(
                    mz,
                    apexRt,
                    options.PeakFindingTolerance,
                    options.MaxMissedScansAllowed,
                    options.MaxRtRange);

                if (adjacentXic.Count < options.MinAdjacentXicPeakCount)
                {
                    break;
                }

                DoubleRange adjacentRtRange = new(adjacentXic.Min(p => p.RetentionTime), adjacentXic.Max(p => p.RetentionTime));
                if (!HasSufficientRtOverlap(seedRtRange, adjacentRtRange, options.MinRtOverlapFraction))
                {
                    break;
                }

                double adjacentMaxIntensity = adjacentXic.Max(p => p.Intensity);
                if (adjacentMaxIntensity < noiseFloor ||
                    adjacentMaxIntensity > frontierIntensity * (1.0 + options.MaxFrontierIntensityIncreaseFraction))
                {
                    break;
                }

                foreach (IIndexedPeak peak in adjacentXic)
                {
                    grownPeaks.Add(peak);
                }

                frontierIntensity = adjacentMaxIntensity;
            }
        }

        return grownPeaks;
    }

    private static HashSet<IIndexedPeak> GrowByConcentricRing(
        PeakIndexingEngine fineIndex,
        ExtractedIonChromatogram seedXic,
        FloodFillOptions options,
        double noiseFloor)
    {
        HashSet<IIndexedPeak> grownPeaks = new(seedXic.Peaks);
        DoubleRange seedRtRange = new(seedXic.StartRT, seedXic.EndRT);
        double apexRt = seedXic.ApexRT;
        double frontierIntensity = seedXic.ApexPeak.Intensity;

        for (int step = 1; step <= options.MaxMzGrowthSteps; step++)
        {
            double tolerance = step * FineBinWidth;
            List<IIndexedPeak> ringXic = fineIndex.GetXic(
                seedXic.ApexPeak.M,
                apexRt,
                new AbsoluteTolerance(tolerance),
                options.MaxMissedScansAllowed,
                options.MaxRtRange);

            if (ringXic.Count < options.NumPeakThreshold)
            {
                continue;
            }

            List<IIndexedPeak> newPeaks = ringXic.Where(p => !grownPeaks.Contains(p)).ToList();
            if (newPeaks.Count == 0)
            {
                continue;
            }

            DoubleRange ringRtRange = new(newPeaks.Min(p => p.RetentionTime), newPeaks.Max(p => p.RetentionTime));
            if (!HasSufficientRtOverlap(seedRtRange, ringRtRange, options.MinRtOverlapFraction))
            {
                break;
            }

            double ringFrontierIntensity = newPeaks.Max(p => p.Intensity);
            if (ringFrontierIntensity < noiseFloor ||
                ringFrontierIntensity > frontierIntensity * (1.0 + options.MaxFrontierIntensityIncreaseFraction))
            {
                break;
            }

            foreach (IIndexedPeak peak in newPeaks)
            {
                grownPeaks.Add(peak);
            }

            frontierIntensity = ringFrontierIntensity;
        }

        return grownPeaks;
    }

    private static bool HasSufficientRtOverlap(DoubleRange a, DoubleRange b, double minOverlapFraction)
    {
        double overlapMin = Math.Max(a.Minimum, b.Minimum);
        double overlapMax = Math.Min(a.Maximum, b.Maximum);
        double overlap = Math.Max(0, overlapMax - overlapMin);
        if (overlap == 0)
        {
            return false;
        }

        double denominator = Math.Max(Math.Min(a.Width, b.Width), 1e-12);
        return overlap / denominator >= minOverlapFraction;
    }

    private static double ComputeNoiseFloorFromPercentile(PeakIndexingEngine fineIndex, double percentile)
    {
        percentile = Math.Clamp(percentile, 0.0, 1.0);

        List<IndexedMassSpectralPeak>[] indexedPeaks = IndexedPeaksField.GetValue(fineIndex) as List<IndexedMassSpectralPeak>[]
            ?? throw new MzLibException("Peak index is not initialized.");

        List<float> intensities = indexedPeaks
            .Where(bin => bin != null)
            .SelectMany(bin => bin)
            .Select(peak => peak.Intensity)
            .OrderBy(i => i)
            .ToList();

        if (intensities.Count == 0)
        {
            return 0;
        }

        int index = (int)Math.Floor((intensities.Count - 1) * percentile);
        return intensities[index];
    }
}
