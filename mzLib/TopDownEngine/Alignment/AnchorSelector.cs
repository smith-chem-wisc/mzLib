using System;
using System.Collections.Generic;
using System.Linq;
using TopDownEngine.Indexing;

namespace TopDownEngine.Alignment;

public class AnchorSelector
{
    private const int ThickBinsPerDalton = 10;

    public AnchorBin[] SelectTopK(IReadOnlyList<ThickIndexView> thickIndexes, int topK, int minFilesForAnchor, double intensityThreshold)
    {
        if (thickIndexes == null)
        {
            throw new ArgumentNullException(nameof(thickIndexes));
        }

        if (topK <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(topK), "topK must be positive.");
        }

        if (minFilesForAnchor <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(minFilesForAnchor), "minFilesForAnchor must be positive.");
        }

        if (intensityThreshold < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(intensityThreshold), "intensityThreshold must be non-negative.");
        }

        var statsByBin = new Dictionary<int, BinStats>();

        foreach (var thickIndex in thickIndexes)
        {
            if (thickIndex == null)
            {
                continue;
            }

            Dictionary<int, float> maximaByBin = thickIndex.GetBinMaximumIntensities((float)intensityThreshold);
            foreach (var kvp in maximaByBin)
            {
                if (!statsByBin.TryGetValue(kvp.Key, out var stats))
                {
                    stats = new BinStats();
                    statsByBin[kvp.Key] = stats;
                }

                stats.FileCount++;
                stats.IntensitySum += kvp.Value;
            }
        }

        return statsByBin
            .Where(kvp => kvp.Value.FileCount >= minFilesForAnchor)
            .OrderByDescending(kvp => kvp.Value.FileCount)
            .ThenByDescending(kvp => kvp.Value.IntensitySum / kvp.Value.FileCount)
            .ThenBy(kvp => kvp.Key)
            .Take(topK)
            .Select(kvp => new AnchorBin(
                kvp.Key,
                BinIndexToMzCenter(kvp.Key),
                kvp.Value.FileCount,
                kvp.Value.IntensitySum / kvp.Value.FileCount))
            .ToArray();
    }

    private static double BinIndexToMzCenter(int thickBinIndex)
    {
        return (thickBinIndex + 0.5) / ThickBinsPerDalton;
    }

    private sealed class BinStats
    {
        public int FileCount { get; set; }
        public double IntensitySum { get; set; }
    }
}
