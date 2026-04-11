using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace TopDownEngine.Features;

public sealed record RtCluster(
    DoubleRange RtRange,
    IReadOnlyList<FeatureGroup> FeatureGroups);

public sealed class RtClusterer
{
    public List<RtCluster> ClusterByRtProximity(
        IReadOnlyList<FeatureGroup> featureGroups,
        double rtProximityWindow)
    {
        if (featureGroups == null)
        {
            throw new ArgumentNullException(nameof(featureGroups));
        }

        if (double.IsNaN(rtProximityWindow) || double.IsInfinity(rtProximityWindow) || rtProximityWindow <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(rtProximityWindow), "rtProximityWindow must be positive and finite.");
        }

        if (featureGroups.Count == 0)
        {
            return new List<RtCluster>();
        }

        List<FeatureGroup> orderedGroups = featureGroups
            .OrderBy(GetRtCenter)
            .ThenBy(f => f.DonorBox.RtRange.Minimum)
            .ThenBy(f => f.DonorBox.RtRange.Maximum)
            .ThenBy(f => f.DonorBox.SourceFile, StringComparer.Ordinal)
            .ThenBy(f => f.DonorBox.MzRange.Minimum)
            .ToList();

        List<RtCluster> clusters = new();
        List<FeatureGroup> currentClusterMembers = new() { orderedGroups[0] };

        double clusterRtMin = orderedGroups[0].DonorBox.RtRange.Minimum;
        double clusterRtMax = orderedGroups[0].DonorBox.RtRange.Maximum;
        double previousCenter = GetRtCenter(orderedGroups[0]);

        for (int i = 1; i < orderedGroups.Count; i++)
        {
            FeatureGroup candidate = orderedGroups[i];
            double candidateCenter = GetRtCenter(candidate);

            if (candidateCenter - previousCenter <= rtProximityWindow)
            {
                currentClusterMembers.Add(candidate);
                clusterRtMin = Math.Min(clusterRtMin, candidate.DonorBox.RtRange.Minimum);
                clusterRtMax = Math.Max(clusterRtMax, candidate.DonorBox.RtRange.Maximum);
                previousCenter = candidateCenter;
                continue;
            }

            clusters.Add(new RtCluster(new DoubleRange(clusterRtMin, clusterRtMax), currentClusterMembers.ToArray()));

            currentClusterMembers = new List<FeatureGroup> { candidate };
            clusterRtMin = candidate.DonorBox.RtRange.Minimum;
            clusterRtMax = candidate.DonorBox.RtRange.Maximum;
            previousCenter = candidateCenter;
        }

        clusters.Add(new RtCluster(new DoubleRange(clusterRtMin, clusterRtMax), currentClusterMembers.ToArray()));
        return clusters;
    }

    private static double GetRtCenter(FeatureGroup group)
    {
        return 0.5 * (group.DonorBox.RtRange.Minimum + group.DonorBox.RtRange.Maximum);
    }
}
