using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using TopDownEngine.Features;

namespace Test.TopDownEngine.Features;

[TestFixture]
public class RtClustererTests
{
    [Test]
    public void ClusterByRtProximity_ThreeFeaturesWithOneGap_ProducesTwoClusters()
    {
        FeatureGroup featureA = BuildFeatureGroup(centerRt: 10.0, sourceFile: "a.raw");
        FeatureGroup featureB = BuildFeatureGroup(centerRt: 10.2, sourceFile: "b.raw");
        FeatureGroup featureC = BuildFeatureGroup(centerRt: 15.0, sourceFile: "c.raw");

        RtClusterer clusterer = new();
        List<RtCluster> clusters = clusterer.ClusterByRtProximity(
            new[] { featureC, featureA, featureB },
            rtProximityWindow: 0.5);

        Assert.That(clusters, Has.Count.EqualTo(2));
        Assert.That(clusters[0].FeatureGroups, Has.Count.EqualTo(2));
        Assert.That(clusters[1].FeatureGroups, Has.Count.EqualTo(1));
        Assert.That(clusters[0].FeatureGroups, Contains.Item(featureA));
        Assert.That(clusters[0].FeatureGroups, Contains.Item(featureB));
        Assert.That(clusters[1].FeatureGroups[0], Is.EqualTo(featureC));
    }

    private static FeatureGroup BuildFeatureGroup(double centerRt, string sourceFile)
    {
        DoubleRange rtRange = new(centerRt - 0.05, centerRt + 0.05);
        FeatureBox box = new(
            MzRange: new MzRange(500.0, 500.1),
            RtRange: rtRange,
            SeedIntensity: 1000,
            TotalIntensity: 5000,
            SourceFile: sourceFile,
            PeakCount: 10);

        return new FeatureGroup(
            DonorBox: box,
            Matches: new Dictionary<SpectraFileInfo, (int KTrue, double PValue, double Intensity)>(),
            TotalSupport: 1);
    }
}
