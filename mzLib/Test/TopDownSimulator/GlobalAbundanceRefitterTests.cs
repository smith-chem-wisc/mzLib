using System;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class GlobalAbundanceRefitterTests
{
    private static MsDataScan BuildMs1Scan(int oneBased, double rt, double[] mz, double[] intensity) =>
        new MsDataScan(
            massSpectrum: new MzSpectrum(mz, intensity, false),
            oneBasedScanNumber: oneBased,
            msnOrder: 1,
            isCentroid: true,
            polarity: Polarity.Positive,
            retentionTime: rt,
            scanWindowRange: new MzRange(100, 3000),
            scanFilter: "f",
            mzAnalyzer: MZAnalyzerType.Orbitrap,
            totalIonCurrent: intensity.Sum(),
            injectionTime: 1.0,
            noiseData: null,
            nativeId: "scan=" + oneBased);

    [Test]
    public void RefitImprovesResidualAndAbundanceErrorForOverlappingProteoforms()
    {
        const int minCharge = 6;
        const int maxCharge = 11;
        const double sigmaMz = 0.012;

        var trueModels = new[]
        {
            new ProteoformModel(
                MonoisotopicMass: 10000.0,
                Abundance: 1.4e6,
                RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.16, Tau: 0.06),
                ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.2, SigmaZ: 1.1),
                Identifier: "pf1"),
            new ProteoformModel(
                MonoisotopicMass: 10000.2,
                Abundance: 8.5e5,
                RtProfile: new EmgProfile(Mu: 20.03, Sigma: 0.15, Tau: 0.06),
                ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.4, SigmaZ: 1.0),
                Identifier: "pf2"),
        };

        var scanArray = BuildCentroidedScansFromForwardModel(trueModels, minCharge, maxCharge, sigmaMz);
        var index = PeakIndexingEngine.InitializeIndexingEngine(scanArray)!;
        var extractor = new GroundTruthExtractor(index, scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

        var truths = new[]
        {
            extractor.Extract(trueModels[0].MonoisotopicMass, rtCenter: 20.0, rtHalfWidth: 0.8, minCharge, maxCharge),
            extractor.Extract(trueModels[1].MonoisotopicMass, rtCenter: 20.03, rtHalfWidth: 0.8, minCharge, maxCharge),
        };

        var initialFits = new[]
        {
            new FittedProteoform(trueModels[0] with { Abundance = trueModels[0].Abundance * 1.8 }, sigmaMz, WidthFitMode.CentroidedFallback, 0, 0),
            new FittedProteoform(trueModels[1] with { Abundance = trueModels[1].Abundance * 0.3 }, sigmaMz, WidthFitMode.CentroidedFallback, 0, 0),
        };

        var refitter = new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(MaxIterations: 20, ConvergenceTolerance: 1e-5));
        var result = refitter.Refit(initialFits, truths, minCharge, maxCharge, sigmaMz);

        double initialError = RelativeErrorSum(initialFits.Select(f => f.Model.Abundance).ToArray(), trueModels.Select(m => m.Abundance).ToArray());
        double finalError = RelativeErrorSum(result.FittedProteoforms.Select(f => f.Model.Abundance).ToArray(), trueModels.Select(m => m.Abundance).ToArray());

        Assert.That(result.FinalResidualFraction, Is.LessThan(result.InitialResidualFraction));
        Assert.That(finalError, Is.LessThan(initialError));
    }

    [Test]
    public void RefitKeepsAbundanceNonNegative()
    {
        const int minCharge = 6;
        const int maxCharge = 11;
        const double sigmaMz = 0.012;

        var trueModels = new[]
        {
            new ProteoformModel(
                MonoisotopicMass: 12000.0,
                Abundance: 9e5,
                RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.20, Tau: 0.07),
                ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.0, SigmaZ: 1.2),
                Identifier: "pf1"),
            new ProteoformModel(
                MonoisotopicMass: 12000.25,
                Abundance: 7e5,
                RtProfile: new EmgProfile(Mu: 20.05, Sigma: 0.20, Tau: 0.07),
                ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.5, SigmaZ: 1.2),
                Identifier: "pf2"),
        };

        var scanArray = BuildCentroidedScansFromForwardModel(trueModels, minCharge, maxCharge, sigmaMz);
        var index = PeakIndexingEngine.InitializeIndexingEngine(scanArray)!;
        var extractor = new GroundTruthExtractor(index, scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

        var truths = new[]
        {
            extractor.Extract(trueModels[0].MonoisotopicMass, rtCenter: 20.0, rtHalfWidth: 0.8, minCharge, maxCharge),
            extractor.Extract(trueModels[1].MonoisotopicMass, rtCenter: 20.05, rtHalfWidth: 0.8, minCharge, maxCharge),
        };

        var initialFits = new[]
        {
            new FittedProteoform(trueModels[0] with { Abundance = 0 }, sigmaMz, WidthFitMode.CentroidedFallback, 0, 0),
            new FittedProteoform(trueModels[1] with { Abundance = trueModels[1].Abundance * 2.5 }, sigmaMz, WidthFitMode.CentroidedFallback, 0, 0),
        };

        var refitter = new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(MaxIterations: 10, ConvergenceTolerance: 1e-4));
        var result = refitter.Refit(initialFits, truths, minCharge, maxCharge, sigmaMz);

        Assert.That(result.FittedProteoforms.All(f => f.Model.Abundance >= 0), Is.True);
    }

    private static MsDataScan[] BuildCentroidedScansFromForwardModel(
        IReadOnlyList<ProteoformModel> trueModels,
        int minCharge,
        int maxCharge,
        double sigmaMz)
    {
        var fm = new ForwardModel(trueModels, minCharge, maxCharge, sigmaMz);
        var kernels = trueModels.Select(m => new IsotopeEnvelopeKernel(m.MonoisotopicMass)).ToArray();

        var scans = new List<MsDataScan>();
        int oneBased = 1;
        for (double rt = 19.2; rt <= 20.8; rt += 0.05)
        {
            var mzList = new List<double>();
            var intList = new List<double>();

            for (int p = 0; p < trueModels.Count; p++)
            {
                for (int z = minCharge; z <= maxCharge; z++)
                {
                    double[] centroids = kernels[p].CentroidMzs(z);
                    for (int i = 0; i < centroids.Length; i++)
                    {
                        mzList.Add(centroids[i]);
                        intList.Add(fm.Evaluate(rt, centroids[i]));
                    }
                }
            }

            var ordered = mzList.Select((m, idx) => (mz: m, intensity: intList[idx])).OrderBy(t => t.mz).ToList();
            scans.Add(BuildMs1Scan(oneBased++, rt, ordered.Select(t => t.mz).ToArray(), ordered.Select(t => t.intensity).ToArray()));
        }

        return scans.ToArray();
    }

    private static double RelativeErrorSum(IReadOnlyList<double> estimate, IReadOnlyList<double> truth)
    {
        double sum = 0;
        for (int i = 0; i < estimate.Count; i++)
        {
            sum += Math.Abs(estimate[i] - truth[i]) / truth[i];
        }

        return sum;
    }
}
