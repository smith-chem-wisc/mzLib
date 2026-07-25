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
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

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
    public void SparseBasisMatchesTheExactBasisToFarBetterThanFitPrecision()
    {
        // The sparsity threshold drops basis terms below a fraction of the sample's own-model
        // basis. That is the one change in this class that can move fitted numbers, so measure how
        // far it actually moves them rather than assuming it is negligible.
        const int minCharge = 6;
        const int maxCharge = 11;
        const double sigmaMz = 0.012;

        var trueModels = new[]
        {
            new ProteoformModel(10000.0, 1.4e6, new EmgProfile(20.0, 0.16, 0.06),
                new GaussianChargeDistribution(8.2, 1.1), "pf1"),
            new ProteoformModel(10000.2, 8.5e5, new EmgProfile(20.03, 0.15, 0.06),
                new GaussianChargeDistribution(8.4, 1.0), "pf2"),
            new ProteoformModel(10001.1, 3.1e5, new EmgProfile(20.06, 0.17, 0.05),
                new GaussianChargeDistribution(8.0, 1.2), "pf3"),
        };

        var scanArray = BuildCentroidedScansFromForwardModel(trueModels, minCharge, maxCharge, sigmaMz);
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

        var truths = trueModels
            .Select(m => extractor.Extract(m.MonoisotopicMass, m.RtProfile.Mu, 0.8, minCharge, maxCharge))
            .ToArray();

        var initialFits = trueModels
            .Select((m, i) => new FittedProteoform(
                m with { Abundance = m.Abundance * (i % 2 == 0 ? 1.8 : 0.3) },
                sigmaMz, WidthFitMode.CentroidedFallback, 0, 0))
            .ToArray();

        GlobalAbundanceRefitResult RefitWith(double threshold) =>
            new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(
                    MaxIterations: 20,
                    ConvergenceTolerance: 1e-5,
                    BasisSparsityThreshold: threshold))
                .Refit(initialFits, truths, minCharge, maxCharge, sigmaMz);

        var exact = RefitWith(0);            // keep every positive basis term
        var sparse = RefitWith(1e-10);       // the shipped default

        Assert.That(sparse.FittedProteoforms, Has.Length.EqualTo(exact.FittedProteoforms.Length));

        for (int i = 0; i < exact.FittedProteoforms.Length; i++)
        {
            double exactAbundance = exact.FittedProteoforms[i].Model.Abundance;
            double sparseAbundance = sparse.FittedProteoforms[i].Model.Abundance;
            double relative = Math.Abs(sparseAbundance - exactAbundance) / exactAbundance;

            Assert.That(relative, Is.LessThan(1e-6),
                $"Proteoform {i} abundance moved by {relative:E3} relative, which is large enough to matter.");
        }

        Assert.That(sparse.FinalResidualFraction,
            Is.EqualTo(exact.FinalResidualFraction).Within(1e-6).Percent);
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
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

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
