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
    public void ResidualFractionFallsMonotonicallyAndConvergesQuickly()
    {
        // Four heavily overlapping proteoforms: same mass to within ~1 Da, same elution window,
        // same charge envelope. This is where a Jacobi sweep — every coordinate reasoning about
        // totals frozen at the start of the sweep — costs the most, because each update is made
        // in ignorance of the three updates beside it.
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
            new ProteoformModel(10000.6, 9.4e5, new EmgProfile(20.01, 0.16, 0.06),
                new GaussianChargeDistribution(8.3, 1.1), "pf4"),
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

        var result = new GlobalAbundanceRefitter(
                new GlobalAbundanceRefitOptions(MaxIterations: 60, ConvergenceTolerance: 1e-5))
            .Refit(initialFits, truths, minCharge, maxCharge, sigmaMz);

        // Monotone descent is the property Gauss-Seidel buys and Jacobi does not guarantee, so it
        // is the sharpest available check that each coordinate really did see its predecessors.
        Assert.That(result.ResidualFractionByIteration, Has.Count.EqualTo(result.IterationsCompleted));
        Assert.That(result.ResidualFractionByIteration[0], Is.LessThan(result.InitialResidualFraction));
        for (int i = 1; i < result.ResidualFractionByIteration.Count; i++)
        {
            Assert.That(result.ResidualFractionByIteration[i],
                Is.LessThanOrEqualTo(result.ResidualFractionByIteration[i - 1]),
                $"Residual energy rose at iteration {i + 1}, which a Gauss-Seidel sweep cannot do.");
        }

        // Measured on this fixture: Jacobi ran all 60 sweeps without reaching 1e-5; Gauss-Seidel
        // converges in 16. The bound is loose enough to absorb platform arithmetic differences and
        // still far below what the stale-totals sweep could reach.
        Assert.That(result.Converged, Is.True);
        Assert.That(result.IterationsCompleted, Is.LessThanOrEqualTo(30));

        for (int i = 0; i < trueModels.Length; i++)
        {
            double relative = Math.Abs(result.FittedProteoforms[i].Model.Abundance - trueModels[i].Abundance)
                              / trueModels[i].Abundance;
            Assert.That(relative, Is.LessThan(1e-4),
                $"Proteoform {i} did not recover its true abundance.");
        }
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

    [Test]
    public void DeduplicatingIsANoOpWhenNoProteoformsOverlap()
    {
        // With no shared peaks there is nothing to double-count, so the deduplicated metric must
        // equal the per-sample one. That is the sharpest available check that the dedup key
        // identifies real peaks rather than merging things it should not.
        const int minCharge = 6;
        const int maxCharge = 11;
        const double sigmaMz = 0.012;

        var trueModels = new[]
        {
            new ProteoformModel(10000.0, 1.4e6, new EmgProfile(20.0, 0.16, 0.06),
                new GaussianChargeDistribution(8.2, 1.1), "pf1"),
            // 3 kDa apart and more than a minute away: no isotopologue and no scan is shared.
            new ProteoformModel(13000.0, 8.5e5, new EmgProfile(21.4, 0.15, 0.06),
                new GaussianChargeDistribution(8.4, 1.0), "pf2"),
        };

        var scanArray = BuildCentroidedScansFromForwardModel(trueModels, minCharge, maxCharge, sigmaMz);
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truths = trueModels
            .Select(m => extractor.Extract(m.MonoisotopicMass, m.RtProfile.Mu, 0.3, minCharge, maxCharge))
            .ToArray();

        // Abundances deliberately off, so the residual is a real quantity rather than float noise.
        var badFits = trueModels
            .Select(m => new FittedProteoform(m with { Abundance = m.Abundance * 1.4 },
                sigmaMz, WidthFitMode.CentroidedFallback, 0, 0))
            .ToArray();

        var deduplicated = Refit(badFits, truths, minCharge, maxCharge, sigmaMz);
        var perSample = Refit(badFits, truths.Select(StripSourcePeakIdentity).ToArray(), minCharge, maxCharge, sigmaMz);

        Assert.That(deduplicated.InitialResiduals!.DistinctPeaks,
            Is.EqualTo(perSample.InitialResiduals!.DistinctPeaks),
            "No sample should have been merged into another when nothing overlaps.");
        Assert.That(deduplicated.InitialResiduals.UnexplainedFraction,
            Is.EqualTo(perSample.InitialResiduals.UnexplainedFraction).Within(1e-9).Percent);
    }

    [Test]
    public void DeduplicatingChangesTheMetricWhenProteoformsOverlap()
    {
        // A crowded region and a sparse one, fitted to different quality. The crowded proteoform is
        // reported twice with masses 0.02 Da apart — the shape of two search results for the same
        // signal — so both claimants match the same experimental peaks. Counting those peaks once
        // per claimant doubles the weight of the crowded region, and because it happens to be the
        // well-explained one, the reported fraction comes out lower than the truth. That weighting
        // depends on how crowded the data is, which is why the number is not comparable between
        // runs until the double-counting is removed.
        const int minCharge = 6;
        const int maxCharge = 11;
        const double sigmaMz = 0.012;

        var crowded = new ProteoformModel(10000.0, 1.4e6, new EmgProfile(20.0, 0.16, 0.06),
            new GaussianChargeDistribution(8.2, 1.1), "crowded");
        var isolated = new ProteoformModel(13000.0, 1.2e6, new EmgProfile(20.0, 0.16, 0.06),
            new GaussianChargeDistribution(8.2, 1.1), "isolated");

        var scanArray = BuildCentroidedScansFromForwardModel(new[] { crowded, isolated }, minCharge, maxCharge, sigmaMz);
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

        var claimants = new[]
        {
            crowded with { MonoisotopicMass = crowded.MonoisotopicMass - 0.02, Abundance = crowded.Abundance / 2, Identifier = "claim-low" },
            crowded with { MonoisotopicMass = crowded.MonoisotopicMass + 0.02, Abundance = crowded.Abundance / 2, Identifier = "claim-high" },
            // Badly wrong on purpose, so the sparse region carries most of the real residual.
            isolated with { Abundance = isolated.Abundance * 0.3 },
        };

        var truths = claimants
            .Select(m => extractor.Extract(m.MonoisotopicMass, m.RtProfile.Mu, 0.8, minCharge, maxCharge))
            .ToArray();
        var fits = claimants
            .Select(m => new FittedProteoform(m, sigmaMz, WidthFitMode.CentroidedFallback, 0, 0))
            .ToArray();

        var deduplicated = Refit(fits, truths, minCharge, maxCharge, sigmaMz);
        var perSample = Refit(fits, truths.Select(StripSourcePeakIdentity).ToArray(), minCharge, maxCharge, sigmaMz);

        Assert.That(deduplicated.InitialResiduals!.DistinctPeaks,
            Is.LessThan(perSample.InitialResiduals!.DistinctPeaks),
            "Overlapping proteoforms must share source peaks, or there is nothing to deduplicate.");
        Assert.That(deduplicated.InitialResiduals.UnexplainedFraction,
            Is.GreaterThan(perSample.InitialResiduals.UnexplainedFraction * 1.05),
            "Double-counting the well-explained crowded region must have been suppressing the fraction.");
    }

    [Test]
    public void SamplesWithNoMatchingPeakAreReportedSeparately()
    {
        // Predicting signal where the instrument recorded nothing is the opposite failure from
        // failing to explain signal that is there, and has the opposite fix. Folding them into one
        // ratio hides which one is happening.
        const int minCharge = 6;
        const int maxCharge = 11;
        const double sigmaMz = 0.012;

        var present = new ProteoformModel(10000.0, 1.4e6, new EmgProfile(20.0, 0.16, 0.06),
            new GaussianChargeDistribution(8.2, 1.1), "present");
        var scanArray = BuildCentroidedScansFromForwardModel(new[] { present }, minCharge, maxCharge, sigmaMz);
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

        // A proteoform the data does not contain: every sample it claims comes back empty.
        var absent = new ProteoformModel(11234.5, 9e5, new EmgProfile(20.0, 0.16, 0.06),
            new GaussianChargeDistribution(8.2, 1.1), "absent");

        var truths = new[]
        {
            extractor.Extract(present.MonoisotopicMass, 20.0, 0.8, minCharge, maxCharge),
            extractor.Extract(absent.MonoisotopicMass, 20.0, 0.8, minCharge, maxCharge),
        };
        var fits = new[]
        {
            new FittedProteoform(present, sigmaMz, WidthFitMode.CentroidedFallback, 0, 0),
            new FittedProteoform(absent, sigmaMz, WidthFitMode.CentroidedFallback, 0, 0),
        };

        var result = new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(MaxIterations: 20, ConvergenceTolerance: 1e-5))
            .Refit(fits, truths, minCharge, maxCharge, sigmaMz);

        Assert.That(result.InitialResiduals!.UnmatchedSamples, Is.GreaterThan(0),
            "The absent proteoform's samples should have found no experimental peak.");
        Assert.That(result.InitialResiduals.OverpredictedFraction, Is.GreaterThan(0),
            "Predicting a proteoform that is not there must show up as overprediction.");

        // Driving the absent proteoform's abundance to zero is exactly what the refit should do,
        // and it should register in the overprediction term rather than the unexplained one.
        Assert.That(result.FittedProteoforms[1].Model.Abundance, Is.LessThan(absent.Abundance * 0.01));
        Assert.That(result.FinalResiduals!.OverpredictedFraction,
            Is.LessThan(result.InitialResiduals.OverpredictedFraction));
    }

    private static GlobalAbundanceRefitResult Refit(
        FittedProteoform[] fits, ProteoformGroundTruth[] truths, int minCharge, int maxCharge, double sigmaMz) =>
        new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(MaxIterations: 1))
            .Refit(fits, truths, minCharge, maxCharge, sigmaMz);

    /// <summary>
    /// The same truth with no peak identity, which takes the refitter's pre-deduplication path:
    /// every sample is treated as its own peak, exactly as before source identity was recorded.
    /// </summary>
    private static ProteoformGroundTruth StripSourcePeakIdentity(ProteoformGroundTruth truth) =>
        new()
        {
            MonoisotopicMass = truth.MonoisotopicMass,
            RetentionTimeCenter = truth.RetentionTimeCenter,
            MinCharge = truth.MinCharge,
            MaxCharge = truth.MaxCharge,
            ZeroBasedScanIndices = truth.ZeroBasedScanIndices,
            ScanTimes = truth.ScanTimes,
            CentroidMzs = truth.CentroidMzs,
            IsotopologueIntensities = truth.IsotopologueIntensities,
            IsotopologuePeakWindows = truth.IsotopologuePeakWindows,
            // SourcePeakIndices/SourcePeakDeltas deliberately left empty — that is what makes this
            // take the pre-deduplication path.
            ChargeXics = truth.ChargeXics,
            MzWindowHalfWidth = truth.MzWindowHalfWidth,
            ApexChargeIndex = truth.ApexChargeIndex,
            ApexScanIndex = truth.ApexScanIndex,
            ApexScanIndexByCharge = truth.ApexScanIndexByCharge,
        };

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
