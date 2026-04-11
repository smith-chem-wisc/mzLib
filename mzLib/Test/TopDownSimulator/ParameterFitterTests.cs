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
public class ParameterFitterTests
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
    public void FitsRoundTripThroughForwardModel()
    {
        // Build synthetic data from a known ProteoformModel, run the full fitter,
        // and verify the fitted model reproduces the inputs.
        const double mass = 10000.0;
        const int minCharge = 6, maxCharge = 11;
        const double sigmaMz = 0.012;
        var truthModel = new ProteoformModel(
            MonoisotopicMass: mass,
            Abundance: 1.5e6,
            RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15));

        var fm = new ForwardModel(new[] { truthModel }, minCharge, maxCharge, sigmaMz);
        var kernel = new IsotopeEnvelopeKernel(mass);
        int nIso = kernel.IsotopologueCount;

        // 81 scans spaced 0.05 min apart, centered on the RT profile's mean (~μ+τ).
        int nS = 81;
        var scans = new List<MsDataScan>();
        for (int s = 0; s < nS; s++)
        {
            double rt = 18.0 + s * 0.05;
            var mzList = new List<double>();
            var intList = new List<double>();
            for (int z = minCharge; z <= maxCharge; z++)
            {
                double[] cent = kernel.CentroidMzs(z);
                for (int i = 0; i < nIso; i++)
                {
                    mzList.Add(cent[i]);
                    intList.Add(fm.Evaluate(rt, cent[i]));
                }
            }
            var ordered = mzList.Select((m, k) => (m, i: intList[k])).OrderBy(t => t.m).ToList();
            scans.Add(BuildMs1Scan(s + 1, rt,
                ordered.Select(t => t.m).ToArray(),
                ordered.Select(t => t.i).ToArray()));
        }

        var scanArray = scans.ToArray();
        var index = PeakIndexingEngine.InitializeIndexingEngine(scanArray)!;
        var extractor = new GroundTruthExtractor(index, scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(mass, rtCenter: 20.0, rtHalfWidth: 2.0,
            minCharge: minCharge, maxCharge: maxCharge);

        var fit = new ParameterFitter(widthFitter: new EnvelopeWidthFitter(fallbackSigmaMz: sigmaMz)).Fit(truth);

        // Width fell back (centroided input), so σ_m == fallback by construction.
        Assert.That(fit.WidthMode, Is.EqualTo(WidthFitMode.CentroidedFallback));
        Assert.That(fit.SigmaMz, Is.EqualTo(sigmaMz));

        // RT profile recovered.
        Assert.That(fit.Model.RtProfile.Mu, Is.EqualTo(truthModel.RtProfile.Mu).Within(5e-3));
        Assert.That(fit.Model.RtProfile.Sigma, Is.EqualTo(truthModel.RtProfile.Sigma).Within(2e-2));
        Assert.That(fit.Model.RtProfile.Tau, Is.EqualTo(truthModel.RtProfile.Tau).Within(2e-2));

        // Charge distribution recovered.
        Assert.That(fit.Model.ChargeDistribution, Is.TypeOf<GaussianChargeDistribution>());
        var gz = (GaussianChargeDistribution)fit.Model.ChargeDistribution;
        Assert.That(gz.MuZ, Is.EqualTo(truthModel.ChargeDistribution.Evaluate(0) > 0 ? gz.MuZ : 0)); // non-strict
        Assert.That(gz.MuZ, Is.EqualTo(8.3).Within(0.15));
        Assert.That(gz.SigmaZ, Is.EqualTo(1.15).Within(0.2));

        // Abundance recovered within a few percent.
        Assert.That(fit.Model.Abundance, Is.EqualTo(truthModel.Abundance).Within(truthModel.Abundance * 0.05));
    }
}
