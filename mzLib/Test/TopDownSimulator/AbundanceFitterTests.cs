using System;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class AbundanceFitterTests
{
    [Test]
    public void RecoversAbundanceFromSyntheticIntensities()
    {
        const double mass = 10000.0;
        const int minCharge = 7, maxCharge = 10;
        const double sigmaMz = 0.01;
        const double trueAbundance = 4.2e5;

        var rt = new EmgProfile(Mu: 20.0, Sigma: 0.25, Tau: 0.1);
        var charge = new GaussianChargeDistribution(MuZ: 8.5, SigmaZ: 1.1);
        var kernel = new IsotopeEnvelopeKernel(mass);
        int nIso = kernel.IsotopologueCount;
        int nC = maxCharge - minCharge + 1;

        // 21 scans across the RT profile.
        int nS = 21;
        var scanTimes = new double[nS];
        for (int s = 0; s < nS; s++) scanTimes[s] = 19.0 + s * 0.1;

        double peakScale = 1.0 / (sigmaMz * Math.Sqrt(2 * Math.PI));

        var centroidMzs = new double[nC][];
        var intensities = new double[nC][][];
        var windows = new PeakSample[nC][][][];
        var chargeXics = new double[nC][];
        for (int c = 0; c < nC; c++)
        {
            int z = minCharge + c;
            centroidMzs[c] = kernel.CentroidMzs(z);
            intensities[c] = new double[nIso][];
            windows[c] = new PeakSample[nIso][][];
            chargeXics[c] = new double[nS];
            for (int i = 0; i < nIso; i++)
            {
                intensities[c][i] = new double[nS];
                windows[c][i] = new PeakSample[nS][];
                for (int s = 0; s < nS; s++) windows[c][i][s] = Array.Empty<PeakSample>();
            }
        }

        for (int s = 0; s < nS; s++)
        {
            double gt = rt.Evaluate(scanTimes[s]);
            for (int c = 0; c < nC; c++)
            {
                int z = minCharge + c;
                double fz = charge.Evaluate(z);
                for (int i = 0; i < nIso; i++)
                {
                    double val = trueAbundance * gt * fz * kernel.Intensity(i) * peakScale;
                    intensities[c][i][s] = val;
                    chargeXics[c][s] += val;
                }
            }
        }

        var truth = new ProteoformGroundTruth
        {
            MonoisotopicMass = mass,
            RetentionTimeCenter = 20.0,
            MinCharge = minCharge,
            MaxCharge = maxCharge,
            ZeroBasedScanIndices = new int[nS],
            ScanTimes = scanTimes,
            CentroidMzs = centroidMzs,
            IsotopologueIntensities = intensities,
            IsotopologuePeakWindows = windows,
            ChargeXics = chargeXics,
            MzWindowHalfWidth = 0.05,
        };

        var fit = new AbundanceFitter().Fit(truth, sigmaMz, rt, charge);

        Assert.That(fit.Abundance, Is.EqualTo(trueAbundance).Within(1e-6 * trueAbundance));
        Assert.That(fit.Residual, Is.LessThan(1e-3));
        Assert.That(fit.SamplesUsed, Is.GreaterThan(0));
    }
}
