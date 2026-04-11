using System;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class ChargeDistributionFitterTests
{
    private static ProteoformGroundTruth BuildTruth(int minCharge, double[] perChargeApex)
    {
        int nC = perChargeApex.Length;
        // Encode apex as a constant single-scan XIC per charge.
        int nS = 3;
        var chargeXics = new double[nC][];
        for (int c = 0; c < nC; c++)
        {
            chargeXics[c] = new double[nS];
            chargeXics[c][1] = perChargeApex[c];
        }
        var centroidMzs = new double[nC][];
        var intensities = new double[nC][][];
        var windows = new PeakSample[nC][][][];
        for (int c = 0; c < nC; c++)
        {
            centroidMzs[c] = new[] { 500.0 };
            intensities[c] = new double[1][];
            intensities[c][0] = new double[nS];
            windows[c] = new PeakSample[1][][];
            windows[c][0] = new PeakSample[nS][];
            for (int s = 0; s < nS; s++) windows[c][0][s] = Array.Empty<PeakSample>();
        }

        return new ProteoformGroundTruth
        {
            MonoisotopicMass = 10000,
            RetentionTimeCenter = 0,
            MinCharge = minCharge,
            MaxCharge = minCharge + nC - 1,
            ZeroBasedScanIndices = new[] { 0, 1, 2 },
            ScanTimes = new[] { 0.0, 1.0, 2.0 },
            CentroidMzs = centroidMzs,
            IsotopologueIntensities = intensities,
            IsotopologuePeakWindows = windows,
            ChargeXics = chargeXics,
            MzWindowHalfWidth = 0.05,
        };
    }

    [Test]
    public void SymmetricGaussianOnChargeRecoversMuSigma()
    {
        const double trueMu = 9.0, trueSigma = 1.2;
        int minZ = 5, maxZ = 13;
        int nC = maxZ - minZ + 1;
        var apex = new double[nC];
        for (int c = 0; c < nC; c++)
        {
            int z = minZ + c;
            double d = (z - trueMu) / trueSigma;
            apex[c] = Math.Exp(-0.5 * d * d);
        }
        var truth = BuildTruth(minZ, apex);

        var fit = new ChargeDistributionFitter().Fit(truth);

        Assert.That(fit.Distribution.MuZ, Is.EqualTo(trueMu).Within(1e-6));
        Assert.That(fit.Distribution.SigmaZ, Is.EqualTo(trueSigma).Within(5e-2));
        Assert.That(fit.ChargesUsed, Is.EqualTo(nC));
    }

    [Test]
    public void SingleChargeFallsBackOnSigma()
    {
        var truth = BuildTruth(minCharge: 8, perChargeApex: new[] { 0.0, 100.0, 0.0 });
        var fit = new ChargeDistributionFitter(fallbackSigmaZ: 1.5).Fit(truth);

        Assert.That(fit.Distribution.MuZ, Is.EqualTo(9.0).Within(1e-9));
        Assert.That(fit.Distribution.SigmaZ, Is.EqualTo(1.5));
        Assert.That(fit.ChargesUsed, Is.EqualTo(1));
    }
}
