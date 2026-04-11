using System;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class RtProfileFitterTests
{
    /// <summary>Builds a minimal ground truth with only the RT-relevant fields populated.</summary>
    private static ProteoformGroundTruth BuildTruth(double[] scanTimes, double[] totalsPerScan)
    {
        int n = scanTimes.Length;
        var chargeXics = new double[1][];
        chargeXics[0] = totalsPerScan;
        var centroidMzs = new double[1][];
        centroidMzs[0] = new[] { 500.0 };
        var intensities = new double[1][][];
        intensities[0] = new double[1][];
        intensities[0][0] = new double[n];
        var windows = new PeakSample[1][][][];
        windows[0] = new PeakSample[1][][];
        windows[0][0] = new PeakSample[n][];
        for (int i = 0; i < n; i++) windows[0][0][i] = Array.Empty<PeakSample>();
        var scanIdx = new int[n];
        for (int i = 0; i < n; i++) scanIdx[i] = i;

        return new ProteoformGroundTruth
        {
            MonoisotopicMass = 10000,
            RetentionTimeCenter = 0,
            MinCharge = 8,
            MaxCharge = 8,
            ZeroBasedScanIndices = scanIdx,
            ScanTimes = scanTimes,
            CentroidMzs = centroidMzs,
            IsotopologueIntensities = intensities,
            IsotopologuePeakWindows = windows,
            ChargeXics = chargeXics,
            MzWindowHalfWidth = 0.05,
        };
    }

    [Test]
    public void SymmetricGaussianRecoversMuAndSigma()
    {
        const double trueMu = 20.0, trueSigma = 0.3;
        int n = 201;
        double[] times = new double[n];
        double[] y = new double[n];
        double dt = 0.02;
        for (int i = 0; i < n; i++)
        {
            double t = trueMu - (n / 2) * dt + i * dt;
            times[i] = t;
            double d = (t - trueMu) / trueSigma;
            y[i] = Math.Exp(-0.5 * d * d);
        }
        var truth = BuildTruth(times, y);

        var fit = new RtProfileFitter().Fit(truth);

        Assert.That(fit.Profile.Mu, Is.EqualTo(trueMu).Within(5e-3));
        Assert.That(fit.Profile.Sigma, Is.EqualTo(trueSigma).Within(5e-3));
        Assert.That(fit.Profile.Tau, Is.LessThanOrEqualTo(0.02));
    }

    [Test]
    public void EmgRecoversMuSigmaTau()
    {
        const double trueMu = 15.0, trueSigma = 0.25, trueTau = 0.35;
        var source = new EmgProfile(trueMu, trueSigma, trueTau);

        int n = 801;
        double dt = 0.01;
        double[] times = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double t = trueMu - 4.0 + i * dt;
            times[i] = t;
            y[i] = source.Evaluate(t);
        }
        var truth = BuildTruth(times, y);

        var fit = new RtProfileFitter().Fit(truth);

        Assert.That(fit.Profile.Mu, Is.EqualTo(trueMu).Within(5e-3));
        Assert.That(fit.Profile.Sigma, Is.EqualTo(trueSigma).Within(1e-2));
        Assert.That(fit.Profile.Tau, Is.EqualTo(trueTau).Within(1e-2));
    }
}
