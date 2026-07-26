using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class EnvelopeWidthFitterTests
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
    public void CentroidedInputFallsBack()
    {
        // Use the exact same centroided synthesis as GroundTruthExtractorTests —
        // one peak per theoretical centroid.
        const double mass = 10000.0;
        const int minCharge = 8, maxCharge = 9;
        var kernel = new IsotopeEnvelopeKernel(mass);
        int nIso = kernel.IsotopologueCount;

        var scans = new List<MsDataScan>();
        double[] rtScale = { 0, 1, 0 };
        for (int s = 0; s < rtScale.Length; s++)
        {
            var mzList = new List<double>();
            var intList = new List<double>();
            for (int z = minCharge; z <= maxCharge; z++)
            {
                double[] cent = kernel.CentroidMzs(z);
                for (int i = 0; i < nIso; i++)
                {
                    mzList.Add(cent[i]);
                    intList.Add(100 * rtScale[s] * kernel.Intensity(i));
                }
            }
            var ordered = mzList.Select((m, k) => (m, i: intList[k])).OrderBy(t => t.m).ToList();
            scans.Add(BuildMs1Scan(s + 1, 20.0 + s * 0.1,
                ordered.Select(t => t.m).ToArray(),
                ordered.Select(t => t.i).ToArray()));
        }

        var scanArray = scans.ToArray();
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(mass, rtCenter: 20.1, rtHalfWidth: 0.5,
            minCharge: minCharge, maxCharge: maxCharge);

        var fit = new EnvelopeWidthFitter(fallbackSigmaMz: 0.013).Fit(truth);

        Assert.That(fit.Mode, Is.EqualTo(WidthFitMode.CentroidedFallback));
        Assert.That(fit.SigmaMz, Is.EqualTo(0.013));
        Assert.That(fit.PeaksUsed, Is.EqualTo(0));
    }

    [Test]
    public void ProfileInputRecoversSigma()
    {
        // Build profile data: each theoretical centroid is rasterized onto a dense
        // m/z grid with a Gaussian of known width. The fitter should recover it.
        const double mass = 10000.0;
        const int minCharge = 8, maxCharge = 9;
        const double trueSigmaMz = 0.008;
        var kernel = new IsotopeEnvelopeKernel(mass);

        var scanArray = BuildProfileScans(kernel, minCharge, maxCharge, trueSigmaMz);
        // Window half-width must be a bit less than half the centroid spacing to avoid
        // neighboring-isotopologue spillover (at z=8 centroids are ~0.125 m/z apart).
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(mass, rtCenter: 20.1, rtHalfWidth: 0.5,
            minCharge: minCharge, maxCharge: maxCharge);

        var fit = new EnvelopeWidthFitter(fallbackSigmaMz: 0.02).Fit(truth);

        Assert.That(fit.Mode, Is.EqualTo(WidthFitMode.Profile));
        Assert.That(fit.SigmaMz, Is.EqualTo(trueSigmaMz).Within(5e-4));
        Assert.That(fit.PeaksUsed, Is.GreaterThan(10));
    }

    [Test]
    public void ProfileInputRecoversSigmaWhenNeighbouringIsotopologuesShareTheWindow()
    {
        // At z = 20 isotopologues are ~0.050 m/z apart, so the extractor's ±0.05 window around one
        // centroid swallows most of each neighbour's peak. Separating them is the entire job of the
        // midpoint boundaries in EnvelopeWidthFitter, and those boundaries are only meaningful if
        // CentroidMzs is ordered by m/z — index i-1 and i+1 have to be the actual neighbours.
        // Against the intensity-ordered envelope this fixture recovered 0.0178, 4.5x the truth.
        const double mass = 10000.0;
        const int minCharge = 20, maxCharge = 21;
        const double trueSigmaMz = 0.004;
        var kernel = new IsotopeEnvelopeKernel(mass);

        double spacing = kernel.CentroidMzs(minCharge)[1] - kernel.CentroidMzs(minCharge)[0];
        Assert.That(spacing, Is.GreaterThan(0).And.LessThan(2 * 0.05),
            "This fixture is only meaningful if neighbouring isotopologues fall inside the window.");

        var scanArray = BuildProfileScans(kernel, minCharge, maxCharge, trueSigmaMz);
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(mass, rtCenter: 20.1, rtHalfWidth: 0.5,
            minCharge: minCharge, maxCharge: maxCharge);

        var fit = new EnvelopeWidthFitter(fallbackSigmaMz: 0.012).Fit(truth);

        Assert.That(fit.Mode, Is.EqualTo(WidthFitMode.Profile));
        Assert.That(fit.SigmaMz, Is.EqualTo(trueSigmaMz).Within(5e-4));
        Assert.That(fit.PeaksUsed, Is.GreaterThan(10));
    }

    /// <summary>
    /// Three scans of profile data in which every theoretical isotopologue of every charge is a
    /// Gaussian of width <paramref name="trueSigmaMz"/>, sampled to ±5σ at 0.001 m/z.
    /// </summary>
    private static MsDataScan[] BuildProfileScans(
        IsotopeEnvelopeKernel kernel, int minCharge, int maxCharge, double trueSigmaMz)
    {
        const double sampleStep = 0.001;
        double halfSpan = 5.0 * trueSigmaMz;
        int nIso = kernel.IsotopologueCount;

        var scans = new List<MsDataScan>();
        double[] rtScale = { 0.2, 1.0, 0.2 };
        for (int s = 0; s < rtScale.Length; s++)
        {
            var mzList = new List<double>();
            var intList = new List<double>();
            for (int z = minCharge; z <= maxCharge; z++)
            {
                double[] cent = kernel.CentroidMzs(z);
                for (int i = 0; i < nIso; i++)
                {
                    double mu = cent[i];
                    double amplitude = 1000.0 * rtScale[s] * kernel.Intensity(i);
                    for (double off = -halfSpan; off <= halfSpan + 1e-12; off += sampleStep)
                    {
                        double mz = mu + off;
                        double y = amplitude * Math.Exp(-(off * off) / (2 * trueSigmaMz * trueSigmaMz));
                        mzList.Add(mz);
                        intList.Add(y);
                    }
                }
            }
            var ordered = mzList.Select((m, k) => (m, i: intList[k])).OrderBy(t => t.m).ToList();
            scans.Add(BuildMs1Scan(s + 1, 20.0 + s * 0.1,
                ordered.Select(t => t.m).ToArray(),
                ordered.Select(t => t.i).ToArray()));
        }

        return scans.ToArray();
    }
}
