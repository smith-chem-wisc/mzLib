using System.Collections.Generic;
using System.Linq;
using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class GroundTruthExtractorTests
{
    private static MsDataScan BuildMs1Scan(int oneBasedScanNumber, double rt, double[] mz, double[] intensity)
    {
        return new MsDataScan(
            massSpectrum: new MzSpectrum(mz, intensity, false),
            oneBasedScanNumber: oneBasedScanNumber,
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
            nativeId: "scan=" + oneBasedScanNumber);
    }

    private const double Mass = 10000.0;
    private const int MinCharge = 7;
    private const int MaxCharge = 10;
    private const double ApexRt = 20.0;

    /// <summary>Chargewise intensity scales, so per-charge XIC sums are predictable. z = 7,8,9,10.</summary>
    private static readonly double[] ChargeScale = { 100.0, 200.0, 300.0, 150.0 };

    /// <summary>Triangular RT profile peaked at the apex; 11 scans spaced 0.1 min apart.</summary>
    private static readonly double[] RtScale = { 0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 0.75, 0.5, 0.25, 0 };

    private static MsDataScan[] BuildScans()
    {
        var kernel = new IsotopeEnvelopeKernel(Mass);
        int nIso = kernel.IsotopologueCount;
        var scans = new List<MsDataScan>();

        for (int s = 0; s < RtScale.Length; s++)
        {
            double rt = ApexRt - 0.5 + s * 0.1;
            var mzList = new List<double>();
            var intList = new List<double>();

            for (int z = MinCharge; z <= MaxCharge; z++)
            {
                int c = z - MinCharge;
                double[] centroids = kernel.CentroidMzs(z);
                for (int i = 0; i < nIso; i++)
                {
                    mzList.Add(centroids[i]);
                    intList.Add(ChargeScale[c] * RtScale[s] * kernel.Intensity(i));
                }
            }

            // Zip + sort by m/z — MzSpectrum expects monotonic x.
            var ordered = mzList.Select((m, k) => (m, i: intList[k])).OrderBy(t => t.m).ToList();
            scans.Add(BuildMs1Scan(
                oneBasedScanNumber: s + 1,
                rt: rt,
                mz: ordered.Select(t => t.m).ToArray(),
                intensity: ordered.Select(t => t.i).ToArray()));
        }

        return scans.ToArray();
    }

    [Test]
    public void ExtractsCorrectTensorShapeAndIntensities()
    {
        const double mass = Mass;
        const int minCharge = MinCharge;
        const int maxCharge = MaxCharge;
        const double apexRt = ApexRt;

        var kernel = new IsotopeEnvelopeKernel(mass);
        int nIso = kernel.IsotopologueCount;
        var chargeScale = ChargeScale;

        var scanArray = BuildScans();
        var extractor = new GroundTruthExtractor(scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

        var truth = extractor.Extract(
            monoisotopicMass: mass,
            rtCenter: apexRt,
            rtHalfWidth: 1.0,
            minCharge: minCharge,
            maxCharge: maxCharge);

        Assert.That(truth.ChargeCount, Is.EqualTo(4));
        Assert.That(truth.ScanCount, Is.EqualTo(11));
        Assert.That(truth.CentroidMzs.Length, Is.EqualTo(4));
        Assert.That(truth.CentroidMzs[0].Length, Is.EqualTo(nIso));

        // Per-charge XIC apex (scans 4–6, rtScale = 1.0) should equal chargeScale (sum of normalized envelope = 1).
        for (int c = 0; c < 4; c++)
        {
            double apex = truth.ChargeXics[c][5]; // middle scan
            Assert.That(apex, Is.EqualTo(chargeScale[c]).Within(1e-3));
        }

        // RT edge scans have rtScale = 0 → XIC should be zero.
        for (int c = 0; c < 4; c++)
        {
            Assert.That(truth.ChargeXics[c][0], Is.EqualTo(0).Within(1e-9));
            Assert.That(truth.ChargeXics[c][^1], Is.EqualTo(0).Within(1e-9));
        }

        // First isotopologue intensity at apex should equal chargeScale[c] * kernel.Intensity(0).
        for (int c = 0; c < 4; c++)
        {
            double expected = chargeScale[c] * kernel.Intensity(0);
            Assert.That(truth.IsotopologueIntensities[c][0][5], Is.EqualTo(expected).Within(1e-3));
        }

        // Peak windows are collected only at the apex (charge, scan) cell, which is the only one
        // any consumer reads. The apex here is the strongest charge at the middle scan.
        Assert.That(truth.MzWindowHalfWidth, Is.EqualTo(0.05));
        Assert.That(truth.IsotopologuePeakWindows.Length, Is.EqualTo(4));

        // Strongest charge wins; the RT profile plateaus over scans 4-6, and ties resolve to the
        // first, so the apex scan is 4.
        Assert.That(truth.ApexChargeIndex,
            Is.EqualTo(System.Array.IndexOf(chargeScale, chargeScale.Max())));
        Assert.That(truth.ApexScanIndex, Is.EqualTo(4));

        int apexC = truth.ApexChargeIndex;
        var window = truth.IsotopologuePeakWindows[apexC][0][truth.ApexScanIndex];
        Assert.That(window.Length, Is.GreaterThanOrEqualTo(1));

        // There should be a peak in the window whose m/z is ~the centroid and whose intensity
        // equals the scalar cell.
        double apexCentroid = truth.CentroidMzs[apexC][0];
        bool foundMatch = window.Any(p =>
            System.Math.Abs(p.Mz - apexCentroid) < 1e-6
            && System.Math.Abs(p.Intensity - truth.IsotopologueIntensities[apexC][0][truth.ApexScanIndex]) < 1e-3);
        Assert.That(foundMatch, Is.True);

        // Every other cell is allocated but empty rather than null, so consumers can index freely.
        for (int c = 0; c < 4; c++)
        {
            Assert.That(truth.IsotopologuePeakWindows[c][0][0], Is.Not.Null.And.Empty);
            Assert.That(truth.IsotopologuePeakWindows[c][0][^1], Is.Not.Null.And.Empty);
        }
    }

    [Test]
    public void ApexCellMatchesTheOneEnvelopeWidthFitterSelects()
    {
        // The extractor populates peak windows only at the apex it computes, and the width fitter
        // independently recomputes an apex to read them. If those ever diverge the fitter silently
        // sees empty windows and falls back, so pin them together.
        var extractor = new GroundTruthExtractor(BuildScans(), ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(Mass, ApexRt, rtHalfWidth: 1.0, minCharge: MinCharge, maxCharge: MaxCharge);

        int bestC = 0, bestS = 0;
        double bestApex = double.NegativeInfinity;
        for (int c = 0; c < truth.ChargeCount; c++)
        {
            for (int s = 0; s < truth.ScanCount; s++)
            {
                double v = truth.ChargeXics[c][s];
                if (v > bestApex)
                {
                    bestApex = v;
                    bestC = c;
                    bestS = s;
                }
            }
        }

        Assert.That(truth.ApexChargeIndex, Is.EqualTo(bestC));
        Assert.That(truth.ApexScanIndex, Is.EqualTo(bestS));
    }
}
