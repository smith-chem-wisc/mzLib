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

    [Test]
    public void ExtractsCorrectTensorShapeAndIntensities()
    {
        const double mass = 10000.0;
        const int minCharge = 7;
        const int maxCharge = 10;
        const double apexRt = 20.0;

        var kernel = new IsotopeEnvelopeKernel(mass);
        int nIso = kernel.IsotopologueCount;

        // Chargewise intensity scales so the extractor can confirm per-charge XIC sums.
        double[] chargeScale = { 100.0, 200.0, 300.0, 150.0 }; // z=7,8,9,10

        // RT profile: triangle peaked at apex; 11 scans spaced 0.1 min apart.
        double[] rtScale = { 0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 0.75, 0.5, 0.25, 0 };

        var scans = new List<MsDataScan>();
        for (int s = 0; s < rtScale.Length; s++)
        {
            double rt = apexRt - 0.5 + s * 0.1;
            var mzList = new List<double>();
            var intList = new List<double>();

            for (int z = minCharge; z <= maxCharge; z++)
            {
                int c = z - minCharge;
                double[] centroids = kernel.CentroidMzs(z);
                for (int i = 0; i < nIso; i++)
                {
                    mzList.Add(centroids[i]);
                    intList.Add(chargeScale[c] * rtScale[s] * kernel.Intensity(i));
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

        var index = PeakIndexingEngine.InitializeIndexingEngine(scans.ToArray())!;
        var extractor = new GroundTruthExtractor(index, ppmTolerance: 20.0);

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
    }
}
