using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Comparison;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class SimulationAndComparisonTests
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
    public void SimulatorBuildsConsistentMsDataFile()
    {
        var model = new ProteoformModel(
            MonoisotopicMass: 10000.0,
            Abundance: 1.5e6,
            RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15));

        double[] scanTimes = Enumerable.Range(0, 11).Select(i => 19.5 + i * 0.1).ToArray();
        var result = new Simulator().Simulate(new[] { model }, minCharge: 6, maxCharge: 11, sigmaMz: 0.012, scanTimes: scanTimes);

        Assert.That(result.Scans, Has.Length.EqualTo(scanTimes.Length));
        Assert.That(result.DataFile.NumSpectra, Is.EqualTo(scanTimes.Length));
        Assert.That(result.Grid.MzGrid.Length, Is.GreaterThan(10));
        Assert.That(result.Scans[0].MassSpectrum.XArray.Length, Is.EqualTo(result.Grid.MzGrid.Length));
    }

    [Test]
    public void ComparisonMetricsAreNearPerfectForMatchingModel()
    {
        const double mass = 10000.0;
        const int minCharge = 6;
        const int maxCharge = 11;
        const double sigmaMz = 0.012;
        var truthModel = new ProteoformModel(
            MonoisotopicMass: mass,
            Abundance: 1.5e6,
            RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15));

        var fm = new ForwardModel(new[] { truthModel }, minCharge, maxCharge, sigmaMz);
        var kernel = new IsotopeEnvelopeKernel(mass);
        int nIso = kernel.IsotopologueCount;

        var scans = new List<MsDataScan>();
        for (int s = 0; s < 81; s++)
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
            scans.Add(BuildMs1Scan(s + 1, rt, ordered.Select(t => t.m).ToArray(), ordered.Select(t => t.i).ToArray()));
        }

        var scanArray = scans.ToArray();
        var index = PeakIndexingEngine.InitializeIndexingEngine(scanArray)!;
        var extractor = new GroundTruthExtractor(index, scanArray, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(mass, rtCenter: 20.0, rtHalfWidth: 2.0, minCharge: minCharge, maxCharge: maxCharge);

        var report = ComparisonReportBuilder.Create(truth, new[] { truthModel }, sigmaMz);

        Assert.That(report.PerScanSpectralAngles.Min(), Is.GreaterThan(0.999));
        Assert.That(report.PerChargeXicCorrelations.All(p => p.Correlation > 0.999), Is.True);
        Assert.That(report.Residuals.ResidualEnergyFraction, Is.LessThan(1e-12));
    }
}
