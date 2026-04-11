using System.Linq;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Plotly.NET.CSharp;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.FileReadingTests;

[TestFixture]
internal class AnalysisExample
{
    [Test]
    [Explicit("Interactive Plotly demo for a synthetic MS1 spectrum")]
    public void PlotSyntheticSpectrum()
    {
        var simulation = BuildSimulation();
        var scan = simulation.Scans[simulation.Scans.Length / 2];

        Chart.Line<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
            .WithTraceInfo("Synthetic MS1 Spectrum")
            .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
            .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
            .WithSize(Width: 1000, Height: 500)
            .Show();
    }

    [Test]
    [Explicit("Interactive Plotly demo for a synthetic charge XIC")]
    public void PlotSyntheticChargeXic()
    {
        const double mass = 10000.0;
        var simulation = BuildSimulation();
        var index = PeakIndexingEngine.InitializeIndexingEngine(simulation.Scans)!;
        var extractor = new GroundTruthExtractor(index, simulation.Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(mass, rtCenter: 20.0, rtHalfWidth: 1.0, minCharge: 6, maxCharge: 11);

        int chargeOffset = 8 - truth.MinCharge;
        Chart.Line<double, double, string>(truth.ScanTimes, truth.ChargeXics[chargeOffset])
            .WithTraceInfo("Charge 8 XIC")
            .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
            .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
            .WithSize(Width: 1000, Height: 500)
            .Show();
    }

    private static SimulationResult BuildSimulation()
    {
        var model = new ProteoformModel(
            MonoisotopicMass: 10000.0,
            Abundance: 1.5e6,
            RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15));

        double[] scanTimes = Enumerable.Range(0, 31).Select(i => 18.5 + i * 0.1).ToArray();
        return new Simulator().Simulate(new[] { model }, minCharge: 6, maxCharge: 11, sigmaMz: 0.012, scanTimes: scanTimes);
    }
}
