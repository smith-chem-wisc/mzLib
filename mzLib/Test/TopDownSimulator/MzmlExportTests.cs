using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class MzmlExportTests
{
    private const double Mass = 10000.0;
    private const int MinCharge = 6;
    private const int MaxCharge = 11;
    private const double SigmaMz = 0.012;

    private string _outputDirectory;

    private static ProteoformModel BuildModel(double mass = Mass, double abundance = 1.5e6, double rtCenter = 20.0) =>
        new ProteoformModel(
            MonoisotopicMass: mass,
            Abundance: abundance,
            RtProfile: new EmgProfile(Mu: rtCenter, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15),
            Identifier: "P00000|TESTPROTEOFORM");

    private static double[] ScanTimes(int count = 11, double start = 19.5, double step = 0.1) =>
        Enumerable.Range(0, count).Select(i => start + i * step).ToArray();

    [SetUp]
    public void SetUp()
    {
        _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownSimulatorMzmlOutput");
        Directory.CreateDirectory(_outputDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_outputDirectory))
            Directory.Delete(_outputDirectory, recursive: true);
    }

    [Test]
    public void CentroidAxisSitsExactlyOnTheoreticalIsotopologueMzs()
    {
        var model = BuildModel();
        var axis = new GridRasterizer().BuildCentroidMzAxis(new[] { model }, MinCharge, MaxCharge, SigmaMz);

        Assert.That(axis, Is.Ordered.Ascending, "The axis must be sorted for the model's binary search.");
        Assert.That(axis, Is.Unique);

        var kernel = new IsotopeEnvelopeKernel(Mass);
        for (int z = MinCharge; z <= MaxCharge; z++)
        {
            foreach (double centroid in kernel.CentroidMzs(z))
            {
                // Every theoretical isotopologue must be represented exactly, not snapped to a grid.
                Assert.That(axis.Any(mz => Math.Abs(mz - centroid) < 1e-9), Is.True,
                    $"Charge {z} isotopologue at {centroid} is missing from the centroid axis.");
            }
        }
    }

    [Test]
    public void CentroidIntensitiesMatchTheForwardModelExactly()
    {
        var model = BuildModel();
        double[] scanTimes = ScanTimes(count: 5);

        var simulation = new Simulator().SimulateCentroided(new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes);
        var forwardModel = new ForwardModel(new[] { model }, MinCharge, MaxCharge, SigmaMz);

        Assert.That(simulation.Scans, Has.Length.EqualTo(scanTimes.Length));
        Assert.That(simulation.Scans.All(s => s.IsCentroid), Is.True);

        for (int s = 0; s < simulation.Scans.Length; s++)
        {
            var scan = simulation.Scans[s];
            for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
            {
                double expected = forwardModel.Evaluate(scanTimes[s], scan.MassSpectrum.XArray[i]);
                Assert.That(scan.MassSpectrum.YArray[i], Is.EqualTo(expected).Within(1e-9).Percent,
                    $"Scan {s + 1} peak {i} deviates from the forward model.");
            }
        }
    }

    [Test]
    public void CentroidSimulationIsFarSmallerThanProfileSimulation()
    {
        var model = BuildModel();
        double[] scanTimes = ScanTimes();

        var profile = new Simulator().Simulate(new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes);
        var centroid = new Simulator().SimulateCentroided(new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes);

        int profilePoints = profile.Grid.MzGrid.Length;
        int centroidPoints = centroid.Grid.MzGrid.Length;

        Assert.That(centroidPoints, Is.GreaterThan(0));
        Assert.That(centroidPoints, Is.LessThan(profilePoints / 10),
            $"Centroid axis ({centroidPoints}) should be far smaller than the profile grid ({profilePoints}).");
    }

    [Test]
    public void ReducerAppliesGlobalRatherThanPerScanFloor()
    {
        // Two scans three orders of magnitude apart. A per-scan relative floor would keep the
        // faint scan's signal; a global floor must discard it.
        var model = BuildModel(abundance: 1.5e6, rtCenter: 20.0);
        double[] scanTimes = { 20.0, 24.0 };

        var simulation = new Simulator().SimulateCentroided(new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes);
        var reduced = SimulatedScanReducer.Reduce(simulation.Scans,
            new ScanReductionOptions { RelativeIntensityThreshold = 1e-4 });

        Assert.That(reduced[0].MassSpectrum.XArray, Is.Not.Empty, "Apex scan should retain peaks.");
        Assert.That(reduced[1].MassSpectrum.XArray, Is.Empty,
            "A scan four minutes off the elution apex should fall below the global floor.");
    }

    [Test]
    public void ReducerFlagsOutputAsCentroidedAndRecomputesTic()
    {
        var model = BuildModel();
        var simulation = new Simulator().SimulateCentroided(new[] { model }, MinCharge, MaxCharge, SigmaMz, ScanTimes());

        var reduced = SimulatedScanReducer.Reduce(simulation.Scans,
            new ScanReductionOptions { RelativeIntensityThreshold = 1e-3 });

        Assert.That(reduced.All(s => s.IsCentroid), Is.True);
        foreach (var scan in reduced)
            Assert.That(scan.TotalIonCurrent, Is.EqualTo(scan.MassSpectrum.YArray.Sum()).Within(1e-6),
                "TIC should reflect the retained peaks, not the pre-threshold spectrum.");
    }

    [Test]
    public void WriteMzmlProducesCentroidedFileThatReadsBackWithMatchingPeaks()
    {
        var model = BuildModel();
        double[] scanTimes = ScanTimes();
        string path = Path.Combine(_outputDirectory, "simulated.mzML");

        var export = new Simulator().WriteMzml(
            new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes, path);

        Assert.That(File.Exists(path), Is.True, "mzML was not written.");
        Assert.That(export.ScanCount, Is.EqualTo(scanTimes.Length));
        Assert.That(export.PeakCount, Is.GreaterThan(0));

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        var scans = roundTripped.GetAllScansList();

        Assert.That(scans, Has.Count.EqualTo(scanTimes.Length));
        Assert.That(scans.Sum(s => s.MassSpectrum.XArray.Length), Is.EqualTo(export.PeakCount),
            "Peak count reported by the exporter should match what the reader finds.");
        Assert.That(scans.All(s => s.MsnOrder == 1), Is.True, "All simulated scans should be MS1.");
        Assert.That(scans.All(s => s.IsCentroid), Is.True, "Centroid flag should survive the round trip.");

        for (int i = 0; i < scans.Count; i++)
            Assert.That(scans[i].RetentionTime, Is.EqualTo(scanTimes[i]).Within(1e-6),
                $"Retention time drifted for scan {i + 1}.");
    }

    [Test]
    public void RoundTrippedPeakPositionsSurviveAtOrbitrapPrecision()
    {
        var model = BuildModel();
        double[] scanTimes = ScanTimes(count: 3);
        string path = Path.Combine(_outputDirectory, "precision.mzML");

        var simulation = new Simulator().SimulateCentroided(new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes);
        var expected = SimulatedScanReducer.Reduce(simulation.Scans, new ScanReductionOptions());

        new Simulator().WriteMzml(new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes, path);

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        var actual = roundTripped.GetAllScansList();

        for (int s = 0; s < actual.Count; s++)
        {
            var expectedMz = expected[s].MassSpectrum.XArray;
            var actualMz = actual[s].MassSpectrum.XArray;
            Assert.That(actualMz, Has.Length.EqualTo(expectedMz.Length), $"Peak count changed in scan {s + 1}.");

            for (int i = 0; i < expectedMz.Length; i++)
            {
                double ppmError = Math.Abs(actualMz[i] - expectedMz[i]) / expectedMz[i] * 1e6;
                Assert.That(ppmError, Is.LessThan(0.1),
                    $"Scan {s + 1} peak {i} shifted by {ppmError:F4} ppm through the mzML round trip.");
            }
        }
    }

    [Test]
    public void WriteMzmlWritesGroundTruthSidecar()
    {
        var model = BuildModel();
        string path = Path.Combine(_outputDirectory, "simulated-with-truth.mzML");

        var export = new Simulator().WriteMzml(
            new[] { model }, MinCharge, MaxCharge, SigmaMz, ScanTimes(), path);

        Assert.That(export.GroundTruthPath, Is.Not.Null);
        Assert.That(File.Exists(export.GroundTruthPath), Is.True);

        var lines = File.ReadAllLines(export.GroundTruthPath);
        Assert.That(lines, Has.Length.EqualTo(2), "Expected a header plus one proteoform row.");

        var header = lines[0].Split('\t');
        var fields = lines[1].Split('\t');
        Assert.That(header, Does.Contain("MonoisotopicMass"));
        Assert.That(fields[Array.IndexOf(header, "Identifier")], Is.EqualTo("P00000|TESTPROTEOFORM"));
        Assert.That(double.Parse(fields[Array.IndexOf(header, "MonoisotopicMass")]), Is.EqualTo(Mass).Within(1e-9));
        Assert.That(double.Parse(fields[Array.IndexOf(header, "ChargeMu")]), Is.EqualTo(8.3).Within(1e-9));
        Assert.That(double.Parse(fields[Array.IndexOf(header, "RtTau")]), Is.EqualTo(0.08).Within(1e-9));
    }

    [Test]
    public void WriteMzmlCanSkipTheSidecar()
    {
        string path = Path.Combine(_outputDirectory, "no-sidecar.mzML");

        var export = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, SigmaMz, ScanTimes(), path,
            writeGroundTruthSidecar: false);

        Assert.That(export.GroundTruthPath, Is.Null);
        Assert.That(File.Exists(Path.ChangeExtension(path, ".groundtruth.tsv")), Is.False);
    }

    [Test]
    public void CoElutingProteoformsBothAppearInOutput()
    {
        var first = BuildModel(mass: 10000.0, rtCenter: 20.0) with { Identifier = "FIRST" };
        var second = BuildModel(mass: 10500.0, rtCenter: 20.05) with { Identifier = "SECOND" };
        string path = Path.Combine(_outputDirectory, "coeluting.mzML");

        var export = new Simulator().WriteMzml(
            new[] { first, second }, MinCharge, MaxCharge, SigmaMz, ScanTimes(), path);

        Assert.That(export.PeakCount, Is.GreaterThan(0));

        var groundTruth = File.ReadAllLines(export.GroundTruthPath!);
        Assert.That(groundTruth, Has.Length.EqualTo(3), "Both proteoforms should be recorded.");

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        var apex = roundTripped.GetAllScansList().OrderByDescending(s => s.TotalIonCurrent).First();

        var firstKernel = new IsotopeEnvelopeKernel(10000.0);
        var secondKernel = new IsotopeEnvelopeKernel(10500.0);

        Assert.That(NearestPpm(apex.MassSpectrum.XArray, firstKernel.CentroidMzs(8)[0]), Is.LessThan(1.0),
            "Expected a peak at the first proteoform's monoisotopic m/z.");
        Assert.That(NearestPpm(apex.MassSpectrum.XArray, secondKernel.CentroidMzs(8)[0]), Is.LessThan(1.0),
            "Expected a peak at the second proteoform's monoisotopic m/z.");
    }

    private static double NearestPpm(double[] observed, double target) =>
        observed.Length == 0
            ? double.MaxValue
            : observed.Min(mz => Math.Abs(mz - target) / target * 1e6);
}
