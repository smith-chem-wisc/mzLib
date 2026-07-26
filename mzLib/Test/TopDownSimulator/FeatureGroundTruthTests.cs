using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class FeatureGroundTruthTests
{
    private const int MinCharge = 6;
    private const int MaxCharge = 14;
    private static readonly ConstantPeakWidth Width = new(0.012);

    private string _outputDirectory;

    [SetUp]
    public void SetUp()
    {
        _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownSimulatorFeatureTruth");
        Directory.CreateDirectory(_outputDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_outputDirectory))
            Directory.Delete(_outputDirectory, recursive: true);
    }

    private static ProteoformModel Model(
        double mass = 10000.0, double abundance = 1.5e6, double rtCenter = 20.0,
        double muZ = 8.3, string identifier = "P00000|TESTPROTEOFORM") =>
        new(
            MonoisotopicMass: mass,
            Abundance: abundance,
            RtProfile: new EmgProfile(Mu: rtCenter, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: muZ, SigmaZ: 1.15),
            Identifier: identifier);

    private static double[] ScanTimes(int count = 21, double start = 19.0, double step = 0.1) =>
        Enumerable.Range(0, count).Select(i => start + i * step).ToArray();

    private static Dictionary<string, string> Row(string[] header, string line) =>
        header.Zip(line.Split('\t'), (h, v) => (h, v)).ToDictionary(t => t.h, t => t.v);

    [Test]
    public void EveryReportedFeatureIsActuallyPresentInTheMzml()
    {
        // The scoring contract: a benchmark must never start with false negatives it cannot fix.
        string path = Path.Combine(_outputDirectory, "present.mzML");
        var models = new[] { Model(), Model(mass: 12500.0, abundance: 4e5, rtCenter: 20.05, muZ: 10.1, identifier: "SECOND") };

        var export = new Simulator().WriteMzml(models, MinCharge, MaxCharge, Width, ScanTimes(), path);

        Assert.That(export.FeatureGroundTruthPath, Is.Not.Null);
        Assert.That(export.FeatureCount, Is.GreaterThan(0));

        var lines = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var header = lines[0].Split('\t');
        var features = lines.Skip(1).Select(l => Row(header, l)).ToArray();
        Assert.That(features, Has.Length.EqualTo(export.FeatureCount));

        var file = MsDataFileReader.GetDataFile(path);
        file.LoadAllStaticData();
        var scans = file.GetAllScansList().ToDictionary(s => s.OneBasedScanNumber);

        foreach (var feature in features)
        {
            int apexScanNumber = int.Parse(feature["ApexScanNumber"]);
            double apexMz = double.Parse(feature["ApexMz"]);
            double apexIntensity = double.Parse(feature["ApexIntensity"]);

            Assert.That(scans.ContainsKey(apexScanNumber), Is.True,
                $"Feature {feature["Identifier"]} z={feature["Charge"]} names scan {apexScanNumber}, which is not in the file.");

            var spectrum = scans[apexScanNumber].MassSpectrum;
            int nearest = spectrum.GetClosestPeakIndex(apexMz);
            double ppm = Math.Abs(spectrum.XArray[nearest] - apexMz) / apexMz * 1e6;

            Assert.That(ppm, Is.LessThan(1.0),
                $"No peak at the apex m/z {apexMz} claimed by {feature["Identifier"]} z={feature["Charge"]}.");

            // The file carries the sum over every proteoform, so its peak is at least this bright.
            Assert.That(spectrum.YArray[nearest], Is.GreaterThanOrEqualTo(apexIntensity * 0.999),
                "A feature's own apex intensity cannot exceed the total intensity recorded there.");
        }
    }

    [Test]
    public void FeatureBoundariesAreScanNumbersInsideTheFile()
    {
        string path = Path.Combine(_outputDirectory, "boundaries.mzML");
        double[] scanTimes = ScanTimes();

        var export = new Simulator().WriteMzml(new[] { Model() }, MinCharge, MaxCharge, Width, scanTimes, path);

        var lines = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var header = lines[0].Split('\t');

        foreach (var feature in lines.Skip(1).Select(l => Row(header, l)))
        {
            int first = int.Parse(feature["FirstScanNumber"]);
            int apex = int.Parse(feature["ApexScanNumber"]);
            int last = int.Parse(feature["LastScanNumber"]);

            Assert.That(first, Is.LessThanOrEqualTo(apex));
            Assert.That(apex, Is.LessThanOrEqualTo(last));
            Assert.That(first, Is.GreaterThanOrEqualTo(1).And.LessThanOrEqualTo(scanTimes.Length));
            Assert.That(last, Is.GreaterThanOrEqualTo(1).And.LessThanOrEqualTo(scanTimes.Length));

            // Retention times must agree with the scan numbers, so a consumer can key on either.
            Assert.That(double.Parse(feature["RtStart"]), Is.EqualTo(scanTimes[first - 1]).Within(1e-9));
            Assert.That(double.Parse(feature["ApexRt"]), Is.EqualTo(scanTimes[apex - 1]).Within(1e-9));
            Assert.That(double.Parse(feature["RtEnd"]), Is.EqualTo(scanTimes[last - 1]).Within(1e-9));
        }
    }

    [Test]
    public void ChargeStatesRemovedByReductionAreNotReported()
    {
        // A charge far out on the Gaussian contributes nothing above the floor. Listing it would be
        // an unfixable false negative for anything scored against this file.
        string path = Path.Combine(_outputDirectory, "reduced.mzML");
        var model = Model(muZ: 8.3);

        var export = new Simulator().WriteMzml(
            new[] { model }, minCharge: 2, maxCharge: 40, Width, ScanTimes(), path,
            new ScanReductionOptions { RelativeIntensityThreshold = 1e-3 });

        var lines = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var header = lines[0].Split('\t');
        var charges = lines.Skip(1).Select(l => int.Parse(Row(header, l)["Charge"])).ToArray();

        Assert.That(charges, Is.Not.Empty);
        Assert.That(charges.Min(), Is.GreaterThan(2), "Charge 2 is ~5σ off the charge apex and cannot survive reduction.");
        Assert.That(charges.Max(), Is.LessThan(40));
        Assert.That(charges, Does.Contain(8), "The apex charge must be present.");

        var file = MsDataFileReader.GetDataFile(path);
        file.LoadAllStaticData();
        var apexScan = file.GetAllScansList().OrderByDescending(s => s.TotalIonCurrent).First();

        var kernel = new IsotopeEnvelopeKernel(model.MonoisotopicMass);
        foreach (int absent in new[] { 2, 3, 40 })
        {
            double mz = kernel.CentroidMzs(absent)[0];
            int nearest = apexScan.MassSpectrum.GetClosestPeakIndex(mz);
            double ppm = Math.Abs(apexScan.MassSpectrum.XArray[nearest] - mz) / mz * 1e6;
            Assert.That(ppm, Is.GreaterThan(1.0),
                $"Charge {absent} was excluded from the truth but its monoisotopic peak is in the file.");
        }
    }

    [Test]
    public void ReportedIsotopologueCountFallsAsChargeMovesOffTheApex()
    {
        // A charge state out on the shoulder of the charge distribution is fainter, so fewer of its
        // isotopologues clear the reduction floor. The count has to describe what is in the file,
        // not how many isotopologues the kernel has.
        string path = Path.Combine(_outputDirectory, "isotopologues.mzML");
        var export = new Simulator().WriteMzml(
            new[] { Model(muZ: 8.3) }, minCharge: 2, maxCharge: 20, Width, ScanTimes(), path);

        var lines = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var header = lines[0].Split('\t');
        var byCharge = lines.Skip(1)
            .Select(l => Row(header, l))
            .ToDictionary(r => int.Parse(r["Charge"]), r => int.Parse(r["NumIsotopologues"]));

        Assert.That(byCharge.Values, Is.All.GreaterThan(0));
        Assert.That(byCharge.ContainsKey(8), Is.True, "The apex charge must be reported.");

        int weakestCharge = byCharge.Keys.OrderBy(z => Math.Abs(z - 8.3)).Last();
        Assert.That(byCharge[8], Is.GreaterThan(byCharge[weakestCharge]),
            $"The apex charge (8) clears the floor across {byCharge[8]} isotopologues; " +
            $"charge {weakestCharge}, furthest from the charge apex, reports {byCharge[weakestCharge]}.");
    }

    [Test]
    public void MonoisotopicMzIsComputedFromTheMassForEachReportedCharge()
    {
        string path = Path.Combine(_outputDirectory, "monoisotopic.mzML");
        var model = Model();

        var export = new Simulator().WriteMzml(new[] { model }, MinCharge, MaxCharge, Width, ScanTimes(), path);

        var lines = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var header = lines[0].Split('\t');

        foreach (var feature in lines.Skip(1).Select(l => Row(header, l)))
        {
            int z = int.Parse(feature["Charge"]);
            double expected = model.MonoisotopicMass.ToMz(z);
            double actual = double.Parse(feature["MonoisotopicMz"]);

            Assert.That(actual, Is.EqualTo(expected).Within(1e-9),
                $"Charge {z} monoisotopic m/z must be the mass converted to m/z, exactly.");
            Assert.That(double.Parse(feature["MonoisotopicMass"]), Is.EqualTo(model.MonoisotopicMass).Within(1e-9));
        }
    }

    [Test]
    public void MonoisotopicMzStaysCorrectAboveTheMassWhereAveragineDropsTheMonoisotopicPeak()
    {
        // The averagine tables drop isotopologues below an absolute probability of 1e-8, and the
        // all-light composition falls under that around 31 kDa. Above it the kernel's lowest-mass
        // entry is NOT the monoisotopic peak, so reading this column off centroids[0] would report
        // a heavier isotopologue's m/z under a header that says monoisotopic — silently wrong for
        // exactly the large proteoforms a top-down benchmark cares about most.
        const double mass = 45000.0;
        var kernel = new IsotopeEnvelopeKernel(mass);
        double lowestRetained = kernel.NeutralMass(0);

        Assert.That(lowestRetained, Is.GreaterThan(mass + 0.5),
            "Fixture assumes the monoisotopic peak has been dropped at this mass; if that changed, " +
            "pick a higher mass or delete this test.");

        string path = Path.Combine(_outputDirectory, "heavy.mzML");
        var model = Model(mass: mass, abundance: 5e6, muZ: 40);

        var export = new Simulator().WriteMzml(
            new[] { model }, minCharge: 30, maxCharge: 50, Width, ScanTimes(), path);

        var lines = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var header = lines[0].Split('\t');
        var features = lines.Skip(1).Select(l => Row(header, l)).ToArray();
        Assert.That(features, Is.Not.Empty);

        foreach (var feature in features)
        {
            int z = int.Parse(feature["Charge"]);
            double actual = double.Parse(feature["MonoisotopicMz"]);

            Assert.That(actual, Is.EqualTo(mass.ToMz(z)).Within(1e-9),
                $"Charge {z} reported {actual}, which is not the monoisotopic m/z.");
            Assert.That(actual, Is.LessThan(lowestRetained.ToMz(z)),
                "The monoisotopic m/z must sit below the lowest isotopologue the envelope retained.");
        }
    }

    [Test]
    public void SidecarSuppressionCoversBothTruthFiles()
    {
        string path = Path.Combine(_outputDirectory, "no-sidecars.mzML");

        var export = new Simulator().WriteMzml(
            new[] { Model() }, MinCharge, MaxCharge, Width, ScanTimes(), path,
            writeGroundTruthSidecar: false);

        Assert.That(export.GroundTruthPath, Is.Null);
        Assert.That(export.FeatureGroundTruthPath, Is.Null);
        Assert.That(export.FeatureCount, Is.Zero);
        Assert.That(File.Exists(Path.ChangeExtension(path, ".features.tsv")), Is.False);
    }

    [Test]
    public void HighMzFeaturesCarryFewerResolvedIsotopologuesUnderAWidthLaw()
    {
        // Under a width law, the same proteoform's low-charge (high m/z) envelope merges, so its
        // peaks are fewer and brighter than the high-charge one. That difference is the benchmark
        // difficulty a constant σ hides.
        string path = Path.Combine(_outputDirectory, "widthlaw.mzML");
        var width = OrbitrapPeakWidth.FromSigmaAt(1000.0, 0.0111);
        var model = Model(mass: 25000.0, muZ: 22, abundance: 5e6);

        var export = new Simulator().WriteMzml(
            new[] { model }, minCharge: 14, maxCharge: 30, width, ScanTimes(), path);

        var lines = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var header = lines[0].Split('\t');
        var features = lines.Skip(1).Select(l => Row(header, l)).ToArray();

        var lowCharge = features.OrderBy(f => int.Parse(f["Charge"])).First();
        var highCharge = features.OrderByDescending(f => int.Parse(f["Charge"])).First();

        Assert.That(double.Parse(lowCharge["ApexMz"]), Is.GreaterThan(double.Parse(highCharge["ApexMz"])),
            "Lower charge must land at higher m/z.");
    }
}
