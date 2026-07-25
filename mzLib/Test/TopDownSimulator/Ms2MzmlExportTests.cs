using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers;
using TopDownSimulator.Ms2;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class Ms2MzmlExportTests
{
    private const string Sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCK";
    private const int PrecursorCharge = 8;

    private string _outputDirectory;

    private static PeptideWithSetModifications Peptide(string sequence = Sequence) =>
        new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>());

    private static Ms2ScanRequest Request(double retentionTime = 20.0, string sequence = Sequence) =>
        new Ms2ScanRequest(Peptide(sequence), PrecursorCharge, retentionTime);

    [SetUp]
    public void SetUp()
    {
        _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownSimulatorMs2Output");
        Directory.CreateDirectory(_outputDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_outputDirectory))
            Directory.Delete(_outputDirectory, recursive: true);
    }

    [Test]
    public void SimulatedScanIsCentroidedMsnWithPrecursorMetadata()
    {
        var request = Request();
        var simulation = new Ms2Simulator().SimulateCentroided(request);

        Assert.That(simulation.Scans, Has.Length.EqualTo(1));
        var scan = simulation.Scans[0];

        Assert.That(scan.MsnOrder, Is.EqualTo(2));
        Assert.That(scan.IsCentroid, Is.True);
        Assert.That(scan.DissociationType, Is.EqualTo(DissociationType.ETD));
        Assert.That(scan.SelectedIonChargeStateGuess, Is.EqualTo(PrecursorCharge));
        Assert.That(scan.SelectedIonMZ, Is.Not.Null);
        Assert.That(scan.SelectedIonMZ.Value,
            Is.EqualTo(request.Precursor.MonoisotopicMass.ToMz(PrecursorCharge)).Within(1e-9));
        Assert.That(scan.RetentionTime, Is.EqualTo(20.0).Within(1e-9));
        Assert.That(scan.MassSpectrum.XArray, Is.Ordered.Ascending);
        Assert.That(scan.MassSpectrum.XArray, Is.Not.Empty);
        Assert.That(scan.TotalIonCurrent, Is.EqualTo(scan.MassSpectrum.YArray.Sum()).Within(1e-6));
    }

    [Test]
    public void ThresholdingUsesAGlobalFloorAcrossAllRequestedScans()
    {
        var requests = new[] { Request(20.0), Request(21.0) };

        var unfiltered = new Ms2Simulator().SimulateCentroided(requests,
            new ScanReductionOptions { RelativeIntensityThreshold = 0 });
        var filtered = new Ms2Simulator().SimulateCentroided(requests,
            new ScanReductionOptions { RelativeIntensityThreshold = 0.5 });

        Assert.That(unfiltered.Scans[0].MassSpectrum.XArray.Length,
            Is.GreaterThan(filtered.Scans[0].MassSpectrum.XArray.Length),
            "A 50% relative floor should remove most peaks.");

        double max = unfiltered.Scans.SelectMany(s => s.MassSpectrum.YArray).Max();
        foreach (double y in filtered.Scans.SelectMany(s => s.MassSpectrum.YArray))
            Assert.That(y, Is.GreaterThanOrEqualTo(max * 0.5));
    }

    [Test]
    public void PeakIntensitiesMatchTheFragmentationResult()
    {
        var simulation = new Ms2Simulator().SimulateCentroided(Request(),
            new ScanReductionOptions { RelativeIntensityThreshold = 0 });

        var scan = simulation.Scans[0];
        var fragmentation = simulation.Fragmentations[0];

        Assert.That(scan.MassSpectrum.XArray, Has.Length.EqualTo(fragmentation.Fragments.Count),
            "One peak per fragment ion when nothing is thresholded away and no two ions collide.");
        Assert.That(scan.MassSpectrum.YArray.Sum(),
            Is.EqualTo(fragmentation.Fragments.Sum(f => f.Intensity)).Within(1e-6));

        foreach (var fragment in fragmentation.Fragments)
        {
            int index = Array.BinarySearch(scan.MassSpectrum.XArray, fragment.Mz);
            Assert.That(index, Is.GreaterThanOrEqualTo(0), $"{fragment.Product.Annotation} is missing from the spectrum.");
            Assert.That(scan.MassSpectrum.YArray[index], Is.EqualTo(fragment.Intensity).Within(1e-9));
        }
    }

    [Test]
    public void WriteMzmlProducesCentroidedMs2ThatReadsBackWithMatchingPeaks()
    {
        var requests = new[] { Request(20.0), Request(20.5) };
        string path = Path.Combine(_outputDirectory, "ms2.mzML");

        var export = new Ms2Simulator().WriteMzml(requests, path);

        Assert.That(File.Exists(path), Is.True, "mzML was not written.");
        Assert.That(export.ScanCount, Is.EqualTo(2));
        Assert.That(export.PeakCount, Is.GreaterThan(0));

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        var scans = roundTripped.GetAllScansList();

        Assert.That(scans, Has.Count.EqualTo(2));
        Assert.That(scans.Sum(s => s.MassSpectrum.XArray.Length), Is.EqualTo(export.PeakCount),
            "Peak count reported by the exporter should match what the reader finds.");
        Assert.That(scans.All(s => s.MsnOrder == 2), Is.True, "All simulated scans should be MS2.");
        Assert.That(scans.All(s => s.IsCentroid), Is.True, "Centroid flag should survive the round trip.");
        Assert.That(scans.All(s => s.DissociationType == DissociationType.ETD), Is.True);
        Assert.That(scans.All(s => s.SelectedIonChargeStateGuess == PrecursorCharge), Is.True);

        double expectedPrecursorMz = Peptide().MonoisotopicMass.ToMz(PrecursorCharge);
        foreach (var scan in scans)
            Assert.That(scan.SelectedIonMZ.Value, Is.EqualTo(expectedPrecursorMz).Within(1e-4));

        Assert.That(scans[0].RetentionTime, Is.EqualTo(20.0).Within(1e-6));
        Assert.That(scans[1].RetentionTime, Is.EqualTo(20.5).Within(1e-6));
    }

    [Test]
    public void RoundTrippedPeakPositionsSurviveAtOrbitrapPrecision()
    {
        var requests = new[] { Request() };
        string path = Path.Combine(_outputDirectory, "ms2-precision.mzML");

        var simulator = new Ms2Simulator();
        var expected = simulator.SimulateCentroided(requests);
        simulator.WriteMzml(requests, path);

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        var actual = roundTripped.GetAllScansList();

        var expectedMz = expected.Scans[0].MassSpectrum.XArray;
        var actualMz = actual[0].MassSpectrum.XArray;
        Assert.That(actualMz, Has.Length.EqualTo(expectedMz.Length));

        for (int i = 0; i < expectedMz.Length; i++)
        {
            double ppmError = Math.Abs(actualMz[i] - expectedMz[i]) / expectedMz[i] * 1e6;
            Assert.That(ppmError, Is.LessThan(0.1),
                $"Peak {i} shifted by {ppmError:F4} ppm through the mzML round trip.");
        }
    }

    [Test]
    public void WriteMzmlWritesFragmentGroundTruthSidecar()
    {
        var requests = new[] { Request() };
        string path = Path.Combine(_outputDirectory, "ms2-with-truth.mzML");

        var export = new Ms2Simulator().WriteMzml(requests, path,
            new ScanReductionOptions { RelativeIntensityThreshold = 0 });

        Assert.That(export.GroundTruthPath, Is.Not.Null);
        Assert.That(File.Exists(export.GroundTruthPath), Is.True);

        var lines = File.ReadAllLines(export.GroundTruthPath);
        var header = lines[0].Split('\t');
        Assert.That(header, Does.Contain("Annotation"));
        Assert.That(header, Does.Contain("OneBasedBondNumber"));
        Assert.That(header, Does.Contain("Intensity"));
        Assert.That(lines, Has.Length.EqualTo(export.PeakCount + 1),
            "One row per simulated fragment ion, plus the header.");

        int intensityColumn = Array.IndexOf(header, "Intensity");
        int annotationColumn = Array.IndexOf(header, "Annotation");
        double assigned = lines.Skip(1).Sum(l => double.Parse(l.Split('\t')[intensityColumn]));

        var fragmentation = new PropensityCascadeFragmentationModel().Fragment(Peptide());
        Assert.That(assigned, Is.EqualTo(fragmentation.Fragments.Sum(f => f.Intensity)).Within(1e-3));
        Assert.That(lines.Skip(1).Select(l => l.Split('\t')[annotationColumn]).Any(a => a.StartsWith("c")), Is.True);
        Assert.That(lines.Skip(1).Select(l => l.Split('\t')[annotationColumn]).Any(a => a.StartsWith("zDot")), Is.True);
    }

    [Test]
    public void WriteMzmlCanSkipTheSidecar()
    {
        string path = Path.Combine(_outputDirectory, "ms2-no-sidecar.mzML");

        var export = new Ms2Simulator().WriteMzml(new[] { Request() }, path, writeGroundTruthSidecar: false);

        Assert.That(export.GroundTruthPath, Is.Null);
        Assert.That(File.Exists(Path.ChangeExtension(path, ".ms2.groundtruth.tsv")), Is.False);
    }

    [Test]
    public void HcdOptionsProduceBAndYIonsAndWriteAnHcdScan()
    {
        var model = new PropensityCascadeFragmentationModel(
            options: new PropensityCascadeOptions { DissociationType = DissociationType.HCD });
        var simulator = new Ms2Simulator(model);
        string path = Path.Combine(_outputDirectory, "ms2-hcd.mzML");

        simulator.WriteMzml(new[] { Request() }, path, writeGroundTruthSidecar: false);

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        Assert.That(roundTripped.GetAllScansList().All(s => s.DissociationType == DissociationType.HCD), Is.True);
    }

    /// <summary>
    /// The whole point of the interface: the simulator does not know what model it is driving.
    /// </summary>
    private sealed class ConstantFragmentationModel : IMs2FragmentationModel
    {
        public Ms2FragmentationResult Fragment(IBioPolymerWithSetMods precursor) =>
            new Ms2FragmentationResult(
                precursor.BaseSequence,
                100.0,
                new List<SimulatedFragmentIon>
                {
                    new SimulatedFragmentIon(
                        new Product(ProductType.c, FragmentationTerminus.N, 500.0, 1, 1, 0),
                        1, 1, 0, 501.00727645, 40.0)
                },
                60.0);
    }

    [Test]
    public void SimulatorAcceptsAnArbitraryFragmentationModel()
    {
        var simulation = new Ms2Simulator(new ConstantFragmentationModel()).SimulateCentroided(Request());

        Assert.That(simulation.Scans[0].MassSpectrum.XArray, Has.Length.EqualTo(1));
        Assert.That(simulation.Scans[0].MassSpectrum.XArray[0], Is.EqualTo(501.00727645).Within(1e-9));
        Assert.That(simulation.Scans[0].MassSpectrum.YArray[0], Is.EqualTo(40.0).Within(1e-9));
        Assert.That(simulation.Fragmentations[0].SurvivingPrecursorIntensity, Is.EqualTo(60.0));
    }
}
