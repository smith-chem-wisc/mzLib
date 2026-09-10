using NUnit.Framework;
using Readers;
using Readers.MaSSSimulator;
using System.IO;
using System.Linq;

namespace Test.MaSSSimulator;

[TestFixture]
public class MaSSSimulatorSpectrumFileTests
{
    [Test]
    public void ReadsSpectraAndAppliesSeparateTruth()
    {
        var spectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-test.masssim.spectra");
        var truthPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-test.rst");
        File.WriteAllText(spectrumPath, "H\tSNR\t2.8\nS\t7\t7\t501.25\nZ\t2\t1000.5\n100.1 25\n200.2 50\n");
        File.WriteAllText(truthPath, "scan:7 peptide:PEPTIDE\n");

        var file = new MaSSSimulatorSpectrumFile(spectrumPath);
        file.ApplyTruth(truthPath);

        Assert.That(file.Results, Has.Count.EqualTo(1));
        Assert.That(file[0].ScanNumber, Is.EqualTo(7));
        Assert.That(file[0].Charge, Is.EqualTo(2));
        Assert.That(file[0].PrecursorMz, Is.EqualTo(501.25));
        Assert.That(file[0].Peaks, Has.Count.EqualTo(2));
        Assert.That(file[0].PeptideSequence, Is.EqualTo("PEPTIDE"));
    }

    [Test]
    public void WriteResultsRoundTripsWithHeaderAndPrecursorMass()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-test.masssim.spectra");
        var source = new MaSSSimulatorSpectrumFile([
            new MaSSSimulatorSpectrum
            {
                ScanNumber = 7,
                PrecursorMz = 501.25,
                Charge = 2,
                PrecursorMass = 1000.5,
                Peaks = [(100.1, 25.0), (200.2, 50.0)]
            }]);

        source.WriteResults(path);
        var loaded = new MaSSSimulatorSpectrumFile(path);

        var lines = File.ReadAllLines(path);
        Assert.That(lines[0], Is.EqualTo("H\tCreationDate\t"));
        Assert.That(lines, Does.Contain("S\t7\t7\t501.25"));
        Assert.That(lines, Does.Contain("Z\t2\t1000.5"));
        Assert.That(lines, Does.Contain("100.1 25"));
        Assert.That(lines, Does.Contain("200.2 50"));
        Assert.That(loaded.Results, Has.Count.EqualTo(1));
        Assert.That(loaded[0].PrecursorMass, Is.EqualTo(1000.5));
        Assert.That(loaded[0].Peaks, Has.Count.EqualTo(2));
    }

    [Test]
    public void WriteMgfOmitsChargeWhenUnknown()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-test.mgf");
        var file = new MaSSSimulatorSpectrumFile([
            new MaSSSimulatorSpectrum
            {
                ScanNumber = 3,
                PrecursorMz = 400.5,
                Charge = 0,
                Peaks = [(100.1, 12.0)]
            }]);

        file.WriteMgf(path);
        var lines = File.ReadAllLines(path);

        Assert.That(lines, Does.Contain("BEGIN IONS"));
        Assert.That(lines, Does.Contain("PEPMASS=400.5"));
        Assert.That(lines, Does.Not.Contain("CHARGE=0+"));
        Assert.That(lines, Does.Contain("SCANS=3"));
        Assert.That(lines, Does.Contain("100.1 12"));
        Assert.That(lines, Does.Contain("END IONS"));
        Assert.That(lines.Any(line => line.StartsWith("CHARGE=")), Is.False);

        var roundTrip = new Mgf(path);
        Assert.That(roundTrip.GetAllScansList(), Has.Count.EqualTo(1));
    }

    [TestCase("S\t7\t7\n", "Invalid MaSS-Simulator spectrum line")]
    [TestCase("Z\t2\t1000.5\n", "MaSS-Simulator Z line without an S line")]
    [TestCase("100.1 25\n", "Peak line without an S line")]
    [TestCase("S\t7\t7\t501.25\n100.1 25 10\n", "Invalid MaSS-Simulator peak line")]
    [TestCase("S\tabc\t7\t501.25\n", "Invalid integer in MaSS-Simulator line")]
    [TestCase("S\t7\t7\tfoo\n", "Invalid number in MaSS-Simulator line")]
    public void LoadResultsRejectsInvalidInput(string contents, string expectedMessage)
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-invalid.masssim.spectra");
        File.WriteAllText(path, contents);

        var file = new MaSSSimulatorSpectrumFile(path);
        var ex = Assert.Throws<System.FormatException>(() => file.LoadResults());

        Assert.That(ex!.Message, Does.Contain(expectedMessage));
    }
}
