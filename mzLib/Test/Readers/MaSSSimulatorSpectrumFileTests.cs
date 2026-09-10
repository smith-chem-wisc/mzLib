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
    public void WritesCasanovoCompatibleMgf()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-test.mgf");
        var file = new MaSSSimulatorSpectrumFile([
            new MaSSSimulatorSpectrum
            {
                ScanNumber = 3,
                PrecursorMz = 400.5,
                Charge = 2,
                Peaks = [(100.1, 12.0)]
            }]);

        file.WriteMgf(path);
        var lines = File.ReadAllLines(path);

        Assert.That(lines, Does.Contain("BEGIN IONS"));
        Assert.That(lines, Does.Contain("PEPMASS=400.5"));
        Assert.That(lines, Does.Contain("CHARGE=2+"));
        Assert.That(lines, Does.Contain("100.1 12"));
        Assert.That(lines, Does.Contain("END IONS"));

        var roundTrip = new Mgf(path);
        Assert.That(roundTrip.GetAllScansList(), Has.Count.EqualTo(1));
    }
}
