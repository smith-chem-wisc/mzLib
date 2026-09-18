using NUnit.Framework;
using Readers.MaSSSimulator;

namespace Test.MaSSSimulator;

[TestFixture]
public class MaSSSimulatorModelsTests
{
    [Test]
    public void PeptideRecordUsesValueEquality()
    {
        var left = new MaSSSimulatorPeptide("PEPTIDE", "P1");
        var right = new MaSSSimulatorPeptide("PEPTIDE", "P1");
        var different = new MaSSSimulatorPeptide("PEPTIDE", "P2");

        Assert.That(left, Is.EqualTo(right));
        Assert.That(left.GetHashCode(), Is.EqualTo(right.GetHashCode()));
        Assert.That(left, Is.Not.EqualTo(different));
    }

    [Test]
    public void SpectrumDefaultsAreEmptyOrUnset()
    {
        var spectrum = new MaSSSimulatorSpectrum();

        Assert.That(spectrum.ScanNumber, Is.EqualTo(0));
        Assert.That(spectrum.PrecursorMz, Is.EqualTo(0));
        Assert.That(spectrum.Charge, Is.EqualTo(0));
        Assert.That(spectrum.PrecursorMass, Is.Null);
        Assert.That(spectrum.Peaks, Is.Empty);
        Assert.That(spectrum.PeptideSequence, Is.Null);
        Assert.That(spectrum.Accession, Is.Null);
    }
}
