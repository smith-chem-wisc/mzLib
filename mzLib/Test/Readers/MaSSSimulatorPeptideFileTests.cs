using NUnit.Framework;
using Readers.MaSSSimulator;
using System.IO;
using System.Linq;

namespace Test.MaSSSimulator;

[TestFixture]
public class MaSSSimulatorPeptideFileTests
{
    [Test]
    public void WritesAndReadsSimulatorPeptideList()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-test.masssim.peptides");
        var source = new MaSSSimulatorPeptideFile([
            new("PEPTIDE", "P1"),
            new("AGAIN", "P2")]);

        source.WriteResults(path);
        var loaded = new MaSSSimulatorPeptideFile(path);

        Assert.That(loaded.Results.Select(p => p.Sequence), Is.EqualTo(new[] { "PEPTIDE", "AGAIN" }));
        Assert.That(File.ReadLines(path).First(), Is.EqualTo("peptide"));
    }

    [Test]
    public void WritesTruthTableWithStableScanNumbers()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "masssim-test.truth.tsv");
        var source = new MaSSSimulatorPeptideFile([new("PEPTIDE", "P1"), new("AGAIN")]);

        source.WriteTruthTable(path);

        var lines = File.ReadAllLines(path);
        Assert.That(lines[0], Is.EqualTo("scan\tpeptide\taccession"));
        Assert.That(lines[1], Is.EqualTo("1\tPEPTIDE\tP1"));
        Assert.That(lines[2], Is.EqualTo("2\tAGAIN\t"));
    }
}
