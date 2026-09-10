using NUnit.Framework;
using Readers.MaSSSimulator;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

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
    public void FromPeptidesPreservesSequenceAndAccession()
    {
        var digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxPeptideLength: 50, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
        var peptide = new Protein("PEPTIDE", "P1").Digest(digestionParams, new List<Modification>(), new List<Modification>()).Single();

        var file = MaSSSimulatorPeptideFile.FromPeptides(new List<IBioPolymerWithSetMods> { peptide });

        Assert.That(file.Results.Select(p => p.Sequence), Is.EqualTo(new[] { "PEPTIDE" }));
        Assert.That(file.Results.Select(p => p.Accession), Is.EqualTo(new[] { "P1" }));
    }

    [Test]
    public void FromProteinsUsesMzLibDigestion()
    {
        var protein = new Protein("PEPTIDE", "P1");
        var digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxPeptideLength: 50, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

        var file = MaSSSimulatorPeptideFile.FromProteins(new[] { protein }, digestionParams);

        Assert.That(file.Results, Is.Not.Empty);
        Assert.That(file.Results.All(p => !string.IsNullOrWhiteSpace(p.Sequence)), Is.True);
        Assert.That(file.Results.All(p => p.Accession == "P1"), Is.True);
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

    [Test]
    public void FromPeptidesRejectsModifiedSequences()
    {
        var oxidation = Mods.GetModification("Oxidation on M");
        Assert.That(oxidation, Is.Not.Null);
        var modified = new PeptideWithSetModifications($"M[{oxidation!.ModificationType}:{oxidation.IdWithMotif}]");

        var ex = Assert.Throws<NotSupportedException>(() => MaSSSimulatorPeptideFile.FromPeptides(new List<IBioPolymerWithSetMods> { modified }));

        Assert.That(ex!.Message, Does.Contain("unmodified peptides only"));
    }
}
