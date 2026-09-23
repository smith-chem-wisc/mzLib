using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics;

[TestFixture]
public class TestBurkmerDatabaseRegression
{
    private static string DataPath(string fileName) =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", fileName);

    [Test]
    public void ModomicsFastaLoadsAllFiveSequences()
    {
        var rnas = RnaDbLoader.LoadRnaFasta(
            DataPath("Burkmers_ModomicsSequences.fasta"), true, DecoyType.None, false, out var errors);

        Assert.That(rnas, Has.Count.EqualTo(5));
        Assert.That(errors, Is.Empty);
    }

    [Test]
    public void ModomicsFastaSkipsUnrepresentableCodeWithWarning()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "invalid-modomics.fasta");
        File.WriteAllText(path, ">id:bad|Name:bad|SOterm:bad|Type:tRNA|Subtype:Ala|Feature:VGC|Species:standard\nUUCAAGUA\u2603UCCAGGAUAGGCU\n");

        try
        {
            var rnas = RnaDbLoader.LoadRnaFasta(path, true, DecoyType.None, false, out var errors);

            Assert.That(rnas, Is.Empty);
            Assert.That(errors.Any(error => error.Contains("bad")), Is.True);
            Assert.That(errors[0], Does.Contain("bad"));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void ModomicsFastaAndXmlProduceIdenticalFixedTopDownOligos()
    {
        var fastaRna = RnaDbLoader.LoadRnaFasta(
            DataPath("Burkmers_ModomicsSequences.fasta"), true, DecoyType.None, false, out var fastaErrors)
            .Single(rna => rna.Name == "22mer#1-Am");
        var xmlRna = RnaDbLoader.LoadRnaXML(
            DataPath("Burkemers_Modomics.xml"), true, DecoyType.None, false,
            [], [], out _, out var xmlErrors)
            .Single(rna => rna.Accession == "22mer#1-Am");

        Assert.That(fastaErrors, Is.Empty);
        Assert.That(xmlErrors, Is.Empty);
        Assert.That(fastaRna.OneBasedFixedModifications.Keys, Is.EqualTo(xmlRna.OneBasedFixedModifications.Keys));

        var fastaProduct = TopDown(fastaRna).Single();
        var xmlProduct = TopDown(xmlRna).Single();

        Assert.That(fastaProduct.FullSequence, Is.EqualTo(xmlProduct.FullSequence));
        Assert.That(fastaProduct.NumFixedMods, Is.EqualTo(1));
        Assert.That(xmlProduct.NumFixedMods, Is.EqualTo(1));
        Assert.That(fastaProduct.MonoisotopicMass, Is.EqualTo(xmlProduct.MonoisotopicMass).Within(1e-9));
    }

    [Test]
    public void PrimaryFastaIsUnmodifiedAndNormalXmlUsesVariableMods()
    {
        var primaryRna = RnaDbLoader.LoadRnaFasta(
            DataPath("Burkmers_PrimarySequence.fasta"), true, DecoyType.None, false, out var fastaErrors)
            .Single(rna => rna.Name == "22mer#1-Am");
        var normalXmlRna = RnaDbLoader.LoadRnaXML(
            DataPath("Burkemers.xml"), true, DecoyType.None, false,
            [], [], out _, out var xmlErrors)
            .Single(rna => rna.Accession == "22mer#1-Am");

        Assert.That(fastaErrors, Is.Empty);
        Assert.That(xmlErrors, Is.Empty);
        Assert.That(primaryRna.OneBasedFixedModifications, Is.Empty);
        Assert.That(normalXmlRna.OneBasedFixedModifications, Is.Empty);
        Assert.That(normalXmlRna.OneBasedPossibleLocalizedModifications, Does.ContainKey(9));

        Assert.That(TopDown(primaryRna), Has.Count.EqualTo(1));
        Assert.That(TopDown(normalXmlRna), Has.Count.EqualTo(1));
    }

    private static List<OligoWithSetMods> TopDown(RNA rna) =>
        rna.Digest(new RnaDigestionParams("top-down") { MaxMods = 0 }, [], []).ToList();
}
