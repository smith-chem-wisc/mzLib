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
public class TestModomicsSequenceTransformation
{
    private const string ModomicsSequence = "GJACUGCBUCUA#UGAA#CA";
    private static readonly int[] ExpectedFixedPositions = [2, 8, 13, 18];

    [Test]
    public void FastaModomicsTransformationCreatesFixedModsForTopDown()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "modomics-transformation.fasta");
        File.WriteAllText(path, ">id:1|Name:test|SOterm:test|Type:t|Subtype:s|Feature:f|Species:Test\n" + ModomicsSequence + "\n");

        try
        {
            var rna = RnaDbLoader.LoadRnaFasta(
                path, true, DecoyType.None, false, out var errors).Single();

            Assert.That(errors, Is.Empty);
            Assert.That(rna.BaseSequence, Is.EqualTo("GUACUGCCUCUAGUGAAGCA"));
            Assert.That(rna.OneBasedFixedModifications.Keys, Is.EquivalentTo(ExpectedFixedPositions));
            AssertTopDownContainsOnlyAnnotatedMods(rna);
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void XmlModomicsTransformationCreatesFixedModsForTopDown()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "modomics-transformation.xml");
        File.WriteAllText(path, $"<mzLibProteinDb><entry><accession>test</accession><name>test</name><sequence length=\"20\">{ModomicsSequence}</sequence></entry></mzLibProteinDb>");

        try
        {
            var rna = RnaDbLoader.LoadRnaXML(
                path, true, DecoyType.None, false, [], [], out _, out var errors).Single();

            Assert.That(errors, Is.Empty);
            Assert.That(rna.BaseSequence, Is.EqualTo("GUACUGCCUCUAGUGAAGCA"));
            Assert.That(rna.OneBasedFixedModifications.Keys, Is.EquivalentTo(ExpectedFixedPositions));
            AssertTopDownContainsOnlyAnnotatedMods(rna);
        }
        finally
        {
            File.Delete(path);
        }
    }

    private static void AssertTopDownContainsOnlyAnnotatedMods(RNA rna)
    {
        var products = rna.Digest(new RnaDigestionParams("top-down") { MaxMods = 0 }, [], []).ToList();

        Assert.That(products, Has.Count.EqualTo(1));
        Assert.That(products[0].NumFixedMods, Is.EqualTo(ExpectedFixedPositions.Length));
        Assert.That(products[0].AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 3, 9, 14, 19 }));
    }
}
