using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics;

/// <summary>
/// These databases encode the same five oligos as unmodified FASTA, fixed-mod
/// MODOMICS FASTA/XML, and variable-mod mzLib XML. The normal XML top-down
/// products must contain the unmodified and MODOMICS-fixed products as a
/// superset, with equivalent chemistry and fragments.
/// </summary>
[TestFixture]
public class TestBurkmerDatabaseRegression
{
    private sealed record ExpectedCase(string Name, string BaseSequence, int? ModifiedPosition);

    private static readonly ExpectedCase[] Cases =
    [
        new("22mer#1-unmod", "UUCAAGUAAUCCAGGAUAGGCU", null),
        new("22mer#1-Am", "UUCAAGUAAUCCAGGAUAGGCU", 9),
        new("22mer#1-m6A", "UUCAAGUAAUCCAGGAUAGGCU", 9),
        new("22mer#1-m6Am", "UUCAAGUAAUCCAGGAUAGGCU", 9),
        new("22mer#2-Am", "UCCCUGAGACCCUAACUUGUGA", 15),
    ];

    private static string DataPath(string fileName) =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", fileName);

    [Test]
    public void AllFiveRepresentationsLoadWithoutWarnings()
    {
        var primary = LoadFasta("Burkmers_PrimarySequence.fasta", out var primaryErrors);
        var modomicsFasta = LoadFasta("Burkmers_ModomicsSequences.fasta", out var modomicsFastaErrors);
        var normalXml = LoadXml("Burkemers.xml", out var normalXmlErrors);
        var modomicsXml = LoadXml("Burkemers_Modomics.xml", out var modomicsXmlErrors);

        Assert.That(primaryErrors, Is.Empty);
        Assert.That(modomicsFastaErrors, Is.Empty);
        Assert.That(normalXmlErrors, Is.Empty);
        Assert.That(modomicsXmlErrors, Is.Empty);

        AssertCaseRecords(primary, false, false);
        AssertCaseRecords(modomicsFasta, true, true);
        AssertCaseRecords(normalXml, false, false);
        AssertCaseRecords(modomicsXml, true, true);
    }

    [Test]
    public void TopDownProductsHaveExpectedFixedVariableAndUnmodifiedForms()
    {
        var primary = LoadFasta("Burkmers_PrimarySequence.fasta", out _).ToDictionary(rna => rna.Name);
        var modomicsFasta = LoadFasta("Burkmers_ModomicsSequences.fasta", out _).ToDictionary(rna => rna.Name);
        var normalXml = LoadXml("Burkemers.xml", out _).ToDictionary(rna => rna.Accession);
        var modomicsXml = LoadXml("Burkemers_Modomics.xml", out _).ToDictionary(rna => rna.Accession);

        foreach (var expected in Cases)
        {
            var primaryProduct = TopDown(primary[expected.Name], 0).Single();
            var fixedFastaProduct = TopDown(modomicsFasta[expected.Name], 0).Single();
            var fixedXmlProduct = TopDown(modomicsXml[expected.Name], 0).Single();
            var variableProducts = TopDown(normalXml[expected.Name], 1);

            if (expected.ModifiedPosition is null)
                AssertProductEquivalent(primaryProduct, fixedFastaProduct, true);
            AssertProductEquivalent(fixedFastaProduct, fixedXmlProduct, expected.ModifiedPosition is null);
            Assert.That(variableProducts.Any(product => EquivalentProducts(product, primaryProduct)), Is.True,
                expected.Name);
            Assert.That(variableProducts.Any(product => EquivalentProducts(product, fixedFastaProduct)), Is.True,
                expected.Name);

            AssertFragmentsEquivalent(primaryProduct, variableProducts.Single(product => EquivalentProducts(product, primaryProduct)));
            AssertFragmentsEquivalent(fixedFastaProduct, variableProducts.Single(product => EquivalentProducts(product, fixedFastaProduct)));
        }
    }

    [Test]
    public void InvalidModomicsSequenceIsSkippedWithRecordWarning()
    {
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "invalid-modomics.fasta");
        File.WriteAllText(path, ">id:bad|Name:bad|SOterm:bad|Type:tRNA|Subtype:Ala|Feature:VGC|Species:standard\nUUCAAGUA\u2603UCCAGGAUAGGCU\n");

        try
        {
            var rnas = RnaDbLoader.LoadRnaFasta(path, true, DecoyType.None, false, out var errors);
            Assert.That(rnas, Is.Empty);
            Assert.That(errors.Any(error => error.Contains("bad")), Is.True);
        }
        finally
        {
            File.Delete(path);
        }
    }

    private static List<RNA> LoadFasta(string fileName, out List<string> errors) =>
        RnaDbLoader.LoadRnaFasta(DataPath(fileName), true, DecoyType.None, false, out errors);

    private static List<RNA> LoadXml(string fileName, out List<string> errors) =>
        RnaDbLoader.LoadRnaXML(DataPath(fileName), true, DecoyType.None, false,
            [], [], out _, out errors);

    private static void AssertCaseRecords(IEnumerable<RNA> rnas, bool expectFixed, bool expectNoPossibleMods)
    {
        var records = rnas.ToDictionary(rna => rna.Name ?? rna.Accession);
        Assert.That(records.Keys, Is.EquivalentTo(Cases.Select(expected => expected.Name)));

        foreach (var expected in Cases)
        {
            var rna = records[expected.Name];
            Assert.That(rna.BaseSequence, Is.EqualTo(expected.BaseSequence), expected.Name);
            Assert.That(rna.OneBasedFixedModifications.Count, Is.EqualTo(expectFixed && expected.ModifiedPosition.HasValue ? 1 : 0), expected.Name);
            if (expected.ModifiedPosition.HasValue && expectFixed)
                Assert.That(rna.OneBasedFixedModifications.Keys, Does.Contain(expected.ModifiedPosition.Value), expected.Name);
            if (expectNoPossibleMods)
                Assert.That(rna.OneBasedPossibleLocalizedModifications, Is.Empty, expected.Name);
        }
    }

    private static List<OligoWithSetMods> TopDown(RNA rna, int maxMods) =>
        rna.Digest(new RnaDigestionParams("top-down") { MaxMods = maxMods }, [], []).ToList();

    private static bool EquivalentProducts(OligoWithSetMods left, OligoWithSetMods right) =>
        left.BaseSequence == right.BaseSequence
        && left.AllModsOneIsNterminus.Keys.SequenceEqual(right.AllModsOneIsNterminus.Keys)
        && left.AllModsOneIsNterminus.All(entry =>
            right.AllModsOneIsNterminus[entry.Key].MonoisotopicMass == entry.Value.MonoisotopicMass);

    private static void AssertProductEquivalent(OligoWithSetMods expected, OligoWithSetMods actual, bool noMods)
    {
        Assert.That(actual.BaseSequence, Is.EqualTo(expected.BaseSequence));
        Assert.That(actual.AllModsOneIsNterminus.Keys, Is.EquivalentTo(expected.AllModsOneIsNterminus.Keys));
        Assert.That(actual.MonoisotopicMass, Is.EqualTo(expected.MonoisotopicMass).Within(1e-9));
        if (noMods)
            Assert.That(actual.AllModsOneIsNterminus, Is.Empty);
    }

    private static void AssertFragmentsEquivalent(OligoWithSetMods expected, OligoWithSetMods actual)
    {
        foreach (var dissociationType in new[] { DissociationType.CID, DissociationType.HCD })
        {
            var expectedProducts = new List<Product>();
            var actualProducts = new List<Product>();
            expected.Fragment(dissociationType, FragmentationTerminus.Both, expectedProducts);
            actual.Fragment(dissociationType, FragmentationTerminus.Both, actualProducts);

            Assert.That(actualProducts.Count, Is.EqualTo(expectedProducts.Count), dissociationType.ToString());
            var expectedSorted = expectedProducts.OrderBy(product => product.ProductType).ThenBy(product => product.Terminus).ThenBy(product => product.FragmentNumber).ToList();
            var actualSorted = actualProducts.OrderBy(product => product.ProductType).ThenBy(product => product.Terminus).ThenBy(product => product.FragmentNumber).ToList();
            for (var i = 0; i < expectedSorted.Count; i++)
            {
                Assert.That(actualSorted[i].ProductType, Is.EqualTo(expectedSorted[i].ProductType));
                Assert.That(actualSorted[i].Terminus, Is.EqualTo(expectedSorted[i].Terminus));
                Assert.That(actualSorted[i].FragmentNumber, Is.EqualTo(expectedSorted[i].FragmentNumber));
                Assert.That(actualSorted[i].NeutralMass, Is.EqualTo(expectedSorted[i].NeutralMass).Within(1e-9));
            }
        }
    }
}
