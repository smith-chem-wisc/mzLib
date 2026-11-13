using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// <see cref="ProteinDbRetriever.UniprotProteomeCatalog"/> keeps the columns of UniProt's proteome catalogue that
/// <see cref="ProteinDbRetriever.UniprotProteomesList"/> has always discarded, most importantly the organism's NCBI
/// taxonomy ID, which is the only link from a proteome to a taxon.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class UniProtProteomeCatalogTests
{
    private string _storageDirectory;

    [SetUp]
    public void SetUp()
    {
        _storageDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "UniProtProteomeCatalogTests_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_storageDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_storageDirectory))
            Directory.Delete(_storageDirectory, recursive: true);
    }

    private string Write(string name, string content)
    {
        string path = Path.Combine(_storageDirectory, name);
        File.WriteAllText(path, content);
        return path;
    }

    private static string Fixture(string name) =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", name);

    [Test]
    [TestCase("b_availableUniProtProteomes.txt")]
    [TestCase("b_availableUniProtProteomes.txt.gz")]
    [TestCase("b_availableUniProtProteomes.zip")]
    public void ReadsEveryColumnOfADownloadedCatalogue(string fileName)
    {
        IReadOnlyDictionary<string, UniProtProteomeEntry> catalog = ProteinDbRetriever.UniprotProteomeCatalog(Fixture(fileName));

        UniProtProteomeEntry human = catalog["UP000005640"];
        Assert.That(human.ProteomeId, Is.EqualTo("UP000005640"));
        Assert.That(human.Organism, Is.EqualTo("Homo sapiens (Human)"));
        Assert.That(human.TaxonId, Is.EqualTo("9606"));
        Assert.That(human.ProteinCount, Is.EqualTo(75776));
        Assert.That(human.Busco, Is.Null, "the human row's BUSCO cell is blank in this catalogue");
        Assert.That(human.Cpd, Is.EqualTo("Outlier"));

        Assert.That(catalog["UP000002311"].TaxonId, Is.EqualTo("559292"));
        Assert.That(catalog.Keys, Has.None.EqualTo("Proteome ID"), "the header is not a proteome");
    }

    /// <summary>The two readers must agree on which proteomes a file holds and what each organism is called.</summary>
    [Test]
    public void AgreesWithUniprotProteomesListOnEveryProteome()
    {
        string path = Fixture("b_availableUniProtProteomes.txt.gz");

        IReadOnlyDictionary<string, UniProtProteomeEntry> catalog = ProteinDbRetriever.UniprotProteomeCatalog(path);
        Dictionary<string, string> list = ProteinDbRetriever.UniprotProteomesList(path);

        // UniprotProteomesList has always returned the header line as an entry; that is not changed here.
        var listed = list.Where(pair => pair.Key != "Proteome ID").ToDictionary(pair => pair.Key, pair => pair.Value);
        Assert.That(catalog.Count, Is.EqualTo(listed.Count));
        Assert.That(catalog.ToDictionary(pair => pair.Key, pair => pair.Value.Organism), Is.EquivalentTo(listed));
    }

    /// <summary>
    /// UniProt's current default header is "Proteome Id / Organism / Organism Id / Protein count": different
    /// capitalisation from the older "Proteome ID", and no BUSCO or CPD columns (observed 2026-11-12).
    /// </summary>
    [Test]
    public void ReadsTheCurrentHeaderWithoutCompletenessColumns()
    {
        string path = Write("proteomes.txt",
            "Proteome Id\tOrganism\tOrganism Id\tProtein count\n" +
            "UP000000625\tEscherichia coli (strain K12)\t83333\t4403\n");

        UniProtProteomeEntry entry = ProteinDbRetriever.UniprotProteomeCatalog(path).Values.Single();

        Assert.That(entry, Is.EqualTo(new UniProtProteomeEntry("UP000000625", "Escherichia coli (strain K12)", "83333", 4403, null, null)));
    }

    [Test]
    public void FindsColumnsByNameWhateverTheirOrder()
    {
        string path = Write("proteomes.txt",
            "Protein count\tOrganism Id\tCPD\tOrganism\tProteome Id\n" +
            "4403\t83333\tStandard\tEscherichia coli (strain K12)\tUP000000625\n");

        UniProtProteomeEntry entry = ProteinDbRetriever.UniprotProteomeCatalog(path).Values.Single();

        Assert.That(entry, Is.EqualTo(new UniProtProteomeEntry("UP000000625", "Escherichia coli (strain K12)", "83333", 4403, null, "Standard")));
    }

    [Test]
    public void ABlankOrMissingCellReadsAsNull()
    {
        string path = Write("proteomes.txt",
            "Proteome Id\tOrganism\tOrganism Id\tProtein count\n" +
            "UP000000001\tNo taxon here\t\tnot-a-number\n" +
            "UP000000002\tShort line\n");

        IReadOnlyDictionary<string, UniProtProteomeEntry> catalog = ProteinDbRetriever.UniprotProteomeCatalog(path);

        Assert.That(catalog["UP000000001"].TaxonId, Is.Null);
        Assert.That(catalog["UP000000001"].ProteinCount, Is.Null);
        Assert.That(catalog["UP000000002"], Is.EqualTo(new UniProtProteomeEntry("UP000000002", "Short line", null, null, null, null)));
    }

    [Test]
    public void SkipsLinesWithNoTabAndLetsTheLastRepeatWin()
    {
        string path = Write("proteomes.txt",
            "Proteome Id\tOrganism\tOrganism Id\n" +
            "UP000000001\tFirst\t1\n" +
            "\n" +
            "truncated\n" +
            "UP000000001\tSecond\t2\n");

        IReadOnlyDictionary<string, UniProtProteomeEntry> catalog = ProteinDbRetriever.UniprotProteomeCatalog(path);

        Assert.That(catalog, Has.Count.EqualTo(1));
        Assert.That(catalog["UP000000001"].Organism, Is.EqualTo("Second"));
        Assert.That(catalog["UP000000001"].TaxonId, Is.EqualTo("2"));
    }

    /// <summary>Without a header there is no telling which column is which, so the reader refuses rather than guessing.</summary>
    [Test]
    public void AFileWithoutAHeader_ThrowsMzLibException()
    {
        string path = Write("proteomes.txt", "UP000005640\tHomo sapiens (Human)\t9606\t75776\n");

        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.UniprotProteomeCatalog(path));
        Assert.That(ex.Message, Does.Contain("proteomes.txt").And.Contain("Proteome ID"));
    }

    [Test]
    public void AnEmptyFile_ThrowsMzLibException()
    {
        string path = Write("proteomes.txt", "");

        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.UniprotProteomeCatalog(path));
        Assert.That(ex.Message, Does.Contain("no header line"));
    }

    [Test]
    [TestCase(null)]
    [TestCase("  ")]
    public void BlankPath_ThrowsArgumentException(string path)
    {
        Assert.Throws<ArgumentException>(() => ProteinDbRetriever.UniprotProteomeCatalog(path));
    }

    [Test]
    public void MissingFile_ThrowsFileNotFound()
    {
        var ex = Assert.Throws<FileNotFoundException>(() => ProteinDbRetriever.UniprotProteomeCatalog("badFilePath"));
        Assert.That(ex.FileName, Is.EqualTo("badFilePath"));
    }

    [Test]
    public void UnsupportedExtension_ThrowsArgumentException()
    {
        string path = Write("proteomes.fasta", "Proteome Id\tOrganism\nUP000005640\tHomo sapiens (Human)\n");

        var ex = Assert.Throws<ArgumentException>(() => ProteinDbRetriever.UniprotProteomeCatalog(path));
        Assert.That(ex.Message, Does.Contain(".fasta"));
    }
}
