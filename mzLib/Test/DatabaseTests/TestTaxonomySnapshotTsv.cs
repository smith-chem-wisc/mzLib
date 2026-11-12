using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Taxonomy;

namespace Test.DatabaseTests;

/// <summary>
/// The two tables a taxonomy is stored as (<see cref="TaxonomySnapshotTsv"/>), and the
/// <see cref="ProteomeTaxonMap"/> read from the second.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestTaxonomySnapshotTsv
{
    private const string NodesHeader = "taxon_id\tparent_taxon_id\trank\tscientific_name\tcommon_name\thidden\n";
    private const string ProteomesHeader = "proteome_id\ttaxon_id\torganism\tprotein_count\tproteome_type\n";

    private static readonly TaxonomyProvenance Provenance =
        new("https://rest.uniprot.org/proteomes/search", "2026_04", "reference:true");

    private string _directory;

    [SetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.WorkDirectory, "TestTaxonomySnapshotTsv_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_directory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_directory))
            Directory.Delete(_directory, recursive: true);
    }

    private string Write(string name, string content)
    {
        string path = Path.Combine(_directory, name);
        File.WriteAllText(path, content);
        return path;
    }

    private static TaxonomyTree Tree() => TaxonomyTree.FromNodes(TestTaxonomyTree.Nodes(), Provenance);

    private static ProteomeTaxonMap Proteomes() => ProteomeTaxonMap.FromProteomes(new[]
    {
        new ProteomeTaxon("UP000005640", 9606, "Homo sapiens", 263886, "Reference proteome"),
        new ProteomeTaxon("UP000000625", 83333, "Escherichia coli (strain K12)", 4403, "Reference and representative proteome"),
        new ProteomeTaxon("UP000000001", 83333, "Escherichia coli (strain K12)", null, null),
    }, Provenance);

    // ---- nodes ---------------------------------------------------------------

    [Test]
    public void Nodes_WriteAHeaderThatSaysWhereTheyCameFrom()
    {
        var text = new StringWriter();
        TaxonomySnapshotTsv.WriteNodes(text, Tree());
        string[] lines = text.ToString().Split('\n');

        Assert.That(lines[0], Is.EqualTo("#!taxonomy_nodes_format 1"));
        Assert.That(lines, Has.Some.StartsWith("#!mzlib_version "));
        Assert.That(lines, Has.Some.EqualTo("#!source https://rest.uniprot.org/proteomes/search"));
        Assert.That(lines, Has.Some.EqualTo("#!source_release 2026_04"));
        Assert.That(lines, Has.Some.EqualTo("#!source_query reference:true"));
        Assert.That(lines, Has.Some.EqualTo(NodesHeader.TrimEnd('\n')));
        Assert.That(lines, Has.Some.EqualTo("2\t131567\tdomain\tBacteria\tbacteria\tfalse"));
        Assert.That(lines, Has.Some.EqualTo("131567\t\tno rank\tcellular organisms\t\ttrue"), "a root has an empty parent");
        Assert.That(text.ToString(), Does.Not.Contain("\r"));
    }

    [Test]
    [TestCase("nodes.tsv")]
    [TestCase("nodes.tsv.gz")]
    public void Nodes_RoundTripThroughAFile(string fileName)
    {
        TaxonomyTree written = Tree();
        string path = Path.Combine(_directory, fileName);

        TaxonomySnapshotTsv.WriteNodes(path, written);
        TaxonomyTree read = TaxonomyTree.Load(path);

        Assert.That(read.Nodes, Is.EqualTo(written.Nodes));
        Assert.That(read.Provenance, Is.EqualTo(Provenance));
        Assert.That(read.AncestorAtRank(83333, TaxonomicRank.Family).Taxon!.ScientificName, Is.EqualTo("Enterobacteriaceae"));
    }

    [Test]
    public void Nodes_WriteTheSameBytesEveryTime()
    {
        var first = new StringWriter();
        var second = new StringWriter();
        TaxonomySnapshotTsv.WriteNodes(first, Tree());
        TaxonomySnapshotTsv.WriteNodes(second, TaxonomyTree.FromNodes(TestTaxonomyTree.Nodes().AsEnumerable().Reverse(), Provenance));

        Assert.That(second.ToString(), Is.EqualTo(first.ToString()), "row order does not depend on the order the taxa were given in");
    }

    /// <summary>A file a user writes by hand needs only the format line, a source line and the columns.</summary>
    [Test]
    public void Nodes_AHandWrittenOverlayAddsToATree()
    {
        string path = Write("my-taxa.tsv",
            "#!taxonomy_nodes_format 1\n" +
            "#!source sequenced in-house\n" +
            NodesHeader +
            "9000001\t561\tspecies\tEscherichia novel\t\tfalse\n" +
            "9000002\t9000001\tstrain\tIn-house isolate 7\t\tfalse\n");

        TaxonomyTree overlaid = Tree().WithOverlay(path);

        Assert.That(overlaid.Count, Is.EqualTo(16));
        Assert.That(overlaid.AncestorAtRank(9000002, TaxonomicRank.Genus).Taxon!.ScientificName, Is.EqualTo("Escherichia"));
    }

    [Test]
    public void Nodes_AnOverlayFileIsNotATreeOnItsOwn()
    {
        string path = Write("my-taxa.tsv",
            "#!taxonomy_nodes_format 1\n#!source sequenced in-house\n" + NodesHeader +
            "9000001\t561\tspecies\tEscherichia novel\t\tfalse\n");

        var ex = Assert.Throws<InvalidDataException>(() => TaxonomyTree.Load(path));
        Assert.That(ex.Message, Does.Contain("my-taxa.tsv").And.Contain("561"));
    }

    [Test]
    public void Nodes_AnOverlayThatRedefinesATaxonIsRefused()
    {
        string path = Write("my-taxa.tsv",
            "#!taxonomy_nodes_format 1\n#!source sequenced in-house\n" + NodesHeader +
            "561\t91347\tgenus\tEscherichia\t\tfalse\n");

        Assert.Throws<InvalidDataException>(() => Tree().WithOverlay(path));
    }

    private static IEnumerable<TestCaseData> BadNodeFiles()
    {
        const string head = "#!taxonomy_nodes_format 1\n#!source test\n" + NodesHeader;
        yield return new TestCaseData("#!taxonomy_nodes_format 2\n#!source test\n" + NodesHeader, "taxonomy_nodes_format").SetName("Nodes_Refuse_AnotherFormatVersion");
        yield return new TestCaseData("#!proteome_taxa_format 1\n#!source test\n" + NodesHeader, "taxonomy_nodes_format").SetName("Nodes_Refuse_AProteomeTable");
        yield return new TestCaseData("#!taxonomy_nodes_format 1\n" + NodesHeader, "source").SetName("Nodes_Refuse_NoSource");
        yield return new TestCaseData("#!taxonomy_nodes_format 1\n#!source test\n#!colour blue\n" + NodesHeader, "colour").SetName("Nodes_Refuse_AnUnknownHeaderKey");
        yield return new TestCaseData("#!taxonomy_nodes_format 1\n#!source test\ntaxon_id\tparent\n", "column header").SetName("Nodes_Refuse_OtherColumns");
        yield return new TestCaseData(head + "2\t\tdomain\tBacteria\n", "line 4").SetName("Nodes_Refuse_AShortRow");
        yield return new TestCaseData(head + "two\t\tdomain\tBacteria\t\tfalse\n", "taxon_id 'two'").SetName("Nodes_Refuse_ATaxonIdThatIsNotANumber");
        yield return new TestCaseData(head + "0\t\tdomain\tBacteria\t\tfalse\n", "taxon_id '0'").SetName("Nodes_Refuse_ATaxonIdOfZero");
        yield return new TestCaseData(head + "2\tx\tdomain\tBacteria\t\tfalse\n", "parent_taxon_id 'x'").SetName("Nodes_Refuse_AParentThatIsNotANumber");
        yield return new TestCaseData(head + "2\t\tdomain\tBacteria\t\tfalse\n2\t\tdomain\tBacteria\t\tfalse\n", "appears twice").SetName("Nodes_Refuse_ARepeatedTaxon");
        yield return new TestCaseData(head + "2\t2\tdomain\tBacteria\t\tfalse\n", "itself").SetName("Nodes_Refuse_ATaxonThatIsItsOwnParent");
        yield return new TestCaseData(head + "2\t\tdomain\t\t\tfalse\n", "scientific_name").SetName("Nodes_Refuse_AnEmptyName");
        yield return new TestCaseData(head + "2\t\tdomain\tBacteria\t\t\n", "hidden").SetName("Nodes_Refuse_AnEmptyHiddenFlag");
        yield return new TestCaseData(head + "2\t\tdomain\tBacteria\t\tyes\n", "hidden").SetName("Nodes_Refuse_AHiddenFlagThatIsNotTrueOrFalse");
        yield return new TestCaseData(head + "2\t3\tdomain\tBacteria\t\tfalse\n", "3").SetName("Nodes_Refuse_AParentThatIsNotInTheFile");
    }

    [TestCaseSource(nameof(BadNodeFiles))]
    public void Nodes_RefuseRatherThanGuess(string content, string expectedInMessage)
    {
        string path = Write("nodes.tsv", content);

        var ex = Assert.Throws<InvalidDataException>(() => TaxonomyTree.Load(path));
        Assert.That(ex.Message, Does.Contain("nodes.tsv").And.Contain(expectedInMessage));
    }

    [Test]
    public void Nodes_AMissingFileIsFileNotFound()
    {
        Assert.Throws<FileNotFoundException>(() => TaxonomyTree.Load(Path.Combine(_directory, "absent.tsv")));
        Assert.Throws<FileNotFoundException>(() => Tree().WithOverlay(Path.Combine(_directory, "absent.tsv")));
    }

    [Test]
    public void Nodes_ANameWithATabCannotBeWritten()
    {
        TaxonomyTree tree = TaxonomyTree.FromNodes(
            new[] { new TaxonomyNode(2, null, "domain", "Bac\tteria", null, false) }, Provenance);

        Assert.Throws<ArgumentException>(() => TaxonomySnapshotTsv.WriteNodes(new StringWriter(), tree));
    }

    // ---- proteomes -----------------------------------------------------------

    [Test]
    public void Proteomes_WriteInIdentifierOrderWithTheirProvenance()
    {
        var text = new StringWriter();
        TaxonomySnapshotTsv.WriteProteomes(text, Proteomes());
        string[] lines = text.ToString().Split('\n');

        Assert.That(lines[0], Is.EqualTo("#!proteome_taxa_format 1"));
        Assert.That(lines, Has.Some.EqualTo("#!source_release 2026_04"));
        int header = Array.IndexOf(lines, ProteomesHeader.TrimEnd('\n'));
        Assert.That(header, Is.GreaterThan(0));
        Assert.That(lines.Skip(header + 1).Take(3), Is.EqualTo(new[]
        {
            "UP000000001\t83333\tEscherichia coli (strain K12)\t\t",
            "UP000000625\t83333\tEscherichia coli (strain K12)\t4403\tReference and representative proteome",
            "UP000005640\t9606\tHomo sapiens\t263886\tReference proteome",
        }));
    }

    [Test]
    [TestCase("proteomes.tsv")]
    [TestCase("proteomes.tsv.gz")]
    public void Proteomes_RoundTripThroughAFile(string fileName)
    {
        ProteomeTaxonMap written = Proteomes();
        string path = Path.Combine(_directory, fileName);

        TaxonomySnapshotTsv.WriteProteomes(path, written);
        ProteomeTaxonMap read = ProteomeTaxonMap.Load(path);

        Assert.That(read.Proteomes, Is.EqualTo(written.Proteomes));
        Assert.That(read.Provenance, Is.EqualTo(Provenance));
    }

    private static IEnumerable<TestCaseData> BadProteomeFiles()
    {
        const string head = "#!proteome_taxa_format 1\n#!source test\n" + ProteomesHeader;
        yield return new TestCaseData("#!taxonomy_nodes_format 1\n#!source test\n" + ProteomesHeader, "proteome_taxa_format").SetName("Proteomes_Refuse_ANodeTable");
        yield return new TestCaseData(head + "\t9606\tHomo sapiens\t\t\n", "proteome_id").SetName("Proteomes_Refuse_AnEmptyIdentifier");
        yield return new TestCaseData(head + "UP1\t9606\tA\t\t\nUP1\t9606\tA\t\t\n", "proteome_id 'UP1'").SetName("Proteomes_Refuse_ARepeatedProteome");
        yield return new TestCaseData(head + "UP1\thuman\tHomo sapiens\t\t\n", "taxon_id 'human'").SetName("Proteomes_Refuse_ATaxonThatIsNotANumber");
        yield return new TestCaseData(head + "UP1\t9606\tHomo sapiens\tmany\t\n", "protein_count 'many'").SetName("Proteomes_Refuse_ACountThatIsNotANumber");
    }

    [TestCaseSource(nameof(BadProteomeFiles))]
    public void Proteomes_RefuseRatherThanGuess(string content, string expectedInMessage)
    {
        string path = Write("proteomes.tsv", content);

        var ex = Assert.Throws<InvalidDataException>(() => ProteomeTaxonMap.Load(path));
        Assert.That(ex.Message, Does.Contain("proteomes.tsv").And.Contain(expectedInMessage));
    }

    // ---- the map -------------------------------------------------------------

    [Test]
    public void Map_FindsAProteomesTaxonAndATaxonsProteomes()
    {
        ProteomeTaxonMap map = Proteomes();

        Assert.That(map.Count, Is.EqualTo(3));
        Assert.That(map.TryGetTaxonId("UP000005640", out long human), Is.True);
        Assert.That(human, Is.EqualTo(9606));
        Assert.That(map.TryGetProteome("UP000000625", out ProteomeTaxon k12), Is.True);
        Assert.That(k12.ProteinCount, Is.EqualTo(4403));
        Assert.That(map.TryGetTaxonId("up000005640", out _), Is.False, "identifiers are compared exactly");
        Assert.That(map.TryGetTaxonId("UP999999999", out _), Is.False);
        Assert.That(map.TryGetTaxonId(null, out _), Is.False);

        Assert.That(map.ProteomesOf(83333).Select(p => p.ProteomeId), Is.EqualTo(new[] { "UP000000001", "UP000000625" }));
        Assert.That(map.ProteomesOf(1), Is.Empty);
    }

    [Test]
    public void Map_RefusesARepeatedOrEmptyProteome()
    {
        var one = new ProteomeTaxon("UP1", 9606, "Homo sapiens", null, null);

        Assert.Throws<InvalidDataException>(() => ProteomeTaxonMap.FromProteomes(new[] { one, one }, Provenance));
        Assert.Throws<InvalidDataException>(() => ProteomeTaxonMap.FromProteomes(new[] { one with { ProteomeId = "" } }, Provenance));
    }

    /// <summary>The downloaded catalogue covers every proteome UniProt lists; rows it cannot place are named, not dropped silently.</summary>
    [Test]
    public void Map_FromTheDownloadedCatalogue_NamesTheRowsItCannotPlace()
    {
        var catalog = new Dictionary<string, UniProtProteomeEntry>
        {
            ["UP000005640"] = new("UP000005640", "Homo sapiens (Human)", "9606", 75776, null, "Outlier"),
            ["UP000000002"] = new("UP000000002", "No taxon", null, 10, null, null),
            ["UP000000001"] = new("UP000000001", "Bad taxon", "nine", 10, null, null),
        };

        ProteomeTaxonMap map = ProteomeTaxonMap.FromCatalog(catalog, "availableUniProtProteomes.txt.gz", out IReadOnlyList<string> skipped);

        Assert.That(map.Proteomes.Single(), Is.EqualTo(new ProteomeTaxon("UP000005640", 9606, "Homo sapiens (Human)", 75776, null)));
        Assert.That(skipped, Is.EqualTo(new[] { "UP000000001", "UP000000002" }));
        Assert.That(map.Provenance.Source, Is.EqualTo("availableUniProtProteomes.txt.gz"));
    }

    [Test]
    public void Map_FromTheCatalogueFixture_PlacesHuman()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "b_availableUniProtProteomes.txt.gz");

        ProteomeTaxonMap map = ProteomeTaxonMap.FromCatalog(ProteinDbRetriever.UniprotProteomeCatalog(path), path, out _);

        Assert.That(map.TryGetTaxonId("UP000005640", out long human), Is.True);
        Assert.That(human, Is.EqualTo(9606));
    }
}
