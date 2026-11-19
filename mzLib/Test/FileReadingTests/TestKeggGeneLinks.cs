using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases.Pathways;
using static Test.FileReadingTests.TestReactomeMembership;

namespace Test.FileReadingTests
{
    /// <summary>
    /// A database's KEGG gene links (<see cref="KeggGeneLinksTsv"/>): the KEGG GENES ids UniProt distributes, verbatim,
    /// for link-outs. No status rows, no KEGG content beyond the id.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestKeggGeneLinks
    {
        private const string Head =
            "#!kegg_gene_links_format 1\n" +
            "#!annotation_db_sha256 0123abcd\n" +
            "accession\tkegg_gene_id\n";

        private string _dir;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "KeggGeneLinks_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private string WriteFile(Action<TextWriter> write)
        {
            string path = Path.Combine(_dir, "t4.tsv");
            using (var writer = new StreamWriter(path)) write(writer);
            return path;
        }

        [Test]
        public void Links_ReadFromEveryEntryWithAKeggId()
        {
            var links = KeggGeneLinksTsv.Links(Fixture());

            Assert.That(links, Is.EqualTo(new[]
            {
                new KeggGeneLinkRow("P11111", "hsa:1111"),
                new KeggGeneLinkRow("P11111_A20T", "hsa:1111"),
                new KeggGeneLinkRow("P22222", "hsa:2222"),
                new KeggGeneLinkRow("P22222", "hsa:2223"),
                new KeggGeneLinkRow("P66666", "hsa:6666"),
                new KeggGeneLinkRow("Q11111", "mmu:1111"),
            }), "ordinal order; the duplicate P66666 once; the decoy skipped");
        }

        [Test]
        public void Links_HaveNoStatusRows()
        {
            var links = KeggGeneLinksTsv.Links(Fixture());

            Assert.That(links.Select(l => l.Accession), Does.Not.Contain("P55555").And.Not.Contain("Q22222"));
            Assert.That(links.All(l => !string.IsNullOrEmpty(l.KeggGeneId)));
            Assert.That(KeggGeneLinksTsv.Columns, Is.EqualTo(new[] { "accession", "kegg_gene_id" }));
        }

        [Test]
        public void Links_SkipNullsAndRefuseANullList()
        {
            Assert.That(KeggGeneLinksTsv.Links(new Protein[] { null }), Is.Empty);
            Assert.Throws<ArgumentNullException>(() => KeggGeneLinksTsv.Links(null));
        }

        [Test]
        public void WriteThenRead_RoundTrips()
        {
            var links = KeggGeneLinksTsv.Links(Fixture());

            var table = KeggGeneLinksTsv.Read(WriteFile(w => KeggGeneLinksTsv.Write(w, links, DbSha)));

            Assert.That(table.Rows, Is.EqualTo(links));
            Assert.That(table.Header["annotation_db_sha256"], Is.EqualTo(DbSha));
            Assert.That(table.Header["kegg_gene_links_format"], Is.EqualTo("1"));
            Assert.That(table.Header.ContainsKey("mzlib_version"), Is.True);
        }

        [TestCase(null)]
        [TestCase("")]
        public void Writer_RefusesAMissingSha_BeforeWritingAnything(string sha)
        {
            using var writer = new StringWriter();

            Assert.Throws<ArgumentException>(() => KeggGeneLinksTsv.Write(writer, KeggGeneLinksTsv.Links(Fixture()), sha));
            Assert.That(writer.ToString(), Is.Empty);
        }

        [TestCase("P1", "hsa:1", "P1", "hsa:1")]
        [TestCase("", "hsa:1", null, null)]
        [TestCase("P1", "", null, null)]
        [TestCase("P1\t", "hsa:1", null, null)]
        [TestCase("P1", "hsa:\n1", null, null)]
        public void Writer_RefusesABadRow_BeforeWritingAnything(string accession, string id, string secondAccession, string secondId)
        {
            var rows = new[] { new KeggGeneLinkRow("P0", "hsa:0"), new KeggGeneLinkRow(accession, id) }.ToList();
            if (secondAccession != null) rows.Add(new KeggGeneLinkRow(secondAccession, secondId));
            using var writer = new StringWriter();

            Assert.Throws<ArgumentException>(() => KeggGeneLinksTsv.Write(writer, rows, DbSha));
            Assert.That(writer.ToString(), Is.Empty);
        }

        [Test]
        public void Writer_RejectsNullArguments()
        {
            Assert.Throws<ArgumentNullException>(() => KeggGeneLinksTsv.Write(null, Array.Empty<KeggGeneLinkRow>(), DbSha));
            Assert.Throws<ArgumentNullException>(() => KeggGeneLinksTsv.Write(new StringWriter(), null, DbSha));
        }

        [TestCase("P1\thsa:1\nP1\thsa:1\n", "appears twice")]
        [TestCase("\thsa:1\n", "accession is empty")]
        [TestCase("P1\t\n", "kegg_gene_id is empty")]
        [TestCase("P1\n", "expected 2 cells")]
        [TestCase("P1\thsa:1\textra\n", "expected 2 cells")]
        public void Readers_Refuse(string rows, string fragment)
        {
            string path = WriteFile(w => w.Write(Head + rows));

            var exception = Assert.Throws<InvalidDataException>(() => KeggGeneLinksTsv.Read(path));
            Assert.That(exception.Message, Does.Contain(fragment));
        }

        [Test]
        public void Readers_RefuseAMissingSha()
        {
            string path = WriteFile(w => w.Write("#!kegg_gene_links_format 1\naccession\tkegg_gene_id\n"));

            var exception = Assert.Throws<InvalidDataException>(() => KeggGeneLinksTsv.Read(path));
            Assert.That(exception.Message, Does.Contain("'annotation_db_sha256' is missing"));
        }

        [Test]
        public void Readers_MissingFile_Throws()
        {
            Assert.Throws<FileNotFoundException>(() => KeggGeneLinksTsv.Read(Path.Combine(_dir, "absent.tsv")));
        }
    }
}
