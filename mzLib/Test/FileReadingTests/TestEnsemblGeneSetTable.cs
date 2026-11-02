using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Security.Cryptography;
using NUnit.Framework;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Ensembl;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The compact gene table: a GTF's gene rows and provenance, written by EnsemblGeneSetWriter and read
    /// back by EnsemblGeneSetReader. A resolution needs only the gene rows, and the human release 116 GTF
    /// is 141 MB compressed (4.66 GB unzipped) against about 0.5 MB for its table.
    ///
    /// The table has one promise: a set read from it is the set read from the GTF, provenance included,
    /// so a resolution counted against either is identical row for row.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestEnsemblGeneSetTable
    {
        private const string Header =
            "#!genome-build GRCh38.p14\n#!genome-version GRCh38\n#!genebuild-last-updated 2025-11\n";

        // Gene 1 has every field; gene 2 has no gene_name; gene 3 has no gene_version. The last row is
        // GAPDH, so humanGAPDH.xml resolves against the fixture.
        private const string Rows =
            "1\thavana\tgene\t100\t900\t.\t-\t.\tgene_id \"ENSG00000000001\"; gene_version \"3\"; gene_name \"GENEA\"; gene_biotype \"protein_coding\";\n" +
            "1\thavana\ttranscript\t100\t900\t.\t-\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; gene_biotype \"protein_coding\";\n" +
            "2\thavana\tgene\t5\t50\t.\t+\t.\tgene_id \"ENSG00000000002\"; gene_version \"1\"; gene_biotype \"lncRNA\";\n" +
            "X\tensembl\tgene\t7\t70\t.\t+\t.\tgene_id \"ENSG00000000003\"; gene_name \"GENEC\"; gene_biotype \"IG_V_gene\";\n" +
            "12\tensembl_havana\tgene\t6534512\t6538374\t.\t+\t.\tgene_id \"ENSG00000111640\"; gene_version \"15\"; gene_name \"GAPDH\"; gene_biotype \"protein_coding\";\n";

        private const string TableColumns = "gene_id\tgene_version\tgene_biotype\tgene_name\tseq_region";

        private string _dir;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "EnsemblGeneSetTable_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private string WritePlain(string name, string text)
        {
            string path = Path.Combine(_dir, name);
            File.WriteAllText(path, text);
            return path;
        }

        private EnsemblGeneSet Gtf(string name = "Homo_sapiens.GRCh38.116.gtf", string rows = Rows) =>
            EnsemblGeneSet.LoadGtf(WritePlain(name, Header + rows));

        private EnsemblGeneSet RoundTrip(EnsemblGeneSet set, string tableName)
        {
            string path = Path.Combine(_dir, tableName);
            EnsemblGeneSetWriter.Write(path, set);
            return EnsemblGeneSetReader.Load(path);
        }

        private static string ValidTable(string sha = "abc", string extraHeader = "", string columns = TableColumns, string rows = "ENSG00000000001\t3\tprotein_coding\tGENEA\t1\n") =>
            "#!ensembl-gene-set-format 1\n#!source-file x.116.gtf.gz\n" +
            (sha == null ? "" : $"#!source-sha256 {sha}\n") + extraHeader + columns + "\n" + rows;

        [TestCase("genes.tsv")]
        [TestCase("genes.tsv.gz")]
        public void RoundTrip_GivesTheSameGenesAndTheGtfsProvenance(string tableName)
        {
            var gtf = Gtf();
            var table = RoundTrip(gtf, tableName);

            Assert.That(table.Genes, Is.EqualTo(gtf.Genes), "every field of every gene, in order");
            Assert.That((table.SourceFileName, table.SourceSha256, table.Release, table.GenomeBuild, table.GenebuildLastUpdated),
                Is.EqualTo((gtf.SourceFileName, gtf.SourceSha256, gtf.Release, gtf.GenomeBuild, gtf.GenebuildLastUpdated)),
                "the provenance is the GTF's, not the table's -- that is what keys a resolution");
            Assert.That(table.SourceFileName, Is.EqualTo("Homo_sapiens.GRCh38.116.gtf"));
        }

        [Test]
        public void RoundTrip_AbsentVersionAndSymbolStayAbsent_NotEmptyOrZero()
        {
            var table = RoundTrip(Gtf(), "genes.tsv");

            table.TryGetGene("ENSG00000000002", out var noName);
            table.TryGetGene("ENSG00000000003", out var noVersion);
            Assert.That(noName.Symbol, Is.Null);
            Assert.That(noVersion.Version, Is.Null);
        }

        [Test]
        public void RoundTrip_AReleaseTheGtfNameDidNotCarryStaysNull()
        {
            var table = RoundTrip(Gtf("genes.gtf"), "genes.tsv");

            Assert.That(table.Release, Is.Null, "an absent header line, not a guessed release");
            Assert.That(File.ReadAllText(Path.Combine(_dir, "genes.tsv")), Does.Not.Contain("#!release"));
        }

        [Test]
        public void Resolution_IsIdenticalAgainstTheTableAndTheGtf()
        {
            var gtf = Gtf();
            var table = RoundTrip(gtf, "genes.tsv.gz");
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "humanGAPDH.xml"),
                true, DecoyType.None, null, false, null, out _);

            string Resolve(EnsemblGeneSet set)
            {
                var output = new StringWriter();
                GeneResolutionTsv.Write(output, new EnsemblGeneResolver(set).ResolveAll(proteins, "sha").ToList());
                return output.ToString();
            }

            string fromGtf = Resolve(gtf);
            Assert.That(fromGtf, Does.Contain("ENSG00000111640"), "the fixture must actually resolve something");
            Assert.That(Resolve(table), Is.EqualTo(fromGtf));
        }

        [Test]
        public void Write_IsDeterministic_WithLfLineEndings()
        {
            var gtf = Gtf();
            string a = Path.Combine(_dir, "a.tsv");
            string b = Path.Combine(_dir, "b.tsv");
            EnsemblGeneSetWriter.Write(a, gtf);
            EnsemblGeneSetWriter.Write(b, gtf);

            Assert.That(File.ReadAllBytes(b), Is.EqualTo(File.ReadAllBytes(a)));
            string text = File.ReadAllText(a);
            Assert.That(text, Does.Not.Contain("\r"), "the same bytes on every platform");
            Assert.That(text, Does.StartWith("#!ensembl-gene-set-format 1\n#!source-file Homo_sapiens.GRCh38.116.gtf\n"));
            Assert.That(text, Does.Contain("\n" + TableColumns + "\nENSG00000000001\t3\tprotein_coding\tGENEA\t1\n"));
            Assert.That(text, Does.Contain("\nENSG00000000002\t1\tlncRNA\t\t2\n"), "no symbol is an empty cell");
        }

        [Test]
        public void Write_AGzPathIsCompressed()
        {
            string path = Path.Combine(_dir, "genes.tsv.gz");
            EnsemblGeneSetWriter.Write(path, Gtf());

            using var gz = new GZipStream(File.OpenRead(path), CompressionMode.Decompress);
            using var reader = new StreamReader(gz);
            Assert.That(reader.ReadLine(), Is.EqualTo("#!ensembl-gene-set-format 1"));
        }

        [Test]
        public void Write_RefusesAValueThatWouldShiftTheRow()
        {
            // The GTF's attribute column is its last, so a tab inside gene_name survives LoadGtf.
            var gtf = Gtf(rows: "1\thavana\tgene\t1\t9\t.\t+\t.\tgene_id \"ENSG00000000009\"; gene_name \"A\tB\"; gene_biotype \"lncRNA\";\n");

            var ex = Assert.Throws<ArgumentException>(() => EnsemblGeneSetWriter.Write(Path.Combine(_dir, "t.tsv"), gtf));
            Assert.That(ex.Message, Does.Contain("gene_name of ENSG00000000009"));
        }

        [Test]
        public void Load_AValidTableReads()
        {
            var set = EnsemblGeneSetReader.Load(WritePlain("t.tsv", ValidTable()));

            Assert.That((set.Count, set.SourceSha256, set.Release), Is.EqualTo((1, "abc", (string)null)));
        }

        private static readonly object[] Malformed =
        {
            new object[] { "#!ensembl-gene-set-format 2\n#!source-file x\n#!source-sha256 a\n" + TableColumns + "\n", "format 1" },
            new object[] { "#!source-file x\n#!ensembl-gene-set-format 1\n#!source-sha256 a\n" + TableColumns + "\n", "first line" },
            new object[] { ValidTable(sha: null), "source-sha256 is required" },
            new object[] { ValidTable(extraHeader: "#!assembly GRCh38\n"), "unknown header key 'assembly'" },
            new object[] { ValidTable(extraHeader: "#!release 116\n#!release 117\n"), "repeated" },
            new object[] { ValidTable(columns: "gene_id\tgene_biotype\tgene_version\tgene_name\tseq_region"), "expected the columns" },
            new object[] { ValidTable(rows: "ENSG00000000001\t3\tprotein_coding\tGENEA\n"), "4 cells, expected 5" },
            new object[] { ValidTable(rows: "ENSG00000000001\tv3\tprotein_coding\tGENEA\t1\n"), "not a number" },
            new object[] { ValidTable(rows: "\t3\tprotein_coding\tGENEA\t1\n"), "required" },
            new object[] { ValidTable(rows: "ENSG00000000001\t3\tprotein_coding\tA\t1\nENSG00000000001\t4\tprotein_coding\tB\t1\n"), "ENSG00000000001 is repeated" },
            new object[] { "", "format 1" },
        };

        [TestCaseSource(nameof(Malformed))]
        public void Load_RefusesRatherThanGuesses(string text, string message)
        {
            var ex = Assert.Throws<InvalidDataException>(() => EnsemblGeneSetReader.Load(WritePlain("bad.tsv", text)));
            Assert.That(ex.Message, Does.Contain(message));
        }

        [Test]
        public void Load_MissingFileThrows()
        {
            Assert.Throws<FileNotFoundException>(() => EnsemblGeneSetReader.Load(Path.Combine(_dir, "absent.tsv")));
        }
    }
}
