using System;
using System.Collections;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Security.Cryptography;
using NUnit.Framework;
using UsefulProteomicsDatabases.Ensembl;

namespace Test.FileReadingTests
{
    /// <summary>
    /// EnsemblGeneSet is the release-pinned gene set a resolution is counted against. It answers
    /// one question the protein database cannot: is this gene id on the assembly we are counting?
    /// Ensembl's cross-references also name ALT-haplotype and patch genes, and counting those as
    /// separate genes made one locus look like twenty (measured on human release 116: a 6.99%
    /// multi-gene rate that is 0.36% against the primary-assembly GTF).
    ///
    /// So the set records what it was read from -- file name, sha256, release, genome build -- and
    /// reads gene rows only.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestEnsemblGeneSet
    {
        private const string Header =
            "#!genome-build GRCh38.p14\n#!genome-version GRCh38\n#!genome-date 2013-12\n" +
            "#!genome-build-accession GCA_000001405.29\n#!genebuild-last-updated 2025-11\n";

        private const string Rows =
            "1\thavana\tgene\t100\t900\t.\t-\t.\tgene_id \"ENSG00000000001\"; gene_version \"3\"; gene_name \"GENEA\"; gene_source \"havana\"; gene_biotype \"protein_coding\";\n" +
            "1\thavana\ttranscript\t100\t900\t.\t-\t.\tgene_id \"ENSG00000000001\"; gene_version \"3\"; transcript_id \"ENST00000000001\"; gene_biotype \"protein_coding\";\n" +
            "1\thavana\texon\t100\t300\t.\t-\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\";\n" +
            "2\thavana\tgene\t5\t50\t.\t+\t.\tgene_id \"ENSG00000000002\"; gene_version \"1\"; gene_source \"havana\"; gene_biotype \"lncRNA\";\n" +
            "X\tensembl\tgene\t7\t70\t.\t+\t.\tgene_id \"ENSG00000000003\"; gene_name \"GENEC\"; gene_biotype \"IG_V_gene\";\n";

        private string _dir;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "EnsemblGeneSet_" + Guid.NewGuid().ToString("N"));
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

        private string WriteGz(string name, string text)
        {
            string path = Path.Combine(_dir, name);
            using var file = File.Create(path);
            using var gz = new GZipStream(file, CompressionLevel.Optimal);
            using var writer = new StreamWriter(gz);
            writer.Write(text);
            return path;
        }

        [Test]
        public void LoadGtf_ReadsGeneRowsOnly()
        {
            var set = EnsemblGeneSet.LoadGtf(WritePlain("Homo_sapiens.GRCh38.116.gtf", Header + Rows));

            Assert.That(set.Count, Is.EqualTo(3), "transcript and exon rows are not genes");
            Assert.That(set.Contains("ENSG00000000001"), Is.True);
            Assert.That(set.Contains("ENST00000000001"), Is.False);
        }

        [Test]
        public void LoadGtf_KeepsBiotypeVersionAndSymbol_SymbolNullWhenAbsent()
        {
            var set = EnsemblGeneSet.LoadGtf(WritePlain("Homo_sapiens.GRCh38.116.gtf", Header + Rows));

            Assert.That(set.TryGetGene("ENSG00000000001", out var a), Is.True);
            Assert.That((a.Biotype, a.Version, a.Symbol, a.SeqRegion), Is.EqualTo(("protein_coding", (int?)3, "GENEA", "1")));

            set.TryGetGene("ENSG00000000002", out var b);
            Assert.That(b.Symbol, Is.Null, "no gene_name is a missing symbol, not an empty one");

            set.TryGetGene("ENSG00000000003", out var c);
            Assert.That(c.Version, Is.Null, "no gene_version is absent, not 0");
            Assert.That(c.Biotype, Is.EqualTo("IG_V_gene"));
        }

        [Test]
        public void LoadGtf_RecordsProvenance()
        {
            string path = WriteGz("Homo_sapiens.GRCh38.116.gtf.gz", Header + Rows);
            var set = EnsemblGeneSet.LoadGtf(path);

            Assert.That(set.SourceFileName, Is.EqualTo("Homo_sapiens.GRCh38.116.gtf.gz"));
            Assert.That(set.Release, Is.EqualTo("116"));
            Assert.That(set.GenomeBuild, Is.EqualTo("GRCh38.p14"));
            Assert.That(set.GenebuildLastUpdated, Is.EqualTo("2025-11"));

            string expected = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();
            Assert.That(set.SourceSha256, Is.EqualTo(expected),
                "the hash of the file as published -- Ensembl's own CHECKSUMS are over the .gz");
        }

        [Test]
        public void LoadGtf_ReleaseIsNullWhenTheFileNameDoesNotCarryOne()
        {
            var set = EnsemblGeneSet.LoadGtf(WritePlain("genes.gtf", Header + Rows));

            Assert.That(set.Release, Is.Null, "not guessed");
        }

        [Test]
        public void LoadGtf_CompressedAndPlainGiveTheSameGenes()
        {
            var plain = EnsemblGeneSet.LoadGtf(WritePlain("a.gtf", Header + Rows));
            var gz = EnsemblGeneSet.LoadGtf(WriteGz("b.gtf.gz", Header + Rows));

            Assert.That(gz.GeneIds, Is.EqualTo(plain.GeneIds));
        }

        [Test]
        public void LoadGtf_GeneRowWithoutAGeneIdIsAnError_NotASkip()
        {
            string path = WritePlain("bad.gtf", Header + "1\thavana\tgene\t1\t9\t.\t+\t.\tgene_biotype \"lncRNA\";\n");

            var ex = Assert.Throws<InvalidDataException>(() => EnsemblGeneSet.LoadGtf(path));
            Assert.That(ex.Message, Does.Contain("line 6"));
        }

        [Test]
        public void LoadGtf_MissingFileThrows()
        {
            Assert.Throws<FileNotFoundException>(() => EnsemblGeneSet.LoadGtf(Path.Combine(_dir, "absent.gtf")));
        }

        [Test]
        public void Contains_IsOnStableIdsOnly()
        {
            var set = EnsemblGeneSet.LoadGtf(WritePlain("a.gtf", Header + Rows));

            Assert.That(set.Contains("ENSG00000000001.3"), Is.False, "callers strip the version; the set does not guess");
            Assert.That(set.GeneIds, Is.Ordered.Using((IComparer)StringComparer.Ordinal));
        }
    }
}
