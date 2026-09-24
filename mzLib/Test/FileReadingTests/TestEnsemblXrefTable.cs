using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Security.Cryptography;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases.Ensembl;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Ensembl's own accession -> gene cross-reference dump, used as a second opinion beside the links a
    /// UniProt XML carries. The two disagree for a reason, not by accident: UniProt links an accession to
    /// every Ensembl transcript that encodes it, including readthrough genes (MSH5 and MSH5-SAPCD1) and
    /// identical paralogs, while Ensembl's xref assigns it where Ensembl's own mapping puts it. Measured
    /// on a reviewed human proteome: 350 multi-gene accessions by the XML, 69 by the xref. Neither is
    /// dropped; every gene row says whether the xref agrees.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestEnsemblXrefTable
    {
        private const string Header =
            "gene_stable_id\ttranscript_stable_id\tprotein_stable_id\txref\tdb_name\tinfo_type\tsource_identity\txref_identity\tlinkage_type\n";

        private static string Row(string gene, string transcript, string xref, string db, string info) =>
            $"{gene}\t{transcript}\tENSP0\t{xref}\t{db}\t{info}\t-\t-\t-\n";

        private string _dir;
        private EnsemblGeneSet _genes;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "EnsemblXref_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
            string gtf = Path.Combine(_dir, "Homo_sapiens.GRCh38.116.gtf");
            File.WriteAllText(gtf, string.Concat(new[] { "ENSG00000000001", "ENSG00000000002", "ENSG00000000003" }.Select(g =>
                $"1\thavana\tgene\t1\t9\t.\t+\t.\tgene_id \"{g}\"; gene_version \"1\"; gene_biotype \"protein_coding\";\n")));
            _genes = EnsemblGeneSet.LoadGtf(gtf);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
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

        private EnsemblXrefTable StandardTable() => EnsemblXrefTable.Load(WriteGz("Homo_sapiens.GRCh38.116.uniprot.tsv.gz",
            Header +
            Row("ENSG00000000001", "ENST1", "P11111", "Uniprot/SWISSPROT", "SEQUENCE_MATCH") +
            Row("ENSG00000000001", "ENST2", "P11111", "Uniprot/SWISSPROT", "DIRECT") +
            Row("ENSG00000000001", "ENST1", "P11111-2", "Uniprot_isoform", "DIRECT") +
            Row("ENSG00000000003", "ENST3", "P33333", "Uniprot/SPTREMBL", "SEQUENCE_MATCH")));

        private static DatabaseReference Transcript(string transcriptId, string versionedGeneId) =>
            new("Ensembl", transcriptId, new List<Tuple<string, string>> { new("gene ID", versionedGeneId) });

        private static Protein Entry(string accession, params DatabaseReference[] references) =>
            new("PEPTIDEK", accession, databaseReferences: references.ToList());

        // ---- the table ----

        [Test]
        public void Load_StrongestLinkPerPair_AndProvenance()
        {
            var table = StandardTable();

            Assert.That(table.TryGetLink("P11111", "ENSG00000000001", out string info), Is.True);
            Assert.That(info, Is.EqualTo("DIRECT"), "DIRECT is an assertion and outranks an inference for the same pair");
            Assert.That(table.TryGetLink("P33333", "ENSG00000000003", out info) && info == "SEQUENCE_MATCH", Is.True);
            Assert.That(table.TryGetLink("P11111", "ENSG00000000002", out _), Is.False);
            Assert.That(table.ContainsAccession("P11111-2"), Is.True, "isoform accessions are rows of their own");

            Assert.That(table.Release, Is.EqualTo("116"));
            Assert.That(table.SourceSha256, Is.EqualTo(Convert.ToHexString(SHA256.HashData(
                File.ReadAllBytes(Path.Combine(_dir, "Homo_sapiens.GRCh38.116.uniprot.tsv.gz")))).ToLowerInvariant()));
        }

        [Test]
        public void Load_RefusesAnUnexpectedLayout()
        {
            string path = WriteGz("x.uniprot.tsv.gz", "gene\txref\nENSG1\tP11111\n");

            var ex = Assert.Throws<InvalidDataException>(() => EnsemblXrefTable.Load(path));
            Assert.That(ex.Message, Does.Contain("columns"));
        }

        [Test]
        public void Load_MissingFileThrows()
        {
            Assert.Throws<FileNotFoundException>(() => EnsemblXrefTable.Load(Path.Combine(_dir, "absent.tsv.gz")));
        }

        // ---- the column it feeds ----

        [Test]
        public void GeneRows_SayWhetherTheXrefAgrees_AndNeitherSourceIsDropped()
        {
            var resolver = new EnsemblGeneResolver(_genes, StandardTable());
            var rows = resolver.Resolve(Entry("P11111",
                Transcript("ENST1.1", "ENSG00000000001.1"),
                Transcript("ENST9.1", "ENSG00000000002.1")), "sha");

            Assert.That(rows.Select(r => r.Outcome).Distinct(), Is.EqualTo(new[] { GeneResolutionOutcome.MultiGene }),
                "the XML's links are reported as they are");
            var byGene = rows.ToDictionary(r => r.GeneId);
            Assert.That(byGene["ENSG00000000001"].EnsemblXrefAgrees, Is.True);
            Assert.That(byGene["ENSG00000000001"].EnsemblXrefInfoType, Is.EqualTo("DIRECT"));
            Assert.That(byGene["ENSG00000000002"].EnsemblXrefAgrees, Is.False, "e.g. a readthrough gene the xref does not assign");
            Assert.That(byGene["ENSG00000000002"].EnsemblXrefInfoType, Is.Null);
            Assert.That(rows.Select(r => r.EnsemblXrefSha256).Distinct().Single(), Is.EqualTo(resolver.Xrefs.SourceSha256));
        }

        [Test]
        public void WithoutAnXrefTable_AgreementIsUnknown_NotFalse()
        {
            var r = new EnsemblGeneResolver(_genes).Resolve(Entry("P11111", Transcript("ENST1.1", "ENSG00000000001.1")), "sha").Single();

            Assert.That(r.EnsemblXrefAgrees, Is.Null);
            Assert.That(r.EnsemblXrefSha256, Is.Null);
        }

        [Test]
        public void RowsWithoutAGene_CarryNoAgreement()
        {
            var r = new EnsemblGeneResolver(_genes, StandardTable()).Resolve(Entry("P99999"), "sha").Single();

            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.NotInSource));
            Assert.That(r.EnsemblXrefAgrees, Is.Null, "there is no gene to agree about");
        }

        [Test]
        public void IsoformAccession_ExactXrefRowFirst_ThenItsEntry()
        {
            var resolver = new EnsemblGeneResolver(_genes, StandardTable());

            var exact = resolver.Resolve(Entry("P11111-2", Transcript("ENST1.1", "ENSG00000000001.1")), "sha").Single();
            Assert.That(exact.EnsemblXrefAgrees, Is.True);

            var viaEntry = resolver.Resolve(Entry("P11111-5", Transcript("ENST1.1", "ENSG00000000001.1")), "sha").Single();
            Assert.That(viaEntry.EnsemblXrefAgrees, Is.True, "P11111-5 has no xref row of its own; its entry does");
        }

        [Test]
        public void Tsv_CarriesTheAgreementColumns()
        {
            var resolver = new EnsemblGeneResolver(_genes, StandardTable());
            var output = new StringWriter();
            GeneResolutionTsv.Write(output, resolver.Resolve(Entry("P11111", Transcript("ENST1.1", "ENSG00000000001.1")), "sha"));

            var lines = output.ToString().Replace("\r", "").TrimEnd('\n').Split('\n');
            var cells = lines[0].Split('\t').Zip(lines[1].Split('\t')).ToDictionary(p => p.First, p => p.Second);
            Assert.That(cells["ensembl_xref_agrees"], Is.EqualTo("true"));
            Assert.That(cells["ensembl_xref_info_type"], Is.EqualTo("DIRECT"));
        }

        // ---- genes only the xref links: rows of their own, never a changed outcome ----

        private EnsemblXrefTable XrefOnlyTable() => EnsemblXrefTable.Load(WriteGz("Rattus_norvegicus.GRCr8.116.uniprot.tsv.gz",
            Header +
            Row("ENSG00000000001", "ENST1", "P44444", "Uniprot/SWISSPROT", "DIRECT") +
            Row("ENSG00000000003", "ENST3", "P44444", "Uniprot/SWISSPROT", "SEQUENCE_MATCH") +
            Row("ENSG00000000003", "ENST3", "P55555", "Uniprot/SWISSPROT", "DIRECT") +
            Row("ENSG00000000099", "ENST99", "P66666", "Uniprot/SWISSPROT", "DIRECT") +
            Row("ENSG00000000002", "ENST2", "P77777", "Uniprot/SWISSPROT", "DIRECT")));

        private static Protein[] XrefOnlyProteins() => new[]
        {
            Entry("P44444", Transcript("ENST1.1", "ENSG00000000001.1")),          // XML one gene, xref adds a second
            Entry("P55555", Transcript("ENST55.1", "ENSG00000055001.1")),         // XML only off the set, xref rescues
            Entry("P66666"),                                                     // xref gene is off the set too
            Entry("P88888"),                                                     // neither source knows it
            new Protein("PEPTIDEK", "P77777", isContaminant: true),              // xref would resolve it
        };

        [Test]
        public void XrefOnlyGene_IsItsOwnRow_AndKeepsTheSearchDatabasesOutcome()
        {
            var rows = new EnsemblGeneResolver(_genes, XrefOnlyTable()).Resolve(XrefOnlyProteins()[0], "sha");

            Assert.That(rows.Select(r => (r.GeneId, r.Source)), Is.EqualTo(new[]
            {
                ("ENSG00000000001", EnsemblGeneResolver.SearchDatabaseSource),
                ("ENSG00000000003", EnsemblGeneResolver.EnsemblXrefSource),
            }));
            Assert.That(rows.Select(r => (r.Outcome, r.GeneCount)).Distinct(), Is.EqualTo(new[] { (GeneResolutionOutcome.Resolved, 1) }),
                "the outcome reports the search database's links; the xref adds an answer, it does not change one");
            var added = rows[1];
            Assert.That((added.EnsemblXrefAgrees, added.EnsemblXrefInfoType, added.VersionedGeneId),
                Is.EqualTo(((bool?)true, "SEQUENCE_MATCH", "ENSG00000000003.1")), "the version is the gene set's");
        }

        [Test]
        public void OffPrimaryOnlyEntry_KeepsItsOutcomeRow_AndGainsTheXrefGene()
        {
            var rows = new EnsemblGeneResolver(_genes, XrefOnlyTable()).Resolve(XrefOnlyProteins()[1], "sha");

            Assert.That(rows.Select(r => (r.Outcome, r.GeneId, r.Source, r.OffPrimaryGenes)), Is.EqualTo(new[]
            {
                (GeneResolutionOutcome.OffPrimaryOnly, (string)null, EnsemblGeneResolver.SearchDatabaseSource, 1),
                (GeneResolutionOutcome.OffPrimaryOnly, "ENSG00000000003", EnsemblGeneResolver.EnsemblXrefSource, 1),
            }));
        }

        [Test]
        public void XrefGeneOutsideTheGeneSet_IsNotAdded()
        {
            var r = new EnsemblGeneResolver(_genes, XrefOnlyTable()).Resolve(XrefOnlyProteins()[2], "sha").Single();

            Assert.That((r.Outcome, r.GeneId), Is.EqualTo((GeneResolutionOutcome.NotInSource, (string)null)),
                "the xref's off-set links are held to the same gene set as the XML's");
        }

        [Test]
        public void Contaminant_GainsNoXrefGene()
        {
            var r = new EnsemblGeneResolver(_genes, XrefOnlyTable()).Resolve(XrefOnlyProteins()[4], "sha").Single();

            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.ContaminantNotMapped));
        }

        [Test]
        public void FilteringToTheSearchDatabaseSource_GivesTheSearchDatabasesViewUnchanged()
        {
            var proteins = XrefOnlyProteins();
            var withXref = new EnsemblGeneResolver(_genes, XrefOnlyTable()).ResolveAll(proteins, "sha")
                .Where(r => r.Source == EnsemblGeneResolver.SearchDatabaseSource);
            var without = new EnsemblGeneResolver(_genes).ResolveAll(proteins, "sha");

            Assert.That(withXref.Select(r => (r.Accession, r.Outcome, r.GeneCount, r.GeneId, r.VersionedGeneId, r.OffPrimaryGenes)),
                Is.EqualTo(without.Select(r => (r.Accession, r.Outcome, r.GeneCount, r.GeneId, r.VersionedGeneId, r.OffPrimaryGenes))));
        }

        [Test]
        public void FilteringToAgreement_GivesEveryXrefGeneInTheSet()
        {
            var table = XrefOnlyTable();
            var proteins = XrefOnlyProteins().Where(p => !p.IsContaminant).ToList();
            var rows = new EnsemblGeneResolver(_genes, table).ResolveAll(proteins, "sha").ToList();

            foreach (var p in proteins)
            {
                var agreeing = rows.Where(r => r.Accession == p.Accession && r.EnsemblXrefAgrees == true).Select(r => r.GeneId);
                var xrefInSet = table.GenesFor(p.Accession).Select(g => g.GeneId).Where(_genes.Contains);
                Assert.That(agreeing, Is.EquivalentTo(xrefInSet), p.Accession);
            }
        }
    }
}
