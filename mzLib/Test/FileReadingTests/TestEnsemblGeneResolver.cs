using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Omics.BioPolymer;
using Proteomics;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Ensembl;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Accession -> stable gene resolution against a release-pinned gene set, reading the gene links
    /// the search database itself carries (Protein.EnsemblGeneReferences).
    ///
    /// Each test pins a distinction that would otherwise collapse into a confidently wrong answer:
    /// a multi-gene accession reduced to a pick, an ALT haplotype counted as a second gene, a
    /// contaminant mapped, "the source has nothing" confused with "not an accession at all".
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestEnsemblGeneResolver
    {
        private const string SearchDbSha = "760984e8d402ade6b1105b811532bdd4041e33e66bb5dc204402d4d6a7be8838";

        private string _dir;
        private EnsemblGeneResolver _resolver;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "EnsemblResolver_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
            string gtf = Path.Combine(_dir, "Homo_sapiens.GRCh38.116.gtf");
            File.WriteAllText(gtf,
                "#!genome-build GRCh38.p14\n" +
                Gene("ENSG00000000001", "GENEA", "protein_coding") +
                Gene("ENSG00000000002", "GENEB", "protein_coding") +
                Gene("ENSG00000111640", "GAPDH", "protein_coding") +
                "1\thavana\tgene\t1\t9\t.\t+\t.\tgene_id \"ENSG00000000004\"; gene_biotype \"lncRNA\";\n");
            _resolver = new EnsemblGeneResolver(EnsemblGeneSet.LoadGtf(gtf));
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private static string Gene(string id, string name, string biotype) =>
            $"1\thavana\tgene\t1\t9\t.\t+\t.\tgene_id \"{id}\"; gene_version \"7\"; gene_name \"{name}\"; gene_biotype \"{biotype}\";\n";

        private static DatabaseReference Transcript(string transcriptId, string versionedGeneId) =>
            new("Ensembl", transcriptId, new List<Tuple<string, string>> { new("gene ID", versionedGeneId) });

        private static Protein Entry(string accession, string uniProtGeneName = null, bool isContaminant = false,
            bool isDecoy = false, params DatabaseReference[] references) =>
            new("PEPTIDEK", accession,
                geneNames: uniProtGeneName == null ? null : new List<Tuple<string, string>> { new("primary", uniProtGeneName) },
                isDecoy: isDecoy, isContaminant: isContaminant, databaseReferences: references.ToList());

        private IReadOnlyList<GeneResolution> Resolve(Protein p) => _resolver.Resolve(p, SearchDbSha);

        [Test]
        public void OneGeneOnTheAssembly_IsResolved_WithTheReleaseSymbolAndBothIdForms()
        {
            var rows = Resolve(Entry("P11111", "UNIPROTNAME", references: new[]
            {
                Transcript("ENST00000000001.1", "ENSG00000000001.7"),
                Transcript("ENST00000000002.1", "ENSG00000000001.7"),
            }));

            var r = rows.Single();
            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.Resolved));
            Assert.That((r.GeneId, r.VersionedGeneId, r.GeneCount), Is.EqualTo(("ENSG00000000001", "ENSG00000000001.7", 1)));
            Assert.That(r.GeneSymbol, Is.EqualTo("GENEA"), "the display label comes from the pinned release");
            Assert.That(r.GeneBiotype, Is.EqualTo("protein_coding"));
            Assert.That(r.UniProtGeneName, Is.EqualTo("UNIPROTNAME"), "what the search database called it, alongside");
        }

        [Test]
        public void MultiGene_IsOneRowPerGene_NeverAPick()
        {
            var rows = Resolve(Entry("P11111", references: new[]
            {
                Transcript("ENST00000000009.1", "ENSG00000000002.7"),
                Transcript("ENST00000000001.1", "ENSG00000000001.7"),
            }));

            Assert.That(rows.Select(r => r.GeneId), Is.EqualTo(new[] { "ENSG00000000001", "ENSG00000000002" }));
            Assert.That(rows.Select(r => r.Outcome).Distinct(), Is.EqualTo(new[] { GeneResolutionOutcome.MultiGene }));
            Assert.That(rows.Select(r => r.GeneCount).Distinct(), Is.EqualTo(new[] { 2 }));
        }

        [Test]
        public void GenesOffTheAssembly_DoNotCountAsExtraGenes_ButAreCountedVisibly()
        {
            var rows = Resolve(Entry("P11111", references: new[]
            {
                Transcript("ENST00000000001.1", "ENSG00000000001.7"),
                Transcript("ENST00000000002.1", "ENSG00000999991.1"),
                Transcript("ENST00000000003.1", "ENSG00000999992.1"),
            }));

            var r = rows.Single();
            Assert.That((r.Outcome, r.GeneId), Is.EqualTo((GeneResolutionOutcome.Resolved, "ENSG00000000001")));
            Assert.That(r.OffPrimaryGenes, Is.EqualTo(2), "dropped, but visibly");
        }

        [Test]
        public void OnlyOffAssemblyGenes_IsItsOwnOutcome_NotNotInSource()
        {
            var r = Resolve(Entry("P11111", references: Transcript("ENST00000000002.1", "ENSG00000999991.1"))).Single();

            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.OffPrimaryOnly));
            Assert.That(r.GeneId, Is.Null);
            Assert.That(r.OffPrimaryGenes, Is.EqualTo(1));
        }

        [Test]
        public void NoGeneLink_IsNotInSource_AndCarriesTheUniProtGeneNameAsALabel()
        {
            // REQ-AGING-8: for entries no source gives a gene id (HERV-K and others), the pinned
            // database's own gene name is the only gene-level fact there is. It travels as a label;
            // the outcome still says no gene id exists.
            var r = Resolve(Entry("P61571", "ERVK-21")).Single();

            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.NotInSource));
            Assert.That(r.GeneId, Is.Null);
            Assert.That(r.GeneSymbol, Is.Null, "no release symbol without a release gene");
            Assert.That(r.UniProtGeneName, Is.EqualTo("ERVK-21"));
        }

        [Test]
        public void Contaminant_IsNeverMapped_EvenWhenItWouldResolve()
        {
            var r = Resolve(Entry("P02768", "ALB", isContaminant: true,
                references: Transcript("ENST00000000001.1", "ENSG00000000001.7"))).Single();

            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.ContaminantNotMapped));
            Assert.That(r.GeneId, Is.Null);
        }

        [Test]
        public void DecoyAccession_IsUnrecognized_NotNotInSource()
        {
            var r = Resolve(Entry("DECOY_P11111", isDecoy: true)).Single();

            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.UnrecognizedAccession));
            Assert.That(r.Namespace, Is.EqualTo(AccessionNamespace.Unrecognized));
        }

        [Test]
        public void IsoformAccession_KeepsItsVerbatimAccessionAndEntry()
        {
            var r = Resolve(Entry("P11111-2", references: Transcript("ENST00000000001.1", "ENSG00000000001.7"))).Single();

            Assert.That((r.Accession, r.EntryAccession, r.Isoform), Is.EqualTo(("P11111-2", "P11111", (int?)2)));
        }

        private static Protein VariantOf(Protein consensus, string variantAccession, List<DatabaseReference> references = null) =>
            new("PEPTIDEK", variantAccession, geneNames: consensus.GeneNames?.ToList(), databaseReferences: references,
                appliedSequenceVariations: new List<SequenceVariation> { new(1, 1, "P", "S", "S70N") },
                nonVariantProtein: consensus);

        [Test]
        public void SequenceVariant_OfAnEntryWithNoGene_IsNotInSource_NotUnrecognized()
        {
            // LoadProteinXML applies UniProt's sequence variants by default and names each proteoform
            // "{accession}_{change}". That is not an accession grammar, but it is a known entry.
            var consensus = Entry("A0A087X1C5", "CYP2D7");
            var r = Resolve(VariantOf(consensus, "A0A087X1C5_S70N")).Single();

            Assert.That(r.Outcome, Is.EqualTo(GeneResolutionOutcome.NotInSource));
            Assert.That((r.Accession, r.EntryAccession, r.Namespace),
                Is.EqualTo(("A0A087X1C5_S70N", "A0A087X1C5", AccessionNamespace.UniProt)),
                "verbatim accession kept; entry and grammar come from the consensus entry");
            Assert.That(r.UniProtGeneName, Is.EqualTo("CYP2D7"));
        }

        [Test]
        public void SequenceVariant_ResolvesThroughItsConsensusEntrysGeneLinks()
        {
            var consensus = Entry("P11111", references: Transcript("ENST00000000001.1", "ENSG00000000001.7"));
            var r = Resolve(VariantOf(consensus, "P11111_S70N", references: new List<DatabaseReference>())).Single();

            Assert.That((r.Outcome, r.GeneId, r.EntryAccession),
                Is.EqualTo((GeneResolutionOutcome.Resolved, "ENSG00000000001", "P11111")));
        }

        [Test]
        public void EveryRow_CarriesItsProvenance()
        {
            var rows = new[]
            {
                Resolve(Entry("P11111", references: Transcript("ENST00000000001.1", "ENSG00000000001.7"))).Single(),
                Resolve(Entry("P22222")).Single(),
            };

            foreach (var r in rows)
            {
                Assert.That(r.SearchDatabaseSha256, Is.EqualTo(SearchDbSha));
                Assert.That(r.Source, Is.EqualTo(EnsemblGeneResolver.SearchDatabaseSource));
                Assert.That(r.GeneSetRelease, Is.EqualTo("116"));
                Assert.That(r.GeneSetSha256, Is.EqualTo(_resolver.GeneSet.SourceSha256));
            }
        }

        [Test]
        public void ResolveAll_GivesEveryProteinAtLeastOneRow()
        {
            var proteins = new[]
            {
                Entry("P11111", references: new[] { Transcript("ENST1.1", "ENSG00000000001.7"), Transcript("ENST2.1", "ENSG00000000002.7") }),
                Entry("P22222"),
                Entry("DECOY_P33333", isDecoy: true),
            };

            var rows = _resolver.ResolveAll(proteins, SearchDbSha).ToList();

            Assert.That(rows.Select(r => r.Accession).Distinct(), Is.EquivalentTo(new[] { "P11111", "P22222", "DECOY_P33333" }));
            Assert.That(rows.Count, Is.EqualTo(4));
        }

        [Test]
        public void RealUniProtXml_Gapdh_Resolves()
        {
            var gapdh = ProteinDbLoader
                .LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "humanGAPDH.xml"),
                    true, DecoyType.None, null, false, null, out _)
                .First(p => p.Accession == "P04406");

            var r = Resolve(gapdh).Single();

            Assert.That((r.Outcome, r.GeneId, r.VersionedGeneId, r.GeneSymbol),
                Is.EqualTo((GeneResolutionOutcome.Resolved, "ENSG00000111640", "ENSG00000111640.15", "GAPDH")));
        }

        [Test]
        public void WriteTsv_OneRowPerResolution_NullsAsEmptyCells_OutcomeAlwaysPresent()
        {
            var rows = new[]
            {
                Resolve(Entry("P11111", references: new[] { Transcript("ENST1.1", "ENSG00000000001.7"), Transcript("ENST2.1", "ENSG00000000002.7") })),
                Resolve(Entry("P22222", "NAMEONLY")),
            }.SelectMany(r => r).ToList();

            var output = new StringWriter();
            GeneResolutionTsv.Write(output, rows);
            var lines = output.ToString().TrimEnd('\n', '\r').Split('\n').Select(l => l.TrimEnd('\r')).ToList();

            var header = lines[0].Split('\t');
            Assert.That(header, Does.Contain("accession").And.Contain("outcome").And.Contain("gene_id")
                .And.Contain("versioned_gene_id").And.Contain("search_database_sha256"));
            Assert.That(lines.Count, Is.EqualTo(1 + 3), "one row per (accession, gene); no '|'-joined cells");

            var notInSource = lines.Skip(1).Select(l => header.Zip(l.Split('\t')).ToDictionary(p => p.First, p => p.Second))
                .Single(c => c["accession"] == "P22222");
            Assert.That(notInSource["outcome"], Is.EqualTo("not_in_source"));
            Assert.That(notInSource["gene_id"], Is.EqualTo(""));
            Assert.That(notInSource["uniprot_gene_name"], Is.EqualTo("NAMEONLY"));
        }
    }
}
