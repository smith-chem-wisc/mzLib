using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases.GeneOntology;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The two files a GO annotation run writes: the (group, term) rows with a "#!" header whose five
    /// counters summarise the run, and one term-to-category table per consumer map, joined to the rows on
    /// (go_id, go_release). These are data-interchange files: their column names and header keys are a
    /// contract with whoever ingests them, and the ingester recounts the header from the rows.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestGoAnnotationTsv
    {
        private static readonly string OntologyPath =
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Ontologies", "go-trimmed.obo");

        private const string DbSha = "0123abcd";
        private static GeneOntologyGraph _go;
        private string _dir;

        [OneTimeSetUp]
        public void LoadOntology() => _go = GeneOntologyGraph.Load(OntologyPath);

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "GoAnnotationTsv_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private static DatabaseReference Go(string id, params string[] evidence) =>
            new("GO", id, evidence.Select(e => Tuple.Create("evidence", e)));

        private static Protein P(string accession, params DatabaseReference[] refs) =>
            new("PEPTIDEK", accession, databaseReferences: refs.ToList());

        private static GoAnnotationGroup Group(string members, double q, bool contaminant = false) =>
            new(members, members.Split('|'), false, contaminant, q);

        private static GoGroupAnnotator Annotator() => new(_go,
            new[] { P("P1", Go("GO:0005634", "ECO:0000314")), P("P2", Go("GO:0005634", "ECO:0000501")), P("P3") }, DbSha);

        private static string WriteAnnotation(IEnumerable<GoAnnotationRow> rows, string sourceSha = null)
        {
            using var writer = new StringWriter();
            GoAnnotationTsv.Write(writer, rows, _go, DbSha, sourceSha);
            return writer.ToString();
        }

        private static Dictionary<string, string> HeaderOf(string text) => text.Split('\n')
            .Where(l => l.StartsWith("#!", StringComparison.Ordinal))
            .Select(l => l.Substring(2).Split(' ', 2))
            .ToDictionary(p => p[0], p => p[1]);

        private static List<Dictionary<string, string>> RowsOf(string text)
        {
            var lines = text.Split('\n').Where(l => l.Length > 0 && !l.StartsWith("#!", StringComparison.Ordinal)).ToList();
            string[] columns = lines[0].Split('\t');
            return lines.Skip(1).Select(l => l.Split('\t'))
                .Select(cells => columns.Zip(cells).ToDictionary(p => p.First, p => p.Second)).ToList();
        }

        [Test]
        public void Header_FormatVersionFirst_ThenProvenance()
        {
            string text = WriteAnnotation(Annotator().Annotate(Group("P1", 0.001)), sourceSha: "feedbeef");
            string[] lines = text.Split('\n');

            Assert.That(lines[0], Is.EqualTo("#!go_annotation_format 1"));
            var header = HeaderOf(text);
            Assert.That(header["go_release"], Is.EqualTo("releases/2026-07-26"));
            Assert.That(header["go_obo_sha256"], Is.EqualTo(_go.SourceSha256));
            Assert.That(header["annotation_db_sha256"], Is.EqualTo(DbSha));
            Assert.That(header["source_file_sha256"], Is.EqualTo("feedbeef"));
            Assert.That(header["counter_q_value_max"], Is.EqualTo("0.01"));
        }

        [Test]
        public void Header_NoSourceFile_OmitsThatLine()
        {
            Assert.That(HeaderOf(WriteAnnotation(Annotator().Annotate(Group("P1", 0.001)))).ContainsKey("source_file_sha256"), Is.False);
        }

        [Test]
        public void Header_HasExactlyFiveCounters()
        {
            var keys = HeaderOf(WriteAnnotation(Annotator().Annotate(Group("P1", 0.001)))).Keys
                .Where(k => k.StartsWith("n_", StringComparison.Ordinal) || k.StartsWith("status_", StringComparison.Ordinal));

            Assert.That(keys, Is.EquivalentTo(new[]
            {
                "n_multi_member_groups", "status_annotated", "status_no_go_terms", "status_no_entry", "status_contaminant"
            }));
        }

        [Test]
        public void Header_CountersMatchRowCounts_AtQValue001()
        {
            var groups = new[]
            {
                Group("P1|P2", 0.001),                       // annotated, two members
                Group("P3", 0.01),                           // no GO terms, exactly at the threshold: counted
                Group("Q9", 0.005),                          // no entry
                Group("Q8", 0.002, contaminant: true),       // contaminant
                Group("P1|P3", 0.2),                         // annotated, multi-member, but above 0.01: not counted
            };
            string text = WriteAnnotation(Annotator().AnnotateAll(groups));
            var header = HeaderOf(text);

            Assert.That(header["n_multi_member_groups"], Is.EqualTo("1"));
            Assert.That(header["status_annotated"], Is.EqualTo("1"));
            Assert.That(header["status_no_go_terms"], Is.EqualTo("1"));
            Assert.That(header["status_no_entry"], Is.EqualTo("1"));
            Assert.That(header["status_contaminant"], Is.EqualTo("1"));

            // What an ingester does: recount from the rows, one vote per group, at the stated threshold.
            var perGroup = RowsOf(text)
                .Where(r => double.Parse(r["q_value"], System.Globalization.CultureInfo.InvariantCulture) <= 0.01)
                .GroupBy(r => r["protein_group"]).Select(g => g.First()).ToList();
            foreach (var status in new[] { "annotated", "no_go_terms", "no_entry", "contaminant" })
            {
                Assert.That(perGroup.Count(r => r["annotation_status"] == status).ToString(), Is.EqualTo(header["status_" + status]));
            }
            Assert.That(perGroup.Count(r => r["n_members"] != "1").ToString(), Is.EqualTo(header["n_multi_member_groups"]));
        }

        [Test]
        public void Rows_IncludeGroupsAboveTheCounterThreshold()
        {
            // Every non-decoy group is in the file whatever its q-value; the consumer filters (go D29).
            var rows = RowsOf(WriteAnnotation(Annotator().AnnotateAll(new[] { Group("P1", 0.001), Group("P3", 0.9) })));

            Assert.That(rows.Select(r => r["protein_group"]).Distinct(), Is.EqualTo(new[] { "P1", "P3" }));
        }

        [Test]
        public void Rows_ColumnsAndCellFormats()
        {
            string text = WriteAnnotation(Annotator().AnnotateAll(new[] { Group("P1|P2", 0.001), Group("P3", 0.0042) }));
            string columnLine = text.Split('\n').First(l => !l.StartsWith("#!", StringComparison.Ordinal));

            Assert.That(columnLine, Is.EqualTo(string.Join('\t', "protein_group", "accession_used", "go_id", "go_name", "aspect",
                "evidence", "inherited", "propagated", "n_members", "n_with", "annotation_status", "q_value",
                "go_release", "go_obo_sha256", "annotation_db_sha256")));

            var rows = RowsOf(text);
            var nucleus = rows.Single(r => r["go_id"] == "GO:0005634");
            Assert.That(nucleus["accession_used"], Is.EqualTo("P1;P2"));
            Assert.That(nucleus["go_name"], Is.EqualTo("nucleus"));
            Assert.That(nucleus["aspect"], Is.EqualTo("cellular_component"));
            Assert.That(nucleus["evidence"], Is.EqualTo("ECO:0000314;ECO:0000501"));
            Assert.That(nucleus["inherited"], Is.EqualTo("false"));
            Assert.That(nucleus["propagated"], Is.EqualTo("false"));
            Assert.That(nucleus["n_members"], Is.EqualTo("2"));
            Assert.That(nucleus["n_with"], Is.EqualTo("2"));
            Assert.That(nucleus["annotation_status"], Is.EqualTo("annotated"));
            Assert.That(nucleus["q_value"], Is.EqualTo("0.001"));
            Assert.That(rows.Single(r => r["go_id"] == "GO:0043226")["propagated"], Is.EqualTo("true"));

            // A term-less row keeps its width: empty cells, never a missing column.
            var termless = rows.Single(r => r["protein_group"] == "P3");
            Assert.That(termless["go_id"], Is.Empty);
            Assert.That(termless["accession_used"], Is.Empty);
            Assert.That(termless["inherited"], Is.Empty);
            Assert.That(termless["propagated"], Is.Empty);
            Assert.That(termless["n_with"], Is.EqualTo("0"));
            Assert.That(termless["annotation_status"], Is.EqualTo("no_go_terms"));
            Assert.That(termless["q_value"], Is.EqualTo("0.0042"));
        }

        [Test]
        public void Rows_GroupsKeepInputOrder_TermsWithinAGroupAreOrdinal()
        {
            var rows = RowsOf(WriteAnnotation(Annotator().AnnotateAll(new[] { Group("P3", 0.001), Group("P1", 0.001) })));

            Assert.That(rows.Select(r => r["protein_group"]).Distinct(), Is.EqualTo(new[] { "P3", "P1" }));
            var ids = rows.Where(r => r["protein_group"] == "P1").Select(r => r["go_id"]).ToList();
            Assert.That(ids, Is.EqualTo(ids.OrderBy(i => i, StringComparer.Ordinal).ToList()));
        }

        [Test]
        public void Output_UsesLineFeedOnly()
        {
            Assert.That(WriteAnnotation(Annotator().Annotate(Group("P1", 0.001))), Does.Not.Contain("\r"));
        }

        [Test]
        public void Writer_RejectsRowsFromAnotherRelease()
        {
            var row = Annotator().Annotate(Group("P1", 0.001))[0] with { GoRelease = "releases/1999-01-01" };

            var ex = Assert.Throws<ArgumentException>(() => WriteAnnotation(new[] { row }));
            Assert.That(ex.Message, Does.Contain("releases/1999-01-01"));
        }

        [Test]
        public void Writer_RejectsRowsFromAnotherGoObo_EvenWithTheSameReleaseName()
        {
            // Two files can claim the same data-version; the hash is what says they are the same ontology.
            var row = Annotator().Annotate(Group("P1", 0.001))[0] with { GoOboSha256 = "different" };

            Assert.Throws<ArgumentException>(() => WriteAnnotation(new[] { row }));
        }

        [Test]
        public void Writer_RejectsSemicolonInsideAnEvidenceCode()
        {
            var row = Annotator().Annotate(Group("P1", 0.001))[0] with { Evidence = new[] { "ECO:1;ECO:2" } };

            Assert.Throws<ArgumentException>(() => WriteAnnotation(new[] { row }));
        }

        [Test]
        public void Writer_RejectsTabInAHeaderValue()
        {
            var rows = Annotator().Annotate(Group("P1", 0.001));

            Assert.Throws<ArgumentException>(() => WriteAnnotation(rows, sourceSha: "ab\tcd"));
        }

        [Test]
        public void Rows_AspectIsWrittenAsGosOwnNamespaceName()
        {
            var annotator = new GoGroupAnnotator(_go, new[] { P("P1", Go("GO:0004365"), Go("GO:0044238")) }, DbSha);
            var rows = RowsOf(WriteAnnotation(annotator.Annotate(Group("P1", 0.001))));

            Assert.That(rows.Single(r => r["go_id"] == "GO:0004365")["aspect"], Is.EqualTo("molecular_function"));
            Assert.That(rows.Single(r => r["go_id"] == "GO:0044238")["aspect"], Is.EqualTo("biological_process"));

            var unknown = Annotator().Annotate(Group("P1", 0.001))[0] with { Aspect = GoAspect.Unknown };
            Assert.That(RowsOf(WriteAnnotation(new[] { unknown }))[0]["aspect"], Is.EqualTo("unknown"));
        }

        [Test]
        public void Rows_NullSet_IsAnEmptyCell()
        {
            var row = Annotator().Annotate(Group("P3", 0.001))[0] with { AccessionUsed = null, Evidence = null };

            var written = RowsOf(WriteAnnotation(new[] { row }))[0];
            Assert.That(written["accession_used"], Is.Empty);
            Assert.That(written["evidence"], Is.Empty);
        }

        [Test]
        public void Writer_RejectsRowsFromAnotherAnnotationDatabase()
        {
            var row = Annotator().Annotate(Group("P1", 0.001))[0] with { AnnotationDbSha256 = "other" };

            Assert.Throws<ArgumentException>(() => WriteAnnotation(new[] { row }));
        }

        [TestCase("P1\tX")]
        [TestCase("P1\nX")]
        [TestCase("\tP1")]
        public void Writer_RejectsTabOrNewlineInValue(string groupName)
        {
            var row = Annotator().Annotate(Group("P1", 0.001))[0] with { ProteinGroup = groupName };

            Assert.Throws<ArgumentException>(() => WriteAnnotation(new[] { row }));
        }

        [Test]
        public void Writer_RejectsSemicolonInsideASetMember()
        {
            // accession_used and evidence are ';'-joined; a member containing ';' would read as two.
            var row = Annotator().Annotate(Group("P1", 0.001))[0] with { AccessionUsed = new[] { "P1;P2" } };

            Assert.Throws<ArgumentException>(() => WriteAnnotation(new[] { row }));
        }

        [Test]
        public void Writer_NoRows_WritesHeaderAndZeroCounters()
        {
            string text = WriteAnnotation(Array.Empty<GoAnnotationRow>());

            Assert.That(HeaderOf(text)["status_annotated"], Is.EqualTo("0"));
            Assert.That(RowsOf(text), Is.Empty);
        }

        [Test]
        public void Writer_NullArguments_Throw()
        {
            using var writer = new StringWriter();
            Assert.Throws<ArgumentNullException>(() => GoAnnotationTsv.Write(null, Array.Empty<GoAnnotationRow>(), _go, DbSha));
            Assert.Throws<ArgumentNullException>(() => GoAnnotationTsv.Write(writer, null, _go, DbSha));
            Assert.Throws<ArgumentNullException>(() => GoAnnotationTsv.Write(writer, Array.Empty<GoAnnotationRow>(), null, DbSha));
        }

        // ---- the per-map term table ----

        private GoCategoryResolver Resolver(string rows, string name = "test_map", string version = "3")
        {
            string path = Path.Combine(_dir, "map.tsv");
            File.WriteAllText(path, $"#!category_map_format 1\n#!map_name {name}\n#!map_version {version}\n" +
                                    "category\tsubcategory\tanchor_go_id\n" + rows);
            return new GoCategoryResolver(GoCategoryMap.Load(path), _go);
        }

        private static string WriteCategories(GoCategoryResolver resolver, IEnumerable<GoAnnotationRow> rows)
        {
            using var writer = new StringWriter();
            GoCategoryTsv.Write(writer, resolver, rows);
            return writer.ToString();
        }

        [Test]
        public void CategoryTable_HeaderCarriesTheJoinKeyAndTheMap()
        {
            var resolver = Resolver("nuc\t\tGO:0005634\n");
            string text = WriteCategories(resolver, Annotator().Annotate(Group("P1", 0.001)));
            string[] lines = text.Split('\n');

            Assert.That(lines[0], Is.EqualTo("#!go_category_format 1"));
            var header = HeaderOf(text);
            Assert.That(header["go_release"], Is.EqualTo("releases/2026-07-26"));
            Assert.That(header["go_obo_sha256"], Is.EqualTo(_go.SourceSha256));
            Assert.That(header["category_map"], Is.EqualTo("test_map 3 " + resolver.Map.SourceSha256));
            Assert.That(lines.First(l => !l.StartsWith("#!", StringComparison.Ordinal)), Is.EqualTo("go_id\tcategory\tsubcategory"));
            Assert.That(text, Does.Not.Contain("\r"));
        }

        [Test]
        public void CategoryTable_CoversEveryAnnotatedTermUnderAnAnchor_AndNothingElse()
        {
            // P1 carries nucleus directly; its 7 ancestors come by propagation. Only terms at or below an anchor
            // get a row: nucleus, organelle and organelle's three subtypes on the way down to nucleus. The CC
            // root, "intracellular anatomical structure" and "cellular anatomical structure" are above every
            // anchor and get none.
            var resolver = Resolver("nuc\t\tGO:0005634\norg\t\tGO:0043226\n");
            var rows = RowsOf(WriteCategories(resolver, Annotator().Annotate(Group("P1", 0.001))));

            Assert.That(rows.Select(r => r["go_id"] + "/" + r["category"]), Is.EqualTo(new[]
            {
                "GO:0005634/nuc", "GO:0005634/org", "GO:0043226/org", "GO:0043227/org", "GO:0043229/org", "GO:0043231/org"
            }));
            Assert.That(rows.All(r => r["subcategory"] == ""), Is.True);
        }

        [Test]
        public void CategoryTable_TermUnderTwoAnchors_IsTwoRows_SubcategoryQualified()
        {
            var resolver = Resolver("org\tnuclear\tGO:0005634\nmembrane_bound\t\tGO:0043227\n");
            var rows = RowsOf(WriteCategories(resolver, Annotator().Annotate(Group("P1", 0.001))))
                .Where(r => r["go_id"] == "GO:0005634").ToList();

            Assert.That(rows.Select(r => r["category"] + "/" + r["subcategory"]),
                Is.EqualTo(new[] { "membrane_bound/", "org/org:nuclear" }));
        }

        [Test]
        public void CategoryTable_TermlessRowsContributeNothing()
        {
            var rows = RowsOf(WriteCategories(Resolver("nuc\t\tGO:0005634\n"), Annotator().Annotate(Group("P3", 0.001))));

            Assert.That(rows, Is.Empty);
        }

        [Test]
        public void CategoryTable_GoReleaseMatchesTheAnnotationFile()
        {
            // The join key is (go_id, go_release). A table computed against one release and joined to rows from
            // another would silently mis-assign, so the writer refuses.
            var row = Annotator().Annotate(Group("P1", 0.001))[0] with { GoRelease = "releases/1999-01-01" };

            Assert.Throws<ArgumentException>(() => WriteCategories(Resolver("nuc\t\tGO:0005634\n"), new[] { row }));
        }

        [Test]
        public void CategoryTable_NullArguments_Throw()
        {
            using var writer = new StringWriter();
            var resolver = Resolver("nuc\t\tGO:0005634\n");
            Assert.Throws<ArgumentNullException>(() => GoCategoryTsv.Write(null, resolver, Array.Empty<GoAnnotationRow>()));
            Assert.Throws<ArgumentNullException>(() => GoCategoryTsv.Write(writer, null, Array.Empty<GoAnnotationRow>()));
            Assert.Throws<ArgumentNullException>(() => GoCategoryTsv.Write(writer, resolver, null));
        }

        [TestCase("two words", "3")]
        [TestCase("m", "v 1")]
        public void CategoryMap_WhitespaceInNameOrVersion_Throws(string name, string version)
        {
            // The table header writes "#!category_map <name> <version> <sha256>", split on spaces.
            Assert.Throws<InvalidDataException>(() => Resolver("nuc\t\tGO:0005634\n", name, version));
        }
    }
}
