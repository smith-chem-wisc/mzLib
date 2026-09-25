using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using NUnit.Framework;
using UsefulProteomicsDatabases.GeneOntology;

namespace Test.FileReadingTests
{
    /// <summary>
    /// A category map is a consumer's file, never mzLib's: it names categories (organelles, complexes,
    /// anything) by anchor GO terms, and a term belongs to a category when the anchor is the term itself
    /// or one of its ancestors. mzLib owns the format and the rule; the consumer owns the rows. These
    /// tests use a deliberately neutral map over the go-trimmed.obo fixture, not any consumer's real one.
    ///
    /// Fixture terms used: GO:0005743 mitochondrial inner membrane, whose ancestors include GO:0005739
    /// mitochondrion and GO:0016020 membrane; GO:0005634 nucleus; GO:0004365 GAPDH activity (MF).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestGoCategoryMap
    {
        private static readonly string OntologyPath =
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Ontologies", "go-trimmed.obo");

        private const string Header = "#!category_map_format 1\n#!map_name test_compartments\n#!map_version 3\n" +
                                      "category\tsubcategory\tanchor_go_id\n";

        private const string Rows =
            "mito\t\tGO:0005739\n" +
            "mito\tinner\tGO:0005743\n" +
            "nuc\t\tGO:0005634\n" +
            "membrane\t\tGO:0016020\n";

        private static GeneOntologyGraph _go;
        private string _dir;

        [OneTimeSetUp]
        public void LoadOntology() => _go = GeneOntologyGraph.Load(OntologyPath);

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "GoCategoryMap_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private string Write(string text, string name = "map.tsv")
        {
            string path = Path.Combine(_dir, name);
            File.WriteAllText(path, text);
            return path;
        }

        private GoCategoryResolver Resolver(string rows = Rows) => new(GoCategoryMap.Load(Write(Header + rows)), _go);

        private static string[] Labels(GoCategoryResolver resolver, string goId) =>
            resolver.Categorize(goId).Select(c => c.Category + "/" + (c.Subcategory ?? "")).ToArray();

        [Test]
        public void CategoryMap_LoadsNameVersionAndSha256()
        {
            string path = Write(Header + Rows);
            string expected;
            using (var stream = File.OpenRead(path))
            {
                expected = Convert.ToHexString(SHA256.HashData(stream)).ToLowerInvariant();
            }

            var map = GoCategoryMap.Load(path);

            Assert.That(map.MapName, Is.EqualTo("test_compartments"));
            Assert.That(map.MapVersion, Is.EqualTo("3"));
            Assert.That(map.SourceSha256, Is.EqualTo(expected));
            Assert.That(map.SourceFileName, Is.EqualTo("map.tsv"));
            Assert.That(map.Anchors, Has.Count.EqualTo(4));
            Assert.That(map.Anchors[1], Is.EqualTo(new GoCategoryAnchor("mito", "inner", "GO:0005743")));
            Assert.That(map.Anchors[0].Subcategory, Is.Null, "an empty cell is no subcategory, not an empty one");
        }

        [TestCase("#!category_map_format 2\n#!map_name m\n#!map_version 1\ncategory\tsubcategory\tanchor_go_id\n", "format")]
        [TestCase("#!map_name m\n#!map_version 1\ncategory\tsubcategory\tanchor_go_id\n", "format")]
        [TestCase("#!category_map_format 1\n#!map_version 1\ncategory\tsubcategory\tanchor_go_id\n", "map_name")]
        [TestCase("#!category_map_format 1\n#!map_name m\ncategory\tsubcategory\tanchor_go_id\n", "map_version")]
        [TestCase("#!category_map_format 1\n#!map_name m\n#!map_name n\n#!map_version 1\ncategory\tsubcategory\tanchor_go_id\n", "map_name")]
        [TestCase("#!category_map_format 1\n#!map_name m\n#!map_version 1\n#!colour blue\ncategory\tsubcategory\tanchor_go_id\n", "colour")]
        [TestCase("#!category_map_format 1\n#!map_name m\n#!map_version 1\ncategory\tanchor_go_id\tsubcategory\n", "header")]
        public void CategoryMap_BadHeader_Throws(string text, string mentioned)
        {
            var ex = Assert.Throws<InvalidDataException>(() => GoCategoryMap.Load(Write(text)));
            Assert.That(ex.Message, Does.Contain(mentioned));
        }

        [TestCase("mito\tGO:0005739\n", "line 5")]
        [TestCase("mito\t\tGO:0005739\textra\n", "line 5")]
        [TestCase("\t\tGO:0005739\n", "line 5")]
        [TestCase("mito\t\t\n", "line 5")]
        [TestCase("mito\t\tmitochondrion\n", "line 5")]
        [TestCase("mito:x\t\tGO:0005739\n", "line 5")]
        [TestCase("mito\tinner:x\tGO:0005743\n", "line 5")]
        [TestCase("mito\t\tGO:0005739\nmito\t\tGO:0005739\n", "line 6")]
        public void CategoryMap_BadRow_ThrowsWithLineNumber(string rows, string line)
        {
            var ex = Assert.Throws<InvalidDataException>(() => GoCategoryMap.Load(Write(Header + rows)));
            Assert.That(ex.Message, Does.Contain("map.tsv " + line));
        }

        [TestCase("#!category_map_format 1\n#!map_name m\n#!map_version\ncategory\tsubcategory\tanchor_go_id\n", "no value")]
        [TestCase("#!category_map_format 1\n#!map_name m\n#!map_version 1\n", "no column header")]
        [TestCase("#!map_name m\n#!map_version 1\ncategory\tsubcategory\tanchor_go_id\n", "(missing)")]
        public void CategoryMap_IncompleteHeader_SaysWhatIsMissing(string text, string mentioned)
        {
            var ex = Assert.Throws<InvalidDataException>(() => GoCategoryMap.Load(Write(text)));
            Assert.That(ex.Message, Does.Contain(mentioned));
        }

        [Test]
        public void CategoryMap_BlankLines_AreIgnored()
        {
            var map = GoCategoryMap.Load(Write(Header.Replace("#!map_version 3\n", "#!map_version 3\n\n") + "\n" + Rows + "\n"));

            Assert.That(map.Anchors, Has.Count.EqualTo(4));
        }

        [Test]
        public void CategoryMap_TwoAnchorsForOneSubcategory_AreAllowed()
        {
            // A category is often reached by two unrelated terms; the rows then share a label.
            var map = GoCategoryMap.Load(Write(Header + "nuc\tenvelope\tGO:0005634\nnuc\tenvelope\tGO:0016020\n"));

            Assert.That(map.Anchors, Has.Count.EqualTo(2));
        }

        [Test]
        public void CategoryMap_MissingFile_ThrowsFileNotFound()
        {
            Assert.Throws<FileNotFoundException>(() => GoCategoryMap.Load(Path.Combine(_dir, "absent.tsv")));
        }

        [Test]
        public void Category_TakesMostSpecificAnchor_SelfOrAncestor()
        {
            var resolver = Resolver();

            // The anchor itself counts (self), and the more specific of two anchors in one category wins:
            // inner membrane is under both mito anchors, and reports only mito/inner.
            Assert.That(Labels(resolver, "GO:0005739"), Is.EquivalentTo(new[] { "mito/" }));
            Assert.That(Labels(resolver, "GO:0005743"), Does.Contain("mito/mito:inner").And.Not.Contain("mito/"));
        }

        [Test]
        public void Category_IsSetValued_TermUnderTwoAnchors()
        {
            // Inner membrane is a mitochondrial part AND a membrane. Both are reported: "most specific"
            // applies within a category, never across categories.
            Assert.That(Labels(Resolver(), "GO:0005743"), Is.EquivalentTo(new[] { "mito/mito:inner", "membrane/" }));
        }

        [Test]
        public void Subcategory_MembersAreSelfQualifying_SoSetColumnsNeverPairByPosition()
        {
            var categories = Resolver().Categorize("GO:0005743");

            // "inner" on its own could belong to any category; "mito:inner" cannot be misread when two
            // set-valued columns are laid side by side.
            Assert.That(categories.Where(c => c.Subcategory != null).Select(c => c.Subcategory), Is.EqualTo(new[] { "mito:inner" }));
        }

        [Test]
        public void TermOutsideEveryAnchor_HasNoCategory()
        {
            // Empty by construction: GAPDH activity is under no anchor. Nothing about MF is special-cased --
            // a consumer's map may anchor any aspect.
            Assert.That(Resolver().Categorize("GO:0004365"), Is.Empty);
        }

        [Test]
        public void Category_AnchorOnAnotherAspect_Works()
        {
            var resolver = Resolver(Rows + "redox\t\tGO:0016491\n");

            Assert.That(Labels(resolver, "GO:0004365"), Is.EquivalentTo(new[] { "redox/" }));
        }

        [Test]
        public void Category_IsDeterministicallyOrdered()
        {
            var labels = Labels(Resolver(), "GO:0005743");

            Assert.That(labels, Is.EqualTo(new[] { "membrane/", "mito/mito:inner" }), "ordinal by category, then subcategory");
        }

        [Test]
        public void Category_TwoUnrelatedSubcategoriesOfOneCategory_AreBothKept_InOrdinalOrder()
        {
            // Inner membrane lies under mitochondrial envelope and under organelle inner membrane, and neither
            // of those is above the other, so both subcategories survive; the cellular_component root is above
            // both and is dropped as less specific. Rows are deliberately written out of order.
            var resolver = Resolver("x\tzeta\tGO:0005740\nx\talpha\tGO:0019866\nx\t\tGO:0005575\n");

            Assert.That(Labels(resolver, "GO:0005743"), Is.EqualTo(new[] { "x/x:alpha", "x/x:zeta" }));
        }

        [Test]
        public void Category_UnknownAnchorInAnUnversionedRelease_SaysSo()
        {
            string obo = Path.Combine(_dir, "unversioned.obo");
            File.WriteAllText(obo, "format-version: 1.2\n\n[Term]\nid: GO:0000001\nname: a\nis_obsolete: true\n");
            var go = GeneOntologyGraph.Load(obo);

            var missing = Assert.Throws<InvalidDataException>(() =>
                new GoCategoryResolver(GoCategoryMap.Load(Write(Header + "gone\t\tGO:9999999\n", "a.tsv")), go));
            var obsolete = Assert.Throws<InvalidDataException>(() =>
                new GoCategoryResolver(GoCategoryMap.Load(Write(Header + "old\t\tGO:0000001\n", "b.tsv")), go));

            Assert.That(missing.Message, Does.Contain("(unversioned)").And.Contain("a.tsv"));
            Assert.That(obsolete.Message, Does.Contain("(unversioned)").And.Contain("b.tsv"));
        }

        [Test]
        public void Category_AltIdTermResolvesLikeItsPrimary()
        {
            var resolver = Resolver();

            Assert.That(Labels(resolver, "GO:0016021"), Is.EqualTo(Labels(resolver, "GO:0016020")));
        }

        [Test]
        public void Category_UnknownTerm_Throws()
        {
            Assert.Throws<ArgumentException>(() => Resolver().Categorize("GO:9999999"));
        }

        [Test]
        public void AnchorNotInOntology_Throws()
        {
            // A map citing a term absent from the loaded release is a mismatch between map and release. It
            // must be loud: silently, the category would just never be assigned.
            var map = GoCategoryMap.Load(Write(Header + "gone\t\tGO:9999999\n"));

            var ex = Assert.Throws<InvalidDataException>(() => new GoCategoryResolver(map, _go));
            Assert.That(ex.Message, Does.Contain("GO:9999999").And.Contain("releases/2026-07-26"));
        }

        [Test]
        public void AnchorOnObsoleteTerm_Throws()
        {
            // An obsolete term is no one's ancestor, so an anchor on it could only ever match itself.
            var map = GoCategoryMap.Load(Write(Header + "old\t\tGO:0006082\n"));

            var ex = Assert.Throws<InvalidDataException>(() => new GoCategoryResolver(map, _go));
            Assert.That(ex.Message, Does.Contain("GO:0006082").And.Contain("obsolete"));
        }

        [Test]
        public void AnchorOnAltId_IsResolvedToItsPrimary()
        {
            var resolver = Resolver("membrane\t\tGO:0016021\n");

            Assert.That(Labels(resolver, "GO:0005743"), Is.EquivalentTo(new[] { "membrane/" }));
        }

        [Test]
        public void Resolver_CarriesItsProvenance()
        {
            var resolver = Resolver();

            Assert.That(resolver.Map.MapName, Is.EqualTo("test_compartments"));
            Assert.That(resolver.Ontology.Release, Is.EqualTo("releases/2026-07-26"));
        }
    }
}
