using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using UsefulProteomicsDatabases.Pathways;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The three tables that depend on the Reactome release alone, written once per release, as GoOntologyTsv does
    /// for go.obo: pathways (id, name, species, human equivalent), hierarchy (parent, child) and ancestors (every
    /// pathway with each ancestor AND itself, so a join of direct memberships through it keeps the direct pathway).
    /// Each file carries its format, the mzLib build that wrote it, the release and both source files' sha256s in
    /// its header; the readers refuse anything else, with the line number, because a table read wrongly is joined
    /// wrongly.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestReactomeReleaseTsv
    {
        private string _tempDirectory;

        [SetUp]
        public void SetUp()
        {
            _tempDirectory = Path.Combine(Path.GetTempPath(), "TestReactomeReleaseTsv_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_tempDirectory);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_tempDirectory))
            {
                Directory.Delete(_tempDirectory, true);
            }
        }

        private string WriteFile(string name, Action<TextWriter> write)
        {
            string path = Path.Combine(_tempDirectory, name);
            using (var writer = new StreamWriter(path))
            {
                write(writer);
            }
            return path;
        }

        private string WriteText(string text)
        {
            string path = Path.Combine(_tempDirectory, "table.tsv");
            File.WriteAllText(path, text);
            return path;
        }

        private const string Shas = "#!reactome_release 96\n#!reactome_pathways_sha256 aa\n#!reactome_relation_sha256 bb\n";

        [Test]
        public void Pathways_WriteThenRead_RoundTrips()
        {
            var release = TestReactomeRelease.LoadFixture();

            var table = ReactomeReleaseTsv.ReadPathways(WriteFile("pathways.tsv", w => ReactomeReleaseTsv.WritePathways(w, release)));

            Assert.That(table.Rows, Is.EqualTo(release.Pathways));
            Assert.That(table.Rows.Single(r => r.PathwayId == "R-GGA-351470").HumanEquivalentId, Is.Null, "an empty cell reads back as null");
        }

        [Test]
        public void Hierarchy_WriteThenRead_RoundTrips()
        {
            var release = TestReactomeRelease.LoadFixture();

            var table = ReactomeReleaseTsv.ReadHierarchy(WriteFile("hierarchy.tsv", w => ReactomeReleaseTsv.WriteHierarchy(w, release)));

            Assert.That(table.Rows, Is.EqualTo(release.Hierarchy));
            Assert.That(table.Rows.Count, Is.EqualTo(19));
        }

        [Test]
        public void Ancestors_WriteThenRead_RoundTripsAndIncludeSelf()
        {
            var release = TestReactomeRelease.LoadFixture();

            var table = ReactomeReleaseTsv.ReadAncestors(WriteFile("ancestors.tsv", w => ReactomeReleaseTsv.WriteAncestors(w, release)));

            Assert.That(table.Rows.Count, Is.EqualTo(63), "25 self pairs plus 38 proper ancestors");
            Assert.That(table.Rows, Does.Contain(new ReactomeAncestorRow("R-HSA-77289", "R-HSA-77289")));
            Assert.That(table.Rows, Does.Contain(new ReactomeAncestorRow("R-HSA-77289", "R-HSA-1430728")));
            Assert.That(table.Rows.Where(r => r.PathwayId == "R-HSA-354192").Select(r => r.AncestorPathwayId),
                Is.EquivalentTo(new[] { "R-HSA-354192", "R-HSA-162582", "R-HSA-76009", "R-HSA-76002", "R-HSA-109582" }));
            foreach (var pathway in release.Pathways)
            {
                Assert.That(table.Rows.Where(r => r.PathwayId == pathway.PathwayId).Select(r => r.AncestorPathwayId),
                    Is.EquivalentTo(release.Ancestors(pathway.PathwayId).Append(pathway.PathwayId)));
            }
        }

        [Test]
        public void Header_CarriesFormatProducerShasAndRelease()
        {
            var release = TestReactomeRelease.LoadFixture();
            string path = WriteFile("pathways.tsv", w => ReactomeReleaseTsv.WritePathways(w, release));

            var header = ReactomeReleaseTsv.ReadPathways(path).Header;

            Assert.That(File.ReadLines(path).First(), Is.EqualTo("#!reactome_pathways_format 1"));
            Assert.That(header["reactome_release"], Is.EqualTo("96"));
            Assert.That(header["reactome_pathways_sha256"], Is.EqualTo(release.PathwaysSha256));
            Assert.That(header["reactome_relation_sha256"], Is.EqualTo(release.RelationSha256));
            Assert.That(header.ContainsKey("mzlib_version"), Is.True);
            Assert.That(header.ContainsKey("mzlib_release"), Is.True);
            Assert.That(File.ReadLines(WriteFile("hierarchy.tsv", w => ReactomeReleaseTsv.WriteHierarchy(w, release))).First(), Is.EqualTo("#!reactome_hierarchy_format 1"));
            Assert.That(File.ReadLines(WriteFile("ancestors.tsv", w => ReactomeReleaseTsv.WriteAncestors(w, release))).First(), Is.EqualTo("#!reactome_ancestors_format 1"));
        }

        [TestCase("pathways", "#!go_terms_format 1\n" + Shas + "pathway_id\tname\tspecies_name\thuman_equivalent_id\n", "must be '#!reactome_pathways_format'")]
        [TestCase("pathways", "#!reactome_pathways_format 1\n#!reactome_release 96\n#!reactome_pathways_sha256 aa\npathway_id\tname\tspecies_name\thuman_equivalent_id\n", "reactome_relation_sha256")]
        [TestCase("pathways", "#!reactome_pathways_format 1\n" + Shas + "pathway_id\tname\tspecies_name\thuman_equivalent_id\nR-HSA-1\tA\tHomo sapiens\tR-HSA-1\nR-HSA-1\tA\tHomo sapiens\tR-HSA-1\n", "repeats")]
        [TestCase("pathways", "#!reactome_pathways_format 1\n" + Shas + "pathway_id\tname\tspecies_name\thuman_equivalent_id\n\tA\tHomo sapiens\t\n", "empty")]
        [TestCase("pathways", "#!reactome_pathways_format 2\n" + Shas + "pathway_id\tname\tspecies_name\thuman_equivalent_id\n", "must be 1")]
        [TestCase("hierarchy", "#!reactome_hierarchy_format 1\n" + Shas + "parent_id\tchild_id\nR-HSA-1\tR-HSA-2\nR-HSA-1\tR-HSA-2\n", "appears twice")]
        [TestCase("hierarchy", "#!reactome_hierarchy_format 1\n" + Shas + "parent_id\tchild_id\nR-HSA-1\n", "expected 2 cells")]
        [TestCase("hierarchy", "#!reactome_hierarchy_format 1\n" + Shas + "parent_id\tchild\n", "column header")]
        [TestCase("ancestors", "#!reactome_ancestors_format 1\n" + Shas + "pathway_id\tancestor_pathway_id\nR-HSA-2\tR-HSA-1\n", "has no (R-HSA-2, R-HSA-2) pair")]
        [TestCase("ancestors", "#!reactome_ancestors_format 1\n" + Shas + "pathway_id\tancestor_pathway_id\nR-HSA-2\t\n", "empty")]
        public void Readers_Refuse(string table, string text, string expected)
        {
            string path = WriteText(text);
            Func<object> read = table switch
            {
                "pathways" => () => ReactomeReleaseTsv.ReadPathways(path),
                "hierarchy" => () => ReactomeReleaseTsv.ReadHierarchy(path),
                _ => () => ReactomeReleaseTsv.ReadAncestors(path),
            };

            var error = Assert.Throws<InvalidDataException>(() => read());

            Assert.That(error.Message, Does.Contain(expected));
        }

        [Test]
        public void Readers_MissingFile_Throws()
        {
            string absent = Path.Combine(_tempDirectory, "absent.tsv");

            Assert.Throws<FileNotFoundException>(() => ReactomeReleaseTsv.ReadPathways(absent));
            Assert.Throws<FileNotFoundException>(() => ReactomeReleaseTsv.ReadHierarchy(absent));
            Assert.Throws<FileNotFoundException>(() => ReactomeReleaseTsv.ReadAncestors(absent));
        }

        [Test]
        public void Writers_RejectNullArguments()
        {
            var release = TestReactomeRelease.LoadFixture();

            Assert.Throws<ArgumentNullException>(() => ReactomeReleaseTsv.WritePathways(null, release));
            Assert.Throws<ArgumentNullException>(() => ReactomeReleaseTsv.WriteHierarchy(TextWriter.Null, null));
            Assert.Throws<ArgumentNullException>(() => ReactomeReleaseTsv.WriteAncestors(TextWriter.Null, null));
        }
    }
}
