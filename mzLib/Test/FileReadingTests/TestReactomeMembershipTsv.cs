using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using UsefulProteomicsDatabases.Pathways;
using static Test.FileReadingTests.TestReactomeMembership;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The stored form of a database's Reactome membership (<see cref="ReactomeMembershipTsv"/>): a strict headed table
    /// that round-trips, names its database, definition and species check, and refuses anything else.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestReactomeMembershipTsv
    {
        private const string Head =
            "#!reactome_membership_format 1\n" +
            "#!definition kegg:DEF-REACTOME-MEMBERSHIP v1\n" +
            "#!annotation_db_sha256 0123abcd\n" +
            "accession\torganism\tpathway_id\tmolecule\tsequence_accession\tinherited\tstatus\n";

        private string _dir;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "ReactomeMembership_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private string WriteFile(string name, Action<TextWriter> write)
        {
            string path = Path.Combine(_dir, name);
            using (var writer = new StreamWriter(path)) write(writer);
            return path;
        }

        private string WriteText(string text) => WriteFile("t1.tsv", w => w.Write(text));

        [Test]
        public void WriteThenRead_RoundTrips()
        {
            var annotation = new ReactomeAccessionAnnotation(Fixture(), DbSha);

            var table = ReactomeMembershipTsv.Read(WriteFile("t1.tsv", w => ReactomeMembershipTsv.Write(w, annotation)));

            Assert.That(table.Rows, Is.EqualTo(annotation.Rows()));
            Assert.That(table.SourceFileName, Is.EqualTo("t1.tsv"));
        }

        [Test]
        public void Header_CarriesShaDefinitionTaxonCodesAndProducer()
        {
            var annotation = new ReactomeAccessionAnnotation(Fixture(), DbSha);

            var table = ReactomeMembershipTsv.Read(WriteFile("t1.tsv", w => ReactomeMembershipTsv.Write(w, annotation)));

            Assert.That(table.Header["reactome_membership_format"], Is.EqualTo("1"));
            Assert.That(table.Header["definition"], Is.EqualTo("kegg:DEF-REACTOME-MEMBERSHIP v1"));
            Assert.That(table.Header["annotation_db_sha256"], Is.EqualTo(DbSha));
            Assert.That(table.Header["species_codes"], Is.EqualTo("NCBITaxon:9606=HSA"));
            Assert.That(table.Header["unvalidated_organisms"], Is.EqualTo("NCBITaxon:10116=3;NCBITaxon:559292=1;NCBITaxon:9031=12"));
            Assert.That(table.Header["dropped_foreign_references"], Is.EqualTo("1"));
            Assert.That(table.Header.ContainsKey("mzlib_version"), Is.True);
            Assert.That(table.Header.ContainsKey("mzlib_release"), Is.True);
        }

        [Test]
        public void Header_EmptyCountsAndMapsAreOmitted()
        {
            var annotation = new ReactomeAccessionAnnotation(new[] { P("Q11111", Mouse) }, DbSha);

            var table = ReactomeMembershipTsv.Read(WriteFile("t1.tsv", w => ReactomeMembershipTsv.Write(w, annotation)));

            Assert.That(table.Header.Keys, Does.Not.Contain("species_codes").And.Not.Contain("unvalidated_organisms")
                .And.Not.Contain("dropped_foreign_references"));
            Assert.That(table.Rows.Single().Status, Is.EqualTo(ReactomeMembershipStatus.SpeciesNotCovered));
        }

        [Test]
        public void Header_IsTheFormatFirstAndEndsLinesWithLf()
        {
            var annotation = new ReactomeAccessionAnnotation(Fixture(), DbSha);
            using var writer = new StringWriter();

            ReactomeMembershipTsv.Write(writer, annotation);

            string text = writer.ToString();
            Assert.That(text, Does.StartWith("#!reactome_membership_format 1\n"));
            Assert.That(text, Does.Not.Contain("\r"));
            Assert.That(text, Does.Contain("\nP55555\tNCBITaxon:9606\t\t\t\t\tno_reference\n"));
            Assert.That(text, Does.Contain("\nP22222\tNCBITaxon:9606\tR-HSA-2001\tP22222-2\tP22222\tfalse\tannotated\n"));
        }

        [TestCase(ReactomeMembershipStatus.Annotated, "annotated")]
        [TestCase(ReactomeMembershipStatus.NoReference, "no_reference")]
        [TestCase(ReactomeMembershipStatus.SpeciesNotCovered, "species_not_covered")]
        public void StatusName_IsLowercaseSnake(ReactomeMembershipStatus status, string name)
        {
            Assert.That(ReactomeMembershipTsv.StatusName(status), Is.EqualTo(name));
        }

        [Test]
        public void StatusName_RefusesAnUndefinedValue()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ReactomeMembershipTsv.StatusName((ReactomeMembershipStatus)99));
        }

        [Test]
        public void Writer_RejectsNullArguments()
        {
            var annotation = new ReactomeAccessionAnnotation(Fixture(), DbSha);

            Assert.Throws<ArgumentNullException>(() => ReactomeMembershipTsv.Write(null, annotation));
            Assert.Throws<ArgumentNullException>(() => ReactomeMembershipTsv.Write(new StringWriter(), null));
        }

        [Test]
        public void Writer_RefusesASeparatorInAValue()
        {
            var annotation = new ReactomeAccessionAnnotation(new[] { P("P1\t1", Human, Reactome("R-HSA-1")) }, DbSha);

            Assert.Throws<ArgumentException>(() => ReactomeMembershipTsv.Write(new StringWriter(), annotation));
        }

        [Test]
        public void Writer_RefusesAnOrganismItsReaderWouldRefuse_BeforeWritingAnything()
        {
            // A taxon id that is not all digits becomes NCBITaxon:9606a, which Read refuses. P1 is valid and sorts
            // first, so an empty output shows every row was checked before the first line was written.
            var annotation = new ReactomeAccessionAnnotation(new[]
            {
                P("P1", Human, Reactome("R-HSA-1")),
                P("P2", "9606a", Reactome("R-HSA-2")),
            }, DbSha);
            using var writer = new StringWriter();

            var exception = Assert.Throws<ArgumentException>(() => ReactomeMembershipTsv.Write(writer, annotation));
            Assert.That(exception.Message, Does.Contain("NCBITaxon:9606a"));
            Assert.That(writer.ToString(), Is.Empty);
        }

        [Test]
        public void Writer_RefusesAnEmptyAccession()
        {
            var annotation = new ReactomeAccessionAnnotation(new[] { P("", Human, Reactome("R-HSA-1")) }, DbSha);
            using var writer = new StringWriter();

            var exception = Assert.Throws<ArgumentException>(() => ReactomeMembershipTsv.Write(writer, annotation));
            Assert.That(exception.Message, Does.Contain("accession is empty"));
            Assert.That(writer.ToString(), Is.Empty);
        }

        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\tannotated_ish\n", "status 'annotated_ish'")]
        [TestCase("P1\tNCBITaxon:9606\t\t\t\t\tannotated\n", "carries no pathway")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\tno_reference\n", "carries a pathway")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\tannotated\nP1\tNCBITaxon:9606\t\t\t\t\tno_reference\n", "pathway-less row and another row")]
        [TestCase("P1\tNCBITaxon:9606\t\t\t\t\tno_reference\nP1\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\tannotated\n", "pathway-less row and another row")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\tP1-2\tQ9\tfalse\tannotated\n", "sequence_accession 'Q9'")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\tP1\tfalse\tannotated\n", "sequence_accession 'P1'")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\tP1-2\t\tfalse\tannotated\n", "sequence_accession ''")]
        [TestCase("P1-2\tNCBITaxon:9606\tR-HSA-1\tP1-3\tP1-3\ttrue\tannotated\n", "inherited row is isoform-scoped")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\t\t\tannotated\n", "inherited is empty")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\t\tyes\tannotated\n", "inherited 'yes'")]
        [TestCase("P1\tNCBITaxon:9606\t\t\t\tfalse\tno_reference\n", "pathway-less row of 'P1' carries a value")]
        [TestCase("P1\tNCBITaxon:9606\t\tP1-2\t\t\tno_reference\n", "pathway-less row of 'P1' carries a value")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\tannotated\nP1\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\tannotated\n", "appears twice")]
        [TestCase("\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\tannotated\n", "accession is empty")]
        [TestCase("P1\t9606\tR-HSA-1\t\t\tfalse\tannotated\n", "organism '9606'")]
        [TestCase("P1\tNCBITaxon:9606\tR-HSA-1\t\t\tfalse\n", "expected 7 cells")]
        public void Readers_Refuse(string rows, string fragment)
        {
            string path = WriteText(Head + rows);

            var exception = Assert.Throws<InvalidDataException>(() => ReactomeMembershipTsv.Read(path));
            Assert.That(exception.Message, Does.Contain(fragment));
        }

        [TestCase("#!reactome_membership_format 2\n#!definition kegg:DEF-REACTOME-MEMBERSHIP v1\n#!annotation_db_sha256 x\n", "must be 1")]
        [TestCase("#!reactome_membership_format 1\n#!annotation_db_sha256 x\n", "'definition' is missing")]
        [TestCase("#!reactome_membership_format 1\n#!definition kegg:DEF-REACTOME-MEMBERSHIP v1\n", "'annotation_db_sha256' is missing")]
        [TestCase("#!reactome_membership_format 1\n#!definition kegg:DEF-REACTOME-MEMBERSHIP v2\n#!annotation_db_sha256 x\n", "definition")]
        [TestCase("#!reactome_membership_format 1\n#!definition kegg:DEF-REACTOME-MEMBERSHIP v1\n#!annotation_db_sha256 x\n#!organism NCBITaxon:9606\n", "unknown header key 'organism'")]
        public void Readers_RefuseAHeader(string header, string fragment)
        {
            string path = WriteText(header + string.Join('\t', ReactomeMembershipTsv.Columns) + "\n");

            var exception = Assert.Throws<InvalidDataException>(() => ReactomeMembershipTsv.Read(path));
            Assert.That(exception.Message, Does.Contain(fragment));
        }

        [Test]
        public void Readers_AcceptAHandWrittenTable()
        {
            string path = WriteText(Head +
                "P1\tNCBITaxon:9606\tR-HSA-1\tP1-2\tP1\tfalse\tannotated\n" +
                "P1\tNCBITaxon:9606\tR-HSA-1\tP1-3\tP1-3\tfalse\tannotated\n" +
                "P1-4\tNCBITaxon:9606\tR-HSA-2\t\t\ttrue\tannotated\n" +
                "P2\t\t\t\t\t\tspecies_not_covered\n");

            var rows = ReactomeMembershipTsv.Read(path).Rows;

            Assert.That(rows, Is.EqualTo(new[]
            {
                new ReactomeMembershipRow("P1", "NCBITaxon:9606", "R-HSA-1", "P1-2", "P1", false, ReactomeMembershipStatus.Annotated),
                new ReactomeMembershipRow("P1", "NCBITaxon:9606", "R-HSA-1", "P1-3", "P1-3", false, ReactomeMembershipStatus.Annotated),
                new ReactomeMembershipRow("P1-4", "NCBITaxon:9606", "R-HSA-2", null, null, true, ReactomeMembershipStatus.Annotated),
                new ReactomeMembershipRow("P2", null, null, null, null, null, ReactomeMembershipStatus.SpeciesNotCovered),
            }));
        }

        [Test]
        public void Readers_MissingFile_Throws()
        {
            Assert.Throws<FileNotFoundException>(() => ReactomeMembershipTsv.Read(Path.Combine(_dir, "absent.tsv")));
        }
    }
}
