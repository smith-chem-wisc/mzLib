using NUnit.Framework;
using Proteomics;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases.GeneOntology;

namespace Test.FileReadingTests.InternalFileReading
{
    /// <summary>
    /// A stored MetaMorpheus protein-group table, turned into the group descriptors the GO annotator takes.
    /// The fixture is #1347's six PXD036557 rows: four targets, one contaminant (P02769) and one decoy.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestProteinGroupFromTsvExtensions
    {
        private static string FixturePath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\ExternalFileTypes\MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv");

        private static string OntologyPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "DataFiles", "Ontologies", "go-trimmed.obo");

        private string _outputDirectory = "";

        [OneTimeSetUp]
        public void SetUp()
        {
            _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinGroupFromTsvExtensions");
            Directory.CreateDirectory(_outputDirectory);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            if (Directory.Exists(_outputDirectory))
                Directory.Delete(_outputDirectory, true);
        }

        private string WriteTable(string name, params string[] rows)
        {
            string path = Path.Combine(_outputDirectory, name);
            File.WriteAllText(path, "Protein Accession\tProtein Decoy/Contaminant/Target\tProtein QValue\n" +
                                    string.Join("\n", rows) + "\n");
            return path;
        }

        [Test]
        public void ToGoAnnotationGroup_CarriesNameMembersFlagsAndQValue()
        {
            var row = new ProteinGroupFromTsvFile(FixturePath).Single(r => r.ProteinGroupName == "P0C0S5|Q71UI9");

            var group = row.ToGoAnnotationGroup();

            Assert.That(group.Name, Is.EqualTo("P0C0S5|Q71UI9"));
            Assert.That(group.MemberAccessions, Is.EqualTo(new[] { "P0C0S5", "Q71UI9" }));
            Assert.That(group.IsDecoy, Is.False);
            Assert.That(group.IsContaminant, Is.False);
            Assert.That(group.QValue, Is.EqualTo(row.QValue));
        }

        [Test]
        public void ToGoAnnotationGroup_Null_Throws()
        {
            Assert.Throws<ArgumentNullException>(() => ((ProteinGroupFromTsv)null!).ToGoAnnotationGroup());
            Assert.Throws<ArgumentNullException>(() => ((IEnumerable<ProteinGroupFromTsv>)null!).ToGoAnnotationGroups().ToList());
        }

        [Test]
        public void Adapter_Pxd036557Fixture_EveryNonDecoyGroupAnnotated()
        {
            var groups = new ProteinGroupFromTsvFile(FixturePath).ToGoAnnotationGroups().ToList();

            Assert.That(groups.Select(g => g.Name),
                Is.EqualTo(new[] { "P68363", "P05141", "P0C0S5|Q71UI9", "P02769", "P63104" }));
            Assert.That(groups.Single(g => g.IsContaminant).Name, Is.EqualTo("P02769"));

            // One member carries a term, the rest have no entry: every group must still come out with a row.
            var go = GeneOntologyGraph.Load(OntologyPath);
            var proteins = new[]
            {
                new Protein("PEPTIDEK", "P68363",
                    databaseReferences: new List<DatabaseReference> { new("GO", "GO:0005634", new[] { Tuple.Create("evidence", "ECO:0000314") }) })
            };
            var rows = new GoGroupAnnotator(go, proteins, "0123abcd").AnnotateAll(groups).ToList();

            Assert.That(rows.Select(r => r.ProteinGroup).Distinct(), Is.EquivalentTo(groups.Select(g => g.Name)));
            Assert.That(rows.Where(r => r.ProteinGroup == "P68363").Select(r => r.Status).Distinct(),
                Is.EqualTo(new[] { GoAnnotationStatus.Annotated }));
            Assert.That(rows.Single(r => r.ProteinGroup == "P02769").Status, Is.EqualTo(GoAnnotationStatus.Contaminant));
        }

        /// <summary>
        /// The reader's IsDecoy is "the label contains D", so an entrapment decoy (ED) is a decoy and never
        /// reaches the annotator, which refuses decoys. An entrapment target (ET) is not a decoy and is kept:
        /// the descriptor has no entrapment field, so it is annotated like any other target.
        /// </summary>
        [Test]
        public void Adapter_EntrapmentDecoy_IsRejected()
        {
            string path = WriteTable("Entrapment_AllProteinGroups.tsv",
                "P1\tT\t0.001", "P2\tC\t0.001", "DECOY_P3\tD\t0.5", "P4\tED\t0.5", "P5\tET\t0.002");

            var groups = new ProteinGroupFromTsvFile(path).ToGoAnnotationGroups().ToList();

            Assert.That(groups.Select(g => g.Name), Is.EqualTo(new[] { "P1", "P2", "P5" }));
            Assert.That(groups.Select(g => g.IsContaminant), Is.EqualTo(new[] { false, true, false }));
            Assert.That(groups.All(g => !g.IsDecoy));
        }

        [Test]
        public void ToGoAnnotationGroup_KeepsTheDecoyFlag_SoTheAnnotatorCanRefuseIt()
        {
            var decoy = new ProteinGroupFromTsvFile(FixturePath).Single(r => r.IsDecoy).ToGoAnnotationGroup();
            var annotator = new GoGroupAnnotator(GeneOntologyGraph.Load(OntologyPath), Array.Empty<Protein>(), "0123abcd");

            Assert.That(decoy.IsDecoy, Is.True);
            Assert.Throws<ArgumentException>(() => annotator.Annotate(decoy).ToList());
        }

        [Test]
        public void AQuantFreeAllProteinGroupsTable_IsReadThroughTheConstructor()
        {
            string path = WriteTable("AllProteinGroups.tsv", "P1|P2\tT\t0.004");

            var group = new ProteinGroupFromTsvFile(path).ToGoAnnotationGroups().Single();

            Assert.That(group.MemberAccessions, Is.EqualTo(new[] { "P1", "P2" }));
            Assert.That(group.QValue, Is.EqualTo(0.004));
        }
    }
}
