using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Security.Cryptography;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.GeneOntology;

namespace Test.FileReadingTests
{
    /// <summary>
    /// GeneOntologyGraph is go.obo read with its edges kept. Protein.GoTerms says which terms a protein
    /// was annotated with; this says what those terms imply. A protein annotated to "mitochondrial inner
    /// membrane" is in the mitochondrion only because the ontology says the one is part_of the other --
    /// drop part_of and that fact is silently lost, which is why most of these tests exist.
    ///
    /// The fixture is a 29-term slice of the real go.obo (release 2026-07-26), closed under is_a and
    /// part_of so every ancestor is present. See DataFiles\Ontologies\PROVENANCE.md for how it was cut.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestGeneOntologyGraph
    {
        private static readonly string FixturePath =
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Ontologies", "go-trimmed.obo");

        private const string MitochondrialInnerMembrane = "GO:0005743";
        private const string Mitochondrion = "GO:0005739";

        private string _dir;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "GeneOntologyGraph_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private string Write(string name, string text)
        {
            string path = Path.Combine(_dir, name);
            File.WriteAllText(path, text);
            return path;
        }

        private static string Term(string id, params string[] tags) =>
            "[Term]\nid: " + id + "\nname: " + id + "\nnamespace: cellular_component\n" +
            string.Concat(tags.Select(t => t + "\n")) + "\n";

        [Test]
        public void Obo_ParsesTermsAndVersion()
        {
            var go = GeneOntologyGraph.Load(FixturePath);

            // The data-version stamp is what makes a pin auditable, so it must never come back empty.
            Assert.That(go.Release, Is.EqualTo("releases/2026-07-26"));
            Assert.That(go.SourceFileName, Is.EqualTo("go-trimmed.obo"));
            Assert.That(go.Count, Is.EqualTo(29));

            Assert.That(go.TryGetTerm(MitochondrialInnerMembrane, out var term), Is.True);
            Assert.That(term.Name, Is.EqualTo("mitochondrial inner membrane"));
            Assert.That(term.Aspect, Is.EqualTo(GoAspect.CellularComponent));
            Assert.That(go.TryGetTerm("GO:0004365", out var gapdh), Is.True);
            Assert.That(gapdh.Aspect, Is.EqualTo(GoAspect.MolecularFunction));
            Assert.That(go.TryGetTerm("GO:0008152", out var metabolism), Is.True);
            Assert.That(metabolism.Aspect, Is.EqualTo(GoAspect.BiologicalProcess));
        }

        [Test]
        public void Obo_ReportsSha256OfTheBytesRead()
        {
            string expected;
            using (var stream = File.OpenRead(FixturePath))
            {
                expected = Convert.ToHexString(SHA256.HashData(stream)).ToLowerInvariant();
            }

            Assert.That(GeneOntologyGraph.Load(FixturePath).SourceSha256, Is.EqualTo(expected));
        }

        [Test]
        public void Obo_AltIdRemapsToPrimary()
        {
            var go = GeneOntologyGraph.Load(FixturePath);

            // GO:0016021 "integral component of membrane" was merged into GO:0016020 "membrane". A UniProt
            // entry older than the merge still cites the old id, and must still land on the term.
            Assert.That(go.TryGetTerm("GO:0016021", out var term), Is.True);
            Assert.That(term.Id, Is.EqualTo("GO:0016020"));
            Assert.That(term.AltIds, Does.Contain("GO:0016021"));
            Assert.That(go.Ancestors("GO:0016021"), Is.EquivalentTo(go.Ancestors("GO:0016020")));
        }

        [Test]
        public void Obo_ObsoleteTermsAreFlagged_NotDropped()
        {
            var go = GeneOntologyGraph.Load(FixturePath);

            // A file citing an obsolete term must still be readable, and the reader needs to know where
            // GO moved it. Dropping the term would make the citation indistinguishable from a typo.
            Assert.That(go.TryGetTerm("GO:0006082", out var term), Is.True);
            Assert.That(term.IsObsolete, Is.True);
            Assert.That(term.ReplacedBy, Is.EqualTo(new[] { "GO:0008152" }));
            Assert.That(go.Ancestors("GO:0006082"), Is.Empty, "an obsolete term has no place in the DAG");
        }

        [Test]
        public void Obo_TypedefStanzasAreNotTerms()
        {
            var go = GeneOntologyGraph.Load(FixturePath);

            Assert.That(go.TryGetTerm("part_of", out _), Is.False);
            Assert.That(go.TryGetTerm("has_part", out _), Is.False);
            Assert.That(go.TermIds.All(id => id.StartsWith("GO:", StringComparison.Ordinal)), Is.True);
        }

        [Test]
        public void Obo_MalformedTagLine_ThrowsWithLineNumber()
        {
            string path = Write("bad.obo", "format-version: 1.2\ndata-version: releases/x\n\n" +
                                           "[Term]\nid: GO:0000001\nthis line has no tag separator\n");

            var ex = Assert.Throws<InvalidDataException>(() => GeneOntologyGraph.Load(path));
            Assert.That(ex.Message, Does.Contain("bad.obo line 6"));
        }

        [Test]
        public void Obo_TermWithoutId_Throws()
        {
            string path = Write("noid.obo", "data-version: x\n\n[Term]\nname: orphan\nnamespace: cellular_component\n");

            Assert.Throws<InvalidDataException>(() => GeneOntologyGraph.Load(path));
        }

        [Test]
        public void Obo_DuplicateTermId_Throws()
        {
            string path = Write("dup.obo", "data-version: x\n\n" + Term("GO:0000001") + Term("GO:0000001"));

            var ex = Assert.Throws<InvalidDataException>(() => GeneOntologyGraph.Load(path));
            Assert.That(ex.Message, Does.Contain("GO:0000001"));
        }

        [Test]
        public void Obo_EdgeToAnAbsentTerm_Throws()
        {
            // An is_a or part_of pointing outside the file would make every ancestor query above it
            // silently short. A truncated download looks exactly like this.
            string path = Write("dangling.obo", "data-version: x\n\n" + Term("GO:0000001", "is_a: GO:9999999 ! gone"));

            var ex = Assert.Throws<InvalidDataException>(() => GeneOntologyGraph.Load(path));
            Assert.That(ex.Message, Does.Contain("GO:9999999"));
        }

        [Test]
        public void Obo_MissingFile_ThrowsFileNotFound()
        {
            Assert.Throws<FileNotFoundException>(() => GeneOntologyGraph.Load(Path.Combine(_dir, "absent.obo")));
        }

        [Test]
        public void Ancestors_FollowsIsA_AndPartOf()
        {
            var go = GeneOntologyGraph.Load(FixturePath);
            var ancestors = go.Ancestors(MitochondrialInnerMembrane);

            // Mitochondrion is reachable ONLY through part_of (inner membrane is_a organelle inner membrane /
            // mitochondrial membrane; mitochondrial membrane part_of mitochondrial envelope part_of
            // mitochondrion). An is_a-only walk finds 6 ancestors; the real answer is 15.
            Assert.That(ancestors, Does.Contain(Mitochondrion));
            Assert.That(ancestors, Does.Contain("GO:0005575"), "reaches the cellular_component root");
            Assert.That(ancestors, Has.Count.EqualTo(15));
            Assert.That(ancestors, Does.Not.Contain(MitochondrialInnerMembrane), "a term is not its own ancestor");
        }

        [Test]
        public void Ancestors_TermWithTwoParents_ReturnsBoth()
        {
            var go = GeneOntologyGraph.Load(FixturePath);

            // GO is a DAG, not a tree: inner membrane has two is_a parents. This is the fact that rules out
            // a single "distance to the term" column (go D22) -- there is no one path to measure.
            Assert.That(go.TryGetTerm(MitochondrialInnerMembrane, out var term), Is.True);
            Assert.That(term.IsAParents, Is.EquivalentTo(new[] { "GO:0019866", "GO:0031966" }));
            Assert.That(go.Ancestors(MitochondrialInnerMembrane), Does.Contain("GO:0019866").And.Contain("GO:0031966"));
        }

        [Test]
        public void Ancestors_TerminatesOnCycle()
        {
            string path = Write("cycle.obo", "data-version: x\n\n" +
                                             Term("GO:0000001", "is_a: GO:0000002") +
                                             Term("GO:0000002", "relationship: part_of GO:0000001 ! back"));
            var go = GeneOntologyGraph.Load(path);

            Assert.That(go.Ancestors("GO:0000001"), Is.EquivalentTo(new[] { "GO:0000002" }));
        }

        [Test]
        public void Ancestors_OtherRelationshipsAreIgnored()
        {
            // has_part points DOWN the hierarchy (a whole has a part), and regulates is not containment at
            // all. Following either would make an annotation imply terms it does not.
            string path = Write("rels.obo", "data-version: x\n\n" +
                                            Term("GO:0000001", "relationship: has_part GO:0000002", "relationship: regulates GO:0000003") +
                                            Term("GO:0000002") + Term("GO:0000003"));

            Assert.That(GeneOntologyGraph.Load(path).Ancestors("GO:0000001"), Is.Empty);
        }

        [Test]
        public void Ancestors_UnknownId_Throws()
        {
            var go = GeneOntologyGraph.Load(FixturePath);

            // A UniProt entry can cite a term newer than the pinned release. That has to be loud: an empty
            // ancestor set would read as "a root", not as "not in this release".
            Assert.Throws<ArgumentException>(() => go.Ancestors("GO:9999999"));
        }

        [Test]
        public void LoadGeneOntology_ExistingFile_DoesNotDownload()
        {
            string path = Path.Combine(_dir, "go.obo");
            File.Copy(FixturePath, path);

            var go = Loaders.LoadGeneOntology(path);

            Assert.That(go.Count, Is.EqualTo(29));
        }

        [Test]
        public void UpdateGeneOntology_StaleTempFromACrashedRun_DoesNotBlockTheUpdate()
        {
            string path = Path.Combine(_dir, "go.obo");
            File.WriteAllText(path + ".temp", "half a download from a run that crashed");
            using var client = new HttpClient(new FixedResponseHandler(File.ReadAllText(FixturePath)));

            Loaders.UpdateGeneOntology(path, client);

            Assert.That(File.Exists(path + ".temp"), Is.False);
            Assert.That(GeneOntologyGraph.Load(path).Count, Is.EqualTo(29));
        }

        [Test]
        public void UpdateGeneOntology_ChangedFile_KeepsTheOldOneAsBackup()
        {
            string path = Path.Combine(_dir, "go.obo");
            File.WriteAllText(path, "data-version: releases/old\n");
            using var client = new HttpClient(new FixedResponseHandler(File.ReadAllText(FixturePath)));

            Loaders.UpdateGeneOntology(path, client);

            Assert.That(GeneOntologyGraph.Load(path).Release, Is.EqualTo("releases/2026-07-26"));
            Assert.That(Directory.GetFiles(_dir, "go.obo*"), Has.Length.EqualTo(2), "the replaced release is kept");
        }

        private sealed class FixedResponseHandler : HttpMessageHandler
        {
            private readonly string _body;
            public FixedResponseHandler(string body) => _body = body;

            protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken) =>
                Task.FromResult(new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent(_body) });
        }
    }
}
