using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Gene Ontology terms are already parsed and stored -- ProteinXmlEntry keeps every
    /// &lt;dbReference&gt; generically, so a UniProt XML's GO annotations arrive on
    /// Protein.DatabaseReferences whether or not anything asks for them. Nothing READS them as GO.
    ///
    /// Protein.GoTerms is a derived view over that list, following NcbiTaxonomyId exactly: a
    /// computed property and a const on Protein, no new constructor parameter and no stored field,
    /// so a Protein read from any source exposes the same thing.
    ///
    /// Three traps these tests pin, all of them real in the data rather than hypothetical:
    ///  - the flat list holds far more than GO (humanGAPDH.xml: 386 dbReferences, 39 of them GO,
    ///    86 PubMed), and the parser applies no depth guard, so filtering by Type is mandatory;
    ///  - UniProt repeats one GO id with different evidence codes, and ProteinDbLoader's duplicate
    ///    -entry merge unions the references rather than collapsing them, so the accessor must
    ///    dedupe by id and union the evidence;
    ///  - properties are a bag of (type, value) tuples in no guaranteed order, so they must be
    ///    matched by type and never by position.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestGeneOntology
    {
        private static string Data(params string[] parts) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, Path.Combine(parts));

        private static Protein LoadGapdh() => ProteinDbLoader
            .LoadProteinXML(Data("DatabaseTests", "humanGAPDH.xml"), true, DecoyType.None, null, false, null, out _)
            .First(p => p.Accession == "P04406");

        /// <summary>
        /// Builds a Protein carrying exactly the supplied references, so a test can express a
        /// shape the fixture does not happen to contain.
        /// </summary>
        private static Protein WithReferences(params DatabaseReference[] references) =>
            new Protein("PEPTIDEK", "P00001", databaseReferences: references.ToList());

        private static DatabaseReference Go(string id, params (string Type, string Value)[] properties) =>
            new DatabaseReference("GO", id,
                properties.Select(p => new Tuple<string, string>(p.Type, p.Value)).ToList());

        [Test]
        public void GoTerms_ParsedFromUniProtXml()
        {
            var protein = LoadGapdh();

            Assert.That(protein.GoTerms.Count, Is.EqualTo(39), "humanGAPDH.xml carries 39 GO dbReferences");

            var cytoplasm = protein.GoTerms.Single(t => t.Id == "GO:0005737");
            Assert.That(cytoplasm.Name, Is.EqualTo("cytoplasm"), "the aspect prefix is not part of the name");
            Assert.That(cytoplasm.Aspect, Is.EqualTo(GoAspect.CellularComponent));
            Assert.That(cytoplasm.EvidenceCodes, Is.EquivalentTo(new[] { "ECO:0000314" }));
            Assert.That(cytoplasm.Projects, Is.EquivalentTo(new[] { "UniProtKB" }));
        }

        [Test]
        public void GoTerms_AspectReadFromTheTermPrefix()
        {
            var protein = LoadGapdh();

            // C:/F:/P: prefixes on the term property are the only place the aspect appears.
            Assert.That(protein.GoTerms.Count(t => t.Aspect == GoAspect.CellularComponent), Is.EqualTo(14));
            Assert.That(protein.GoTerms.Count(t => t.Aspect == GoAspect.MolecularFunction), Is.EqualTo(8));
            Assert.That(protein.GoTerms.Count(t => t.Aspect == GoAspect.BiologicalProcess), Is.EqualTo(17));
            Assert.That(protein.GoTerms.Any(t => t.Aspect == GoAspect.Unknown), Is.False);
        }

        [Test]
        public void GoTerms_UnprefixedTermIsUnknownAspectAndKeepsItsWholeName()
        {
            // Absence of a prefix must not silently eat the first character of the name, and must
            // not be guessed at -- an organelle claim built on a guessed aspect is worse than none.
            var protein = WithReferences(Go("GO:0000001", ("term", "mitochondrion inheritance")));

            var term = protein.GoTerms.Single();
            Assert.That(term.Aspect, Is.EqualTo(GoAspect.Unknown));
            Assert.That(term.Name, Is.EqualTo("mitochondrion inheritance"));
        }

        [Test]
        public void GoTerms_ExcludeNonGoDatabaseReferences()
        {
            var protein = LoadGapdh();

            // The flat list also holds PubMed, DOI, EMBL, ChEBI and organism references, and the
            // parser applies no depth guard, so nested references arrive here too.
            Assert.That(protein.DatabaseReferences.Count, Is.GreaterThan(300));
            Assert.That(protein.DatabaseReferences.Any(r => r.Type == "PubMed"), Is.True,
                "if this fixture ever loses its PubMed refs the test stops proving anything");
            Assert.That(protein.GoTerms.Select(t => t.Id), Is.All.StartWith("GO:"));
            Assert.That(protein.GoTerms.Count, Is.LessThan(protein.DatabaseReferences.Count));
        }

        [Test]
        public void GoTerms_DedupedById_EvidenceAndProjectsUnioned()
        {
            // The same GO id recurs with different ECO codes, and ProteinDbLoader's duplicate-entry
            // merge unions the reference lists rather than collapsing them.
            var protein = WithReferences(
                Go("GO:0005737", ("term", "C:cytoplasm"), ("evidence", "ECO:0000314"), ("project", "UniProtKB")),
                Go("GO:0005737", ("term", "C:cytoplasm"), ("evidence", "ECO:0007005"), ("project", "HPA")),
                Go("GO:0005829", ("term", "C:cytosol"), ("evidence", "ECO:0000314"), ("project", "MGI")));

            Assert.That(protein.GoTerms.Count, Is.EqualTo(2), "one entry per distinct GO id");

            var cytoplasm = protein.GoTerms.Single(t => t.Id == "GO:0005737");
            Assert.That(cytoplasm.EvidenceCodes, Is.EquivalentTo(new[] { "ECO:0000314", "ECO:0007005" }));
            Assert.That(cytoplasm.Projects, Is.EquivalentTo(new[] { "UniProtKB", "HPA" }));
        }

        [Test]
        public void GoTerms_RepeatedEvidenceCodeIsNotCountedTwice()
        {
            var protein = WithReferences(
                Go("GO:0005737", ("term", "C:cytoplasm"), ("evidence", "ECO:0000314"), ("project", "UniProtKB")),
                Go("GO:0005737", ("term", "C:cytoplasm"), ("evidence", "ECO:0000314"), ("project", "UniProtKB")));

            var term = protein.GoTerms.Single();
            Assert.That(term.EvidenceCodes, Is.EquivalentTo(new[] { "ECO:0000314" }));
            Assert.That(term.Projects, Is.EquivalentTo(new[] { "UniProtKB" }));
        }

        [Test]
        public void GoTerms_PropertiesMatchedByType_NeverByPosition()
        {
            // Nothing guarantees property order; ProteinDbWriter re-sorts them on write
            // (evidence, project, term), which is already a different order from the XML.
            var protein = WithReferences(
                Go("GO:0005739", ("project", "UniProtKB"), ("evidence", "ECO:0000314"), ("term", "C:mitochondrion")));

            var term = protein.GoTerms.Single();
            Assert.That(term.Name, Is.EqualTo("mitochondrion"));
            Assert.That(term.Aspect, Is.EqualTo(GoAspect.CellularComponent));
            Assert.That(term.EvidenceCodes, Is.EquivalentTo(new[] { "ECO:0000314" }));
        }

        [Test]
        public void GoTerms_ReferenceWithNoPropertiesStillYieldsTheTerm()
        {
            // The id is the annotation. A missing term/evidence/project must not drop the row,
            // because a dropped row is indistinguishable from a protein that was never annotated.
            var protein = WithReferences(new DatabaseReference("GO", "GO:0005737", null));

            var term = protein.GoTerms.Single();
            Assert.That(term.Id, Is.EqualTo("GO:0005737"));
            Assert.That(term.Name, Is.Empty);
            Assert.That(term.Aspect, Is.EqualTo(GoAspect.Unknown));
            Assert.That(term.EvidenceCodes, Is.Empty);
        }

        [Test]
        public void Decoys_CarryNoGoTerms()
        {
            // Extends the binding Decoys_DoNotInheritUnrelatedDatabaseReferences: only the taxonomy
            // travels onto a decoy, so enrichment backgrounds stay target-only.
            var proteins = ProteinDbLoader.LoadProteinXML(Data("DatabaseTests", "humanGAPDH.xml"),
                true, DecoyType.Reverse, null, false, null, out _);

            var target = proteins.First(p => !p.IsDecoy);
            var decoy = proteins.First(p => p.IsDecoy);

            Assert.That(target.GoTerms, Is.Not.Empty);
            Assert.That(decoy.GoTerms, Is.Empty, "a decoy sequence is not annotated in GO");
        }

        [Test]
        public void GoTerms_EmptyRatherThanNullWhenNothingIsAnnotated()
        {
            var protein = new Protein("PEPTIDEK", "P00001");
            Assert.That(protein.GoTerms, Is.Not.Null);
            Assert.That(protein.GoTerms, Is.Empty);
        }

        [Test]
        public void GeneOntologyDatabaseReferenceType_IsReExportedFromTheLoader()
        {
            // Same shape as NcbiTaxonomyDatabaseReferenceType: callers that already reference the
            // loader should not have to reach into Proteomics for the string.
            Assert.That(ProteinDbLoader.GeneOntologyDatabaseReferenceType,
                Is.EqualTo(Protein.GeneOntologyDatabaseReferenceType));
            Assert.That(Protein.GeneOntologyDatabaseReferenceType, Is.EqualTo("GO"));
        }
    }
}
