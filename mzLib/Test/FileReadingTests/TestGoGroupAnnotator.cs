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
    /// GoGroupAnnotator turns protein groups into one row per (group, GO term) without collapsing the
    /// group: every member's terms are kept, and each row says which members carry the term
    /// (accession_used) and how many (n_with of n_members). A group with no term still gets a row, and
    /// annotation_status says why.
    ///
    /// Ontology: the go-trimmed.obo fixture. GO:0005743 mitochondrial inner membrane lies below GO:0005739
    /// mitochondrion (via part_of) and below GO:0043226 organelle; GO:0005634 nucleus lies below
    /// GO:0043226 too. GO:0004365 is GAPDH activity (MF).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestGoGroupAnnotator
    {
        private static readonly string OntologyPath =
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Ontologies", "go-trimmed.obo");

        private const string InnerMembrane = "GO:0005743";
        private const string Mitochondrion = "GO:0005739";
        private const string Nucleus = "GO:0005634";
        private const string Organelle = "GO:0043226";
        private const string Gapdh = "GO:0004365";
        private const string DbSha = "0123abcd";

        private static GeneOntologyGraph _go;

        [OneTimeSetUp]
        public void LoadOntology() => _go = GeneOntologyGraph.Load(OntologyPath);

        private static DatabaseReference Go(string id, params string[] evidence) =>
            new("GO", id, new[] { Tuple.Create("term", "C:whatever UniProt wrote") }
                .Concat(evidence.Select(e => Tuple.Create("evidence", e))));

        private static Protein P(string accession, params DatabaseReference[] refs) =>
            new("PEPTIDEK", accession, databaseReferences: refs.ToList());

        private static GoAnnotationGroup Group(string members, bool contaminant = false, bool decoy = false, double q = 0.001) =>
            new(members, members.Split('|'), decoy, contaminant, q);

        private static GoGroupAnnotator Annotator(params Protein[] proteins) => new(_go, proteins, DbSha);

        private static GoAnnotationRow Row(IEnumerable<GoAnnotationRow> rows, string goId) => rows.Single(r => r.GoId == goId);

        [Test]
        public void Rows_OneRowPerTermCarriedByAnyMember()
        {
            // The worked example from go's thread 008: all three members are mitochondrial, one is nuclear.
            var annotator = Annotator(
                P("P1", Go(Mitochondrion, "ECO:0000314")),
                P("P2", Go(Mitochondrion, "ECO:0000314")),
                P("P3", Go(Mitochondrion, "ECO:0000314"), Go(Nucleus, "ECO:0000501")));

            var rows = annotator.Annotate(Group("P1|P2|P3"));

            Assert.That(Row(rows, Mitochondrion).NWith, Is.EqualTo(3));
            Assert.That(Row(rows, Nucleus).NWith, Is.EqualTo(1));
            Assert.That(rows.All(r => r.NMembers == 3), Is.True);
            Assert.That(rows.Select(r => r.GoId).Distinct().Count(), Is.EqualTo(rows.Count), "one row per term");
        }

        [Test]
        public void Rows_ReconstructViews_UnionAndConsensus_NoLeadingView()
        {
            var rows = Annotator(
                    P("P1", Go(Mitochondrion)),
                    P("P2", Go(Mitochondrion), Go(Nucleus)))
                .Annotate(Group("P1|P2"));

            // Union = every row. Consensus = n_with == n_members. There is no "leading" view to reconstruct:
            // MetaMorpheus orders members alphabetically and never picks one (go D25).
            var consensus = rows.Where(r => r.NWith == r.NMembers).Select(r => r.GoId).ToList();
            Assert.That(consensus, Does.Contain(Mitochondrion).And.Not.Contain(Nucleus));
            Assert.That(rows.Select(r => r.GoId), Does.Contain(Nucleus));
        }

        [Test]
        public void AccessionUsed_IsSetValued_ExactlyTheCarryingMembers()
        {
            var rows = Annotator(
                    P("P1", Go(InnerMembrane)),
                    P("P2", Go(Nucleus)),
                    P("P3"))
                .Annotate(Group("P1|P2|P3"));

            Assert.That(Row(rows, Mitochondrion).AccessionUsed, Is.EqualTo(new[] { "P1" }));
            Assert.That(Row(rows, Organelle).AccessionUsed, Is.EqualTo(new[] { "P1", "P2" }));
            foreach (var row in rows)
            {
                Assert.That(row.NWith, Is.EqualTo(row.AccessionUsed.Count), "n_with is |accession_used|, always");
            }
        }

        [Test]
        public void Propagation_AncestorRowsCarryUnionOfDescendantEvidence()
        {
            // Thread 008's example: direct inner membrane {314} and direct nucleus {501} produce a propagated
            // organelle row carrying both codes.
            var rows = Annotator(
                    P("P1", Go(InnerMembrane, "ECO:0000314")),
                    P("P2", Go(Nucleus, "ECO:0000501")))
                .Annotate(Group("P1|P2"));

            var organelle = Row(rows, Organelle);
            Assert.That(organelle.Propagated, Is.True);
            Assert.That(organelle.Evidence, Is.EqualTo(new[] { "ECO:0000314", "ECO:0000501" }));

            // So a downstream "drop rows backed only by ECO:0000501" filter still bites after propagation:
            // mitochondrion came only from the 314-backed inner membrane, nucleus only from 501.
            Assert.That(Row(rows, Mitochondrion).Evidence, Is.EqualTo(new[] { "ECO:0000314" }));
            Assert.That(Row(rows, Nucleus).Evidence, Is.EqualTo(new[] { "ECO:0000501" }));
            Assert.That(Row(rows, Nucleus).Propagated, Is.False);
        }

        [Test]
        public void Propagated_IsFalse_WhenAnyCarryingMemberHasTheTermDirectly()
        {
            var rows = Annotator(
                    P("P1", Go(Mitochondrion, "ECO:0000269")),
                    P("P2", Go(InnerMembrane, "ECO:0000314")))
                .Annotate(Group("P1|P2"));

            var mito = Row(rows, Mitochondrion);
            Assert.That(mito.Propagated, Is.False, "P1 carries it directly");
            Assert.That(mito.AccessionUsed, Is.EqualTo(new[] { "P1", "P2" }));
            Assert.That(mito.Evidence, Is.EqualTo(new[] { "ECO:0000269", "ECO:0000314" }));
            Assert.That(mito.Inherited, Is.False);
        }

        [Test]
        public void Rows_NameAndAspectComeFromTheOntology_NotFromUniProtsText()
        {
            var row = Row(Annotator(P("P1", Go(Gapdh))).Annotate(Group("P1")), Gapdh);

            Assert.That(row.GoName, Is.EqualTo("glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity"));
            Assert.That(row.Aspect, Is.EqualTo(GoAspect.MolecularFunction));
        }

        [Test]
        public void Rows_AltIdCitedByUniProt_LandsOnThePrimaryTerm()
        {
            // GO:0016021 was merged into GO:0016020. A member citing the old id and a member citing the new one
            // carry the same term.
            var rows = Annotator(P("P1", Go("GO:0016021", "ECO:0000255")), P("P2", Go("GO:0016020", "ECO:0000314")))
                .Annotate(Group("P1|P2"));

            Assert.That(rows.Any(r => r.GoId == "GO:0016021"), Is.False);
            Assert.That(Row(rows, "GO:0016020").AccessionUsed, Is.EqualTo(new[] { "P1", "P2" }));
        }

        [Test]
        public void Rows_ObsoleteTermCitedByUniProt_IsKept_WithoutPropagation()
        {
            var rows = Annotator(P("P1", Go("GO:0006082"))).Annotate(Group("P1"));

            Assert.That(rows.Select(r => r.GoId), Is.EqualTo(new[] { "GO:0006082" }));
            Assert.That(rows[0].Status, Is.EqualTo(GoAnnotationStatus.Annotated));
        }

        [Test]
        public void Annotator_TermAbsentFromTheRelease_Throws_NamingEveryMissingId()
        {
            // A UniProt XML newer than the pinned go.obo. Refuse rather than silently drop the term or treat it
            // as a root: the fix is a matching release, and the message has to say which ids.
            var ex = Assert.Throws<InvalidDataException>(() =>
                Annotator(P("P1", Go("GO:9999998")), P("P2", Go("GO:9999999"), Go(Nucleus))));

            Assert.That(ex.Message, Does.Contain("GO:9999998").And.Contain("GO:9999999").And.Contain("releases/2026-07-26"));
        }

        [Test]
        public void EveryNonDecoyGroup_GetsAtLeastOneRow_NoGoTerms()
        {
            var rows = Annotator(P("P1")).Annotate(Group("P1"));

            Assert.That(rows, Has.Count.EqualTo(1));
            Assert.That(rows[0].Status, Is.EqualTo(GoAnnotationStatus.NoGoTerms));
            Assert.That(rows[0].GoId, Is.Null);
            Assert.That(rows[0].GoName, Is.Null);
            Assert.That(rows[0].Aspect, Is.Null);
            Assert.That(rows[0].Propagated, Is.Null, "no term, so neither true nor false");
            Assert.That(rows[0].Inherited, Is.Null);
            Assert.That(rows[0].AccessionUsed, Is.Empty);
            Assert.That(rows[0].Evidence, Is.Empty);
            Assert.That(rows[0].NWith, Is.EqualTo(0));
        }

        [Test]
        public void MemberAbsentFromTheAnnotationDatabase_ReadsNoEntry()
        {
            var rows = Annotator(P("P1")).Annotate(Group("Q9"));

            Assert.That(rows.Single().Status, Is.EqualTo(GoAnnotationStatus.NoEntry));
        }

        [Test]
        public void AnnotationStatus_Precedence_NoEntryBeatsNoGoTerms()
        {
            // One member found without terms, one not found at all: the more specific absence wins (go D19).
            var rows = Annotator(P("P1")).Annotate(Group("P1|Q9"));

            Assert.That(rows.Single().Status, Is.EqualTo(GoAnnotationStatus.NoEntry));
        }

        [Test]
        public void AnnotationStatus_Precedence_ContaminantBeatsNoEntry()
        {
            var rows = Annotator().Annotate(Group("Q9", contaminant: true));

            Assert.That(rows.Single().Status, Is.EqualTo(GoAnnotationStatus.Contaminant));
        }

        [Test]
        public void AnnotatedContaminant_ReadsAnnotated_NotContaminant()
        {
            // The invariant consumers build against: a non-empty go_id always means annotated. So an annotated
            // keratin group reads annotated, and "contaminant" does not catch every contaminant group (go D19).
            var rows = Annotator(P("P1", Go(Mitochondrion))).Annotate(Group("P1", contaminant: true));

            Assert.That(rows.All(r => r.GoId != null && r.Status == GoAnnotationStatus.Annotated), Is.True);
        }

        [Test]
        public void ContaminantGroups_AreIncluded()
        {
            Assert.That(Annotator(P("P1", Go(Nucleus))).Annotate(Group("P1", contaminant: true)), Is.Not.Empty);
        }

        [Test]
        public void DecoyGroup_IsRejectedAtTheBoundary()
        {
            // Decoys carry no GO (go D3). A decoy group would read no_entry and pass for a real absence.
            Assert.Throws<ArgumentException>(() => Annotator(P("P1")).Annotate(Group("DECOY_P1", decoy: true)));
        }

        [Test]
        public void DecoyProteinsInTheAnnotationSet_AreIgnored()
        {
            var decoy = new Protein("PEPTIDEK", "P1", isDecoy: true, databaseReferences: new List<DatabaseReference> { Go(Nucleus) });

            Assert.That(new GoGroupAnnotator(_go, new[] { decoy }, DbSha).Annotate(Group("P1")).Single().Status,
                Is.EqualTo(GoAnnotationStatus.NoEntry));
        }

        [Test]
        public void SameAccessionInTwoDatabases_TermsAreUnioned()
        {
            // A keratin is in the target database and in the contaminant database. Neither copy is "the" entry.
            var annotator = Annotator(P("P1", Go(Nucleus, "ECO:0000314")), P("P1", Go(Nucleus, "ECO:0000501"), Go(Gapdh)));
            var rows = annotator.Annotate(Group("P1"));

            Assert.That(Row(rows, Nucleus).Evidence, Is.EqualTo(new[] { "ECO:0000314", "ECO:0000501" }));
            Assert.That(rows.Select(r => r.GoId), Does.Contain(Gapdh));
        }

        [Test]
        public void EveryRow_CarriesTheGroupAndProvenance()
        {
            var rows = Annotator(P("P1", Go(InnerMembrane)), P("P2")).Annotate(Group("P1|P2", q: 0.0042));

            foreach (var row in rows)
            {
                Assert.That(row.ProteinGroup, Is.EqualTo("P1|P2"));
                Assert.That(row.QValue, Is.EqualTo(0.0042));
                Assert.That(row.GoRelease, Is.EqualTo("releases/2026-07-26"));
                Assert.That(row.GoOboSha256, Is.EqualTo(_go.SourceSha256));
                Assert.That(row.AnnotationDbSha256, Is.EqualTo(DbSha));
            }
        }

        [Test]
        public void Rows_AreOrderedByGoId()
        {
            var ids = Annotator(P("P1", Go(InnerMembrane), Go(Gapdh))).Annotate(Group("P1")).Select(r => r.GoId).ToList();

            Assert.That(ids, Is.EqualTo(ids.OrderBy(i => i, StringComparer.Ordinal).ToList()));
            Assert.That(ids, Has.Count.EqualTo(1 + 15 + 1 + 6), "inner membrane + its 15 ancestors + GAPDH + its 6");
        }

        [Test]
        public void AnnotateAll_KeepsInputOrder_AndEveryGroup()
        {
            var rows = Annotator(P("P1", Go(Nucleus))).AnnotateAll(new[] { Group("Q9"), Group("P1") }).ToList();

            Assert.That(rows.First().ProteinGroup, Is.EqualTo("Q9"));
            Assert.That(rows.Select(r => r.ProteinGroup).Distinct(), Is.EqualTo(new[] { "Q9", "P1" }));
        }

        [Test]
        public void Isoform_InheritsItsEntrysTerms_FlaggedInherited()
        {
            // A FASTA search with isoforms reports P04406-2; UniProt XML has only the entry P04406. The isoform
            // takes the entry's terms, marked inherited, because an isoform can differ precisely in where it is.
            var rows = Annotator(P("P04406", Go(Nucleus, "ECO:0000314"))).Annotate(Group("P04406-2"));

            var nucleus = Row(rows, Nucleus);
            Assert.That(nucleus.Inherited, Is.True);
            Assert.That(nucleus.AccessionUsed, Is.EqualTo(new[] { "P04406-2" }), "the member as the search named it");
            Assert.That(nucleus.Evidence, Is.EqualTo(new[] { "ECO:0000314" }));
            Assert.That(Row(rows, Organelle).Inherited, Is.True, "propagated rows inherit the flag");
        }

        [Test]
        public void Inherited_IsFalse_WhenAnyCarryingMemberHasItsOwnEntry()
        {
            var rows = Annotator(P("P04406", Go(Nucleus))).Annotate(Group("P04406|P04406-2"));

            var nucleus = Row(rows, Nucleus);
            Assert.That(nucleus.Inherited, Is.False);
            Assert.That(nucleus.AccessionUsed, Is.EqualTo(new[] { "P04406", "P04406-2" }));
        }

        [Test]
        public void Isoform_WhoseEntryIsAbsent_ReadsNoEntry()
        {
            Assert.That(Annotator(P("P04406", Go(Nucleus))).Annotate(Group("Q13148-3")).Single().Status,
                Is.EqualTo(GoAnnotationStatus.NoEntry));
        }

        [Test]
        public void Isoform_WhoseEntryHasNoTerms_ReadsNoGoTerms()
        {
            // The entry was found, so the absence is of terms, not of an entry.
            Assert.That(Annotator(P("P04406")).Annotate(Group("P04406-2")).Single().Status,
                Is.EqualTo(GoAnnotationStatus.NoGoTerms));
        }

        [Test]
        public void Isoform_WithItsOwnEntry_UsesItAndIsNotInherited()
        {
            var rows = Annotator(P("P04406", Go(Nucleus)), P("P04406-2", Go(Gapdh))).Annotate(Group("P04406-2"));

            Assert.That(rows.Any(r => r.GoId == Nucleus), Is.False, "its own entry wins; the parent's is not merged in");
            Assert.That(Row(rows, Gapdh).Inherited, Is.False);
        }

        [Test]
        public void AccessionOutsideUniProtsGrammar_NeverInherits()
        {
            // ProteinAccession parses and never repairs: a hyphenated name that is not a UniProt isoform
            // is not stripped to something that happens to exist.
            var rows = Annotator(P("contam", Go(Nucleus))).Annotate(Group("contam-2"));

            Assert.That(rows.Single().Status, Is.EqualTo(GoAnnotationStatus.NoEntry));
        }

        [Test]
        public void NullArguments_ThrowArgumentNull_NotNullReference()
        {
            Assert.Throws<ArgumentNullException>(() => new GoGroupAnnotator(null, Array.Empty<Protein>(), DbSha));
            Assert.Throws<ArgumentNullException>(() => new GoGroupAnnotator(_go, null, DbSha));
            Assert.Throws<ArgumentNullException>(() => Annotator().Annotate(null));
        }

        [Test]
        public void Group_WithNoMembers_Throws()
        {
            Assert.Throws<ArgumentException>(() => Annotator().Annotate(new GoAnnotationGroup("", Array.Empty<string>(), false, false, 0)));
        }
    }
}
