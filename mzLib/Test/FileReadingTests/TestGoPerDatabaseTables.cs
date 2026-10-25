using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases.Annotation;
using UsefulProteomicsDatabases.GeneOntology;

namespace Test.FileReadingTests
{
    /// <summary>
    /// GO annotation stored once per database, not once per search. Four tables replace a search's group x term
    /// file:
    /// <list type="bullet">
    /// <item>(a) a database's direct annotations (<see cref="GoAccessionTermsTsv"/>), one per database and go.obo;</item>
    /// <item>the release's terms and closure (<see cref="GoOntologyTsv"/>), one per go.obo;</item>
    /// <item>(b) a map's categories over the whole release (<see cref="GoCategoryTsv"/>), one per go.obo and map;</item>
    /// <item>(c) a search's member resolution (<see cref="MemberResolutionTsv"/>), one per search, holding no GO.</item>
    /// </list>
    /// A catalog view joins them. The test that matters most is <see cref="Join_OfTheStoredTables_EqualsTheAnnotator"/>:
    /// the join, computed from the tables as read back from disk, gives exactly GoGroupAnnotator's rows.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestGoPerDatabaseTables
    {
        private static readonly string OntologyPath =
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Ontologies", "go-trimmed.obo");

        private const string InnerMembrane = "GO:0005743";
        private const string Mitochondrion = "GO:0005739";
        private const string Nucleus = "GO:0005634";
        private const string Organelle = "GO:0043226";
        private const string Gapdh = "GO:0004365";
        private const string Membrane = "GO:0016020";
        private const string MembraneAltId = "GO:0016021";
        private const string ObsoleteTerm = "GO:0006082";
        private const string DbSha = "0123abcd";

        private static GeneOntologyGraph _go;
        private string _dir;

        [OneTimeSetUp]
        public void LoadOntology() => _go = GeneOntologyGraph.Load(OntologyPath);

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "GoPerDatabase_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private static DatabaseReference Go(string id, params string[] evidence) =>
            new("GO", id, new[] { Tuple.Create("term", "C:whatever UniProt wrote") }
                .Concat(evidence.Select(e => Tuple.Create("evidence", e))));

        private static Protein P(string accession, params DatabaseReference[] refs) =>
            new("PEPTIDEK", accession, databaseReferences: refs.ToList());

        private static GoAnnotationGroup Group(string members, bool contaminant = false, double q = 0.001, string name = null) =>
            new(name ?? members, members.Split('|'), false, contaminant, q);

        // ---- (c) the member resolver: presence in the database only, no GO ----

        [Test]
        public void Resolver_MemberInTheDatabase_IsItself()
        {
            var r = new MemberResolver(new[] { P("P04406") }).Resolve("P04406");

            Assert.That(r, Is.EqualTo(new MemberResolution("P04406", "P04406", MemberResolutionKind.Self, false, "P04406")));
        }

        // The sequence the member is: what an annotation naming one isoform (a Reactome molecule) joins on.
        [TestCase("P04406", "P04406", "P04406", "P04406")]
        [TestCase("P04406_A20T", "P04406", "P04406", "P04406")]
        [TestCase("P04406_A20T", "P04406|P04406_A20T", "P04406_A20T", "P04406")]
        [TestCase("P04406-2", "P04406", "P04406", "P04406-2")]
        [TestCase("P04406-2_A20T", "P04406", "P04406", "P04406-2")]
        [TestCase("P04406-2_A20T", "P04406|P04406-2", "P04406-2", "P04406-2")]
        [TestCase("P04406-2_A20T_G31", "P04406", "P04406", "P04406-2")]
        [TestCase("CON_P02768", "CON_P02768", "CON_P02768", "CON_P02768")]
        public void Resolver_SequenceAccession_IsTheMemberWithoutItsVariantSuffix(string member, string database,
            string annotatedAs, string sequence)
        {
            var r = new MemberResolver(database.Split('|').Select(a => P(a))).Resolve(member);

            Assert.That(r.AnnotatedAs, Is.EqualTo(annotatedAs));
            Assert.That(r.SequenceAccession, Is.EqualTo(sequence));
            Assert.That(MemberResolver.SequenceAccessionOf(member), Is.EqualTo(sequence));
        }

        [Test]
        public void Resolver_SequenceAccession_NullWhenNoEntry_AndTheStaticRefusesEmpty()
        {
            Assert.That(new MemberResolver(new[] { P("P04406") }).Resolve("Random_P04406_f1").SequenceAccession, Is.Null);
            Assert.That(MemberResolver.SequenceAccessionOf("Random_P04406_f1"), Is.EqualTo("Random_P04406_f1"));
            Assert.Throws<ArgumentException>(() => MemberResolver.SequenceAccessionOf(""));
        }

        [TestCase("P04406-2", "P04406")]
        [TestCase("P04406_A20T", "P04406")]
        [TestCase("P04406_A20T_G31", "P04406")]
        [TestCase("P04406-2_A20T", "P04406")]
        public void Resolver_IsoformOrVariantAbsentFromTheDatabase_TakesItsEntry(string member, string entry)
        {
            var r = new MemberResolver(new[] { P("P04406") }).Resolve(member);

            Assert.That(r.Kind, Is.EqualTo(MemberResolutionKind.Inherited));
            Assert.That(r.AnnotatedAs, Is.EqualTo(entry));
            Assert.That(r.Accession, Is.EqualTo(member), "the member as the search named it");
        }

        [Test]
        public void Resolver_VariantOfAnIsoformInTheDatabase_TakesTheIsoform()
        {
            var r = new MemberResolver(new[] { P("P04406"), P("P04406-2") }).Resolve("P04406-2_A20T");

            Assert.That(r.AnnotatedAs, Is.EqualTo("P04406-2"));
            Assert.That(r.Kind, Is.EqualTo(MemberResolutionKind.Inherited));
        }

        [Test]
        public void Resolver_DependsOnPresenceOnly_NotOnWhatAnEntryIsAnnotatedWith()
        {
            // An isoform in the database is itself even with no GO: whether it then takes its entry's terms is the
            // database annotation's rule (D37), so the same resolution serves GO, pathways or anything else.
            var r = new MemberResolver(new[] { P("P04406", Go(Nucleus)), P("P04406-2") }).Resolve("P04406-2");

            Assert.That(r.Kind, Is.EqualTo(MemberResolutionKind.Self));
        }

        [TestCase("Q13148-3")]
        [TestCase("Random_P04406_f1")]
        [TestCase("Random_foreign_P04406")]
        [TestCase("NP_000001_A20T")]
        [TestCase("contam-2")]
        public void Resolver_NoEntryByTheGrammar_ReadsNoEntry(string member)
        {
            var r = new MemberResolver(new[] { P("P04406"), P("contam"), P("NP_000001") }).Resolve(member);

            Assert.That(r.Kind, Is.EqualTo(MemberResolutionKind.NoEntry));
            Assert.That(r.AnnotatedAs, Is.Null);
        }

        [Test]
        public void Resolver_Entrapment_ByAccessionOrByTheDatabasesMark()
        {
            var flagged = new Protein("PEPTIDEK", "X1", isEntrapment: true);
            var resolver = new MemberResolver(new[] { flagged, P("P1") });

            Assert.That(resolver.Resolve("X1").IsEntrapment, Is.True);
            Assert.That(resolver.Resolve("Random_P1_f0").IsEntrapment, Is.True);
            Assert.That(resolver.Resolve("P1").IsEntrapment, Is.False);
        }

        [Test]
        public void Resolver_DecoysAreNotEntries()
        {
            var decoy = new Protein("PEPTIDEK", "P1", isDecoy: true);

            Assert.That(new MemberResolver(new[] { decoy }).Resolve("P1").Kind, Is.EqualTo(MemberResolutionKind.NoEntry));
        }

        [Test]
        public void Resolver_ResolveAll_OneRowPerDistinctMember_InOrdinalOrder()
        {
            var rows = new MemberResolver(new[] { P("P2") }).ResolveAll(new[] { "P2", "P1", "P2", "", null });

            Assert.That(rows.Select(r => r.Accession), Is.EqualTo(new[] { "P1", "P2" }));
        }

        [Test]
        public void Resolver_BadArguments_Throw()
        {
            Assert.Throws<ArgumentNullException>(() => new MemberResolver(null));
            Assert.Throws<ArgumentException>(() => new MemberResolver(Array.Empty<Protein>()).Resolve(""));
            Assert.Throws<ArgumentNullException>(() => new MemberResolver(Array.Empty<Protein>()).ResolveAll(null));
        }

        // ---- (a) a database's direct annotations ----

        [Test]
        public void Accession_DirectTermsOnly_EveryEntryHasARow()
        {
            var rows = new GoAccessionAnnotation(_go, new[] { P("P1", Go(InnerMembrane, "ECO:0000314")), P("P2") }, DbSha).Rows();

            Assert.That(rows.Where(r => r.Accession == "P1").Select(r => r.GoId), Is.EqualTo(new[] { InnerMembrane }),
                "no ancestors: they are the closure's");
            var p2 = rows.Single(r => r.Accession == "P2");
            Assert.That(p2.GoId, Is.Null, "in the database with no GO is a row, not an absence");
            Assert.That(p2.Inherited, Is.Null);
        }

        [Test]
        public void Accession_IsoformInTheDatabaseWithNoGo_TakesItsEntrysTerms_Inherited()
        {
            var rows = new GoAccessionAnnotation(_go, new[] { P("P04406", Go(Nucleus, "ECO:0000314")), P("P04406-2") }, DbSha).Rows();

            var isoform = rows.Single(r => r.Accession == "P04406-2");
            Assert.That(isoform.GoId, Is.EqualTo(Nucleus));
            Assert.That(isoform.Inherited, Is.True);
            Assert.That(isoform.Evidence, Is.EqualTo(new[] { "ECO:0000314" }));
            Assert.That(rows.Single(r => r.Accession == "P04406").Inherited, Is.False);
        }

        [Test]
        public void Accession_IsoformWithItsOwnGo_KeepsIt()
        {
            var rows = new GoAccessionAnnotation(_go, new[] { P("P04406", Go(Nucleus)), P("P04406-2", Go(Gapdh)) }, DbSha).Rows();

            Assert.That(rows.Where(r => r.Accession == "P04406-2").Select(r => (r.GoId, r.Inherited)),
                Is.EqualTo(new[] { (Gapdh, (bool?)false) }));
        }

        [Test]
        public void Accession_AltIdLandsOnThePrimary_ObsoleteIsKept_DuplicatesAreUnioned()
        {
            var rows = new GoAccessionAnnotation(_go, new[]
            {
                P("P1", Go(MembraneAltId, "ECO:0000255")),
                P("P2", Go(ObsoleteTerm)),
                P("P3", Go(Nucleus, "ECO:0000314")),
                P("P3", Go(Nucleus, "ECO:0000501")),
            }, DbSha).Rows();

            Assert.That(rows.Single(r => r.Accession == "P1").GoId, Is.EqualTo(Membrane));
            Assert.That(rows.Single(r => r.Accession == "P2").GoId, Is.EqualTo(ObsoleteTerm));
            Assert.That(rows.Single(r => r.Accession == "P3").Evidence, Is.EqualTo(new[] { "ECO:0000314", "ECO:0000501" }));
        }

        [Test]
        public void Accession_UnknownId_RefusedUnlessSkipped()
        {
            var proteins = new[] { P("P1", Go("GO:9999999"), Go(Nucleus)) };

            Assert.Throws<InvalidDataException>(() => new GoAccessionAnnotation(_go, proteins, DbSha));
            var skipped = new GoAccessionAnnotation(_go, proteins, DbSha, skipUnknownGoIds: true);
            Assert.That(skipped.UnresolvedGoIds, Is.EqualTo(new[] { "GO:9999999" }));
            Assert.That(skipped.Rows().Select(r => r.GoId), Is.EqualTo(new[] { Nucleus }));
        }

        // ---- round trips ----

        private string WriteFile(string name, Action<TextWriter> write)
        {
            string path = Path.Combine(_dir, name);
            using (var writer = new StreamWriter(path)) write(writer);
            return path;
        }

        [Test]
        public void AccessionTerms_RoundTrip()
        {
            var annotation = new GoAccessionAnnotation(_go,
                new[] { P("P04406", Go(Nucleus, "ECO:0000314", "ECO:0000501")), P("P04406-2"), P("P3") }, DbSha);

            var table = GoAccessionTermsTsv.Read(WriteFile("a.tsv", w => GoAccessionTermsTsv.Write(w, annotation)));

            Assert.That(table.Header["go_obo_sha256"], Is.EqualTo(_go.SourceSha256));
            Assert.That(table.Header["go_release"], Is.EqualTo("releases/2026-07-26"));
            Assert.That(table.Header["annotation_db_sha256"], Is.EqualTo(DbSha));
            Assert.That(table.Header.ContainsKey("unresolved_go_ids"), Is.False);
            Assert.That(table.Rows.Select(Key), Is.EqualTo(annotation.Rows().Select(Key)));
        }

        private static string Key(GoAccessionTermRow r) =>
            $"{r.Accession}|{r.GoId}|{r.Aspect}|{string.Join(",", r.Evidence)}|{r.Inherited}";

        [Test]
        public void AccessionTerms_HeaderCountsSkippedIds()
        {
            var annotation = new GoAccessionAnnotation(_go, new[] { P("P1", Go("GO:9999999")) }, DbSha, skipUnknownGoIds: true);

            var table = GoAccessionTermsTsv.Read(WriteFile("a.tsv", w => GoAccessionTermsTsv.Write(w, annotation)));

            Assert.That(table.Header["unresolved_go_ids"], Is.EqualTo("1"));
        }

        [Test]
        public void MemberResolution_RoundTrip()
        {
            var rows = new MemberResolver(new[] { P("P04406") })
                .ResolveAll(new[] { "P04406", "P04406-2", "Random_P04406_f0", "Q9ZZZ9" });

            var table = MemberResolutionTsv.Read(WriteFile("c.tsv", w => MemberResolutionTsv.Write(w, rows, DbSha, "feed")));

            Assert.That(table.Rows, Is.EqualTo(rows));
            Assert.That(table.Header["annotation_db_sha256"], Is.EqualTo(DbSha));
            Assert.That(table.Header["source_file_sha256"], Is.EqualTo("feed"));
        }

        [Test]
        public void MemberResolution_Writer_RefusesAContradictionOrARepeat()
        {
            using var writer = new StringWriter();
            Assert.Throws<ArgumentException>(() => MemberResolutionTsv.Write(writer,
                new[] { new MemberResolution("P1", "P2", MemberResolutionKind.Self, false, "P1") }, DbSha));
            Assert.Throws<ArgumentException>(() => MemberResolutionTsv.Write(writer,
                new[] { new MemberResolution("P1", "P1", MemberResolutionKind.Self, false, "P1"),
                        new MemberResolution("P1", "P1", MemberResolutionKind.Self, false, "P1") }, DbSha));
            Assert.Throws<ArgumentException>(() => MemberResolutionTsv.Write(writer,
                new[] { new MemberResolution("P1-2", "P1", MemberResolutionKind.Inherited, false, "P1") }, DbSha));
            Assert.Throws<ArgumentException>(() => MemberResolutionTsv.Write(writer, Array.Empty<MemberResolution>(), ""));
        }

        [Test]
        public void Terms_RoundTrip_ObsoleteIncluded()
        {
            var table = GoOntologyTsv.ReadTerms(WriteFile("terms.tsv", w => GoOntologyTsv.WriteTerms(w, _go)));

            Assert.That(table.Rows, Has.Count.EqualTo(_go.Count));
            Assert.That(table.Rows.Single(r => r.GoId == ObsoleteTerm).IsObsolete, Is.True);
            Assert.That(table.Rows.Single(r => r.GoId == Gapdh).GoName,
                Is.EqualTo("glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity"));
            Assert.That(table.Rows.Any(r => r.GoId == MembraneAltId), Is.False, "alternative ids are not terms");
        }

        [Test]
        public void Closure_EveryTermWithItselfAndItsAncestors_ObsoleteOnlyItself()
        {
            var table = GoOntologyTsv.ReadClosure(WriteFile("closure.tsv", w => GoOntologyTsv.WriteClosure(w, _go)));
            var up = table.Rows.ToLookup(r => r.GoId, r => r.AncestorGoId);

            foreach (string id in _go.TermIds)
            {
                Assert.That(up[id], Is.EquivalentTo(_go.Ancestors(id).Append(id)), id);
            }
            Assert.That(up[ObsoleteTerm], Is.EqualTo(new[] { ObsoleteTerm }));
            Assert.That(up[InnerMembrane], Does.Contain(Mitochondrion).And.Contain(Organelle));
        }

        [Test]
        public void Closure_Reader_RefusesATermWithoutItsSelfPair()
        {
            string path = Path.Combine(_dir, "closure.tsv");
            File.WriteAllText(path, $"#!go_closure_format 1\n#!go_obo_sha256 x\ngo_id\tancestor_go_id\n{InnerMembrane}\t{Mitochondrion}\n");

            var ex = Assert.Throws<InvalidDataException>(() => GoOntologyTsv.ReadClosure(path));
            Assert.That(ex.Message, Does.Contain(InnerMembrane));
        }

        [Test]
        public void Categories_OverTheWholeRelease_AgreeWithTheTableForAnySearchsTerms()
        {
            string mapPath = Path.Combine(_dir, "map.tsv");
            File.WriteAllText(mapPath, "#!category_map_format 1\n#!map_name organelle\n#!map_version 1\n" +
                                       "category\tsubcategory\tanchor_go_id\n" +
                                       $"mitochondrion\t\t{Mitochondrion}\nmitochondrion\tinner_membrane\t{InnerMembrane}\nnucleus\t\t{Nucleus}\n");
            var resolver = new GoCategoryResolver(GoCategoryMap.Load(mapPath), _go);
            var groupRows = new GoGroupAnnotator(_go, new[] { P("P1", Go(InnerMembrane)), P("P2", Go(Nucleus)) }, DbSha)
                .AnnotateAll(new[] { Group("P1|P2") });

            var whole = GoCategoryTsv.Read(WriteFile("b.tsv", w => GoCategoryTsv.Write(w, resolver)));
            var forSearch = GoCategoryTsv.Read(WriteFile("b_search.tsv", w => GoCategoryTsv.Write(w, resolver, groupRows)));

            var searchTerms = groupRows.Select(r => r.GoId).ToHashSet();
            Assert.That(whole.Rows.Where(r => searchTerms.Contains(r.GoId)), Is.EqualTo(forSearch.Rows));
            Assert.That(whole.Rows.Single(r => r.GoId == InnerMembrane).Subcategory, Is.EqualTo("mitochondrion:inner_membrane"));
            Assert.That(whole.Header["category_map"], Does.StartWith("organelle 1 "));
        }

        // ---- the readers refuse what they cannot read correctly ----

        [TestCase("#!go_accession_terms_format 2\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\n", "must be 1")]
        [TestCase("#!go_obo_sha256 x\n#!go_accession_terms_format 1\n", "first header line")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\n#!bogus 1\n", "unknown header key")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\naccession\tgo_id\taspect\tevidence\tinherited\n", "annotation_db_sha256")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\n", "column header")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\tGO:0005634\n", "cells")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\t\tcellular_component\t\t\n", "term-less")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\t\t\t\t\nP1\tGO:0005634\tcellular_component\t\tfalse\n", "term-less row and another")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\tGO:0005634\tplace\t\tfalse\n", "not a GO namespace")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\tGO:0005634\tcellular_component\t\tmaybe\n", "true or false")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!go_obo_sha256 x\n", "appears twice")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\n", "no column header")]
        [TestCase("accession\tgo_id\taspect\tevidence\tinherited\n", "(missing)")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\n\tGO:0005634\tcellular_component\t\tfalse\n", "accession is empty")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\tGO:0005634\tcellular_component\t\tfalse\nP1\t\t\t\t\n", "term-less row and another")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\tGO:0005634\t\t\tfalse\n", "lacks its aspect")]
        [TestCase("#!go_accession_terms_format 1\n#!go_obo_sha256 x\n#!annotation_db_sha256 y\naccession\tgo_id\taspect\tevidence\tinherited\nP1\tGO:0005634\tcellular_component\t\tfalse\n\nP1\tGO:0005634\tcellular_component\t\ttrue\n", "appears twice")]
        public void AccessionTerms_Reader_Refuses(string text, string expected)
        {
            string path = Path.Combine(_dir, "a.tsv");
            File.WriteAllText(path, text);

            var ex = Assert.Throws<InvalidDataException>(() => GoAccessionTermsTsv.Read(path));
            Assert.That(ex.Message, Does.Contain(expected));
        }

        [TestCase("P1\tP2\tself\tfalse\tP1\n", "contradicts resolution")]
        [TestCase("P1\t\tinherited\tfalse\tP1\n", "contradicts resolution")]
        [TestCase("P1\tP1\tno_entry\tfalse\t\n", "contradicts resolution")]
        [TestCase("P1\tP1\tsomewhere\tfalse\tP1\n", "is not self")]
        [TestCase("P1\tP1\tself\t\tP1\n", "entrapment is empty")]
        [TestCase("P1\tP1\tself\tfalse\tP1\nP1\tP1\tself\tfalse\tP1\n", "appears twice")]
        [TestCase("P1-2\tP1\tinherited\tfalse\tP1\n", "sequence_accession 'P1' contradicts")]
        [TestCase("P1_A20T\tP1\tinherited\tfalse\t\n", "sequence_accession '' contradicts")]
        [TestCase("X1\t\tno_entry\tfalse\tX1\n", "sequence_accession 'X1' contradicts")]
        [TestCase("P1\tP1\tself\tfalse\n", "cells")]
        [TestCase("\t\tno_entry\tfalse\t\n", "accession is empty")]
        public void MemberResolution_Reader_Refuses(string rows, string expected)
        {
            string path = Path.Combine(_dir, "c.tsv");
            File.WriteAllText(path, "#!member_resolution_format 1\n#!annotation_db_sha256 y\naccession\tannotated_as\tresolution\tentrapment\tsequence_accession\n" + rows);

            var ex = Assert.Throws<InvalidDataException>(() => MemberResolutionTsv.Read(path));
            Assert.That(ex.Message, Does.Contain(expected));
        }

        [TestCase("go_id\tgo_name\taspect\tobsolete\n\tx\tcellular_component\tfalse\n", "empty or repeats")]
        [TestCase("go_id\tgo_name\taspect\tobsolete\nGO:1\tx\tcellular_component\tfalse\nGO:1\tx\tcellular_component\tfalse\n", "empty or repeats")]
        public void Terms_Reader_Refuses(string body, string expected)
        {
            string path = Path.Combine(_dir, "terms.tsv");
            File.WriteAllText(path, "#!go_terms_format 1\n#!go_obo_sha256 x\n" + body);

            var ex = Assert.Throws<InvalidDataException>(() => GoOntologyTsv.ReadTerms(path));
            Assert.That(ex.Message, Does.Contain(expected));
        }

        [TestCase("GO:1\t\n", "closure cell is empty")]
        [TestCase("GO:1\tGO:1\nGO:1\tGO:1\n", "appears twice")]
        public void Closure_Reader_Refuses(string rows, string expected)
        {
            string path = Path.Combine(_dir, "closure.tsv");
            File.WriteAllText(path, "#!go_closure_format 1\n#!go_obo_sha256 x\ngo_id\tancestor_go_id\n" + rows);

            var ex = Assert.Throws<InvalidDataException>(() => GoOntologyTsv.ReadClosure(path));
            Assert.That(ex.Message, Does.Contain(expected));
        }

        [TestCase("\tnucleus\t\n", "go_id or category is empty")]
        [TestCase("GO:1\tnucleus\tinner_membrane\n", "is not qualified")]
        [TestCase("GO:1\tnucleus\t\nGO:1\tnucleus\t\n", "repeats an earlier row")]
        public void Categories_Reader_Refuses(string rows, string expected)
        {
            string path = Path.Combine(_dir, "b.tsv");
            File.WriteAllText(path, "#!go_category_format 1\n#!go_obo_sha256 x\n#!category_map organelle 1 y\ngo_id\tcategory\tsubcategory\n" + rows);

            var ex = Assert.Throws<InvalidDataException>(() => GoCategoryTsv.Read(path));
            Assert.That(ex.Message, Does.Contain(expected));
        }

        [Test]
        public void AccessionTerms_Writer_RefusesAnEvidenceCodeHoldingTheSetSeparator()
        {
            var annotation = new GoAccessionAnnotation(_go, new[] { P("P1", Go(Nucleus, "ECO:1;2")) }, DbSha);

            using var writer = new StringWriter();
            Assert.Throws<ArgumentException>(() => GoAccessionTermsTsv.Write(writer, annotation));
        }

        [Test]
        public void Resolver_Contains_IgnoresNullProteinsAndNullAccessions()
        {
            var resolver = new MemberResolver(new[] { P("P1"), null });

            Assert.That(resolver.Contains("P1"), Is.True);
            Assert.That(resolver.Contains("P2"), Is.False);
            Assert.That(resolver.Contains(null), Is.False);
            Assert.That(resolver.IsEntrapment(null), Is.False);
        }

        [Test]
        public void Accession_TryGetTerms_UnknownOrNullAccession_IsFalse()
        {
            var annotation = new GoAccessionAnnotation(_go, new[] { P("P1", Go(Nucleus)) }, DbSha);

            Assert.That(annotation.TryGetTerms("P2", out _, out _), Is.False);
            Assert.That(annotation.TryGetTerms(null, out _, out _), Is.False);
        }

        [Test]
        public void Readers_MissingFile_Throws()
        {
            string path = Path.Combine(_dir, "absent.tsv");
            Assert.Throws<FileNotFoundException>(() => GoAccessionTermsTsv.Read(path));
            Assert.Throws<FileNotFoundException>(() => MemberResolutionTsv.Read(path));
            Assert.Throws<FileNotFoundException>(() => GoOntologyTsv.ReadTerms(path));
            Assert.Throws<FileNotFoundException>(() => GoOntologyTsv.ReadClosure(path));
            Assert.Throws<FileNotFoundException>(() => GoCategoryTsv.Read(path));
        }

        // ---- the join a catalog view makes equals the annotator ----

        /// <summary>
        /// A database and groups that exercise every rule of go:DEF-GROUP-GO-ANNOTATION v1: a member carrying only
        /// some of its group's terms; an isoform absent from the database, and one present with no GO (D37); a
        /// variant, a variant of an isoform, and a variant of an isoform with no GO; entrapment by accession and by
        /// the database's mark; contaminants, annotated and not; an obsolete id, an unknown id (skipped) and an
        /// alternative id; a duplicate entry; and each annotation status.
        /// </summary>
        internal static (Protein[] Database, GoAnnotationGroup[] Groups) Fixture()
        {
            var database = new[]
            {
                P("P11111", Go(Mitochondrion, "ECO:0000314"), Go(Nucleus, "ECO:0000501")),
                P("P22222", Go(InnerMembrane, "ECO:0000314")),
                P("P33333"),
                P("P04406", Go(Nucleus, "ECO:0000314"), Go(Gapdh, "ECO:0000269")),
                P("P04406-2"),
                P("Q13148", Go(InnerMembrane, "ECO:0000250")),
                P("Q13148-3", Go(Gapdh, "ECO:0000314")),
                P("P55555", Go(MembraneAltId, "ECO:0000255")),
                P("P66666", Go(ObsoleteTerm, "ECO:0000304")),
                P("P77777", Go("GO:9999999", "ECO:0000314"), Go(Nucleus, "ECO:0000501")),
                P("P88888", Go(Nucleus, "ECO:0000314")),
                P("P88888", Go(Nucleus, "ECO:0000501"), Go(Gapdh, "ECO:0000501")),
                P("P35908", Go(Nucleus, "ECO:0000501")),
                new Protein("PEPTIDEK", "X11111", isEntrapment: true, databaseReferences: new List<DatabaseReference> { Go(Mitochondrion) }),
                P("Random_foreign_Q99999", Go(Nucleus, "ECO:0000314")),
                new Protein("PEPTIDEK", "P99999", isDecoy: true, databaseReferences: new List<DatabaseReference> { Go(Nucleus) }),
            };
            var groups = new[]
            {
                Group("P11111|P22222|P33333", q: 0.004),
                Group("P04406|P04406-2"),
                Group("P04406-2"),
                Group("P04406-3"),
                Group("P04406_A20T"),
                Group("Q13148-3_A20T"),
                Group("Q13148-2_A20T"),
                Group("P04406-2_A20T"),
                Group("P55555|P66666"),
                Group("P77777"),
                Group("P88888", q: 0.02),
                Group("P33333"),
                Group("Q9ZZZ9"),
                Group("P33333|Q9ZZZ9"),
                Group("Q9ZZZ9", contaminant: true, name: "Q9ZZZ9 (contaminant)"),
                Group("P35908", contaminant: true),
                Group("X11111|Random_P11111_f0|Random_foreign_Q99999"),
                Group("NP_000001|P99999"),
            };
            return (database, groups);
        }

        /// <summary>
        /// The catalog view, written from the stored tables alone, with only equality joins and counts: each group
        /// member to its resolution (c), its annotated_as to the database's direct rows (a), each direct term to its
        /// closure pairs, and names and aspects from the term table. A member's evidence for a term is the union of
        /// the evidence of its direct terms at or below it; it is direct when it is annotated to the term itself;
        /// inherited when its resolution or its entry's row says so.
        /// </summary>
        internal static List<GoAnnotationRow> JoinView(IEnumerable<GoAnnotationGroup> groups,
            AnnotationTable<MemberResolution> c, AnnotationTable<GoAccessionTermRow> a,
            AnnotationTable<GoClosureRow> closure, AnnotationTable<GoTermRow> terms)
        {
            var resolution = c.Rows.ToDictionary(r => r.Accession, StringComparer.Ordinal);
            var direct = a.Rows.Where(r => r.GoId != null).ToLookup(r => r.Accession, StringComparer.Ordinal);
            var up = closure.Rows.ToLookup(r => r.GoId, r => r.AncestorGoId, StringComparer.Ordinal);
            var term = terms.Rows.ToDictionary(r => r.GoId, StringComparer.Ordinal);
            a.Header.TryGetValue("go_release", out var release);
            string obo = a.Header["go_obo_sha256"];
            string db = a.Header["annotation_db_sha256"];

            var output = new List<GoAnnotationRow>();
            foreach (var group in groups)
            {
                var members = group.MemberAccessions.Where(m => !string.IsNullOrEmpty(m)).Distinct(StringComparer.Ordinal).ToList();
                var entrapment = members.Where(m => resolution[m].IsEntrapment).OrderBy(m => m, StringComparer.Ordinal).ToList();
                var used = new SortedDictionary<string, SortedSet<string>>(StringComparer.Ordinal);
                var directMembers = new Dictionary<string, SortedSet<string>>(StringComparer.Ordinal);
                var inheritedMembers = new Dictionary<string, SortedSet<string>>(StringComparer.Ordinal);
                var evidence = new Dictionary<string, SortedSet<string>>(StringComparer.Ordinal);
                var evidenceByMember = new Dictionary<string, SortedDictionary<string, SortedSet<string>>>(StringComparer.Ordinal);
                bool anyNoEntry = false;

                foreach (string member in members)
                {
                    var r = resolution[member];
                    if (r.Kind == MemberResolutionKind.NoEntry)
                    {
                        anyNoEntry = true;
                        continue;
                    }
                    foreach (var d in direct[r.AnnotatedAs])
                    {
                        bool inherited = r.Kind == MemberResolutionKind.Inherited || d.Inherited == true;
                        foreach (string t in up[d.GoId])
                        {
                            Add(used, t, member);
                            if (t == d.GoId) Add(directMembers, t, member);
                            if (inherited) Add(inheritedMembers, t, member);
                            foreach (string e in d.Evidence) Add(evidence, t, e);
                            if (!evidenceByMember.TryGetValue(t, out var byMember))
                                evidenceByMember[t] = byMember = new SortedDictionary<string, SortedSet<string>>(StringComparer.Ordinal);
                            if (!byMember.TryGetValue(member, out var codes)) byMember[member] = codes = new SortedSet<string>(StringComparer.Ordinal);
                            codes.UnionWith(d.Evidence);
                        }
                    }
                }

                if (used.Count == 0)
                {
                    var status = group.IsContaminant ? GoAnnotationStatus.Contaminant
                        : anyNoEntry ? GoAnnotationStatus.NoEntry : GoAnnotationStatus.NoGoTerms;
                    output.Add(new GoAnnotationRow(group.Name, Array.Empty<string>(), null, null, null, Array.Empty<string>(),
                        null, null, members.Count, 0, status, group.QValue, release, obo, db, Array.Empty<string>(),
                        Array.Empty<string>(), new Dictionary<string, IReadOnlyList<string>>(), entrapment));
                    continue;
                }
                foreach (var (t, carriers) in used)
                {
                    var d = directMembers.GetValueOrDefault(t) ?? new SortedSet<string>();
                    var i = inheritedMembers.GetValueOrDefault(t) ?? new SortedSet<string>();
                    output.Add(new GoAnnotationRow(group.Name, carriers.ToList(), t, term[t].GoName, term[t].Aspect,
                        evidence.GetValueOrDefault(t)?.ToList() ?? new List<string>(),
                        Inherited: i.Count == carriers.Count, Propagated: d.Count == 0,
                        members.Count, carriers.Count, GoAnnotationStatus.Annotated, group.QValue, release, obo, db,
                        AccessionDirect: d.ToList(), AccessionInherited: i.ToList(),
                        EvidenceByMember: evidenceByMember[t].ToDictionary(m => m.Key, m => (IReadOnlyList<string>)m.Value.ToList()),
                        EntrapmentMembers: entrapment));
                }
            }
            return output;
        }

        private static void Add<T>(IDictionary<string, SortedSet<T>> sets, string key, T value)
        {
            if (!sets.TryGetValue(key, out var set)) sets[key] = set = new SortedSet<T>();
            set.Add(value);
        }

        private static string Serialize(IEnumerable<GoAnnotationRow> rows, IReadOnlyCollection<string> unresolved)
        {
            using var writer = new StringWriter();
            GoAnnotationTsv.Write(writer, rows, _go, DbSha, null, unresolved);
            return writer.ToString();
        }

        [Test]
        public void Join_OfTheStoredTables_EqualsTheAnnotator()
        {
            var (database, groups) = Fixture();
            var annotation = new GoAccessionAnnotation(_go, database, DbSha, skipUnknownGoIds: true);
            var resolver = new MemberResolver(database);
            var members = resolver.ResolveAll(groups.SelectMany(g => g.MemberAccessions));

            // Through the files, as a store holds them: nothing is passed in memory.
            var a = GoAccessionTermsTsv.Read(WriteFile("a.tsv", w => GoAccessionTermsTsv.Write(w, annotation)));
            var c = MemberResolutionTsv.Read(WriteFile("c.tsv", w => MemberResolutionTsv.Write(w, members, DbSha)));
            var closure = GoOntologyTsv.ReadClosure(WriteFile("closure.tsv", w => GoOntologyTsv.WriteClosure(w, _go)));
            var terms = GoOntologyTsv.ReadTerms(WriteFile("terms.tsv", w => GoOntologyTsv.WriteTerms(w, _go)));

            var expected = new GoGroupAnnotator(annotation, resolver).AnnotateAll(groups);
            var joined = JoinView(groups, c, a, closure, terms);

            Assert.That(Serialize(joined, annotation.UnresolvedGoIds), Is.EqualTo(Serialize(expected, annotation.UnresolvedGoIds)));
            Assert.That(expected.Select(r => r.Status).Distinct(), Is.EquivalentTo(Enum.GetValues<GoAnnotationStatus>()),
                "the fixture reaches every status");
        }

        [Test]
        public void Join_TheAnnotatorBuiltFromProteins_EqualsTheOneBuiltFromItsParts()
        {
            var (database, groups) = Fixture();

            var fromProteins = new GoGroupAnnotator(_go, database, DbSha, skipUnknownGoIds: true).AnnotateAll(groups);
            var fromParts = new GoGroupAnnotator(new GoAccessionAnnotation(_go, database, DbSha, skipUnknownGoIds: true),
                new MemberResolver(database)).AnnotateAll(groups);

            Assert.That(Serialize(fromParts, new[] { "GO:9999999" }), Is.EqualTo(Serialize(fromProteins, new[] { "GO:9999999" })));
        }

        [Test]
        public void Annotator_FromParts_NullArguments_Throw()
        {
            Assert.Throws<ArgumentNullException>(() => new GoGroupAnnotator(null, new MemberResolver(Array.Empty<Protein>())));
            Assert.Throws<ArgumentNullException>(() => new GoGroupAnnotator(new GoAccessionAnnotation(_go, Array.Empty<Protein>(), DbSha), null));
        }
    }
}
