using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases.Annotation;
using UsefulProteomicsDatabases.Pathways;
using static Test.FileReadingTests.TestReactomeMembership;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The acceptance test of kegg:DEF-REACTOME-MEMBERSHIP v1. The stored tables -- a search's member resolution (c),
    /// the database's Reactome membership (T1) and KEGG gene links (T4) -- are written to disk and read back, joined by
    /// equality only as a catalog view would, and compared with a hand-written expectation. Reactome has no production
    /// group annotator, so the expectation is written by hand to be independent of the code under test.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestReactomePerDatabaseTables
    {
        private static readonly string[] MemberLadder = { "no_entry", "species_not_covered", "no_reference" };

        private string _dir;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "ReactomePerDatabase_" + Guid.NewGuid().ToString("N"));
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

        internal sealed record Group(string Name, string[] Members, bool IsContaminant);

        internal sealed record MemberAnswer(string Status, string[] Pathways);

        internal sealed record GroupAnswer(string Status, int Members, IReadOnlyDictionary<string, string[]> Carriers, string[] Entrapment);

        private static Group G(string name, string members, bool contaminant = false) => new(name, members.Split('|'), contaminant);

        /// <summary>
        /// Groups over <see cref="TestReactomeMembership.Fixture"/>: each status on the ladder; a contaminant group whose
        /// member is a target entry, and one with no pathway; a group whose only annotated member is foreign-proteome
        /// entrapment; members absent from the database (a collision counter, an entrapment copy, a contaminant); an
        /// isoform the database holds as its own entry and one it does not; the displayed isoform written as its own id
        /// (the X-1 corner); variants of isoforms; and organisms kept unfiltered.
        /// </summary>
        internal static Group[] Groups() => new[]
        {
            G("G1", "P11111|P11111_A20T"),
            G("G2", "P22222|P22222-3|P22222-2"),
            G("G3", "P11111", contaminant: true),
            G("G4", "CON_P02768", contaminant: true),
            G("G5", "Random_foreign_Q99999|P55555"),
            G("G6", "P11111_2|P55555"),
            G("G7", "P11111_2|Q11111"),
            G("G8", "Q11111|P55555"),
            G("G9", "P55555"),
            G("G10", "P11111|Random_P11111_f0|CON_P02768"),
            G("G11", "P22222-3_A20T|P22222-5|P22222-4"),
            G("G12", "Q33333|Q44441"),
            G("G13", "P33333|P66666"),
            G("G14", "CON_P02768"),
            G("G15", "P77777|P77777-2|P77777-3"),
            G("G16", "P77777-3"),
        };

        private static readonly string[] None = Array.Empty<string>();

        /// <summary>Each member's pathways and status, by hand from the definition.</summary>
        private static readonly Dictionary<string, MemberAnswer> ExpectedMembers = new(StringComparer.Ordinal)
        {
            ["P11111"] = new("annotated", new[] { "R-HSA-1001", "R-HSA-1002" }),
            ["P11111_A20T"] = new("annotated", new[] { "R-HSA-1001", "R-HSA-1002" }),
            ["P11111_2"] = new("no_entry", None), // a collision counter names another entry, never P11111's
            ["P22222"] = new("annotated", new[] { "R-HSA-2001", "R-HSA-2003" }), // the displayed isoform's row and the whole entry
            ["P22222-2"] = new("annotated", new[] { "R-HSA-2003" }), // the X-1 corner: whole-entry rows only, not repaired
            ["P22222-3"] = new("annotated", new[] { "R-HSA-2002", "R-HSA-2003" }), // inherited whole entry + its isoform row
            ["P22222-3_A20T"] = new("annotated", new[] { "R-HSA-2002", "R-HSA-2003" }),
            ["P22222-4"] = new("annotated", new[] { "R-HSA-2004" }), // its own reference only
            ["P22222-5"] = new("annotated", new[] { "R-HSA-2003" }), // absent isoform: the entry's whole-entry rows
            ["P33333"] = new("annotated", new[] { "R-HSA-3001" }), // R-MMU-3002 dropped
            ["P55555"] = new("no_reference", None),
            ["P66666"] = new("annotated", new[] { "R-HSA-6001", "R-HSA-6002" }),
            ["P77777"] = new("no_reference", None), // D20: the entry's only reference is P77777-2's, and P77777 displays no isoform
            ["P77777-2"] = new("annotated", new[] { "R-HSA-7701" }),
            ["P77777-3"] = new("no_reference", None), // D20: not "annotated" with nothing
            ["Q11111"] = new("species_not_covered", None),
            ["Q33333"] = new("annotated", new[] { "R-HSA-4003" }), // rat is unvalidated: kept unfiltered
            ["Q44441"] = new("annotated", new[] { "R-GGA-5001", "R-GGA-5002", "R-GGA-5003" }),
            ["Random_foreign_Q99999"] = new("annotated", new[] { "R-SCE-7001" }),
            ["Random_P11111_f0"] = new("no_entry", None), // entrapment never borrows its target's rows
            ["CON_P02768"] = new("no_entry", None),
        };

        private static Dictionary<string, string[]> C(params (string Pathway, string Carriers)[] pairs) =>
            pairs.ToDictionary(p => p.Pathway, p => p.Carriers.Split('|'), StringComparer.Ordinal);

        /// <summary>Each group's answer, by hand from the definition's group rules.</summary>
        private static readonly Dictionary<string, GroupAnswer> ExpectedGroups = new(StringComparer.Ordinal)
        {
            ["G1"] = new("annotated", 2, C(("R-HSA-1001", "P11111|P11111_A20T"), ("R-HSA-1002", "P11111|P11111_A20T")), None),
            ["G2"] = new("annotated", 3, C(("R-HSA-2001", "P22222"), ("R-HSA-2002", "P22222-3"), ("R-HSA-2003", "P22222|P22222-2|P22222-3")), None),
            ["G3"] = new("annotated", 1, C(("R-HSA-1001", "P11111"), ("R-HSA-1002", "P11111")), None),
            ["G4"] = new("contaminant", 1, C(), None),
            ["G5"] = new("annotated", 2, C(("R-SCE-7001", "Random_foreign_Q99999")), new[] { "Random_foreign_Q99999" }),
            ["G6"] = new("no_entry", 2, C(), None),
            ["G7"] = new("no_entry", 2, C(), None),
            ["G8"] = new("species_not_covered", 2, C(), None),
            ["G9"] = new("no_reference", 1, C(), None),
            ["G10"] = new("annotated", 3, C(("R-HSA-1001", "P11111"), ("R-HSA-1002", "P11111")), new[] { "Random_P11111_f0" }),
            ["G11"] = new("annotated", 3, C(("R-HSA-2002", "P22222-3_A20T"), ("R-HSA-2003", "P22222-3_A20T|P22222-5"), ("R-HSA-2004", "P22222-4")), None),
            ["G12"] = new("annotated", 2, C(("R-GGA-5001", "Q44441"), ("R-GGA-5002", "Q44441"), ("R-GGA-5003", "Q44441"), ("R-HSA-4003", "Q33333")), None),
            ["G13"] = new("annotated", 2, C(("R-HSA-3001", "P33333"), ("R-HSA-6001", "P66666"), ("R-HSA-6002", "P66666")), None),
            ["G14"] = new("no_entry", 1, C(), None), // the same member as G4 without the group flag: contaminant is the flag's
            ["G15"] = new("annotated", 3, C(("R-HSA-7701", "P77777-2")), None),
            ["G16"] = new("no_reference", 1, C(), None),
        };

        /// <summary>The stored tables of one database and one search, through files only.</summary>
        private (AnnotationTable<MemberResolution> C, AnnotationTable<ReactomeMembershipRow> T1, AnnotationTable<KeggGeneLinkRow> T4)
            StoredTables(Protein[] database, string tag = "")
        {
            var members = new MemberResolver(database).ResolveAll(Groups().SelectMany(g => g.Members));
            var annotation = new ReactomeAccessionAnnotation(database, DbSha);
            var c = MemberResolutionTsv.Read(WriteFile($"c{tag}.tsv", w => MemberResolutionTsv.Write(w, members, DbSha)));
            var t1 = ReactomeMembershipTsv.Read(WriteFile($"t1{tag}.tsv", w => ReactomeMembershipTsv.Write(w, annotation)));
            var t4 = KeggGeneLinksTsv.Read(WriteFile($"t4{tag}.tsv", w => KeggGeneLinksTsv.Write(w, KeggGeneLinksTsv.Links(database), DbSha)));
            return (c, t1, t4);
        }

        /// <summary>
        /// The catalog view's member join, from the stored tables with equality joins only: whole-entry rows on
        /// annotated_as, isoform rows on sequence_accession, DISTINCT. no_entry when (c) says so; annotated when any row
        /// carries a pathway; otherwise the status of the annotated_as entry, except that an annotated entry gives
        /// no_reference (D20: its references are scoped to other isoforms).
        /// </summary>
        internal static Dictionary<string, MemberAnswer> JoinMembers(AnnotationTable<MemberResolution> c, AnnotationTable<ReactomeMembershipRow> t1)
        {
            var whole = t1.Rows.Where(r => r.SequenceAccession == null).ToLookup(r => r.Accession, StringComparer.Ordinal);
            var isoform = t1.Rows.Where(r => r.SequenceAccession != null).ToLookup(r => r.SequenceAccession, StringComparer.Ordinal);
            var entryStatus = t1.Rows.GroupBy(r => r.Accession, StringComparer.Ordinal)
                .ToDictionary(g => g.Key, g => ReactomeMembershipTsv.StatusName(g.First().Status), StringComparer.Ordinal);

            var answers = new Dictionary<string, MemberAnswer>(StringComparer.Ordinal);
            foreach (var r in c.Rows)
            {
                if (r.Kind == MemberResolutionKind.NoEntry)
                {
                    answers[r.Accession] = new MemberAnswer("no_entry", None);
                    continue;
                }
                var pathways = whole[r.AnnotatedAs].Concat(isoform[r.SequenceAccession])
                    .Where(row => row.PathwayId != null)
                    .Select(row => row.PathwayId)
                    .Distinct(StringComparer.Ordinal)
                    .OrderBy(p => p, StringComparer.Ordinal)
                    .ToArray();
                string entry = entryStatus[r.AnnotatedAs];
                answers[r.Accession] = pathways.Length > 0
                    ? new MemberAnswer("annotated", pathways)
                    : new MemberAnswer(entry == "annotated" ? "no_reference" : entry, None);
            }
            return answers;
        }

        /// <summary>
        /// The catalog view's group step: one member set for status, union and counts. annotated when any member carries a
        /// pathway; contaminant from the group's flag; otherwise the highest member status on the ladder. Entrapment is a
        /// label.
        /// </summary>
        internal static GroupAnswer JoinGroup(Group group, IReadOnlyDictionary<string, MemberAnswer> members,
            AnnotationTable<MemberResolution> c)
        {
            var distinct = group.Members.Distinct(StringComparer.Ordinal).ToList();
            var carriers = distinct
                .SelectMany(m => members[m].Pathways.Select(p => (Pathway: p, Member: m)))
                .GroupBy(x => x.Pathway, StringComparer.Ordinal)
                .OrderBy(x => x.Key, StringComparer.Ordinal)
                .ToDictionary(x => x.Key, x => x.Select(y => y.Member).OrderBy(m => m, StringComparer.Ordinal).ToArray(), StringComparer.Ordinal);
            string status = carriers.Count > 0 ? "annotated"
                : group.IsContaminant ? "contaminant"
                : MemberLadder.First(s => distinct.Any(m => members[m].Status == s));
            var resolution = c.Rows.ToDictionary(r => r.Accession, StringComparer.Ordinal);
            var entrapment = distinct.Where(m => resolution[m].IsEntrapment).OrderBy(m => m, StringComparer.Ordinal).ToArray();
            return new GroupAnswer(status, distinct.Count, carriers, entrapment);
        }

        private static string Show(MemberAnswer a) => $"{a.Status}[{string.Join(",", a.Pathways)}]";

        private static string Show(GroupAnswer a) =>
            $"{a.Status} n={a.Members} {{{string.Join("; ", a.Carriers.Select(kv => $"{kv.Key}:{string.Join(",", kv.Value)}"))}}} entrapment=[{string.Join(",", a.Entrapment)}]";

        private Dictionary<string, GroupAnswer> JoinAll(out Dictionary<string, MemberAnswer> members)
        {
            var (c, t1, _) = StoredTables(Fixture());
            members = JoinMembers(c, t1);
            var m = members;
            return Groups().ToDictionary(g => g.Name, g => JoinGroup(g, m, c), StringComparer.Ordinal);
        }

        [Test]
        public void Join_OfTheStoredTables_EqualsTheHandWrittenExpectation()
        {
            var groups = JoinAll(out var members);

            Assert.That(members.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key} {Show(kv.Value)}"),
                Is.EqualTo(ExpectedMembers.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key} {Show(kv.Value)}")));
            Assert.That(groups.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key} {Show(kv.Value)}"),
                Is.EqualTo(ExpectedGroups.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key} {Show(kv.Value)}")));
        }

        [Test]
        public void Join_KeggLinksFollowAnnotatedAs()
        {
            var (c, _, t4) = StoredTables(Fixture());
            var links = t4.Rows.ToLookup(r => r.Accession, r => r.KeggGeneId, StringComparer.Ordinal);
            var resolution = c.Rows.ToDictionary(r => r.Accession, StringComparer.Ordinal);

            string[] LinksOf(string member) => resolution[member].AnnotatedAs is { } entry ? links[entry].ToArray() : None;

            Assert.That(LinksOf("P22222-5"), Is.EqualTo(new[] { "hsa:2222", "hsa:2223" }));
            Assert.That(LinksOf("Q11111"), Is.EqualTo(new[] { "mmu:1111" }));
            Assert.That(LinksOf("P55555"), Is.Empty);
            Assert.That(LinksOf("P11111_2"), Is.Empty);
        }

        [Test]
        public void Join_SameMemberSameAnswer_WhetherOrNotTheDatabaseHoldsTheIsoform()
        {
            var (cWith, t1With, _) = StoredTables(Fixture(withIsoformEntryP22222_3: true), "with");
            var (cWithout, t1Without, _) = StoredTables(Fixture(withIsoformEntryP22222_3: false), "without");
            var with = JoinMembers(cWith, t1With);
            var without = JoinMembers(cWithout, t1Without);

            Assert.That(cWith.Rows.Single(r => r.Accession == "P22222-3").Kind, Is.EqualTo(MemberResolutionKind.Self));
            Assert.That(cWithout.Rows.Single(r => r.Accession == "P22222-3").Kind, Is.EqualTo(MemberResolutionKind.Inherited));
            foreach (string member in new[] { "P22222-3", "P22222-3_A20T" })
            {
                Assert.That(Show(without[member]), Is.EqualTo(Show(with[member])), member);
                Assert.That(Show(with[member]), Is.EqualTo(Show(ExpectedMembers[member])), member);
            }
        }

        [Test]
        public void Join_MemberOutsideItsEntrysIsoformScopedRows_ReadsNoReference()
        {
            var (_, t1, _) = StoredTables(Fixture());
            var groups = JoinAll(out var members);

            Assert.That(t1.Rows.Where(r => r.Accession == "P77777").Select(r => r.Status), Is.All.EqualTo(ReactomeMembershipStatus.Annotated));
            Assert.That(Show(members["P77777-3"]), Is.EqualTo("no_reference[]"));
            Assert.That(Show(members["P77777-2"]), Is.EqualTo("annotated[R-HSA-7701]"));
            Assert.That(groups["G16"].Status, Is.EqualTo("no_reference"));
        }

        [Test]
        public void Group_StatusIsTheHighestOnTheLadder()
        {
            var groups = JoinAll(out _);

            Assert.That(groups["G6"].Status, Is.EqualTo("no_entry"), "no_entry > no_reference");
            Assert.That(groups["G7"].Status, Is.EqualTo("no_entry"), "no_entry > species_not_covered");
            Assert.That(groups["G8"].Status, Is.EqualTo("species_not_covered"), "species_not_covered > no_reference");
            Assert.That(groups["G9"].Status, Is.EqualTo("no_reference"));
            Assert.That(groups["G4"].Status, Is.EqualTo("contaminant"), "contaminant > no_entry");
        }

        [Test]
        public void Group_ContaminantComesFromTheGroupFlag_AndAnnotatedStillWins()
        {
            var groups = JoinAll(out _);

            Assert.That(groups["G3"].Status, Is.EqualTo("annotated"));
            Assert.That(groups["G3"].Carriers.Keys, Is.EqualTo(new[] { "R-HSA-1001", "R-HSA-1002" }));
            Assert.That(groups["G4"].Status, Is.EqualTo("contaminant"));
            Assert.That(groups["G14"].Status, Is.EqualTo("no_entry"), "the same CON_ member without the flag");
        }

        [Test]
        public void Group_CountsAndUnionUseEveryMember()
        {
            var g10 = JoinAll(out _)["G10"];

            Assert.That(g10.Members, Is.EqualTo(3), "no_entry, entrapment and contaminant members all count");
            Assert.That(g10.Carriers.Keys, Is.EqualTo(new[] { "R-HSA-1001", "R-HSA-1002" }));
            Assert.That(g10.Carriers.Values.Select(v => v.Length), Is.All.EqualTo(1));
            Assert.That(g10.Entrapment, Is.EqualTo(new[] { "Random_P11111_f0" }));
        }

        [Test]
        public void Group_ForeignEntrapmentOnlyAnnotatedMember_ReadsAnnotatedWithItsPathwaysAndIsLabelled()
        {
            var g5 = JoinAll(out _)["G5"];

            Assert.That(g5.Status, Is.EqualTo("annotated"));
            Assert.That(g5.Carriers.Keys, Is.EqualTo(new[] { "R-SCE-7001" }));
            Assert.That(g5.Entrapment, Is.EqualTo(new[] { "Random_foreign_Q99999" }));
            Assert.That(g5.Members, Is.EqualTo(2));
        }

        [Test]
        public void Fixture_ReachesEveryStatusAndRule()
        {
            var (c, t1, _) = StoredTables(Fixture());
            var groups = JoinAll(out var members);

            Assert.That(t1.Rows.Select(r => r.Status).Distinct(), Is.EquivalentTo(Enum.GetValues<ReactomeMembershipStatus>()));
            Assert.That(t1.Rows.Any(r => r.Inherited == true), "an inherited row");
            Assert.That(t1.Rows.Any(r => r.Molecule != null && r.SequenceAccession != r.Molecule), "a displayed-isoform row");
            Assert.That(t1.Rows.Any(r => r.Molecule != null && r.SequenceAccession == r.Molecule), "another isoform's row");
            Assert.That(t1.Header["dropped_foreign_references"], Is.Not.EqualTo("0"));
            Assert.That(t1.Header["unvalidated_organisms"].Split(';'), Has.Length.EqualTo(3));
            Assert.That(c.Rows.Select(r => r.Kind).Distinct(), Is.EquivalentTo(Enum.GetValues<MemberResolutionKind>()));
            Assert.That(c.Rows.Any(r => r.IsEntrapment && r.Kind == MemberResolutionKind.NoEntry), "an entrapment copy absent");
            Assert.That(c.Rows.Any(r => r.IsEntrapment && r.Kind == MemberResolutionKind.Self), "a foreign entrapment entry");
            Assert.That(members.Values.Select(m => m.Status).Distinct(),
                Is.EquivalentTo(new[] { "annotated", "no_entry", "no_reference", "species_not_covered" }));
            Assert.That(groups.Values.Select(g => g.Status).Distinct(),
                Is.EquivalentTo(new[] { "annotated", "contaminant", "no_entry", "species_not_covered", "no_reference" }));
            Assert.That(ExpectedMembers.Keys, Is.EquivalentTo(Groups().SelectMany(g => g.Members).Distinct()), "every member has an expectation");
        }
    }
}
