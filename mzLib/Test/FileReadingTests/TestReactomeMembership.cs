using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases.Pathways;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Reactome membership per annotation database (kegg:DEF-REACTOME-MEMBERSHIP v1): the rows
    /// <see cref="ReactomeAccessionAnnotation"/> builds from a target set's proteins, rule by rule, over
    /// <see cref="Fixture"/>.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestReactomeMembership
    {
        internal const string DbSha = "0123abcd";
        internal const string Human = "9606";
        internal const string Mouse = "10090";
        internal const string Rat = "10116";
        internal const string Chicken = "9031";
        internal const string Yeast = "559292";

        internal static DatabaseReference Taxon(string id) =>
            new("NCBI Taxonomy", id, new List<Tuple<string, string>>());

        internal static DatabaseReference Reactome(string id, string molecule = null) =>
            new("Reactome", id, new[] { Tuple.Create(ReactomePathway.PathwayNamePropertyType, "name of " + id) }, molecule);

        internal static DatabaseReference Kegg(string id) =>
            new("KEGG", id, new List<Tuple<string, string>>());

        internal static Protein P(string accession, string taxon, params DatabaseReference[] references) =>
            new("PEPTIDEK", accession, databaseReferences: references.Prepend(Taxon(taxon)).ToList());

        internal static Protein Displaying(string accession, string displayedIsoform, string taxon, params DatabaseReference[] references) =>
            new("PEPTIDEK", accession, databaseReferences: references.Prepend(Taxon(taxon)).ToList(),
                uniProtSequenceAttributes: new UniProtSequenceAttributes(8, 0, "", DateTime.MinValue, 1, null,
                    UniProtSequenceAttributes.FragmentType.unspecified, displayedIsoform));

        /// <summary>
        /// A target set that reaches every rule of the definition:
        /// <list type="bullet">
        /// <item>P11111: two whole-entry pathways and one KEGG id; P11111_A20T is an applied variant of it.</item>
        /// <item>P22222 displays isoform P22222-2: R-HSA-2001 is scoped to the displayed isoform, R-HSA-2002 to
        /// P22222-3, and R-HSA-2003 covers the whole entry.</item>
        /// <item>P22222-3: an isoform entry with no reference of its own, so it inherits R-HSA-2003 only (D37).</item>
        /// <item>P22222-4: an isoform entry with its own reference, which it keeps alone.</item>
        /// <item>P33333: human, with one foreign-prefix id (R-MMU-3002), dropped and counted.</item>
        /// <item>P55555: human with no reference: no_reference.</item>
        /// <item>P66666 twice, with different references, which are unioned.</item>
        /// <item>P77777: its only reference is scoped to isoform P77777-2, so a member of another of its isoforms takes
        /// no pathway although the entry reads annotated (D20).</item>
        /// <item>Human has 13 references, 12 of them HSA: a strict majority of at least 10, so its code is HSA.</item>
        /// <item>Q11111, Q22222: mouse, with no Reactome reference at all: species_not_covered.</item>
        /// <item>Q33331-Q33333: rat, 3 references (one R-HSA): below the minimum, kept unfiltered, unvalidated.</item>
        /// <item>Q44441-Q44444: chicken, 12 references split 6 GGA / 6 HSA: a tie, kept unfiltered, unvalidated.</item>
        /// <item>Random_foreign_Q99999: a foreign-proteome entrapment entry (yeast), annotated by its own reference.</item>
        /// <item>A decoy with references, which is skipped.</item>
        /// </list>
        /// </summary>
        internal static Protein[] Fixture(bool withIsoformEntryP22222_3 = true)
        {
            var proteins = new List<Protein>
            {
                P("P11111", Human, Reactome("R-HSA-1001"), Reactome("R-HSA-1002"), Kegg("hsa:1111")),
                P("P11111_A20T", Human, Reactome("R-HSA-1001"), Reactome("R-HSA-1002"), Kegg("hsa:1111")),
                Displaying("P22222", "P22222-2", Human,
                    Reactome("R-HSA-2001", "P22222-2"), Reactome("R-HSA-2002", "P22222-3"), Reactome("R-HSA-2003"),
                    Kegg("hsa:2223"), Kegg("hsa:2222")),
                P("P22222-3", Human),
                P("P22222-4", Human, Reactome("R-HSA-2004")),
                P("P33333", Human, Reactome("R-HSA-3001"), Reactome("R-MMU-3002")),
                P("P55555", Human),
                P("P66666", Human, Reactome("R-HSA-6001"), Kegg("hsa:6666")),
                P("P66666", Human, Reactome("R-HSA-6002"), Kegg("hsa:6666")),
                P("P77777", Human, Reactome("R-HSA-7701", "P77777-2")),
                P("Q11111", Mouse, Kegg("mmu:1111")),
                P("Q22222", Mouse),
                P("Q33331", Rat, Reactome("R-RNO-4001")),
                P("Q33332", Rat, Reactome("R-RNO-4002")),
                P("Q33333", Rat, Reactome("R-HSA-4003")),
                P("Q44441", Chicken, Reactome("R-GGA-5001"), Reactome("R-GGA-5002"), Reactome("R-GGA-5003")),
                P("Q44442", Chicken, Reactome("R-GGA-5004"), Reactome("R-GGA-5005"), Reactome("R-GGA-5006")),
                P("Q44443", Chicken, Reactome("R-HSA-5007"), Reactome("R-HSA-5008"), Reactome("R-HSA-5009")),
                P("Q44444", Chicken, Reactome("R-HSA-5010"), Reactome("R-HSA-5011"), Reactome("R-HSA-5012")),
                P("Random_foreign_Q99999", Yeast, Reactome("R-SCE-7001")),
                new Protein("PEPTIDEK", "DECOY_P11111", isDecoy: true,
                    databaseReferences: new List<DatabaseReference> { Taxon(Human), Reactome("R-HSA-1001"), Kegg("hsa:1111") }),
            };
            if (!withIsoformEntryP22222_3)
            {
                proteins.RemoveAll(p => p.Accession == "P22222-3");
            }
            return proteins.ToArray();
        }

        private static ReactomeAccessionAnnotation Build() => new(Fixture(), DbSha);

        private static ReactomeMembershipRow[] RowsOf(string accession) =>
            Build().Rows().Where(r => r.Accession == accession).ToArray();

        private static string[] PathwaysOf(string accession) =>
            RowsOf(accession).Where(r => r.PathwayId != null).Select(r => r.PathwayId).ToArray();

        [Test]
        public void EveryEntry_HasAtLeastOneRow()
        {
            var annotation = Build();
            var expected = Fixture().Where(p => !p.IsDecoy).Select(p => p.Accession).Distinct().OrderBy(a => a, StringComparer.Ordinal);

            Assert.That(annotation.Accessions, Is.EqualTo(expected));
            Assert.That(annotation.Rows().Select(r => r.Accession).Distinct(), Is.EquivalentTo(expected));
        }

        [Test]
        public void AnnotatedStatus_IffTheRowCarriesAPathway()
        {
            var rows = Build().Rows();

            Assert.That(rows, Is.Not.Empty);
            foreach (var row in rows)
            {
                Assert.That(row.Status == ReactomeMembershipStatus.Annotated, Is.EqualTo(row.PathwayId != null), row.ToString());
            }
            var pathwayless = rows.Where(r => r.PathwayId == null).ToList();
            Assert.That(pathwayless.All(r => r.Molecule == null && r.SequenceAccession == null && r.Inherited == null));
        }

        [Test]
        public void WholeEntryRow_HasNullMoleculeAndSequenceAccession()
        {
            var row = RowsOf("P22222").Single(r => r.PathwayId == "R-HSA-2003");

            Assert.That(row, Is.EqualTo(new ReactomeMembershipRow("P22222", "NCBITaxon:9606", "R-HSA-2003", null, null, false,
                ReactomeMembershipStatus.Annotated)));
        }

        [Test]
        public void DisplayedIsoformRow_SequenceAccessionIsTheEntrysSequenceAccession()
        {
            var row = RowsOf("P22222").Single(r => r.PathwayId == "R-HSA-2001");

            Assert.That(row.Molecule, Is.EqualTo("P22222-2"));
            Assert.That(row.SequenceAccession, Is.EqualTo("P22222"));
            Assert.That(row.Inherited, Is.False);
        }

        [Test]
        public void OtherIsoformRow_SequenceAccessionIsTheMolecule()
        {
            var row = RowsOf("P22222").Single(r => r.PathwayId == "R-HSA-2002");

            Assert.That(row.Molecule, Is.EqualTo("P22222-3"));
            Assert.That(row.SequenceAccession, Is.EqualTo("P22222-3"));
        }

        [Test]
        public void Inheritance_IsoformWithoutReferences_TakesWholeEntryRowsFlaggedInherited()
        {
            Assert.That(RowsOf("P22222-3"), Is.EqualTo(new[]
            {
                new ReactomeMembershipRow("P22222-3", "NCBITaxon:9606", "R-HSA-2003", null, null, true, ReactomeMembershipStatus.Annotated)
            }));
        }

        [Test]
        public void Inheritance_NeverCopiesMoleculeRows()
        {
            var rows = Build().Rows();

            Assert.That(rows.Where(r => r.Inherited == true).All(r => r.Molecule == null && r.SequenceAccession == null));
            Assert.That(RowsOf("P22222-3").Select(r => r.PathwayId), Does.Not.Contain("R-HSA-2001").And.Not.Contain("R-HSA-2002"));
        }

        [Test]
        public void Inheritance_IsoformWithOwnReferences_KeepsOnlyItsOwn()
        {
            Assert.That(RowsOf("P22222-4"), Is.EqualTo(new[]
            {
                new ReactomeMembershipRow("P22222-4", "NCBITaxon:9606", "R-HSA-2004", null, null, false, ReactomeMembershipStatus.Annotated)
            }));
        }

        [Test]
        public void Inheritance_AppliedVariantIsNotAnIsoform_KeepsItsOwnRowsUninherited()
        {
            Assert.That(PathwaysOf("P11111_A20T"), Is.EqualTo(new[] { "R-HSA-1001", "R-HSA-1002" }));
            Assert.That(RowsOf("P11111_A20T").All(r => r.Inherited == false));
        }

        [Test]
        public void Organism_WrittenAsNcbiTaxonCurie()
        {
            Assert.That(RowsOf("P11111").Select(r => r.Organism).Distinct(), Is.EqualTo(new[] { "NCBITaxon:9606" }));
            Assert.That(RowsOf("Q11111").Single().Organism, Is.EqualTo("NCBITaxon:10090"));
            Assert.That(RowsOf("Random_foreign_Q99999").Single().Organism, Is.EqualTo("NCBITaxon:559292"));
        }

        [Test]
        public void Organism_Missing_WrittenEmptyAndTreatedAsOneOrganism()
        {
            var proteins = new[]
            {
                new Protein("PEPTIDEK", "CONTAM_1"),
                new Protein("PEPTIDEK", "CONTAM_2"),
                P("P11111", Human, Reactome("R-HSA-1")),
            };
            var annotation = new ReactomeAccessionAnnotation(proteins, DbSha);
            var rows = annotation.Rows();

            Assert.That(rows.Where(r => r.Accession.StartsWith("CONTAM_")), Is.EqualTo(new[]
            {
                new ReactomeMembershipRow("CONTAM_1", null, null, null, null, null, ReactomeMembershipStatus.SpeciesNotCovered),
                new ReactomeMembershipRow("CONTAM_2", null, null, null, null, null, ReactomeMembershipStatus.SpeciesNotCovered),
            }));
            Assert.That(annotation.UnvalidatedOrganisms.Keys, Is.EqualTo(new[] { "NCBITaxon:9606" }));
        }

        [Test]
        public void Organism_Missing_WithReferences_IsOneOrganismNamedUnknown()
        {
            var proteins = new[]
            {
                new Protein("PEPTIDEK", "X1", databaseReferences: new List<DatabaseReference> { Reactome("R-HSA-1") }),
                new Protein("PEPTIDEK", "X2", databaseReferences: new List<DatabaseReference> { Reactome("R-MMU-2") }),
                new Protein("PEPTIDEK", "X3"),
            };
            var annotation = new ReactomeAccessionAnnotation(proteins, DbSha);

            Assert.That(ReactomeAccessionAnnotation.UnknownOrganism, Is.EqualTo("unknown"));
            Assert.That(annotation.UnvalidatedOrganisms, Is.EqualTo(new Dictionary<string, int> { ["unknown"] = 2 }));
            Assert.That(annotation.Rows().Select(r => (r.Accession, r.Organism, r.PathwayId, r.Status)), Is.EqualTo(new[]
            {
                ("X1", (string)null, "R-HSA-1", ReactomeMembershipStatus.Annotated),
                ("X2", (string)null, "R-MMU-2", ReactomeMembershipStatus.Annotated),
                ("X3", (string)null, (string)null, ReactomeMembershipStatus.NoReference),
            }));
        }

        [Test]
        public void Species_CodeIsTheStrictMajorityPrefix()
        {
            var annotation = Build();

            Assert.That(annotation.SpeciesCodes, Is.EqualTo(new Dictionary<string, string> { ["NCBITaxon:9606"] = "HSA" }));
        }

        [Test]
        public void Species_ForeignPrefixIdsAreDroppedAndCounted()
        {
            var annotation = Build();

            Assert.That(PathwaysOf("P33333"), Is.EqualTo(new[] { "R-HSA-3001" }));
            Assert.That(annotation.DroppedForeignReferences, Is.EqualTo(1));
        }

        [Test]
        public void Species_BelowTenReferences_KeptUnfilteredAndListedUnvalidated()
        {
            var annotation = Build();

            Assert.That(annotation.UnvalidatedOrganisms["NCBITaxon:10116"], Is.EqualTo(3));
            Assert.That(annotation.UnvalidatedOrganisms["NCBITaxon:559292"], Is.EqualTo(1));
            Assert.That(PathwaysOf("Q33333"), Is.EqualTo(new[] { "R-HSA-4003" }), "a foreign-looking id is kept without a code");
            Assert.That(annotation.SpeciesCodes.ContainsKey("NCBITaxon:10116"), Is.False);
        }

        [Test]
        public void Species_NineReferences_IsBelowTheMinimum()
        {
            var proteins = Enumerable.Range(1, 9).Select(i => P($"Q{50000 + i}", Rat, Reactome($"R-RNO-{i}"))).ToList();
            var nine = new ReactomeAccessionAnnotation(proteins, DbSha);
            proteins.Add(P("Q50010", Rat, Reactome("R-RNO-10")));
            var ten = new ReactomeAccessionAnnotation(proteins, DbSha);

            Assert.That(nine.SpeciesCodes, Is.Empty);
            Assert.That(nine.UnvalidatedOrganisms["NCBITaxon:10116"], Is.EqualTo(9));
            Assert.That(ten.SpeciesCodes["NCBITaxon:10116"], Is.EqualTo("RNO"));
            Assert.That(ten.UnvalidatedOrganisms, Is.Empty);
        }

        [Test]
        public void Species_Tie_IsNotAMajority_KeptUnfilteredAndListedUnvalidated()
        {
            var annotation = Build();

            Assert.That(annotation.UnvalidatedOrganisms["NCBITaxon:9031"], Is.EqualTo(12));
            Assert.That(annotation.SpeciesCodes.ContainsKey("NCBITaxon:9031"), Is.False);
            Assert.That(PathwaysOf("Q44443"), Is.EqualTo(new[] { "R-HSA-5007", "R-HSA-5008", "R-HSA-5009" }));
        }

        [Test]
        public void Species_PluralityIsNotAMajority()
        {
            // 5 HSA, 4 GGA, 3 MMU: HSA leads but holds 5 of 12, not more than half.
            var proteins = new[]
            {
                P("Q60001", Chicken, Enumerable.Range(1, 5).Select(i => Reactome($"R-HSA-{i}")).ToArray()),
                P("Q60002", Chicken, Enumerable.Range(1, 4).Select(i => Reactome($"R-GGA-{i}")).ToArray()),
                P("Q60003", Chicken, Enumerable.Range(1, 3).Select(i => Reactome($"R-MMU-{i}")).ToArray()),
            };
            var annotation = new ReactomeAccessionAnnotation(proteins, DbSha);

            Assert.That(annotation.SpeciesCodes, Is.Empty);
            Assert.That(annotation.UnvalidatedOrganisms["NCBITaxon:9031"], Is.EqualTo(12));
            Assert.That(annotation.DroppedForeignReferences, Is.EqualTo(0));
        }

        [Test]
        public void Species_OrganismWithoutAnyReference_IsSpeciesNotCovered()
        {
            Assert.That(RowsOf("Q11111"), Is.EqualTo(new[]
            {
                new ReactomeMembershipRow("Q11111", "NCBITaxon:10090", null, null, null, null, ReactomeMembershipStatus.SpeciesNotCovered)
            }));
            Assert.That(RowsOf("Q22222").Single().Status, Is.EqualTo(ReactomeMembershipStatus.SpeciesNotCovered));
        }

        [Test]
        public void NoReference_InACoveredOrganism()
        {
            Assert.That(RowsOf("P55555"), Is.EqualTo(new[]
            {
                new ReactomeMembershipRow("P55555", "NCBITaxon:9606", null, null, null, null, ReactomeMembershipStatus.NoReference)
            }));
        }

        [Test]
        public void DuplicateAccession_RowsAreUnioned()
        {
            Assert.That(PathwaysOf("P66666"), Is.EqualTo(new[] { "R-HSA-6001", "R-HSA-6002" }));
        }

        [Test]
        public void Decoys_AreSkipped()
        {
            var annotation = Build();

            Assert.That(annotation.Accessions, Does.Not.Contain("DECOY_P11111"));
            Assert.That(new ReactomeAccessionAnnotation(new Protein[] { null, Fixture().Last() }, DbSha).Accessions, Is.Empty);
        }

        [Test]
        public void Rows_AreInOrdinalOrder()
        {
            var rows = Build().Rows();
            var ordered = rows
                .OrderBy(r => r.Accession, StringComparer.Ordinal)
                .ThenBy(r => r.PathwayId ?? "", StringComparer.Ordinal)
                .ThenBy(r => r.Molecule ?? "", StringComparer.Ordinal);

            Assert.That(rows, Is.EqualTo(ordered));
        }

        [TestCase(null)]
        [TestCase("")]
        public void Constructor_RefusesMissingAnnotationDbSha(string sha)
        {
            Assert.Throws<ArgumentException>(() => new ReactomeAccessionAnnotation(Fixture(), sha));
        }

        [Test]
        public void Constructor_RefusesNullProteins()
        {
            Assert.Throws<ArgumentNullException>(() => new ReactomeAccessionAnnotation(null, DbSha));
        }
    }
}
