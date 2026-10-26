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
    /// Reactome pathways and KEGG gene ids are already parsed and stored: ProteinXmlEntry keeps every
    /// &lt;dbReference&gt; generically, and since the isoform-scope change it keeps each one's &lt;molecule&gt; too.
    /// Protein.ReactomePathways and Protein.KeggGeneIds are derived views over that list, following GoTerms and
    /// EnsemblGeneIds: a const and a computed property on Protein, no constructor parameter and no stored field.
    ///
    /// What these pin:
    ///  - filtering by the exact Type, because the flat list also holds GO, PubMed, organism and every other reference;
    ///  - one pathway per Reactome id, with the isoforms of repeats unioned, and a whole-entry reference covering any
    ///    isoform-scoped repeat of the same id (Molecules empty means the whole entry);
    ///  - properties matched by type, never by position, because ProteinDbWriter re-sorts them;
    ///  - a KEGG reference is a bare gene id (hsa:2597) with no pathway: UniProt carries no KEGG pathway membership.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestReactomePathways
    {
        private static string Data(params string[] parts) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, Path.Combine(parts));

        private static List<Protein> Load(string fileName, DecoyType decoyType = DecoyType.None) =>
            ProteinDbLoader.LoadProteinXML(Data("DatabaseTests", fileName), true, decoyType, null, false, null, out _);

        private static Protein LoadGapdh() => Load("humanGAPDH.xml").Single(p => p.Accession == "P04406");

        private static Protein WithReferences(params DatabaseReference[] references) =>
            new Protein("PEPTIDEK", "P00001", databaseReferences: references.ToList());

        private static DatabaseReference Reactome(string id, string name, string molecule = null) =>
            new DatabaseReference("Reactome", id,
                name == null ? new List<Tuple<string, string>>() : new List<Tuple<string, string>> { new Tuple<string, string>("pathway name", name) },
                molecule);

        private static DatabaseReference Kegg(string id) =>
            new DatabaseReference("KEGG", id, new List<Tuple<string, string>>());

        [Test]
        public void ReactomePathways_ReadFromHumanGapdh()
        {
            var pathways = LoadGapdh().ReactomePathways;

            Assert.That(pathways.Select(p => (p.Id, p.Name)), Is.EquivalentTo(new[]
            {
                ("R-HSA-70171", "Glycolysis"),
                ("R-HSA-70263", "Gluconeogenesis"),
            }));
            Assert.That(pathways.All(p => p.Molecules.Count == 0), Is.True, "GAPDH's Reactome references cover the whole entry");
        }

        [Test]
        public void ReactomePathways_EmptyRatherThanNullWhenNothingIsAnnotated()
        {
            Assert.That(new Protein("PEPTIDEK", "P00001").ReactomePathways, Is.Not.Null.And.Empty);
            Assert.That(ReactomePathway.FromDatabaseReferences(null), Is.Not.Null.And.Empty);
        }

        [Test]
        public void ReactomePathways_FilterIsExactType()
        {
            var references = new List<DatabaseReference>
            {
                Reactome("R-HSA-70171", "Glycolysis"),
                new DatabaseReference("reactome", "R-HSA-1", null),
                new DatabaseReference("REACTOME", "R-HSA-2", null),
                new DatabaseReference("Reactome ", "R-HSA-3", null),
                new DatabaseReference("GO", "GO:0005737", null),
                null,
            };

            Assert.That(ReactomePathway.FromDatabaseReferences(references).Select(p => p.Id), Is.EqualTo(new[] { "R-HSA-70171" }));
        }

        [Test]
        public void ReactomePathways_DuplicateIdsMergedAndMoleculesUnioned()
        {
            var protein = WithReferences(
                Reactome("R-HSA-77289", "Mitochondrial Fatty Acid Beta-Oxidation", "P00001-2"),
                Reactome("R-HSA-77289", "Mitochondrial Fatty Acid Beta-Oxidation", "P00001-3"));

            var pathway = protein.ReactomePathways.Single();

            Assert.That(pathway.Id, Is.EqualTo("R-HSA-77289"));
            Assert.That(pathway.Name, Is.EqualTo("Mitochondrial Fatty Acid Beta-Oxidation"));
            Assert.That(pathway.Molecules, Is.EquivalentTo(new[] { "P00001-2", "P00001-3" }));
        }

        [Test]
        public void ReactomePathways_WholeEntryReferenceCoversTheIsoformScopedOne()
        {
            var protein = WithReferences(
                Reactome("R-HSA-77289", "Mitochondrial Fatty Acid Beta-Oxidation", "P00001-2"),
                Reactome("R-HSA-77289", "Mitochondrial Fatty Acid Beta-Oxidation"));

            Assert.That(protein.ReactomePathways.Single().Molecules, Is.Empty,
                "a reference for the whole entry covers every isoform, so the merged pathway is whole-entry");
        }

        [Test]
        public void ReactomePathways_MoleculesReadFromIsoformScopeFixture()
        {
            var protein = Load("isoformScope.xml").Single(p => p.Accession == "Q9ZZZ1");

            var byId = protein.ReactomePathways.ToDictionary(p => p.Id, p => p.Molecules.OrderBy(m => m, StringComparer.Ordinal).ToList());

            Assert.That(byId.Keys, Is.EquivalentTo(new[] { "R-HSA-77289", "R-HSA-70171", "R-HSA-556833" }));
            Assert.That(byId["R-HSA-77289"], Is.EqualTo(new[] { "Q9ZZZ1-1" }));
            Assert.That(byId["R-HSA-556833"], Is.EqualTo(new[] { "Q9ZZZ1-2" }));
            Assert.That(byId["R-HSA-70171"], Is.Empty);
        }

        [Test]
        public void ReactomePathways_MissingNameKeptAsEmpty()
        {
            var pathway = WithReferences(Reactome("R-HSA-70171", null)).ReactomePathways.Single();

            Assert.That(pathway.Id, Is.EqualTo("R-HSA-70171"), "a pathway without a name is still a membership; dropping it would read as 'not annotated'");
            Assert.That(pathway.Name, Is.EqualTo(""));
        }

        [Test]
        public void ReactomePathways_PropertiesMatchedByTypeNotPosition()
        {
            var reference = new DatabaseReference("Reactome", "R-HSA-70171", new List<Tuple<string, string>>
            {
                new Tuple<string, string>("evidence", "ECO:0000304"),
                new Tuple<string, string>("pathway name", "Glycolysis"),
            });

            Assert.That(WithReferences(reference).ReactomePathways.Single().Name, Is.EqualTo("Glycolysis"));
        }

        [Test]
        public void Decoys_CarryNoReactomePathways()
        {
            var decoys = Load("isoformScope.xml", DecoyType.Reverse).Where(p => p.IsDecoy).ToList();

            Assert.That(decoys, Is.Not.Empty);
            Assert.That(decoys.SelectMany(d => d.ReactomePathways), Is.Empty, "only the taxonomy reference travels onto a decoy");
            Assert.That(decoys.SelectMany(d => d.KeggGeneIds), Is.Empty);
        }

        [Test]
        public void KeggGeneIds_ReadFromHumanGapdh()
        {
            var fixture = Load("isoformScope.xml");

            Assert.That(LoadGapdh().KeggGeneIds, Is.EqualTo(new[] { "hsa:2597" }));
            Assert.That(fixture.Single(p => p.Accession == "Q9ZZZ1").KeggGeneIds, Is.EqualTo(new[] { "hsa:9000001" }));
            Assert.That(fixture.Single(p => p.Accession == "Q9ZZZ2").KeggGeneIds, Is.Not.Null.And.Empty);
        }

        [Test]
        public void KeggGeneIds_DistinctAndInOrdinalOrder()
        {
            var protein = WithReferences(Kegg("mmu:2"), Kegg("mmu:10"), Kegg("mmu:2"), new DatabaseReference("kegg", "mmu:3", null));

            Assert.That(protein.KeggGeneIds, Is.EqualTo(new[] { "mmu:10", "mmu:2" }));
        }

        [Test]
        public void PathwayDatabaseReferenceTypes_AreReexportedOnProteinDbLoader()
        {
            Assert.That(ProteinDbLoader.ReactomeDatabaseReferenceType, Is.EqualTo("Reactome"));
            Assert.That(Protein.ReactomeDatabaseReferenceType, Is.EqualTo("Reactome"));
            Assert.That(ProteinDbLoader.KeggDatabaseReferenceType, Is.EqualTo("KEGG"));
            Assert.That(Protein.KeggDatabaseReferenceType, Is.EqualTo("KEGG"));
        }

        [Test]
        public void DisplayedIsoformId_IsTheSequenceAttributesValue()
        {
            var fixture = Load("isoformScope.xml");

            Assert.That(fixture.Single(p => p.Accession == "Q9ZZZ1").DisplayedIsoformId, Is.EqualTo("Q9ZZZ1-2"));
            Assert.That(fixture.Single(p => p.Accession == "Q9ZZZ2").DisplayedIsoformId, Is.Null);
            Assert.That(LoadGapdh().DisplayedIsoformId, Is.EqualTo("P04406-1"));
        }
    }
}
