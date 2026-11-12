using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using UsefulProteomicsDatabases.Taxonomy;

namespace Test.DatabaseTests;

/// <summary>
/// The taxonomy embedded in mzLib (<see cref="TaxonomyTree.LoadPinned"/>, <see cref="ProteomeTaxonMap.LoadPinned"/>).
/// </summary>
/// <remarks>
/// These pin specific lineages on purpose. The snapshot is refreshed by hand (see
/// UsefulProteomicsDatabases/Resources/PROVENANCE.md), and NCBI does move and rename taxa: if a refresh puts one of
/// these organisms in a different family, that has to be a failing test and a decision, not a silent change in
/// every result grouped by taxon afterwards.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestPinnedTaxonomy
{
    private const long Human = 9606, Mouse = 10090, YeastS288c = 559292, EColiK12 = 83333, SarsCov2 = 2697049;

    private static string NameAt(long taxonId, string rank) =>
        TaxonomyTree.LoadPinned().AncestorAtRank(taxonId, rank).Taxon?.ScientificName;

    [Test]
    public void TheSnapshotSaysWhichReleaseItIs()
    {
        TaxonomyTree tree = TaxonomyTree.LoadPinned();
        ProteomeTaxonMap proteomes = ProteomeTaxonMap.LoadPinned();

        Assert.That(tree.Provenance.Source, Is.EqualTo(UniProtTaxonomyClient.ProteomeSearchUrl));
        Assert.That(tree.Provenance.Query, Is.EqualTo(UniProtTaxonomyClient.ReferenceProteomesQuery));
        Assert.That(tree.Provenance.Release, Does.Match(@"^\d{4}_\d{2}$"));
        Assert.That(proteomes.Provenance, Is.EqualTo(tree.Provenance), "the two tables come from one run");
        Assert.That(TaxonomyTree.LoadPinned(), Is.SameAs(tree), "read once");
    }

    [Test]
    public void ItHoldsTheReferenceProteomesAndTheTreeAboveThem()
    {
        TaxonomyTree tree = TaxonomyTree.LoadPinned();
        ProteomeTaxonMap proteomes = ProteomeTaxonMap.LoadPinned();

        Assert.That(proteomes.Count, Is.GreaterThan(30000));
        Assert.That(tree.Count, Is.GreaterThan(60000));
        Assert.That(tree.Roots.Select(r => r.ScientificName), Is.EquivalentTo(new[] { "cellular organisms", "Viruses" }));
        Assert.That(tree.Nodes.Select(n => n.ScientificName), Has.None.Empty);
    }

    /// <summary>
    /// A handful of reference proteomes sit under a taxon UniProt names by number only (see PROVENANCE.md). Those
    /// are the only proteomes the tree cannot place, and there must stay few of them.
    /// </summary>
    [Test]
    public void NearlyEveryProteomeHasAPlaceInTheTree()
    {
        TaxonomyTree tree = TaxonomyTree.LoadPinned();
        var unplaced = ProteomeTaxonMap.LoadPinned().Proteomes.Where(p => !tree.Contains(p.TaxonId)).ToList();

        Assert.That(unplaced.Count, Is.LessThan(100));
        Assert.That(unplaced.Select(p => p.Organism), Is.All.Empty, "an unplaced proteome is one UniProt gave no organism name");
    }

    [Test]
    public void HumanLineage()
    {
        Assert.That(ProteomeTaxonMap.LoadPinned().TryGetTaxonId("UP000005640", out long taxon), Is.True);
        Assert.That(taxon, Is.EqualTo(Human));

        Assert.That(NameAt(Human, TaxonomicRank.Species), Is.EqualTo("Homo sapiens"));
        Assert.That(NameAt(Human, TaxonomicRank.Genus), Is.EqualTo("Homo"));
        Assert.That(NameAt(Human, TaxonomicRank.Family), Is.EqualTo("Hominidae"));
        Assert.That(NameAt(Human, TaxonomicRank.Order), Is.EqualTo("Primates"));
        Assert.That(NameAt(Human, TaxonomicRank.Class), Is.EqualTo("Mammalia"));
        Assert.That(NameAt(Human, TaxonomicRank.Phylum), Is.EqualTo("Chordata"));
        Assert.That(NameAt(Human, TaxonomicRank.Kingdom), Is.EqualTo("Metazoa"));
        Assert.That(NameAt(Human, TaxonomicRank.Domain), Is.EqualTo("Eukaryota"));
    }

    [Test]
    public void MouseSharesAClassWithHumanButNotAnOrder()
    {
        Assert.That(ProteomeTaxonMap.LoadPinned().TryGetTaxonId("UP000000589", out long taxon), Is.True);
        Assert.That(taxon, Is.EqualTo(Mouse));

        Assert.That(NameAt(Mouse, TaxonomicRank.Species), Is.EqualTo("Mus musculus"));
        Assert.That(NameAt(Mouse, TaxonomicRank.Family), Is.EqualTo("Muridae"));
        Assert.That(NameAt(Mouse, TaxonomicRank.Order), Is.EqualTo("Rodentia"));
        Assert.That(NameAt(Mouse, TaxonomicRank.Class), Is.EqualTo(NameAt(Human, TaxonomicRank.Class)));

        var byOrder = TaxonomyTree.LoadPinned().GroupByRank(new long?[] { Human, Mouse }, id => id, TaxonomicRank.Order);
        var byClass = TaxonomyTree.LoadPinned().GroupByRank(new long?[] { Human, Mouse }, id => id, TaxonomicRank.Class);
        Assert.That(byOrder, Has.Count.EqualTo(2));
        Assert.That(byClass, Has.Count.EqualTo(1));
    }

    [Test]
    public void YeastLineage()
    {
        Assert.That(ProteomeTaxonMap.LoadPinned().TryGetTaxonId("UP000002311", out long taxon), Is.True);
        Assert.That(taxon, Is.EqualTo(YeastS288c));

        Assert.That(NameAt(YeastS288c, TaxonomicRank.Species), Is.EqualTo("Saccharomyces cerevisiae"));
        Assert.That(NameAt(YeastS288c, TaxonomicRank.Family), Is.EqualTo("Saccharomycetaceae"));
        Assert.That(NameAt(YeastS288c, TaxonomicRank.Kingdom), Is.EqualTo("Fungi"));
    }

    [Test]
    public void EColiK12Lineage()
    {
        Assert.That(ProteomeTaxonMap.LoadPinned().TryGetTaxonId("UP000000625", out long taxon), Is.True);
        Assert.That(taxon, Is.EqualTo(EColiK12));

        Assert.That(NameAt(EColiK12, TaxonomicRank.Species), Is.EqualTo("Escherichia coli"));
        Assert.That(NameAt(EColiK12, TaxonomicRank.Genus), Is.EqualTo("Escherichia"));
        Assert.That(NameAt(EColiK12, TaxonomicRank.Family), Is.EqualTo("Enterobacteriaceae"));
        Assert.That(NameAt(EColiK12, TaxonomicRank.Order), Is.EqualTo("Enterobacterales"));
        Assert.That(NameAt(EColiK12, TaxonomicRank.Class), Is.EqualTo("Gammaproteobacteria"));
        Assert.That(NameAt(EColiK12, TaxonomicRank.Phylum), Is.EqualTo("Pseudomonadota"));
        Assert.That(NameAt(EColiK12, TaxonomicRank.Domain), Is.EqualTo("Bacteria"));
    }

    [Test]
    public void SarsCov2Lineage()
    {
        Assert.That(ProteomeTaxonMap.LoadPinned().TryGetTaxonId("UP000464024", out long taxon), Is.True);
        Assert.That(taxon, Is.EqualTo(SarsCov2));

        Assert.That(NameAt(SarsCov2, TaxonomicRank.Species), Is.EqualTo("Betacoronavirus pandemicum"));
        Assert.That(NameAt(SarsCov2, TaxonomicRank.Genus), Is.EqualTo("Betacoronavirus"));
        Assert.That(NameAt(SarsCov2, TaxonomicRank.Family), Is.EqualTo("Coronaviridae"));
        Assert.That(NameAt(SarsCov2, TaxonomicRank.Order), Is.EqualTo("Nidovirales"));
        Assert.That(TaxonomyTree.LoadPinned().AncestorAtRank(SarsCov2, TaxonomicRank.Domain).Status,
            Is.EqualTo(TaxonRankStatus.RankNotAssigned), "viruses belong to no domain");
        Assert.That(TaxonomyTree.LoadPinned().AncestorAtRank(SarsCov2, "realm").Taxon!.ScientificName, Is.EqualTo("Riboviria"));
    }

    [Test]
    public void TopLevelGroupsAreTheDomainsAndViruses()
    {
        TaxonomyTree tree = TaxonomyTree.LoadPinned();

        Assert.That(tree.TopLevelGroup(Human).Taxon!.ScientificName, Is.EqualTo("Eukaryota"));
        Assert.That(tree.TopLevelGroup(EColiK12).Taxon!.ScientificName, Is.EqualTo("Bacteria"));
        Assert.That(tree.TopLevelGroup(SarsCov2).Taxon!.ScientificName, Is.EqualTo("Viruses"));

        var groups = ProteomeTaxonMap.LoadPinned().Proteomes
            .Where(p => tree.Contains(p.TaxonId))
            .GroupBy(p => tree.TopLevelGroup(p.TaxonId).Taxon!.ScientificName)
            .ToDictionary(g => g.Key, g => g.Count());
        Assert.That(groups.Keys, Is.EquivalentTo(new[] { "Bacteria", "Archaea", "Eukaryota", "Viruses" }));
        Assert.That(groups["Bacteria"], Is.GreaterThan(10000));
        Assert.That(groups["Viruses"], Is.GreaterThan(10000));
    }

    [Test]
    public void HumanAndMouseMeetAboveTheirOrders()
    {
        TaxonomyTree tree = TaxonomyTree.LoadPinned();

        TaxonomyNode shared = tree.LowestCommonAncestor(Human, Mouse).Taxon!;

        Assert.That(shared.ScientificName, Is.EqualTo("Euarchontoglires"));
        Assert.That(tree.LowestCommonAncestor(Human, SarsCov2).Status, Is.EqualTo(TaxonRankStatus.RankNotAssigned));
    }

    /// <summary>The reason the status exists: thousands of reference viral proteomes have a family and no order.</summary>
    [Test]
    public void ManyVirusesHaveNoOrder()
    {
        TaxonomyTree tree = TaxonomyTree.LoadPinned();
        var placed = ProteomeTaxonMap.LoadPinned().Proteomes.Where(p => tree.Contains(p.TaxonId)).ToList();

        var byOrder = tree.GroupByRank(placed, p => (long?)p.TaxonId, TaxonomicRank.Order);
        TaxonGroup<ProteomeTaxon> noOrder = byOrder.Single(g => g.Status == TaxonRankStatus.RankNotAssigned);

        Assert.That(noOrder.Items.Count, Is.GreaterThan(1000));
        Assert.That(byOrder.Sum(g => g.Items.Count), Is.EqualTo(placed.Count), "no proteome is dropped by the grouping");
    }
}
