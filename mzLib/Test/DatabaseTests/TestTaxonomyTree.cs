using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases.Taxonomy;

namespace Test.DatabaseTests;

/// <summary>
/// <see cref="TaxonomyTree"/> on a small hand-built tree: two enterobacteria that share a family but not a genus, and
/// a virus whose family has no order above it. The bacterial taxa carry their real NCBI identifiers and ranks; the
/// virus branch (identifiers 900001 to 900003) is invented, because the test needs a lineage with a rank missing.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestTaxonomyTree
{
    private const long CellularOrganisms = 131567, Bacteria = 2, Pseudomonadota = 1224, Gammaproteobacteria = 1236,
        Enterobacterales = 91347, Enterobacteriaceae = 543, Escherichia = 561, EColi = 562, EColiK12 = 83333,
        Salmonella = 590, SalmonellaEnterica = 28901, Viruses = 900001, OrphanFamily = 900002, OrphanVirus = 900003;

    private static readonly TaxonomyProvenance Provenance = new("hand-built for TestTaxonomyTree");

    internal static List<TaxonomyNode> Nodes() => new()
    {
        new(CellularOrganisms, null, "no rank", "cellular organisms", null, true),
        new(Bacteria, CellularOrganisms, "domain", "Bacteria", "bacteria", false),
        new(Pseudomonadota, Bacteria, "phylum", "Pseudomonadota", "proteobacteria", false),
        new(Gammaproteobacteria, Pseudomonadota, "class", "Gammaproteobacteria", null, false),
        new(Enterobacterales, Gammaproteobacteria, "order", "Enterobacterales", null, false),
        new(Enterobacteriaceae, Enterobacterales, "family", "Enterobacteriaceae", "enterobacteria", false),
        new(Escherichia, Enterobacteriaceae, "genus", "Escherichia", null, false),
        new(EColi, Escherichia, "species", "Escherichia coli", "E. coli", true),
        new(EColiK12, EColi, "strain", "Escherichia coli (strain K12)", null, false),
        new(Salmonella, Enterobacteriaceae, "genus", "Salmonella", null, false),
        new(SalmonellaEnterica, Salmonella, "species", "Salmonella enterica", null, false),
        new(Viruses, null, "acellular root", "Viruses", null, false),
        new(OrphanFamily, Viruses, "family", "Orphanviridae", null, false),
        new(OrphanVirus, OrphanFamily, "species", "Orphan virus", null, false),
    };

    internal static TaxonomyTree Tree() => TaxonomyTree.FromNodes(Nodes(), Provenance);

    [Test]
    public void HoldsEveryTaxonInIdentifierOrder()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.Count, Is.EqualTo(14));
        Assert.That(tree.Nodes.Select(n => n.TaxonId), Is.Ordered.Ascending);
        Assert.That(tree.Roots.Select(n => n.TaxonId), Is.EqualTo(new[] { CellularOrganisms, Viruses }));
        Assert.That(tree.Provenance, Is.EqualTo(Provenance));
        Assert.That(tree.Contains(EColi), Is.True);
        Assert.That(tree.Contains(1), Is.False);
        Assert.That(tree.TryGetNode(Escherichia, out var genus), Is.True);
        Assert.That(genus.ScientificName, Is.EqualTo("Escherichia"));
        Assert.That(tree.TryGetNode(1, out _), Is.False);
    }

    [Test]
    public void Ancestors_AreNearestFirstAndEndAtARoot()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.Ancestors(EColiK12).Select(n => n.TaxonId), Is.EqualTo(new[]
        {
            EColi, Escherichia, Enterobacteriaceae, Enterobacterales, Gammaproteobacteria, Pseudomonadota, Bacteria, CellularOrganisms
        }));
        Assert.That(tree.Ancestors(CellularOrganisms), Is.Empty, "a root has none");
        Assert.That(tree.Ancestors(1), Is.Empty, "nor does a taxon the tree does not hold");
    }

    [Test]
    public void Children_AreTheTaxaDirectlyBeneath()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.Children(Enterobacteriaceae).Select(n => n.TaxonId), Is.EqualTo(new[] { Escherichia, Salmonella }));
        Assert.That(tree.Children(EColiK12), Is.Empty);
        Assert.That(tree.Children(1), Is.Empty);
    }

    [Test]
    public void AncestorAtRank_FindsTheTaxonItselfOrItsNearestAncestorAtThatRank()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.AncestorAtRank(EColiK12, TaxonomicRank.Species).Taxon!.TaxonId, Is.EqualTo(EColi));
        Assert.That(tree.AncestorAtRank(EColi, TaxonomicRank.Species).Taxon!.TaxonId, Is.EqualTo(EColi),
            "a species is its own species");
        Assert.That(tree.AncestorAtRank(EColiK12, TaxonomicRank.Family).Taxon!.TaxonId, Is.EqualTo(Enterobacteriaceae));
        Assert.That(tree.AncestorAtRank(SalmonellaEnterica, TaxonomicRank.Family).Taxon!.TaxonId, Is.EqualTo(Enterobacteriaceae));
        Assert.That(tree.AncestorAtRank(EColiK12, TaxonomicRank.Family).Status, Is.EqualTo(TaxonRankStatus.Found));
    }

    /// <summary>A hidden taxon is still an ancestor: E. coli is marked hidden in UniProt's lineage of strain K12.</summary>
    [Test]
    public void AncestorAtRank_SeesHiddenTaxa()
    {
        TaxonomyNode species = Tree().AncestorAtRank(EColiK12, TaxonomicRank.Species).Taxon!;

        Assert.That(species.IsHidden, Is.True);
        Assert.That(species.ScientificName, Is.EqualTo("Escherichia coli"));
    }

    /// <summary>NCBI renamed the top level from "superkingdom" to "domain"; a query for either finds a node ranked as the other.</summary>
    [Test]
    [TestCase("domain")]
    [TestCase("superkingdom")]
    [TestCase("Domain")]
    [TestCase(" domain ")]
    public void AncestorAtRank_TreatsSuperkingdomAndDomainAsOneLevel(string rank)
    {
        Assert.That(Tree().AncestorAtRank(EColiK12, rank).Taxon!.TaxonId, Is.EqualTo(Bacteria));
    }

    /// <summary>"This taxon has no order" and "this taxon is not in the tree" are different answers.</summary>
    [Test]
    public void AncestorAtRank_TellsAMissingRankFromAMissingTaxon()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.AncestorAtRank(OrphanVirus, TaxonomicRank.Order),
            Is.EqualTo(new TaxonRankLookup(TaxonRankStatus.RankNotAssigned, null)));
        Assert.That(tree.AncestorAtRank(EColi, TaxonomicRank.Species).Status, Is.EqualTo(TaxonRankStatus.Found));
        Assert.That(tree.AncestorAtRank(Escherichia, TaxonomicRank.Species),
            Is.EqualTo(new TaxonRankLookup(TaxonRankStatus.RankNotAssigned, null)), "a genus has no species");
        Assert.That(tree.AncestorAtRank(1, TaxonomicRank.Order),
            Is.EqualTo(new TaxonRankLookup(TaxonRankStatus.TaxonNotInTree, null)));
    }

    [Test]
    [TestCase(null)]
    [TestCase("")]
    [TestCase("no rank")]
    [TestCase("clade")]
    [TestCase("flavour")]
    public void AncestorAtRank_RefusesARankWithNoPlaceInTheHierarchy(string rank)
    {
        Assert.Throws<ArgumentException>(() => Tree().AncestorAtRank(EColi, rank));
    }

    /// <summary>Viruses have no domain, so no single rank splits a search's databases into bacteria, eukaryotes and viruses.</summary>
    [Test]
    public void TopLevelGroup_IsTheDomainOrFailingThatTheRoot()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.TopLevelGroup(EColiK12).Taxon!.ScientificName, Is.EqualTo("Bacteria"));
        Assert.That(tree.TopLevelGroup(Bacteria).Taxon!.ScientificName, Is.EqualTo("Bacteria"));
        Assert.That(tree.TopLevelGroup(OrphanVirus).Taxon!.ScientificName, Is.EqualTo("Viruses"));
        Assert.That(tree.TopLevelGroup(CellularOrganisms).Taxon!.TaxonId, Is.EqualTo(CellularOrganisms), "a root above the domain is its own group");
        Assert.That(tree.TopLevelGroup(1), Is.EqualTo(new TaxonRankLookup(TaxonRankStatus.TaxonNotInTree, null)));
    }

    [Test]
    public void LowestCommonAncestor_IsTheNearestTaxonTwoOrganismsShare()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.LowestCommonAncestor(EColiK12, SalmonellaEnterica).Taxon!.TaxonId, Is.EqualTo(Enterobacteriaceae));
        Assert.That(tree.LowestCommonAncestor(SalmonellaEnterica, EColiK12).Taxon!.TaxonId, Is.EqualTo(Enterobacteriaceae),
            "the order of the two does not matter");
        Assert.That(tree.LowestCommonAncestor(EColiK12, EColi).Taxon!.TaxonId, Is.EqualTo(EColi), "a taxon counts as its own ancestor");
        Assert.That(tree.LowestCommonAncestor(EColi, EColi).Taxon!.TaxonId, Is.EqualTo(EColi));
        Assert.That(tree.LowestCommonAncestor(EColiK12, Escherichia).Taxon!.TaxonId, Is.EqualTo(Escherichia));
    }

    [Test]
    public void LowestCommonAncestor_TellsNothingSharedFromATaxonThatIsNotThere()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.LowestCommonAncestor(EColi, OrphanVirus),
            Is.EqualTo(new TaxonRankLookup(TaxonRankStatus.RankNotAssigned, null)), "a bacterium and a virus are under different roots");
        Assert.That(tree.LowestCommonAncestor(EColi, 1), Is.EqualTo(new TaxonRankLookup(TaxonRankStatus.TaxonNotInTree, null)));
        Assert.That(tree.LowestCommonAncestor(1, EColi), Is.EqualTo(new TaxonRankLookup(TaxonRankStatus.TaxonNotInTree, null)));
    }

    [Test]
    public void TryGetLineage_GivesTheEightMajorRanks()
    {
        TaxonomyTree tree = Tree();

        Assert.That(tree.TryGetLineage(EColiK12, out TaxonomyLineage lineage), Is.True);
        Assert.That(lineage.Taxon.TaxonId, Is.EqualTo(EColiK12));
        Assert.That(lineage.Domain!.ScientificName, Is.EqualTo("Bacteria"));
        Assert.That(lineage.Kingdom, Is.Null, "this tree gives the bacteria no kingdom");
        Assert.That(lineage.Phylum!.ScientificName, Is.EqualTo("Pseudomonadota"));
        Assert.That(lineage.Class!.ScientificName, Is.EqualTo("Gammaproteobacteria"));
        Assert.That(lineage.Order!.ScientificName, Is.EqualTo("Enterobacterales"));
        Assert.That(lineage.Family!.ScientificName, Is.EqualTo("Enterobacteriaceae"));
        Assert.That(lineage.Genus!.ScientificName, Is.EqualTo("Escherichia"));
        Assert.That(lineage.Species!.ScientificName, Is.EqualTo("Escherichia coli"));

        Assert.That(tree.TryGetLineage(OrphanVirus, out TaxonomyLineage virus), Is.True);
        Assert.That(virus.Order, Is.Null);
        Assert.That(virus.Family!.ScientificName, Is.EqualTo("Orphanviridae"));

        Assert.That(tree.TryGetLineage(1, out _), Is.False);
    }

    [Test]
    public void GroupByRank_PutsEveryItemInExactlyOneGroup()
    {
        var databases = new (string Name, long? Taxon)[]
        {
            ("k12", EColiK12), ("salmonella", SalmonellaEnterica), ("ecoli", EColi), ("orphan", OrphanVirus),
            ("unknown", 1), ("unstated", null), ("also unstated", null),
        };

        IReadOnlyList<TaxonGroup<(string Name, long? Taxon)>> byGenus =
            Tree().GroupByRank(databases, d => d.Taxon, TaxonomicRank.Genus);

        Assert.That(byGenus.Select(g => (g.Status, g.Taxon?.ScientificName, string.Join(",", g.Items.Select(i => i.Name)))),
            Is.EqualTo(new[]
            {
                (TaxonRankStatus.Found, "Escherichia", "k12,ecoli"),
                (TaxonRankStatus.Found, "Salmonella", "salmonella"),
                (TaxonRankStatus.RankNotAssigned, (string)null, "orphan"),
                (TaxonRankStatus.TaxonNotInTree, (string)null, "unknown"),
                (TaxonRankStatus.TaxonNotStated, (string)null, "unstated,also unstated"),
            }));
        Assert.That(byGenus.Sum(g => g.Items.Count), Is.EqualTo(databases.Length));
    }

    [Test]
    public void GroupByRank_AtFamilyBringsTheTwoGeneraTogether()
    {
        var byFamily = Tree().GroupByRank(new long?[] { EColiK12, SalmonellaEnterica, OrphanVirus }, id => id, TaxonomicRank.Family);

        Assert.That(byFamily.Select(g => (g.Taxon!.ScientificName, g.Items.Count)),
            Is.EqualTo(new[] { ("Enterobacteriaceae", 2), ("Orphanviridae", 1) }));
    }

    [Test]
    public void GroupByRank_RefusesARankWithNoPlaceInTheHierarchy()
    {
        Assert.Throws<ArgumentException>(() => Tree().GroupByRank(new long?[] { EColi }, id => id, "clade"));
    }

    // ---- building refuses rather than repairs --------------------------------

    [Test]
    public void FromNodes_RefusesARepeatedIdentifier()
    {
        var nodes = Nodes();
        nodes.Add(new TaxonomyNode(EColi, Escherichia, "species", "Escherichia coli again", null, false));

        var ex = Assert.Throws<InvalidDataException>(() => TaxonomyTree.FromNodes(nodes, Provenance));
        Assert.That(ex.Message, Does.Contain("562").And.Contain("twice"));
    }

    [Test]
    public void FromNodes_RefusesAParentThatIsNotThere()
    {
        var nodes = Nodes();
        nodes.Add(new TaxonomyNode(777, 778, "species", "Lost species", null, false));

        var ex = Assert.Throws<InvalidDataException>(() => TaxonomyTree.FromNodes(nodes, Provenance));
        Assert.That(ex.Message, Does.Contain("777").And.Contain("778"));
    }

    [Test]
    public void FromNodes_RefusesACycle()
    {
        var nodes = new List<TaxonomyNode>
        {
            new(1, 3, "genus", "A", null, false),
            new(2, 1, "species", "B", null, false),
            new(3, 2, "strain", "C", null, false),
            new(4, null, "domain", "Fine", null, false),
        };

        var ex = Assert.Throws<InvalidDataException>(() => TaxonomyTree.FromNodes(nodes, Provenance));
        Assert.That(ex.Message, Does.Contain("cycle"));
    }

    [Test]
    public void FromNodes_AcceptsAnEmptyTree()
    {
        TaxonomyTree tree = TaxonomyTree.FromNodes(Array.Empty<TaxonomyNode>(), Provenance);

        Assert.That(tree.Count, Is.Zero);
        Assert.That(tree.Roots, Is.Empty);
    }

    // ---- overlays ------------------------------------------------------------

    [Test]
    public void WithOverlay_AddsTaxaBeneathExistingOnesAndBeneathEachOther()
    {
        TaxonomyTree tree = Tree();
        var additions = new[]
        {
            new TaxonomyNode(9000002, 9000001, "strain", "In-house isolate 7", null, false),
            new TaxonomyNode(9000001, Escherichia, "species", "Escherichia novel", null, false),
        };

        TaxonomyTree overlaid = tree.WithOverlay(additions);

        Assert.That(overlaid.Count, Is.EqualTo(tree.Count + 2));
        Assert.That(overlaid.AncestorAtRank(9000002, TaxonomicRank.Family).Taxon!.TaxonId, Is.EqualTo(Enterobacteriaceae));
        Assert.That(overlaid.AncestorAtRank(9000002, TaxonomicRank.Species).Taxon!.TaxonId, Is.EqualTo(9000001));
        Assert.That(overlaid.Provenance, Is.EqualTo(tree.Provenance));
        Assert.That(tree.Contains(9000001), Is.False, "the original tree is unchanged");
    }

    [Test]
    public void WithOverlay_AcceptsATaxonIdenticalToOneAlreadyThere()
    {
        TaxonomyTree tree = Tree();

        TaxonomyTree overlaid = tree.WithOverlay(new[] { Nodes().Single(n => n.TaxonId == EColi) });

        Assert.That(overlaid.Count, Is.EqualTo(tree.Count));
    }

    /// <summary>Adding is allowed; redefining is not: a replaced parent would silently move every organism beneath it.</summary>
    [Test]
    public void WithOverlay_RefusesToRedefineATaxon()
    {
        var moved = new TaxonomyNode(Escherichia, Enterobacterales, "genus", "Escherichia", null, false);

        var ex = Assert.Throws<InvalidDataException>(() => Tree().WithOverlay(new[] { moved }));
        Assert.That(ex.Message, Does.Contain("561").And.Contain("redefine"));
    }

    [Test]
    public void WithOverlay_RefusesARepeatAndADanglingParent()
    {
        var twice = new TaxonomyNode(9000001, Escherichia, "species", "Escherichia novel", null, false);
        Assert.Throws<InvalidDataException>(() => Tree().WithOverlay(new[] { twice, twice }));

        var dangling = new TaxonomyNode(9000001, 424242, "species", "Escherichia novel", null, false);
        var ex = Assert.Throws<InvalidDataException>(() => Tree().WithOverlay(new[] { dangling }));
        Assert.That(ex.Message, Does.Contain("424242"));
    }

    // ---- identifiers and ranks -----------------------------------------------

    [Test]
    [TestCase("9606", true, 9606)]
    [TestCase(" 9606 ", true, 9606)]
    [TestCase(null, false, 0)]
    [TestCase("", false, 0)]
    [TestCase("0", false, 0)]
    [TestCase("-5", false, 0)]
    [TestCase("9606.0", false, 0)]
    [TestCase("NCBITaxon:9606", false, 0)]
    public void TryParseTaxonId_ReadsAPositiveWholeNumber(string text, bool expected, long expectedId)
    {
        Assert.That(TaxonomyTree.TryParseTaxonId(text, out long taxonId), Is.EqualTo(expected));
        if (expected)
        {
            Assert.That(taxonId, Is.EqualTo(expectedId));
        }
    }

    [Test]
    public void TaxonomicRank_OrdersRanksFromTheRootDown()
    {
        Assert.That(TaxonomicRank.MajorRanks, Is.EqualTo(new[]
        {
            "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"
        }));
        for (int i = 1; i < TaxonomicRank.MajorRanks.Count; i++)
        {
            Assert.That(TaxonomicRank.Compare(TaxonomicRank.MajorRanks[i - 1], TaxonomicRank.MajorRanks[i]), Is.LessThan(0));
        }
        Assert.That(TaxonomicRank.Compare("species", "family"), Is.GreaterThan(0));
        Assert.That(TaxonomicRank.Compare("Superkingdom", "domain"), Is.Zero);
        Assert.That(TaxonomicRank.Compare("clade", "family"), Is.Null, "an unranked taxon can sit at any depth");
        Assert.That(TaxonomicRank.Compare(null, "family"), Is.Null);
    }

    [Test]
    public void TaxonomicRank_KnowsWhichNamesAreLevels()
    {
        Assert.That(TaxonomicRank.IsRanked("family"), Is.True);
        Assert.That(TaxonomicRank.IsRanked("Subfamily"), Is.True);
        Assert.That(TaxonomicRank.IsRanked("no rank"), Is.False);
        Assert.That(TaxonomicRank.IsRanked("clade"), Is.False);
        Assert.That(TaxonomicRank.IsRanked(""), Is.False);
        Assert.That(TaxonomicRank.IsRanked(null), Is.False);

        Assert.That(TaxonomicRank.AreSameRank("family", "FAMILY"), Is.True);
        Assert.That(TaxonomicRank.AreSameRank("family", "subfamily"), Is.False);
        Assert.That(TaxonomicRank.AreSameRank("clade", "clade"), Is.False, "two unranked taxa are not at one level");
        Assert.That(TaxonomicRank.AreSameRank("strain", "serotype"), Is.False,
            "both sit below species, but they are different kinds of group");
        Assert.That(TaxonomicRank.AreSameRank("strain", "strain"), Is.True);
    }
}
