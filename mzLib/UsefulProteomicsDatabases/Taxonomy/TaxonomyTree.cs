#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases.Taxonomy
{
    /// <summary>What <see cref="TaxonomyTree.AncestorAtRank"/> found.</summary>
    public enum TaxonRankStatus
    {
        /// <summary>The taxon, or one of its ancestors, has the rank asked for.</summary>
        Found,

        /// <summary>
        /// The taxon is in the tree but neither it nor any ancestor has the rank asked for. This is a fact about the
        /// taxonomy, not a gap in the data: many viruses belong to no order, for example.
        /// </summary>
        RankNotAssigned,

        /// <summary>The tree has no taxon with that identifier, so nothing can be said about its lineage.</summary>
        TaxonNotInTree,

        /// <summary>No taxon identifier was given at all. Only <see cref="TaxonomyTree.GroupByRank{T}"/> reports this.</summary>
        TaxonNotStated,
    }

    /// <summary>The answer to "which taxon at this rank does that taxon belong to?".</summary>
    /// <param name="Status">Whether a taxon was found and, when not, why not.</param>
    /// <param name="Taxon">The taxon at the rank asked for; null unless <paramref name="Status"/> is <see cref="TaxonRankStatus.Found"/>.</param>
    public sealed record TaxonRankLookup(TaxonRankStatus Status, TaxonomyNode? Taxon);

    /// <summary>The items <see cref="TaxonomyTree.GroupByRank{T}"/> placed under one taxon, or under one reason for having none.</summary>
    /// <param name="Status">Why these items are together: they share a taxon, or they share a reason for lacking one.</param>
    /// <param name="Taxon">The shared taxon; null unless <paramref name="Status"/> is <see cref="TaxonRankStatus.Found"/>.</param>
    /// <param name="Items">The items, in the order they were given.</param>
    public sealed record TaxonGroup<T>(TaxonRankStatus Status, TaxonomyNode? Taxon, IReadOnlyList<T> Items);

    /// <summary>
    /// A taxonomy held as a tree: every taxon has an identifier, a rank and at most one parent. It answers what a flat
    /// table of kingdom-to-species names cannot: which family an organism belongs to, whether two organisms share a
    /// genus, and what sits beneath a taxon.
    /// </summary>
    /// <remarks>
    /// <para>
    /// A tree is immutable once built. It need not hold the whole NCBI taxonomy, and the ones mzLib builds do not: they
    /// hold the organisms of a set of UniProt proteomes and every ancestor of each, which is all that is needed to roll
    /// those organisms up to any rank. See <see cref="UniProtTaxonomyClient"/> for how one is built and
    /// <see cref="LoadPinned"/> for the one shipped with mzLib.
    /// </para>
    /// <para>
    /// Building refuses rather than repairs: an identifier that appears twice, a parent that is not in the tree, and a
    /// cycle are all errors. A tree that loaded is one whose every ancestor walk terminates at a root.
    /// </para>
    /// </remarks>
    public sealed class TaxonomyTree
    {
        internal const string PinnedNodesResource = "UsefulProteomicsDatabases.Resources.taxonomy_nodes.tsv.gz";

        private static readonly Lazy<TaxonomyTree> Pinned = new(() =>
            TaxonomySnapshotTsv.ReadNodesResource(PinnedNodesResource));

        private readonly Dictionary<long, TaxonomyNode> _nodes;
        private readonly Dictionary<long, List<TaxonomyNode>> _children;
        private readonly List<TaxonomyNode> _ordered;
        private readonly List<TaxonomyNode> _roots;

        private TaxonomyTree(Dictionary<long, TaxonomyNode> nodes, TaxonomyProvenance provenance)
        {
            _nodes = nodes;
            Provenance = provenance;
            _ordered = nodes.Values.OrderBy(n => n.TaxonId).ToList();
            _children = new Dictionary<long, List<TaxonomyNode>>();
            _roots = new List<TaxonomyNode>();
            foreach (TaxonomyNode node in _ordered)
            {
                if (node.ParentTaxonId is long parent)
                {
                    if (!_children.TryGetValue(parent, out var siblings))
                    {
                        _children[parent] = siblings = new List<TaxonomyNode>();
                    }
                    siblings.Add(node);
                }
                else
                {
                    _roots.Add(node);
                }
            }
        }

        /// <summary>
        /// Where the tree's taxa came from. It describes the base tree only: a tree made by
        /// <see cref="WithOverlay(IEnumerable{TaxonomyNode})"/> keeps its base's provenance, so taxa a caller added are
        /// not recorded in it.
        /// </summary>
        public TaxonomyProvenance Provenance { get; }

        /// <summary>How many taxa the tree holds.</summary>
        public int Count => _nodes.Count;

        /// <summary>Every taxon, in ascending order of identifier.</summary>
        public IReadOnlyList<TaxonomyNode> Nodes => _ordered;

        /// <summary>The taxa with no parent, in ascending order of identifier.</summary>
        public IReadOnlyList<TaxonomyNode> Roots => _roots;

        /// <summary>
        /// Builds a tree from its taxa.
        /// </summary>
        /// <exception cref="InvalidDataException">A taxon has no name, an identifier appears twice, a parent is not
        /// among the taxa, or the parent links form a cycle.</exception>
        public static TaxonomyTree FromNodes(IEnumerable<TaxonomyNode> nodes, TaxonomyProvenance provenance)
        {
            ArgumentNullException.ThrowIfNull(nodes);
            ArgumentNullException.ThrowIfNull(provenance);

            var byId = new Dictionary<long, TaxonomyNode>();
            foreach (TaxonomyNode node in nodes)
            {
                if (!byId.TryAdd(node.TaxonId, node))
                {
                    throw new InvalidDataException($"Taxon {Id(node.TaxonId)} appears twice.");
                }
            }
            Validate(byId);
            return new TaxonomyTree(byId, provenance);
        }

        /// <summary>
        /// Reads a tree from a node table written by <see cref="TaxonomySnapshotTsv.WriteNodes(TextWriter, TaxonomyTree)"/>.
        /// A path ending in ".gz" is read as gzip.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file is not a node table, or its taxa do not form a tree.</exception>
        public static TaxonomyTree Load(string nodesPath) => TaxonomySnapshotTsv.ReadNodes(nodesPath);

        /// <summary>
        /// The tree shipped inside mzLib: the organisms of UniProt's reference proteomes and all their ancestors, as
        /// of the release recorded in <see cref="Provenance"/>. Read once, on first use.
        /// </summary>
        /// <remarks>
        /// It is pinned, not downloaded, so that the same search annotated a year apart places an organism in the same
        /// family. See Resources/PROVENANCE.md for its source and how it is refreshed. A taxon it lacks can be added
        /// with <see cref="WithOverlay(string)"/>; a fuller tree can be built with <see cref="UniProtTaxonomyClient"/>.
        /// </remarks>
        public static TaxonomyTree LoadPinned() => Pinned.Value;

        /// <summary>
        /// A new tree holding this tree's taxa and those in <paramref name="additions"/>. This is how a caller adds
        /// organisms the tree lacks: a strain sequenced in-house, a metagenome bin, a taxon newer than the snapshot.
        /// </summary>
        /// <remarks>
        /// Adding is allowed; redefining is not. An addition whose identifier the tree already has must be identical
        /// to the taxon there, or it is refused: a silently replaced parent would move every organism beneath it.
        /// An addition may name a parent in this tree or among the additions. The new tree keeps this tree's
        /// <see cref="Provenance"/>, which therefore does not describe the additions.
        /// </remarks>
        /// <exception cref="InvalidDataException">An addition redefines an existing taxon, repeats, has no name, names
        /// a parent that is in neither set, or closes a cycle.</exception>
        public TaxonomyTree WithOverlay(IEnumerable<TaxonomyNode> additions)
        {
            ArgumentNullException.ThrowIfNull(additions);

            var merged = new Dictionary<long, TaxonomyNode>(_nodes);
            var added = new HashSet<long>();
            foreach (TaxonomyNode node in additions)
            {
                if (_nodes.TryGetValue(node.TaxonId, out var existing))
                {
                    if (existing != node)
                    {
                        throw new InvalidDataException(
                            $"Taxon {Id(node.TaxonId)} ({node.ScientificName}) is already in the tree as " +
                            $"'{existing.ScientificName}' with a different definition. An overlay may add taxa but not redefine them.");
                    }
                    continue;
                }
                if (!added.Add(node.TaxonId))
                {
                    throw new InvalidDataException($"Taxon {Id(node.TaxonId)} appears twice in the overlay.");
                }
                merged[node.TaxonId] = node;
            }
            Validate(merged);
            return new TaxonomyTree(merged, Provenance);
        }

        /// <summary>
        /// A new tree holding this tree's taxa and those in the node table at <paramref name="nodesPath"/>; see
        /// <see cref="WithOverlay(IEnumerable{TaxonomyNode})"/> for the rules.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file is not a node table, or an addition is refused.</exception>
        public TaxonomyTree WithOverlay(string nodesPath) =>
            WithOverlay(TaxonomySnapshotTsv.ReadNodeRows(nodesPath).Rows);

        /// <summary>Whether the tree holds a taxon with this identifier.</summary>
        public bool Contains(long taxonId) => _nodes.ContainsKey(taxonId);

        /// <summary>Finds a taxon by identifier.</summary>
        public bool TryGetNode(long taxonId, out TaxonomyNode node) => _nodes.TryGetValue(taxonId, out node!);

        /// <summary>
        /// The taxa directly beneath <paramref name="taxonId"/>, in ascending order of identifier. Empty for a taxon
        /// with none, and for an identifier the tree does not hold.
        /// </summary>
        public IReadOnlyList<TaxonomyNode> Children(long taxonId) =>
            _children.TryGetValue(taxonId, out var children) ? children : Array.Empty<TaxonomyNode>();

        /// <summary>
        /// Every taxon above <paramref name="taxonId"/>, nearest first and ending at a root. The taxon itself is not
        /// included. Empty for a root, and for an identifier the tree does not hold.
        /// </summary>
        public IReadOnlyList<TaxonomyNode> Ancestors(long taxonId)
        {
            var ancestors = new List<TaxonomyNode>();
            if (!_nodes.TryGetValue(taxonId, out var node))
            {
                return ancestors;
            }
            while (node.ParentTaxonId is long parent)
            {
                node = _nodes[parent];
                ancestors.Add(node);
            }
            return ancestors;
        }

        /// <summary>
        /// The taxon at <paramref name="rank"/> that <paramref name="taxonId"/> belongs to: the taxon itself when it
        /// has that rank, otherwise its nearest ancestor that does.
        /// </summary>
        /// <param name="taxonId">The taxon to place.</param>
        /// <param name="rank">A rank name such as <see cref="TaxonomicRank.Family"/>; compared with
        /// <see cref="TaxonomicRank.AreSameRank"/>.</param>
        /// <returns>
        /// The taxon found, or a status saying why there is none. "Not in the tree" and "has no such rank" are kept
        /// apart because they call for different handling: the first is missing data, the second is an answer.
        /// </returns>
        /// <exception cref="ArgumentException"><paramref name="rank"/> is not a rank with a place in the hierarchy
        /// (it is null, empty, "no rank", "clade" or unknown).</exception>
        public TaxonRankLookup AncestorAtRank(long taxonId, string rank)
        {
            if (!TaxonomicRank.IsRanked(rank))
            {
                throw new ArgumentException($"'{rank}' is not a rank with a place in the hierarchy.", nameof(rank));
            }
            if (!_nodes.TryGetValue(taxonId, out var node))
            {
                return new TaxonRankLookup(TaxonRankStatus.TaxonNotInTree, null);
            }
            while (true)
            {
                if (TaxonomicRank.AreSameRank(node.Rank, rank))
                {
                    return new TaxonRankLookup(TaxonRankStatus.Found, node);
                }
                if (node.ParentTaxonId is not long parent)
                {
                    return new TaxonRankLookup(TaxonRankStatus.RankNotAssigned, null);
                }
                node = _nodes[parent];
            }
        }

        /// <summary>
        /// The broadest named group <paramref name="taxonId"/> belongs to: its domain (Bacteria, Archaea, Eukaryota)
        /// when it has one, otherwise the root of its lineage, which for a virus is "Viruses". This is the level a
        /// search's databases are usually first split at, and no single rank gives it, because viruses have no domain.
        /// </summary>
        /// <returns>The group, or <see cref="TaxonRankStatus.TaxonNotInTree"/> when the tree does not hold the taxon.
        /// Never <see cref="TaxonRankStatus.RankNotAssigned"/>: every taxon in the tree has a root.</returns>
        public TaxonRankLookup TopLevelGroup(long taxonId)
        {
            if (!_nodes.TryGetValue(taxonId, out var node))
            {
                return new TaxonRankLookup(TaxonRankStatus.TaxonNotInTree, null);
            }
            while (true)
            {
                if (TaxonomicRank.AreSameRank(node.Rank, TaxonomicRank.Domain) || node.ParentTaxonId is not long parent)
                {
                    return new TaxonRankLookup(TaxonRankStatus.Found, node);
                }
                node = _nodes[parent];
            }
        }

        /// <summary>
        /// The nearest taxon that is an ancestor of both <paramref name="firstTaxonId"/> and
        /// <paramref name="secondTaxonId"/>, counting each taxon as an ancestor of itself: the genus two species share,
        /// the family two genera share. How far up it sits says how closely two organisms are related.
        /// </summary>
        /// <returns>
        /// The shared taxon; <see cref="TaxonRankStatus.TaxonNotInTree"/> when the tree lacks either taxon; or
        /// <see cref="TaxonRankStatus.RankNotAssigned"/> when both are in the tree but under different roots, as a
        /// bacterium and a virus are, and so share nothing.
        /// </returns>
        public TaxonRankLookup LowestCommonAncestor(long firstTaxonId, long secondTaxonId)
        {
            if (!_nodes.TryGetValue(firstTaxonId, out var first) || !_nodes.TryGetValue(secondTaxonId, out var second))
            {
                return new TaxonRankLookup(TaxonRankStatus.TaxonNotInTree, null);
            }
            var above = new HashSet<long> { first.TaxonId };
            foreach (TaxonomyNode ancestor in Ancestors(firstTaxonId))
            {
                above.Add(ancestor.TaxonId);
            }
            while (true)
            {
                if (above.Contains(second.TaxonId))
                {
                    return new TaxonRankLookup(TaxonRankStatus.Found, second);
                }
                if (second.ParentTaxonId is not long parent)
                {
                    return new TaxonRankLookup(TaxonRankStatus.RankNotAssigned, null);
                }
                second = _nodes[parent];
            }
        }

        /// <summary>
        /// The taxon's place at each of the eight major ranks, in one object. False when the tree does not hold the
        /// taxon.
        /// </summary>
        public bool TryGetLineage(long taxonId, out TaxonomyLineage lineage)
        {
            if (!_nodes.TryGetValue(taxonId, out var taxon))
            {
                lineage = null!;
                return false;
            }
            TaxonomyNode? At(string rank) => AncestorAtRank(taxonId, rank).Taxon;
            lineage = new TaxonomyLineage(taxon,
                At(TaxonomicRank.Domain), At(TaxonomicRank.Kingdom), At(TaxonomicRank.Phylum), At(TaxonomicRank.Class),
                At(TaxonomicRank.Order), At(TaxonomicRank.Family), At(TaxonomicRank.Genus), At(TaxonomicRank.Species));
            return true;
        }

        /// <summary>
        /// Groups items by the taxon each belongs to at <paramref name="rank"/>: organisms by family, databases by
        /// genus, proteins by phylum.
        /// </summary>
        /// <param name="items">The things to group.</param>
        /// <param name="taxonIdOf">The taxon of an item, or null when it states none.</param>
        /// <param name="rank">The rank to group at, e.g. <see cref="TaxonomicRank.Family"/>.</param>
        /// <returns>
        /// One group per taxon found, in ascending order of identifier, followed by at most one group for each reason
        /// an item could not be placed (<see cref="TaxonRankStatus.RankNotAssigned"/>,
        /// <see cref="TaxonRankStatus.TaxonNotInTree"/>, <see cref="TaxonRankStatus.TaxonNotStated"/>, in that order).
        /// Every item is in exactly one group, so nothing is dropped silently.
        /// </returns>
        /// <exception cref="ArgumentException"><paramref name="rank"/> is not a rank with a place in the hierarchy.</exception>
        public IReadOnlyList<TaxonGroup<T>> GroupByRank<T>(IEnumerable<T> items, Func<T, long?> taxonIdOf, string rank)
        {
            ArgumentNullException.ThrowIfNull(items);
            ArgumentNullException.ThrowIfNull(taxonIdOf);
            if (!TaxonomicRank.IsRanked(rank))
            {
                throw new ArgumentException($"'{rank}' is not a rank with a place in the hierarchy.", nameof(rank));
            }

            var placed = new SortedDictionary<long, (TaxonomyNode Taxon, List<T> Items)>();
            var unplaced = new Dictionary<TaxonRankStatus, List<T>>();
            foreach (T item in items)
            {
                TaxonRankLookup lookup = taxonIdOf(item) is long taxonId
                    ? AncestorAtRank(taxonId, rank)
                    : new TaxonRankLookup(TaxonRankStatus.TaxonNotStated, null);

                if (lookup.Taxon is TaxonomyNode taxon)
                {
                    if (!placed.TryGetValue(taxon.TaxonId, out var group))
                    {
                        placed[taxon.TaxonId] = group = (taxon, new List<T>());
                    }
                    group.Items.Add(item);
                }
                else
                {
                    if (!unplaced.TryGetValue(lookup.Status, out var group))
                    {
                        unplaced[lookup.Status] = group = new List<T>();
                    }
                    group.Add(item);
                }
            }

            var groups = placed.Values
                .Select(g => new TaxonGroup<T>(TaxonRankStatus.Found, g.Taxon, g.Items))
                .ToList();
            foreach (TaxonRankStatus status in new[]
                     {
                         TaxonRankStatus.RankNotAssigned, TaxonRankStatus.TaxonNotInTree, TaxonRankStatus.TaxonNotStated
                     })
            {
                if (unplaced.TryGetValue(status, out var group))
                {
                    groups.Add(new TaxonGroup<T>(status, null, group));
                }
            }
            return groups;
        }

        /// <summary>
        /// Reads a taxon identifier written as text, which is how a protein carries one
        /// (<c>Protein.NcbiTaxonomyId</c>) and how UniProt's proteome catalogue lists one. False for null, blank,
        /// anything that is not a whole number, and zero or less.
        /// </summary>
        public static bool TryParseTaxonId(string? text, out long taxonId)
        {
            return long.TryParse(text?.Trim(), NumberStyles.None, CultureInfo.InvariantCulture, out taxonId) && taxonId > 0;
        }

        /// <summary>
        /// Refuses a taxon with no name (the node table cannot hold one), a parent that is absent and a cycle. Each
        /// taxon is walked once.
        /// </summary>
        private static void Validate(Dictionary<long, TaxonomyNode> nodes)
        {
            foreach (TaxonomyNode node in nodes.Values)
            {
                if (string.IsNullOrEmpty(node.ScientificName))
                {
                    throw new InvalidDataException($"Taxon {Id(node.TaxonId)} has no scientific_name.");
                }
                if (node.ParentTaxonId is long parent && !nodes.ContainsKey(parent))
                {
                    throw new InvalidDataException(
                        $"Taxon {Id(node.TaxonId)} ({node.ScientificName}) names parent {Id(parent)}, which is not among the taxa.");
                }
            }

            // A taxon is settled once its walk to a root has been seen to terminate. A walk that meets a taxon
            // it has itself already passed has found a cycle.
            var settled = new HashSet<long>();
            var path = new HashSet<long>();
            foreach (TaxonomyNode start in nodes.Values)
            {
                path.Clear();
                TaxonomyNode node = start;
                while (!settled.Contains(node.TaxonId))
                {
                    if (!path.Add(node.TaxonId))
                    {
                        throw new InvalidDataException(
                            $"The parent links form a cycle through taxon {Id(node.TaxonId)} ({node.ScientificName}).");
                    }
                    if (node.ParentTaxonId is not long parent)
                    {
                        break;
                    }
                    node = nodes[parent];
                }
                settled.UnionWith(path);
            }
        }

        private static string Id(long taxonId) => taxonId.ToString(CultureInfo.InvariantCulture);
    }
}
