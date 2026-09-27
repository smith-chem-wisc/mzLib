using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>
    /// What one gene has in another species. Every gene gets exactly one, so "the source found nothing",
    /// "the source had nothing to compare against" and "the source never looked" never share an empty
    /// result. Compara does not say whether an absent ortholog is absent from biology or from the
    /// release; these classes say how far its inference got instead.
    /// </summary>
    public enum OrthologyStatus
    {
        /// <summary>At least one ortholog in the target species.</summary>
        HasOrtholog,

        /// <summary>
        /// The gene's tree holds a gene of the target species, and Compara called no ortholog to it. The
        /// source had the chance; the closest the data comes to "no ortholog exists".
        /// </summary>
        NoEdgeInSharedTree,

        /// <summary>The gene's tree holds no gene of the target species, so there was nothing to compare.</summary>
        TreeLacksTargetSpecies,

        /// <summary>In the gene set, but in no gene tree. No call was attempted.</summary>
        NotInAnyTree,

        /// <summary>Not in any of the snapshot's gene sets: for a gene id supplied from elsewhere, e.g.
        /// resolved against a different release.</summary>
        NotInGeneSet
    }

    /// <summary>
    /// One gene per requested species, in the order the species were given, where every pair of them is
    /// an ortholog. <see cref="AllOneToOne"/> is true when every one of those orthologies is one-to-one.
    /// </summary>
    public sealed record OrthologTuple(IReadOnlyList<string> GeneIds, bool AllOneToOne);

    /// <summary>
    /// Compara's homologies among a chosen set of species, with the gene sets and gene trees they are
    /// counted against: the content of an orthology store snapshot, before it is written anywhere.
    ///
    /// The unit is the species PAIR, because Compara asserts relationships only pairwise. Each pair
    /// (a, b) with a &lt;= b in ordinal order holds every relationship between them, one row per
    /// relationship, oriented so side A is species a; the pair (s, s) holds s's paralogs. Nothing is
    /// picked, merged or collapsed, and no group is formed by chaining calls: <see cref="SpeciesSet"/>
    /// requires an ortholog between EVERY two members of a tuple, because Compara's calls are not
    /// transitive.
    /// </summary>
    public sealed class OrthologySnapshot
    {
        private readonly Dictionary<(string, string), List<ComparaHomology>> _pairs;
        private readonly Dictionary<string, string> _speciesOfGene;
        private readonly Dictionary<string, HashSet<string>> _speciesInTree;

        /// <summary>gene -> partner species -> (partner gene, homology_type), both directions, orthologs only.</summary>
        private readonly Dictionary<string, Dictionary<string, List<(string Gene, string Type)>>> _orthologs;

        private OrthologySnapshot(string release, IReadOnlyList<string> species,
            IReadOnlyDictionary<string, EnsemblGeneSet> geneSets, ComparaGeneTreeContent geneTrees,
            IReadOnlyList<ComparaHomologyDump> dumps, Dictionary<(string, string), List<ComparaHomology>> pairs,
            Dictionary<string, string> speciesOfGene, int identicalDuplicates)
        {
            Release = release;
            Species = species;
            GeneSets = geneSets;
            GeneTrees = geneTrees;
            Dumps = dumps;
            _pairs = pairs;
            _speciesOfGene = speciesOfGene;
            IdenticalDuplicatesDropped = identicalDuplicates;

            _speciesInTree = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);
            foreach (var (gene, sp) in speciesOfGene)
            {
                if (geneTrees.TryGetMember(gene, out var m))
                {
                    if (!_speciesInTree.TryGetValue(m.TreeId, out var set))
                    {
                        _speciesInTree[m.TreeId] = set = new HashSet<string>(StringComparer.Ordinal);
                    }
                    set.Add(sp);
                }
            }

            _orthologs = new Dictionary<string, Dictionary<string, List<(string, string)>>>(StringComparer.Ordinal);
            foreach (var row in pairs.Values.SelectMany(r => r).Where(r => r.Class == HomologyClass.Ortholog))
            {
                AddOrtholog(row.GeneA, row.SpeciesB, row.GeneB, row.HomologyType);
                AddOrtholog(row.GeneB, row.SpeciesA, row.GeneA, row.HomologyType);
            }
        }

        /// <summary>The Ensembl release every input belongs to.</summary>
        public string Release { get; }

        /// <summary>The species, Ensembl production names (e.g. "homo_sapiens"), in ordinal order.</summary>
        public IReadOnlyList<string> Species { get; }

        public IReadOnlyDictionary<string, EnsemblGeneSet> GeneSets { get; }

        public ComparaGeneTreeContent GeneTrees { get; }

        /// <summary>One homology dump per species, in species order.</summary>
        public IReadOnlyList<ComparaHomologyDump> Dumps { get; }

        /// <summary>
        /// Rows met in both of a pair's dumps and identical in every value. Zero in release 116, where the
        /// dumps partition; counted so that changing is visible.
        /// </summary>
        public int IdenticalDuplicatesDropped { get; }

        /// <summary>Every pair (a, b), a &lt;= b, including (s, s), in ordinal order.</summary>
        public IEnumerable<(string A, string B)> Pairs =>
            Species.SelectMany((a, i) => Species.Skip(i).Select(b => (a, b)));

        /// <summary>
        /// Every relationship between <paramref name="a"/> and <paramref name="b"/>, side A being the
        /// species that sorts first, in order of homology id. Empty when there is none.
        /// </summary>
        public IReadOnlyList<ComparaHomology> Pair(string a, string b)
        {
            RequireSpecies(a, nameof(a));
            RequireSpecies(b, nameof(b));
            var key = string.CompareOrdinal(a, b) <= 0 ? (a, b) : (b, a);
            return _pairs.TryGetValue(key, out var rows) ? rows : Array.Empty<ComparaHomology>();
        }

        /// <summary>The ortholog rows between two different species, oriented so side A is <paramref name="from"/>.</summary>
        public IEnumerable<ComparaHomology> Orthologs(string from, string to)
        {
            RequireDistinct(from, to);
            return Pair(from, to)
                .Where(r => r.Class == HomologyClass.Ortholog)
                .Select(r => r.SpeciesA == from ? r : r.Swapped());
        }

        /// <summary>The species whose gene set holds <paramref name="geneId"/>, or null.</summary>
        public string SpeciesOf(string geneId) =>
            geneId != null && _speciesOfGene.TryGetValue(geneId, out var sp) ? sp : null;

        /// <summary>The status of one gene in <paramref name="target"/>. Any gene id is accepted.</summary>
        public OrthologyStatus StatusOf(string geneId, string target)
        {
            RequireSpecies(target, nameof(target));
            string species = SpeciesOf(geneId);
            if (species == null)
            {
                return OrthologyStatus.NotInGeneSet;
            }
            RequireDistinct(species, target);

            if (_orthologs.TryGetValue(geneId, out var partners) && partners.ContainsKey(target))
            {
                return OrthologyStatus.HasOrtholog;
            }
            if (!GeneTrees.TryGetMember(geneId, out var member))
            {
                return OrthologyStatus.NotInAnyTree;
            }
            return _speciesInTree[member.TreeId].Contains(target)
                ? OrthologyStatus.NoEdgeInSharedTree
                : OrthologyStatus.TreeLacksTargetSpecies;
        }

        /// <summary>
        /// Every gene of <paramref name="from"/>'s gene set with its status in <paramref name="to"/>, in
        /// ordinal order. All biotypes are included: which genes are the denominator is the caller's
        /// choice, and it is not the same as protein_coding.
        /// </summary>
        public IEnumerable<(EnsemblGene Gene, OrthologyStatus Status)> PairStatus(string from, string to)
        {
            RequireDistinct(from, to);
            return GeneSets[from].Genes.Select(g => (g, StatusOf(g.GeneId, to)));
        }

        /// <summary>
        /// The tuples of one gene per species in which every two members are orthologs, genes in the
        /// order <paramref name="species"/> is given. A tuple is never formed by chaining: A~B and B~C
        /// without A~C is not a tuple. Tuples are in ordinal order of their genes.
        /// </summary>
        public IEnumerable<OrthologTuple> SpeciesSet(params string[] species)
        {
            ArgumentNullException.ThrowIfNull(species);
            if (species.Length < 2 || species.Distinct(StringComparer.Ordinal).Count() != species.Length)
            {
                throw new ArgumentException("At least two distinct species are required.", nameof(species));
            }
            foreach (string s in species)
            {
                RequireSpecies(s, nameof(species));
            }

            return Tuples(species.ToArray());
        }

        private IEnumerable<OrthologTuple> Tuples(string[] species)
        {
            foreach (string first in GeneSets[species[0]].GeneIds)
            {
                var partial = new List<(string[] Genes, bool OneToOne)> { (new[] { first }, true) };
                for (int k = 1; k < species.Length && partial.Count > 0; k++)
                {
                    var next = new List<(string[], bool)>();
                    foreach (var (genes, oneToOne) in partial)
                    {
                        foreach (var (candidate, type) in Partners(genes[0], species[k]))
                        {
                            bool all = oneToOne && type == "ortholog_one2one";
                            bool clique = true;
                            for (int j = 1; j < genes.Length && clique; j++)
                            {
                                string edge = EdgeType(genes[j], species[k], candidate);
                                clique = edge != null;
                                all &= edge == "ortholog_one2one";
                            }
                            if (clique)
                            {
                                next.Add((genes.Append(candidate).ToArray(), all));
                            }
                        }
                    }
                    partial = next;
                }

                foreach (var (genes, oneToOne) in partial.OrderBy(t => string.Join("\t", t.Genes), StringComparer.Ordinal))
                {
                    yield return new OrthologTuple(genes, oneToOne);
                }
            }
        }

        /// <summary>
        /// Builds a snapshot, refusing inputs that would make it wrong rather than incomplete-looking.
        /// </summary>
        /// <param name="release">The Ensembl release. Every input that names its release must name this one.</param>
        /// <param name="geneSets">Species (Ensembl production name) -> its primary-assembly gene set.</param>
        /// <param name="geneTrees">Gene-tree membership, unrestricted or restricted to these same gene sets.</param>
        /// <param name="dumps">Exactly one homology dump per species, each loaded for at least these species.
        /// A pair's relationships sit in either of its two dumps, arbitrarily, so a missing dump would
        /// silently lose them.</param>
        /// <exception cref="ArgumentException">The inputs do not fit together: a release differs, a dump is
        /// missing, repeated or loaded for fewer species, or the gene trees were restricted to other gene sets.</exception>
        /// <exception cref="InvalidDataException">The data breaks the store's rules: a gene in two species'
        /// gene sets, a homology naming a gene its species' gene set does not hold, or one homology id with
        /// two different rows.</exception>
        public static OrthologySnapshot Build(string release, IReadOnlyDictionary<string, EnsemblGeneSet> geneSets,
            ComparaGeneTreeContent geneTrees, IEnumerable<ComparaHomologyDump> dumps)
        {
            ArgumentException.ThrowIfNullOrEmpty(release);
            ArgumentNullException.ThrowIfNull(geneSets);
            ArgumentNullException.ThrowIfNull(geneTrees);
            ArgumentNullException.ThrowIfNull(dumps);
            if (geneSets.Count == 0)
            {
                throw new ArgumentException("At least one species is required.", nameof(geneSets));
            }

            var species = geneSets.Keys.OrderBy(s => s, StringComparer.Ordinal).ToList();
            var speciesOfGene = new Dictionary<string, string>(StringComparer.Ordinal);
            foreach (string s in species)
            {
                var set = geneSets[s] ?? throw new ArgumentException($"The gene set of {s} is null.", nameof(geneSets));
                RequireRelease(set.Release, release, set.SourceFileName);
                foreach (string gene in set.GeneIds)
                {
                    if (!speciesOfGene.TryAdd(gene, s))
                    {
                        throw new InvalidDataException($"{gene} is in the gene sets of both {speciesOfGene[gene]} and {s}.");
                    }
                }
            }

            RequireRelease(geneTrees.Release, release, geneTrees.SourceFileName);
            if (geneTrees.RestrictedToGeneSets != null)
            {
                var missing = species.Where(s => !geneTrees.RestrictedToGeneSets.Contains(geneSets[s].SourceSha256)).ToList();
                if (missing.Count > 0)
                {
                    throw new ArgumentException(
                        $"{geneTrees.SourceFileName} was restricted to other gene sets; it may lack genes of {string.Join(", ", missing)}.",
                        nameof(geneTrees));
                }
            }

            var dumpOf = new Dictionary<string, ComparaHomologyDump>(StringComparer.Ordinal);
            foreach (var dump in dumps)
            {
                if (dump == null)
                {
                    throw new ArgumentException("A dump is null.", nameof(dumps));
                }
                RequireRelease(dump.Release, release, dump.SourceFileName);
                if (dump.Genome == null || !geneSets.ContainsKey(dump.Genome))
                {
                    throw new ArgumentException(
                        $"{dump.SourceFileName} is a dump of {dump.Genome ?? "no genome (it is empty)"}, which is not one of the species.",
                        nameof(dumps));
                }
                if (!dumpOf.TryAdd(dump.Genome, dump))
                {
                    throw new ArgumentException($"Two dumps of {dump.Genome} were given.", nameof(dumps));
                }
                var notLoaded = species.Where(s => !dump.Species.Contains(s)).ToList();
                if (notLoaded.Count > 0)
                {
                    throw new ArgumentException(
                        $"{dump.SourceFileName} was loaded without {string.Join(", ", notLoaded)}.", nameof(dumps));
                }
            }
            var noDump = species.Where(s => !dumpOf.ContainsKey(s)).ToList();
            if (noDump.Count > 0)
            {
                throw new ArgumentException(
                    $"No homology dump for {string.Join(", ", noDump)}. Each pair's relationships are split between " +
                    "its two species' dumps, so every species' dump is required.", nameof(dumps));
            }

            var byId = new Dictionary<string, ComparaHomology>(StringComparer.Ordinal);
            int duplicates = 0;
            foreach (string s in species)
            {
                var dump = dumpOf[s];
                foreach (var row in dump.Rows)
                {
                    if (!geneSets.ContainsKey(row.SpeciesA) || !geneSets.ContainsKey(row.SpeciesB))
                    {
                        continue;
                    }
                    RequireGene(row.GeneA, row.SpeciesA, row, geneSets);
                    RequireGene(row.GeneB, row.SpeciesB, row, geneSets);

                    var oriented = string.CompareOrdinal(row.SpeciesA, row.SpeciesB) > 0 ? row.Swapped() : row;
                    if (byId.TryGetValue(row.HomologyId, out var prior))
                    {
                        if (!prior.SameRelationship(oriented))
                        {
                            throw new InvalidDataException(
                                $"Homology {row.HomologyId} has different rows in {prior.SourceFile} and {row.SourceFile}.");
                        }
                        duplicates++;
                        continue;
                    }
                    byId[row.HomologyId] = oriented;
                }
            }

            var pairs = byId.Values
                .GroupBy(r => (r.SpeciesA, r.SpeciesB))
                .ToDictionary(g => g.Key,
                    g => g.OrderBy(r => r.HomologyId.Length).ThenBy(r => r.HomologyId, StringComparer.Ordinal).ToList());

            return new OrthologySnapshot(release, species, geneSets, geneTrees,
                species.Select(s => dumpOf[s]).ToList(), pairs, speciesOfGene, duplicates);
        }

        /// <summary>The status as written in a table, e.g. "tree_lacks_target_species".</summary>
        public static string StatusName(OrthologyStatus status) => status switch
        {
            OrthologyStatus.HasOrtholog => "has_ortholog",
            OrthologyStatus.NoEdgeInSharedTree => "no_edge_in_shared_tree",
            OrthologyStatus.TreeLacksTargetSpecies => "tree_lacks_target_species",
            OrthologyStatus.NotInAnyTree => "not_in_any_tree",
            OrthologyStatus.NotInGeneSet => "not_in_gene_set",
            _ => throw new ArgumentOutOfRangeException(nameof(status), status, null)
        };

        private void AddOrtholog(string gene, string partnerSpecies, string partner, string type)
        {
            if (!_orthologs.TryGetValue(gene, out var bySpecies))
            {
                _orthologs[gene] = bySpecies = new Dictionary<string, List<(string, string)>>(StringComparer.Ordinal);
            }
            if (!bySpecies.TryGetValue(partnerSpecies, out var list))
            {
                bySpecies[partnerSpecies] = list = new List<(string, string)>();
            }
            list.Add((partner, type));
        }

        private IEnumerable<(string Gene, string Type)> Partners(string gene, string species) =>
            _orthologs.TryGetValue(gene, out var bySpecies) && bySpecies.TryGetValue(species, out var list)
                ? list.OrderBy(p => p.Gene, StringComparer.Ordinal)
                : Enumerable.Empty<(string, string)>();

        /// <summary>The homology_type of the ortholog between two genes, or null when there is none.</summary>
        private string EdgeType(string gene, string partnerSpecies, string partner) =>
            Partners(gene, partnerSpecies).Where(p => p.Gene == partner).Select(p => p.Type).FirstOrDefault();

        private void RequireSpecies(string species, string parameter)
        {
            if (species == null || !GeneSets.ContainsKey(species))
            {
                throw new ArgumentException($"'{species}' is not one of the snapshot's species.", parameter);
            }
        }

        private void RequireDistinct(string from, string to)
        {
            RequireSpecies(from, nameof(from));
            RequireSpecies(to, nameof(to));
            if (from == to)
            {
                throw new ArgumentException("Orthology is between two different species; use Pair(s, s) for paralogs.");
            }
        }

        private static void RequireRelease(string found, string expected, string source)
        {
            if (found != null && found != expected)
            {
                throw new ArgumentException($"{source} is from release {found}, not {expected}.");
            }
        }

        private static void RequireGene(string gene, string species, ComparaHomology row,
            IReadOnlyDictionary<string, EnsemblGeneSet> geneSets)
        {
            if (!geneSets[species].Contains(gene))
            {
                throw new InvalidDataException(
                    $"{row.SourceFile}: homology {row.HomologyId} names {gene}, which is not in the {species} gene set " +
                    $"({geneSets[species].SourceFileName}).");
            }
        }
    }
}
