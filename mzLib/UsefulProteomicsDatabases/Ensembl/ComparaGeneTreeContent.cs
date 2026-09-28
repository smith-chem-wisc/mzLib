using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>A gene's place in Compara's gene trees.</summary>
    /// <param name="GeneId">The stable gene id.</param>
    /// <param name="TreeId">Compara's gene-tree stable id, e.g. "ENSGT00940000153582".</param>
    /// <param name="CanonicalProteinId">The protein Compara used for the gene (the row flagged Y).</param>
    public sealed record ComparaGeneTreeMember(string GeneId, string TreeId, string CanonicalProteinId);

    /// <summary>
    /// Compara's gene-tree membership dump (compara/&lt;collection&gt;.GeneTree_content.default.e&lt;release&gt;.txt.gz):
    /// which tree each gene belongs to. It has no header; its four columns are tree id, protein id, gene
    /// id and whether that protein is the gene's canonical one (Y or N).
    ///
    /// A gene tree is the group the source asserts. It is what lets a store say "the tree holding this
    /// gene has no gene of the other species" apart from "it has one, and no ortholog was called", and
    /// it is the only grouping here that is not derived by chaining pairwise calls.
    ///
    /// The reader refuses rather than guesses: a row without four columns, a flag other than Y or N, a
    /// gene in two trees, or a kept gene without exactly one canonical protein each throw. In release
    /// 116 every one of its 4.2 million genes is in one tree with one canonical protein.
    /// </summary>
    public sealed class ComparaGeneTreeContent
    {
        private static readonly Regex ReleaseInFileName =
            new(@"\.e(\d+)\.txt(\.gz)?$", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        private static readonly Regex CollectionInFileName =
            new(@"^([A-Za-z_]+)\.GeneTree_content\.", RegexOptions.Compiled);

        private readonly Dictionary<string, ComparaGeneTreeMember> _members;

        private ComparaGeneTreeContent(Dictionary<string, ComparaGeneTreeMember> members, string sourceFileName,
            string sourceSha256, string release, string collection, IReadOnlyList<string> restrictedTo)
        {
            _members = members;
            SourceFileName = sourceFileName;
            SourceSha256 = sourceSha256;
            Release = release;
            Collection = collection;
            RestrictedToGeneSets = restrictedTo;
        }

        public string SourceFileName { get; }

        /// <summary>Lower-case hex sha256 of the file's bytes as read (the compressed bytes for a .gz).</summary>
        public string SourceSha256 { get; }

        /// <summary>The Ensembl release parsed from the file name, or null when it carries none.</summary>
        public string Release { get; }

        /// <summary>The Compara collection parsed from the file name (e.g. "vertebrates"), or null.</summary>
        public string Collection { get; }

        /// <summary>
        /// The sha256s of the gene sets the members were restricted to, or null when every gene was kept.
        /// A consumer checks this rather than trusting that a missing gene is outside every tree.
        /// </summary>
        public IReadOnlyList<string> RestrictedToGeneSets { get; }

        /// <summary>The number of genes kept.</summary>
        public int Count => _members.Count;

        /// <summary>Every kept member, in ordinal order of gene id.</summary>
        public IEnumerable<ComparaGeneTreeMember> Members =>
            _members.Values.OrderBy(m => m.GeneId, StringComparer.Ordinal);

        public bool TryGetMember(string geneId, out ComparaGeneTreeMember member)
        {
            member = null;
            return geneId != null && _members.TryGetValue(geneId, out member);
        }

        /// <summary>
        /// Reads the dump. With <paramref name="restrictTo"/>, only genes in those gene sets are kept, and
        /// their sha256s are recorded in <see cref="RestrictedToGeneSets"/>; without it, every gene is kept
        /// (millions, for the vertebrate collection). Every row is checked either way.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file breaks one of the rules in the class summary.</exception>
        public static ComparaGeneTreeContent Load(string path, IEnumerable<EnsemblGeneSet> restrictTo = null)
        {
            List<EnsemblGeneSet> sets = restrictTo?.ToList();
            if (sets != null && sets.Any(s => s == null))
            {
                throw new ArgumentException("A gene set is null.", nameof(restrictTo));
            }

            string sha256 = EnsemblFile.Sha256(path, "Compara gene-tree content");
            string name = Path.GetFileName(path);
            var trees = new Dictionary<string, string>(StringComparer.Ordinal);
            var canonical = new Dictionary<string, List<string>>(StringComparer.Ordinal);

            using (var reader = EnsemblFile.OpenText(path))
            {
                int lineNumber = 0;
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    lineNumber++;
                    string[] c = line.Split('\t');
                    if (c.Length != 4)
                    {
                        throw new InvalidDataException($"{name} line {lineNumber}: {c.Length} cells, expected 4.");
                    }

                    string tree = c[0], protein = c[1], gene = c[2], flag = c[3];
                    if (flag != "Y" && flag != "N")
                    {
                        throw new InvalidDataException($"{name} line {lineNumber}: canonical flag '{flag}' is neither Y nor N.");
                    }

                    if (sets != null && !sets.Any(s => s.Contains(gene)))
                    {
                        continue;
                    }

                    if (!trees.TryAdd(gene, tree) && trees[gene] != tree)
                    {
                        throw new InvalidDataException(
                            $"{name} line {lineNumber}: {gene} is in two gene trees, {trees[gene]} and {tree}.");
                    }

                    if (flag == "Y")
                    {
                        if (!canonical.TryGetValue(gene, out var proteins))
                        {
                            canonical[gene] = proteins = new List<string>();
                        }
                        proteins.Add(protein);
                    }
                }
            }

            var members = new Dictionary<string, ComparaGeneTreeMember>(StringComparer.Ordinal);
            foreach (var (gene, tree) in trees)
            {
                int n = canonical.TryGetValue(gene, out var proteins) ? proteins.Count : 0;
                if (n != 1)
                {
                    throw new InvalidDataException($"{name}: {gene} has {n} canonical proteins, expected 1.");
                }
                members[gene] = new ComparaGeneTreeMember(gene, tree, proteins[0]);
            }

            var release = ReleaseInFileName.Match(name);
            var collection = CollectionInFileName.Match(name);
            return new ComparaGeneTreeContent(members, name, sha256,
                release.Success ? release.Groups[1].Value : null,
                collection.Success ? collection.Groups[1].Value : null,
                sets?.Select(s => s.SourceSha256).Distinct(StringComparer.Ordinal)
                    .OrderBy(s => s, StringComparer.Ordinal).ToList());
        }
    }
}
