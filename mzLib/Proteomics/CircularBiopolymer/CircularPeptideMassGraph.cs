using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Proteomics.CircularBiopolymer
{
    public sealed class CircularPeptideMassGraph
    {
        public sealed class AminoAcid
        {
            public char Symbol { get; }
            public double MonoisotopicMass { get; }

            public AminoAcid(char symbol, double monoisotopicMass)
            {
                Symbol = symbol;
                MonoisotopicMass = monoisotopicMass;
            }
        }

        public sealed class CompositionNode
        {
            public string Key { get; }
            public int Length { get; }
            public double Mass { get; }
            public IReadOnlyDictionary<char, int> Counts { get; }

            public HashSet<CompositionNode> Parents { get; } = new();
            public HashSet<CompositionNode> Children { get; } = new();

            internal CompositionNode(string key, int length, double mass, Dictionary<char, int> countsSnapshot)
            {
                Key = key;
                Length = length;
                Mass = mass;
                Counts = countsSnapshot;
            }

            public override string ToString() => $"{Key} (L={Length}, Mass={Mass:F5})";
        }

        private readonly IReadOnlyList<AminoAcid> _alphabet;
        private readonly Dictionary<char, AminoAcid> _aaBySymbol;
        private readonly Dictionary<string, CompositionNode> _nodes = new(StringComparer.Ordinal);

        public CircularPeptideMassGraph(IEnumerable<AminoAcid> alphabet)
        {
            _alphabet = alphabet?.ToList() ?? throw new ArgumentNullException(nameof(alphabet));
            if (_alphabet.Count == 0) throw new ArgumentException("Alphabet must contain at least one residue.", nameof(alphabet));
            _aaBySymbol = _alphabet.ToDictionary(a => a.Symbol, a => a);
        }

        public IReadOnlyCollection<CompositionNode> Nodes => _nodes.Values;

        public void Build(int maxLength)
        {
            if (maxLength < 1) throw new ArgumentOutOfRangeException(nameof(maxLength), "maxLength must be >= 1.");
            var counts = _alphabet.ToDictionary(a => a.Symbol, _ => 0);

            foreach (var aa in _alphabet)
            {
                counts[aa.Symbol] = 1;
                var leaf = GetOrCreate(counts);
                ExpandFrom(leaf, counts, currentLength: 1, maxLength);
                counts[aa.Symbol] = 0;
            }
        }

        private CompositionNode GetOrCreate(Dictionary<char, int> counts)
        {
            string key = CanonicalKey(counts);
            if (_nodes.TryGetValue(key, out var existing))
                return existing;

            int length = counts.Values.Sum();
            double mass = counts.Sum(kv => kv.Value * _aaBySymbol[kv.Key].MonoisotopicMass);
            var snapshot = counts.Where(kv => kv.Value > 0).ToDictionary(kv => kv.Key, kv => kv.Value);

            var node = new CompositionNode(key, length, mass, snapshot);
            _nodes.Add(key, node);
            return node;
        }

        private void ExpandFrom(CompositionNode child, Dictionary<char, int> counts, int currentLength, int maxLength)
        {
            if (currentLength >= maxLength)
                return;

            foreach (var aa in _alphabet)
            {
                counts[aa.Symbol]++;
                var parent = GetOrCreate(counts);
                if (child.Parents.Add(parent))
                    parent.Children.Add(child);
                ExpandFrom(parent, counts, currentLength + 1, maxLength);
                counts[aa.Symbol]--;
            }
        }

        private static string CanonicalKey(Dictionary<char, int> counts)
        {
            return string.Concat(counts
                .Where(kv => kv.Value > 0)
                .OrderBy(kv => kv.Key)
                .Select(kv => $"{kv.Key}{kv.Value}"));
        }

        private sealed class Coverage
        {
            public int Mask;
            public double[] BestErrorPerPeak;
            public Coverage(int k)
            {
                Mask = 0;
                BestErrorPerPeak = Enumerable.Repeat(double.PositiveInfinity, k).ToArray();
            }
        }

        private static int PopCount(int x) => BitOperations.PopCount((uint)x);

        /// <summary>
        /// Fast composition-level predicate from motifs:
        /// requires each residue's count in the composition to be at least
        /// the maximum count it appears in any single motif (lower bound).
        /// Example: motifs ["ACG", "AGG"] ? A?1, C?1, G?2.
        /// </summary>
        public static Func<CompositionNode, bool> MustCoverMotifResidueCounts(params string[] motifs)
        {
            if (motifs == null || motifs.Length == 0)
                return _ => true;

            var required = new Dictionary<char, int>();
            foreach (var motif in motifs.Where(m => !string.IsNullOrWhiteSpace(m)))
            {
                var local = new Dictionary<char, int>();
                foreach (char c in motif)
                {
                    local.TryGetValue(c, out int v);
                    local[c] = v + 1;
                }
                foreach (var kv in local)
                {
                    if (!required.TryGetValue(kv.Key, out int cur) || kv.Value > cur)
                        required[kv.Key] = kv.Value;
                }
            }
            if (required.Count == 0) return _ => true;

            return node =>
            {
                foreach (var kv in required)
                {
                    if (!node.Counts.TryGetValue(kv.Key, out int have) || have < kv.Value)
                        return false;
                }
                return true;
            };
        }

        /// <summary>
        /// Original ranking API (kept for compatibility).
        /// </summary>
        public IReadOnlyList<(CompositionNode node, int explainedCount, double totalError)>
            RankParentsByCoverage(double[] peaks, double toleranceDa, int targetParentLength, int topK = 5)
        {
            return RankParentsByCoverage(peaks, toleranceDa, targetParentLength, topK, candidateFilter: null);
        }

        /// <summary>
        /// Rank parents at target length. Optionally apply a candidate-node predicate to prune the search space
        /// (e.g., motif-based filter). Peaks that do not match any node naturally do not contribute to coverage.
        /// </summary>
        public IReadOnlyList<(CompositionNode node, int explainedCount, double totalError)>
            RankParentsByCoverage(double[] peaks,
                                  double toleranceDa,
                                  int targetParentLength,
                                  int topK,
                                  Func<CompositionNode, bool> candidateFilter)
        {
            if (peaks == null || peaks.Length == 0) throw new ArgumentException("Provide at least one peak.", nameof(peaks));
            if (toleranceDa <= 0) throw new ArgumentOutOfRangeException(nameof(toleranceDa), "Tolerance must be positive.");
            if (!_nodes.Any()) throw new InvalidOperationException("Graph is empty; call Build() first.");

            int K = peaks.Length;
            var coverage = new Dictionary<CompositionNode, Coverage>();

            Coverage GetCov(CompositionNode n)
            {
                if (!coverage.TryGetValue(n, out var c))
                    coverage[n] = c = new Coverage(K);
                return c;
            }

            // Seed: for each peak, match nodes and propagate upward
            foreach (var (peak, idx) in peaks.Select((p, i) => (p, i)))
            {
                foreach (var match in _nodes.Values)
                {
                    double err = Math.Abs(match.Mass - peak);
                    if (err <= toleranceDa)
                        PropagateUp(match, idx, err, GetCov);
                }
            }

            // Assemble and (optionally) filter candidate parents
            var candidates = coverage
                .Where(kv => kv.Key.Length == targetParentLength)
                .Where(kv => candidateFilter == null || candidateFilter(kv.Key))
                .Select(kv =>
                {
                    int explained = PopCount(kv.Value.Mask);
                    double totalErr = 0.0;
                    for (int i = 0; i < K; i++)
                        if (((kv.Value.Mask >> i) & 1) == 1)
                            totalErr += kv.Value.BestErrorPerPeak[i];
                    return (node: kv.Key, explainedCount: explained, totalError: totalErr);
                })
                .OrderByDescending(r => r.explainedCount)
                .ThenBy(r => r.totalError)
                .ThenBy(r => r.node.Mass)
                .Take(topK)
                .ToList();

            return candidates;

            void PropagateUp(CompositionNode node, int peakIndex, double error, Func<CompositionNode, Coverage> covFactory)
            {
                var cov = covFactory(node);
                if (error + 1e-12 >= cov.BestErrorPerPeak[peakIndex]) return;

                cov.BestErrorPerPeak[peakIndex] = error;
                cov.Mask |= (1 << peakIndex);

                foreach (var parent in node.Parents)
                    PropagateUp(parent, peakIndex, error, covFactory);
            }
        }
    }
}