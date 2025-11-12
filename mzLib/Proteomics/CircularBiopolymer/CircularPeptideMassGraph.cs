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

            // Parents are enough for upward traversal (coverage propagation).
            public HashSet<CompositionNode> Parents { get; } = new();

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

        // Mass index for fast peak-to-node matching
        private bool _massIndexDirty = true;
        private double[] _sortedMasses = Array.Empty<double>();
        private CompositionNode[] _sortedNodes = Array.Empty<CompositionNode>();

        public CircularPeptideMassGraph(IEnumerable<AminoAcid> alphabet)
        {
            _alphabet = alphabet?.ToList() ?? throw new ArgumentNullException(nameof(alphabet));
            if (_alphabet.Count == 0) throw new ArgumentException("Alphabet must contain at least one residue.", nameof(alphabet));
            _aaBySymbol = _alphabet.ToDictionary(a => a.Symbol, a => a);
        }

        public IReadOnlyCollection<CompositionNode> Nodes => _nodes.Values;

        /// <summary>
        /// Build nodes for all compositions (with repetition) of length 1..maxLength.
        /// Uses a non-decreasing index generator to avoid duplicates and builds level-by-level
        /// to minimize redundant work. This version constructs each composition exactly once.
        /// </summary>
        public void Build(int maxLength)
        {
            BuildInternal(maxLength, pruneFeasible: null);
        }

        /// <summary>
        /// Build nodes up to maxLength while pruning branches that cannot meet
        /// the provided per-residue lower bounds by target length (feasibility check).
        /// This can drastically reduce the search space when you have motif-derived minima.
        /// </summary>
        public void Build(int maxLength, Dictionary<char, int> lowerBounds)
        {
            // Convert lowerBounds (symbol -> minCount) to array aligned with alphabet order for speed
            int[] required = new int[_alphabet.Count];
            if (lowerBounds != null)
            {
                for (int i = 0; i < _alphabet.Count; i++)
                {
                    var sym = _alphabet[i].Symbol;
                    required[i] = lowerBounds.TryGetValue(sym, out var v) ? v : 0;
                }
            }
            else
            {
                // all zeros => no pruning
                required = null;
            }

            BuildInternal(maxLength, pruneFeasible: required);
        }

        /// <summary>
        /// Build lower bounds from motifs (same aggregation logic as MustCoverMotifResidueCounts).
        /// Useful to prune during graph construction.
        /// </summary>
        public static Dictionary<char, int> MotifLowerBounds(params string[] motifs)
        {
            var required = new Dictionary<char, int>();
            if (motifs == null) return required;

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
            return required;
        }

        private void BuildInternal(int maxLength, int[] pruneFeasible)
        {
            if (maxLength < 1) throw new ArgumentOutOfRangeException(nameof(maxLength), "maxLength must be >= 1.");

            int n = _alphabet.Count;
            int[] counts = new int[n];

            var layers = new List<List<CompositionNode>>(maxLength + 1);
            layers.Add(new List<CompositionNode>()); // dummy for length 0

            for (int L = 1; L <= maxLength; L++)
            {
                var layer = new List<CompositionNode>();
                Compose(0, L, layer); // FIX: write nodes into the current layer
                layers.Add(layer);

                // Connect L -> L+1 parents (only need upward traversal)
                if (L < maxLength)
                {
                    foreach (var node in layer)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            var parentCounts = new Dictionary<char, int>(node.Counts);
                            char sym = _alphabet[i].Symbol;
                            parentCounts.TryGetValue(sym, out int cur);
                            parentCounts[sym] = cur + 1;
                            var parent = GetOrCreate(parentCounts);
                            node.Parents.Add(parent);
                        }
                    }
                }
            }

            _massIndexDirty = true;
            return;

            // Generate all non-decreasing-index compositions of 'remaining' across [start..n-1]
            void Compose(int start, int remaining, List<CompositionNode> targetLayer)
            {
                if (remaining == 0)
                {
                    var snap = new Dictionary<char, int>();
                    double mass = 0;
                    for (int i = 0; i < n; i++)
                    {
                        if (counts[i] > 0)
                        {
                            snap[_alphabet[i].Symbol] = counts[i];
                            mass += counts[i] * _alphabet[i].MonoisotopicMass;
                        }
                    }
                    var node = GetOrCreate(snap, mass);
                    targetLayer.Add(node); // FIX: add to the correct layer
                    return;
                }

                if (pruneFeasible != null)
                {
                    int need = 0;
                    for (int i = 0; i < n; i++)
                    {
                        int deficit = pruneFeasible[i] - counts[i];
                        if (deficit > 0) need += deficit;
                    }
                    if (need > remaining)
                        return;
                }

                for (int i = start; i < n; i++)
                {
                    counts[i]++;
                    Compose(i, remaining - 1, targetLayer);
                    counts[i]--;
                }
            }
        }

        private CompositionNode GetOrCreate(Dictionary<char, int> countsSnapshot, double? knownMass = null)
        {
            string key = CanonicalKey(countsSnapshot);
            if (_nodes.TryGetValue(key, out var existing))
                return existing;

            int length = countsSnapshot.Values.Sum();
            double mass = knownMass ?? countsSnapshot.Sum(kv => kv.Value * _aaBySymbol[kv.Key].MonoisotopicMass);

            var node = new CompositionNode(key, length, mass, new Dictionary<char, int>(countsSnapshot));
            _nodes.Add(key, node);
            _massIndexDirty = true;
            return node;
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
        /// Rank parents at target length. Optionally apply a candidate-node predicate to prune the search space.
        /// Uses a mass index for O(log N + m) peak matching instead of scanning all nodes.
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

            EnsureMassIndex();

            int K = peaks.Length;
            var coverage = new Dictionary<CompositionNode, Coverage>(capacity: Math.Min(_sortedNodes.Length, 8192));

            Coverage GetCov(CompositionNode n)
            {
                if (!coverage.TryGetValue(n, out var c))
                    coverage[n] = c = new Coverage(K);
                return c;
            }

            // Seed: for each peak, match nodes via binary search window and propagate upward
            for (int i = 0; i < K; i++)
            {
                double p = peaks[i];
                foreach (var match in WindowByMass(p, toleranceDa))
                {
                    double err = Math.Abs(match.Mass - p);
                    PropagateUp(match, i, err, GetCov);
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

        private void EnsureMassIndex()
        {
            if (!_massIndexDirty && _sortedMasses.Length == _nodes.Count)
                return;

            // Build sorted arrays of (mass, node) for fast window queries
            var arr = _nodes.Values.ToArray();
            Array.Sort(arr, static (a, b) => a.Mass.CompareTo(b.Mass));
            _sortedNodes = arr;
            _sortedMasses = arr.Select(n => n.Mass).ToArray();
            _massIndexDirty = false;
        }

        private IEnumerable<CompositionNode> WindowByMass(double target, double tol)
        {
            if (_sortedMasses.Length == 0) yield break;

            double loVal = target - tol;
            double hiVal = target + tol;

            int lo = LowerBound(_sortedMasses, loVal);
            for (int i = lo; i < _sortedMasses.Length && _sortedMasses[i] <= hiVal; i++)
                yield return _sortedNodes[i];

            static int LowerBound(double[] a, double value)
            {
                int l = 0, r = a.Length;
                while (l < r)
                {
                    int m = (l + r) >> 1;
                    if (a[m] < value) l = m + 1;
                    else r = m;
                }
                return l;
            }
        }
    }
}