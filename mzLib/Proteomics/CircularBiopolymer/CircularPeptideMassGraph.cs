using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Proteomics.CircularBiopolymer
{
    /// <summary>
    /// Order-agnostic, mass-labeled, compositional graph for circular peptides.
    /// Nodes represent residue compositions (multisets) rather than ordered sequences.
    /// An edge connects a composition of length L to all compositions of length L+1
    /// obtainable by adding exactly one residue (i.e., all possible parents).
    ///
    /// This captures all permutations in a single network. Traversal can then be
    /// performed from matched fragment masses upward to rank likely parents.
    /// </summary>
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

        /// <summary>
        /// Build nodes for all compositions (with repetition) of length 1..maxLength.
        /// Connect each node to its length+1 parents by adding one residue.
        /// </summary>
        public void Build(int maxLength)
        {
            if (maxLength < 1) throw new ArgumentOutOfRangeException(nameof(maxLength), "maxLength must be >= 1.");

            var counts = _alphabet.ToDictionary(a => a.Symbol, _ => 0);

            // Seed all leaves (length 1) and recursively expand upward.
            foreach (var aa in _alphabet)
            {
                counts[aa.Symbol] = 1;
                var leaf = GetOrCreate(counts);
                ExpandFrom(leaf, counts, currentLength: 1, maxLength);
                counts[aa.Symbol] = 0; // backtrack
            }
        }

        private CompositionNode GetOrCreate(Dictionary<char, int> counts)
        {
            string key = CanonicalKey(counts);
            if (_nodes.TryGetValue(key, out var existing))
                return existing;

            int length = counts.Values.Sum();
            double mass = counts.Sum(kv => kv.Value * _aaBySymbol[kv.Key].MonoisotopicMass);

            var snapshot = counts
                .Where(kv => kv.Value > 0)
                .ToDictionary(kv => kv.Key, kv => kv.Value);

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

                counts[aa.Symbol]--; // backtrack
            }
        }

        private static string CanonicalKey(Dictionary<char, int> counts)
        {
            // Sort by residue symbol; omit zero counts; format e.g. A2C1G1
            return string.Concat(counts
                .Where(kv => kv.Value > 0)
                .OrderBy(kv => kv.Key)
                .Select(kv => $"{kv.Key}{kv.Value}"));
        }

        private sealed class Coverage
        {
            public int Mask;                // bitmask: which peaks are explained by any descendant
            public double[] BestErrorPerPeak;

            public Coverage(int k)
            {
                Mask = 0;
                BestErrorPerPeak = Enumerable.Repeat(double.PositiveInfinity, k).ToArray();
            }
        }

        private static int PopCount(int x) => BitOperations.PopCount((uint)x);

        /// <summary>
        /// Given a set of peaks (fragment masses) and a tolerance, rank likely parent compositions at the target length.
        /// Strategy:
        /// - For each peak, find all nodes within tolerance and recursively propagate a coverage bit up to all ancestors.
        /// - Each ancestor collects a mask of matched peaks and best mass error per peak.
        /// Ranking: by explained peak count (desc), then total error (asc), then mass (asc).
        ///
        /// Note: peaks that do not match any node within tolerance can be ignored as noise by setting
        /// ignorePeaksNotInGraph=true (default). When enabled, non-matching peaks do not penalize candidates.
        /// </summary>
        /// <param name="peaks">Observed fragment masses.</param>
        /// <param name="toleranceDa">Absolute tolerance in Daltons for mass matching.</param>
        /// <param name="targetParentLength">Length at which to consider parents (typically the top length you built).</param>
        /// <param name="topK">Max number of results to return.</param>
        /// <param name="ignorePeaksNotInGraph">If true, peaks that do not match any node are treated as noise and ignored.</param>
        public IReadOnlyList<(CompositionNode node, int explainedCount, double totalError)>
            RankParentsByCoverage(double[] peaks, double toleranceDa, int targetParentLength, int topK = 5, bool ignorePeaksNotInGraph = true)
        {
            if (peaks == null || peaks.Length == 0) throw new ArgumentException("Provide at least one peak.", nameof(peaks));
            if (toleranceDa <= 0) throw new ArgumentOutOfRangeException(nameof(toleranceDa), "Tolerance must be positive.");
            if (!_nodes.Any()) throw new InvalidOperationException("Graph is empty; call Build() first.");

            // Optionally filter out peaks that have no match anywhere in the graph (noise)
            var nodeList = _nodes.Values.ToList();
            int K = peaks.Length;
            int[] matchableIdx = Enumerable.Range(0, K)
                .Where(i => nodeList.Any(n => Math.Abs(n.Mass - peaks[i]) <= toleranceDa))
                .ToArray();

            double[] workPeaks = peaks;
            int[] indexMap = null; // maps workPeaks index -> original peak index (unused outside scoring)
            if (ignorePeaksNotInGraph && matchableIdx.Length < K)
            {
                workPeaks = matchableIdx.Select(i => peaks[i]).ToArray();
                indexMap = matchableIdx;
            }

            int W = workPeaks.Length;
            if (W == 0)
            {
                // Nothing to explain (all noise); return empty
                return Array.Empty<(CompositionNode node, int explainedCount, double totalError)>();
            }

            var coverage = new Dictionary<CompositionNode, Coverage>();
            Coverage GetCov(CompositionNode n)
            {
                if (!coverage.TryGetValue(n, out var c))
                    coverage[n] = c = new Coverage(W);
                return c;
            }

            // Seed: for each work-peak, match nodes and propagate upward
            for (int j = 0; j < W; j++)
            {
                double peak = workPeaks[j];
                foreach (var match in nodeList)
                {
                    double err = Math.Abs(match.Mass - peak);
                    if (err <= toleranceDa)
                    {
                        PropagateUp(match, j, err, GetCov);
                    }
                }
            }

            // Assemble candidate parents at the requested length
            var candidates = coverage
                .Where(kv => kv.Key.Length == targetParentLength)
                .Select(kv =>
                {
                    int explained = PopCount(kv.Value.Mask);
                    double totalErr = 0.0;
                    for (int i = 0; i < W; i++)
                    {
                        if (((kv.Value.Mask >> i) & 1) == 1)
                            totalErr += kv.Value.BestErrorPerPeak[i];
                    }
                    return (node: kv.Key, explainedCount: explained, totalError: totalErr);
                })
                .OrderByDescending(r => r.explainedCount)
                .ThenBy(r => r.totalError)
                .ThenBy(r => r.node.Mass)
                .Take(topK)
                .ToList();

            return candidates;

            void PropagateUp(CompositionNode node, int workPeakIndex, double error, Func<CompositionNode, Coverage> covFactory)
            {
                var cov = covFactory(node);
                // Only propagate if we improve the best error for this peak at this node
                if (error + 1e-12 >= cov.BestErrorPerPeak[workPeakIndex])
                    return;

                cov.BestErrorPerPeak[workPeakIndex] = error;
                cov.Mask |= (1 << workPeakIndex);

                foreach (var parent in node.Parents)
                    PropagateUp(parent, workPeakIndex, error, covFactory);
            }
        }
    }
}