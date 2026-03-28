using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Proteomics.CircularBiopolymer
{
    public sealed class CircularPeptideMassGraph
    {
        public const ushort MassScale = 10; // store masses as (Daltons * MassScale), e.g. tenths of a Dalton

        private static ushort Scale(double mono)
        {
            double scaled = Math.Round(mono * MassScale);
            if (scaled < 0 || scaled > ushort.MaxValue)
                throw new ArgumentOutOfRangeException(nameof(mono), $"Scaled mass {scaled} out of ushort range");
            return (ushort)scaled;
        }

        public sealed class AminoAcid
        {
            public char Symbol { get; }
            public ushort ScaledMonoisotopicMass { get; } // tenths of a Dalton

            public AminoAcid(char symbol, double monoisotopicMass)
            {
                Symbol = symbol;
                ScaledMonoisotopicMass = Scale(monoisotopicMass);
            }
        }

        public sealed class CompositionNode
        {
            public string Key { get; }
            public int Length { get; }
            public ushort Mass { get; } // scaled mass (sum of residue scaled masses)
            public IReadOnlyDictionary<char, int> Counts { get; }

            public HashSet<CompositionNode> Parents { get; } = new();

            internal CompositionNode(string key, int length, ushort mass, Dictionary<char, int> countsSnapshot)
            {
                Key = key;
                Length = length;
                Mass = mass;
                Counts = countsSnapshot;
            }

            public override string ToString() => $"{Key} (L={Length}, Mass={Mass / (double)MassScale:F3})";
        }

        private readonly IReadOnlyList<AminoAcid> _alphabet;
        private readonly Dictionary<char, AminoAcid> _aaBySymbol;
        private readonly Dictionary<string, CompositionNode> _nodes = new(StringComparer.Ordinal);

        private bool _massIndexDirty = true;
        private ushort[] _sortedMasses = Array.Empty<ushort>();
        private CompositionNode[] _sortedNodes = Array.Empty<CompositionNode>();

        public CircularPeptideMassGraph(IEnumerable<AminoAcid> alphabet)
        {
            _alphabet = alphabet?.ToList() ?? throw new ArgumentNullException(nameof(alphabet));
            if (_alphabet.Count == 0) throw new ArgumentException("Alphabet must contain at least one residue.", nameof(alphabet));
            _aaBySymbol = _alphabet.ToDictionary(a => a.Symbol, a => a);
        }

        public IReadOnlyCollection<CompositionNode> Nodes => _nodes.Values;

        public void Build(int maxLength) => BuildInternal(maxLength, pruneFeasible: null);

        public void Build(int maxLength, Dictionary<char, int> lowerBounds)
        {
            int[] required = null;
            if (lowerBounds != null)
            {
                required = new int[_alphabet.Count];
                for (int i = 0; i < _alphabet.Count; i++)
                {
                    var sym = _alphabet[i].Symbol;
                    required[i] = lowerBounds.TryGetValue(sym, out var v) ? v : 0;
                }
            }
            BuildInternal(maxLength, pruneFeasible: required);
        }

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

            for (int L = 1; L <= maxLength; L++)
            {
                var layer = new List<CompositionNode>();
                Compose(0, L, layer);
                // Connect to parents for next layer
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

            void Compose(int start, int remaining, List<CompositionNode> targetLayer)
            {
                if (remaining == 0)
                {
                    var snap = new Dictionary<char, int>();
                    int length = 0;
                    int totalScaledMass = 0;
                    for (int i = 0; i < n; i++)
                    {
                        if (counts[i] > 0)
                        {
                            int c = counts[i];
                            snap[_alphabet[i].Symbol] = c;
                            length += c;
                            totalScaledMass += c * _alphabet[i].ScaledMonoisotopicMass;
                        }
                    }
                    if (totalScaledMass > ushort.MaxValue)
                        throw new InvalidOperationException($"Scaled mass {totalScaledMass} exceeds ushort limit");
                    var node = GetOrCreate(snap, (ushort)totalScaledMass);
                    targetLayer.Add(node);
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

        private CompositionNode GetOrCreate(Dictionary<char, int> countsSnapshot, ushort? knownMass = null)
        {
            string key = CanonicalKey(countsSnapshot);
            if (_nodes.TryGetValue(key, out var existing))
                return existing;

            int length = countsSnapshot.Values.Sum();
            int totalScaled = knownMass ?? countsSnapshot.Sum(kv => kv.Value * _aaBySymbol[kv.Key].ScaledMonoisotopicMass);
            if (totalScaled > ushort.MaxValue)
                throw new InvalidOperationException($"Scaled mass {totalScaled} exceeds ushort limit");
            var node = new CompositionNode(key, length, (ushort)totalScaled, new Dictionary<char, int>(countsSnapshot));
            _nodes.Add(key, node);
            _massIndexDirty = true;
            return node;
        }

        private static string CanonicalKey(Dictionary<char, int> counts) =>
            string.Concat(counts.Where(kv => kv.Value > 0).OrderBy(kv => kv.Key).Select(kv => $"{kv.Key}{kv.Value}"));

        private sealed class Coverage
        {
            public int Mask;
            public ushort[] BestErrorPerPeak; // scaled units
            public Coverage(int k)
            {
                Mask = 0;
                BestErrorPerPeak = Enumerable.Repeat<ushort>(ushort.MaxValue, k).ToArray();
            }
        }

        private static int PopCount(int x) => BitOperations.PopCount((uint)x);

        public static Func<CompositionNode, bool> MustCoverMotifResidueCounts(params string[] motifs)
        {
            if (motifs == null || motifs.Length == 0) return _ => true;
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

        public IReadOnlyList<(CompositionNode node, int explainedCount, ushort totalError)>
            RankParentsByCoverage(ushort[] peaks, ushort toleranceScaled, int targetParentLength, int topK = 5) =>
            RankParentsByCoverage(peaks, toleranceScaled, targetParentLength, topK, candidateFilter: null);

        public IReadOnlyList<(CompositionNode node, int explainedCount, ushort totalError)>
            RankParentsByCoverage(ushort[] peaks,
                                  ushort toleranceScaled,
                                  int targetParentLength,
                                  int topK,
                                  Func<CompositionNode, bool> candidateFilter)
        {
            if (peaks == null || peaks.Length == 0) throw new ArgumentException("Provide at least one peak.", nameof(peaks));
            if (toleranceScaled == 0) throw new ArgumentOutOfRangeException(nameof(toleranceScaled), "Tolerance must be > 0 (scaled units).");
            if (!_nodes.Any()) throw new InvalidOperationException("Graph is empty; call Build() first.");

            EnsureMassIndex();

            int K = peaks.Length;
            var coverage = new Dictionary<CompositionNode, Coverage>(Math.Min(_sortedNodes.Length, 8192));
            Coverage GetCov(CompositionNode n)
            {
                if (!coverage.TryGetValue(n, out var c))
                    coverage[n] = c = new Coverage(K);
                return c;
            }

            for (int i = 0; i < K; i++)
            {
                ushort p = peaks[i];
                foreach (var match in WindowByMass(p, toleranceScaled))
                {
                    ushort err = (ushort)Math.Abs(match.Mass - p);
                    PropagateUp(match, i, err, GetCov);
                }
            }

            var candidates = coverage
                .Where(kv => kv.Key.Length == targetParentLength)
                .Where(kv => candidateFilter == null || candidateFilter(kv.Key))
                .Select(kv =>
                {
                    int explained = PopCount(kv.Value.Mask);
                    int totalErr = 0;
                    for (int i = 0; i < K; i++)
                        if (((kv.Value.Mask >> i) & 1) == 1)
                            totalErr += kv.Value.BestErrorPerPeak[i];
                    return (
                        node: kv.Key,
                        explainedCount: explained,
                        totalError: (ushort)Math.Min(totalErr, ushort.MaxValue)
                    );
                })
                .OrderByDescending(r => r.explainedCount)
                .ThenBy(r => r.totalError)
                .ThenBy(r => r.node.Mass)
                .Take(topK)
                .ToList();

            return candidates;

            void PropagateUp(CompositionNode node, int peakIndex, ushort error, Func<CompositionNode, Coverage> covFactory)
            {
                var cov = covFactory(node);
                if (error >= cov.BestErrorPerPeak[peakIndex]) return;
                cov.BestErrorPerPeak[peakIndex] = error;
                cov.Mask |= (1 << peakIndex);
                foreach (var parent in node.Parents)
                    PropagateUp(parent, peakIndex, error, covFactory);
            }
        }

        private void EnsureMassIndex()
        {
            if (!_massIndexDirty && _sortedMasses.Length == _nodes.Count) return;
            var arr = _nodes.Values.ToArray();
            Array.Sort(arr, static (a, b) => a.Mass.CompareTo(b.Mass));
            _sortedNodes = arr;
            _sortedMasses = arr.Select(n => n.Mass).ToArray();
            _massIndexDirty = false;
        }

        private IEnumerable<CompositionNode> WindowByMass(ushort target, ushort tol)
        {
            if (_sortedMasses.Length == 0) yield break;
            int loVal = target - tol < 0 ? 0 : target - tol;
            int hiVal = target + tol;
            int lo = LowerBound(_sortedMasses, (ushort)loVal);
            for (int i = lo; i < _sortedMasses.Length && _sortedMasses[i] <= hiVal; i++)
                yield return _sortedNodes[i];

            static int LowerBound(ushort[] a, ushort value)
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