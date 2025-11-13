using System;
using System.Linq;
using NUnit.Framework;
using Proteomics.CircularBiopolymer;

namespace Test.CircularPeptide
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class CircularPeptideMassGraphTests
    {
        private static ushort S(double da) => (ushort)Math.Round(da * CircularPeptideMassGraph.MassScale);
        private static ushort MinTol(double da)
        {
            // Scale and ensure at least 1 (to avoid zero-tolerance and ambiguous Math.Max overloads)
            var scaled = S(da);
            return scaled == 0 ? (ushort)1 : scaled;
        }

        [Test]
        public void Build_CompositionGraph_ForFourResidues_UptoLen4()
        {
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('V', 99.06841),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
            };

            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 4);

            Assert.That(graph.Nodes.Count, Is.EqualTo(69));
            var leafA = graph.Nodes.Single(n => n.Key == "A1");
            var parentKeys = leafA.Parents.Select(p => p.Key).OrderBy(k => k).ToArray();
            Assert.That(parentKeys, Is.EquivalentTo(new[] { "A1C1", "A1G1", "A1V1", "A2" }));
            var top = graph.Nodes.Single(n => n.Key == "A1C1G1V1");
            Assert.Multiple(() =>
            {
                Assert.That(top.Length, Is.EqualTo(4));
                Assert.That(top.Parents.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void RankParents_FromFourPeaks_PrefersBalancedParent()
        {
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('V', 99.06841),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
            };
            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 4);

            ushort a = S(71.03711), c = S(103.00919), v = S(99.06841), g = S(57.02146);
            var peaks = new[] { a, c, (ushort)(a + g), (ushort)(c + v) };
            ushort tol = MinTol(0.05); // was (ushort)Math.Max(1, S(0.05)) causing ambiguity

            var ranked = graph.RankParentsByCoverage(peaks, tol, 4, topK: 3);
            Assert.That(ranked.Count, Is.GreaterThan(0));
            var best = ranked[0];
            Assert.Multiple(() =>
            {
                Assert.That(best.node.Key, Is.EqualTo("A1C1G1V1"));
                Assert.That(best.explainedCount, Is.EqualTo(4));
                Assert.That(best.totalError, Is.LessThanOrEqualTo((ushort)(4 * tol)));
            });
        }

        [Test]
        public void RankParents_IgnoresNoisePeaks_ByDefault()
        {
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('V', 99.06841),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
            };
            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 4);

            ushort a = S(71.03711), c = S(103.00919), v = S(99.06841), g = S(57.02146);
            var peaks = new[] { a, c, (ushort)(a + g), (ushort)(c + v), S(12.34567) };
            ushort tol = MinTol(0.05);

            var ranked = graph.RankParentsByCoverage(peaks, tol, 4, topK: 3);
            Assert.That(ranked.Count, Is.GreaterThan(0));
            var best = ranked[0];
            Assert.Multiple(() =>
            {
                Assert.That(best.node.Key, Is.EqualTo("A1C1G1V1"));
                Assert.That(best.explainedCount, Is.EqualTo(4));
            });
        }

        [Test]
        public void RecoverFiveDistinctFiveMers_RandomizedSingletonPeaks()
        {
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('D', 115.02694),
                new CircularPeptideMassGraph.AminoAcid('E', 129.04259),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
            };
            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 5);

            string expectedKey = string.Concat(alphabet.OrderBy(a => a.Symbol).Select(a => $"{a.Symbol}1"));
            var rnd = new Random(12345);
            const int trials = 5;
            ushort tol = MinTol(0.05);

            for (int t = 0; t < trials; t++)
            {
                var permuted = alphabet.OrderBy(_ => rnd.Next()).ToArray();
                var peaks = permuted.Select(a => a.ScaledMonoisotopicMass).ToArray();
                var ranked = graph.RankParentsByCoverage(peaks, tol, 5, topK: 1);
                Assert.That(ranked.Count, Is.GreaterThan(0));
                var best = ranked[0];
                Assert.Multiple(() =>
                {
                    Assert.That(best.node.Key, Is.EqualTo(expectedKey));
                    Assert.That(best.explainedCount, Is.EqualTo(5));
                    Assert.That(best.totalError, Is.LessThanOrEqualTo((ushort)(5 * tol)));
                });
            }
        }

        [Test]
        public void RankParents_WithMotifCompositionFilter_PrunesCandidates_Positive()
        {
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('D', 115.02694),
                new CircularPeptideMassGraph.AminoAcid('E', 129.04259),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
            };
            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 5);

            var peaks = alphabet.Select(a => a.ScaledMonoisotopicMass).ToArray();
            var filter = CircularPeptideMassGraph.MustCoverMotifResidueCounts("AG", "CD");
            ushort tol = MinTol(0.05);

            var ranked = graph.RankParentsByCoverage(peaks, tol, 5, topK: 3, candidateFilter: filter);
            Assert.That(ranked.Count, Is.GreaterThan(0));
            Assert.That(ranked[0].node.Key, Is.EqualTo("A1C1D1E1G1"));
        }

        [Test]
        public void RankParents_WithMotifCompositionFilter_RejectsImpossible()
        {
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('D', 115.02694),
                new CircularPeptideMassGraph.AminoAcid('E', 129.04259),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
            };
            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 5);

            var peaks = alphabet.Select(a => a.ScaledMonoisotopicMass).ToArray();
            var filter = CircularPeptideMassGraph.MustCoverMotifResidueCounts("GG", "ACDE");
            ushort tol = MinTol(0.05);

            var ranked = graph.RankParentsByCoverage(peaks, tol, 5, topK: 10, candidateFilter: filter);
            Assert.That(ranked.Count, Is.EqualTo(0));
        }

        [Test, Category("Slow")]
        public void FindSpecificTenMer_WithTwoMotifs_ReducesCandidatesAndFindsParent()
        {
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('D', 115.02694),
                new CircularPeptideMassGraph.AminoAcid('E', 129.04259),
                new CircularPeptideMassGraph.AminoAcid('F', 147.06841),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
                new CircularPeptideMassGraph.AminoAcid('H', 137.05891),
                new CircularPeptideMassGraph.AminoAcid('K', 128.09496),
                new CircularPeptideMassGraph.AminoAcid('M', 131.04049),
                new CircularPeptideMassGraph.AminoAcid('Y', 163.06333),
            };
            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 10);

            string expectedKey = string.Concat(alphabet.OrderBy(a => a.Symbol).Select(a => $"{a.Symbol}1"));
            var peaks = alphabet.Select(a => a.ScaledMonoisotopicMass).ToArray();
            var motifFilter = CircularPeptideMassGraph.MustCoverMotifResidueCounts("ACDG", "EFHM");
            ushort tol = MinTol(0.05);

            var ranked = graph.RankParentsByCoverage(peaks, tol, 10, topK: 3, candidateFilter: motifFilter);
            Assert.That(ranked.Count, Is.GreaterThan(0));
            var best = ranked[0];
            Assert.Multiple(() =>
            {
                Assert.That(best.node.Key, Is.EqualTo(expectedKey));
                Assert.That(best.explainedCount, Is.EqualTo(10));
                Assert.That(best.totalError, Is.LessThanOrEqualTo((ushort)(10 * tol)));
            });
        }
    }
}