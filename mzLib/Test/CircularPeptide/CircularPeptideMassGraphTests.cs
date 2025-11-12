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

            Assert.That(graph.Nodes.Count, Is.EqualTo(69), "Unexpected node count for compositions up to length 4.");

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

            double a = 71.03711, c = 103.00919, v = 99.06841, g = 57.02146;
            var peaks = new[] { a, c, a + g, c + v };

            var ranked = graph.RankParentsByCoverage(peaks, toleranceDa: 0.01, targetParentLength: 4, topK: 3);

            Assert.That(ranked.Count, Is.GreaterThan(0), "No ranked parents found.");
            var best = ranked[0];

            Assert.Multiple(() =>
            {
                Assert.That(best.node.Key, Is.EqualTo("A1C1G1V1"), "Expected balanced composition to rank first.");
                Assert.That(best.explainedCount, Is.EqualTo(4));
                Assert.That(best.totalError, Is.LessThanOrEqualTo(4 * 0.01 + 1e-9));
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

            double a = 71.03711, c = 103.00919, v = 99.06841, g = 57.02146;
            var peaksWithNoise = new[] { a, c, a + g, c + v, 12.34567 };

            var ranked = graph.RankParentsByCoverage(peaksWithNoise, toleranceDa: 0.01, targetParentLength: 4, topK: 3);
            Assert.That(ranked.Count, Is.GreaterThan(0));

            var best = ranked[0];
            Assert.Multiple(() =>
            {
                Assert.That(best.node.Key, Is.EqualTo("A1C1G1V1"));
                Assert.That(best.explainedCount, Is.EqualTo(4));
            });
        }

        // Faster randomized test: 5 trials of 5-mers with 5 distinct residues
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
            const double tolDa = 0.01;

            for (int t = 0; t < trials; t++)
            {
                var permuted = alphabet.OrderBy(_ => rnd.Next()).ToArray();
                var peaks = permuted.Select(a => a.MonoisotopicMass).ToArray();

                var ranked = graph.RankParentsByCoverage(peaks, toleranceDa: tolDa, targetParentLength: 5, topK: 1);
                Assert.That(ranked.Count, Is.GreaterThan(0), "No parent candidates ranked.");
                var best = ranked[0];

                Assert.Multiple(() =>
                {
                    Assert.That(best.node.Key, Is.EqualTo(expectedKey));
                    Assert.That(best.explainedCount, Is.EqualTo(5));
                    Assert.That(best.totalError, Is.LessThanOrEqualTo(5 * tolDa + 1e-9));
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

            // Peaks from balanced A1C1D1E1G1
            var peaks = alphabet.Select(a => a.MonoisotopicMass).ToArray();

            // Motifs require A,G,C,D to be present (E optional). Balanced passes the filter.
            var filter = CircularPeptideMassGraph.MustCoverMotifResidueCounts("AG", "CD");

            var ranked = graph.RankParentsByCoverage(peaks, toleranceDa: 0.01, targetParentLength: 5, topK: 3, candidateFilter: filter);
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

            // Balanced peaks A1C1D1E1G1
            var peaks = alphabet.Select(a => a.MonoisotopicMass).ToArray();

            // Require G>=2 AND the presence of A,C,D,E. At length=5 this is impossible, so zero candidates should pass.
            var filter = CircularPeptideMassGraph.MustCoverMotifResidueCounts("GG", "ACDE");

            var ranked = graph.RankParentsByCoverage(peaks, toleranceDa: 0.01, targetParentLength: 5, topK: 10, candidateFilter: filter);
            Assert.That(ranked.Count, Is.EqualTo(0), "No candidate should pass a G>=2 + ACDE requirement for a balanced 5-mer.");
        }
        [Test, Category("Slow")]
        public void FindSpecificTenMer_WithTwoMotifs_ReducesCandidatesAndFindsParent()
        {
            // Alphabet of 10 distinct residues with unique monoisotopic masses
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

            // IMPORTANT: Build the full graph (do NOT prune at build time).
            // Build-time pruning with motif lower-bounds removes early-length nodes (e.g., length-1 leaves).
            // Since we seed coverage from singleton peaks, we must keep length-1 nodes available.
            graph.Build(maxLength: 10);

            // Target is the 10-mer containing each residue exactly once.
            string expectedKey = string.Concat(alphabet
                .OrderBy(a => a.Symbol)
                .Select(a => $"{a.Symbol}1"));

            // Provide singleton peaks (one per residue).
            var peaks = alphabet.Select(a => a.MonoisotopicMass).ToArray();

            // Apply motif pruning at ranking time (cheap and keeps leaves intact).
            // Two separate 4-mers lower-bound the composition space.
            var motifFilter = CircularPeptideMassGraph.MustCoverMotifResidueCounts("ACDG", "EFHM");

            const double tolDa = 0.01;

            var ranked = graph.RankParentsByCoverage(
                peaks,
                toleranceDa: tolDa,
                targetParentLength: 10,
                topK: 3,
                candidateFilter: motifFilter);

            Assert.That(ranked.Count, Is.GreaterThan(0), "No candidates after motif pruning.");
            var best = ranked[0];

            Assert.Multiple(() =>
            {
                Assert.That(best.node.Key, Is.EqualTo(expectedKey));
                Assert.That(best.explainedCount, Is.EqualTo(10));
                Assert.That(best.totalError, Is.LessThanOrEqualTo(10 * tolDa + 1e-9));
            });
        }
    }
}