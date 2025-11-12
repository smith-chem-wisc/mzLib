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
            // Example alphabet (monoisotopic masses in Da)
            var alphabet = new[]
            {
                new CircularPeptideMassGraph.AminoAcid('A', 71.03711),
                new CircularPeptideMassGraph.AminoAcid('C', 103.00919),
                new CircularPeptideMassGraph.AminoAcid('V', 99.06841),
                new CircularPeptideMassGraph.AminoAcid('G', 57.02146),
            };

            var graph = new CircularPeptideMassGraph(alphabet);
            graph.Build(maxLength: 4);

            // Number of compositions with repetition for n=4 residues:
            // Sum over L=1..4 of C(n+L-1, L) = C(4,1)+C(5,2)+C(6,3)+C(7,4) = 4+10+20+35 = 69
            Assert.That(graph.Nodes.Count, Is.EqualTo(69), "Unexpected node count for compositions up to length 4.");

            // Leaf 'A1' must have 4 parents at length 2: A2, A1C1, A1G1, A1V1
            var leafA = graph.Nodes.Single(n => n.Key == "A1");
            var parentKeys = leafA.Parents.Select(p => p.Key).OrderBy(k => k).ToArray();
            Assert.That(parentKeys, Is.EquivalentTo(new[] { "A1C1", "A1G1", "A1V1", "A2" }));

            // Balanced top node 'A1C1G1V1' exists and has no parents (top layer)
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

            // Construct four peaks derived from the balanced composition {A,C,G,V}
            // e.g., fragments: A, C, A+G, C+V
            double a = 71.03711;
            double c = 103.00919;
            double v = 99.06841;
            double g = 57.02146;

            var peaks = new[] { a, c, a + g, c + v };

            var ranked = graph.RankParentsByCoverage(peaks, toleranceDa: 0.01, targetParentLength: 4, topK: 3);

            Assert.That(ranked.Count, Is.GreaterThan(0), "No ranked parents found.");
            var best = ranked[0];

            Assert.Multiple(() =>
            {
                Assert.That(best.node.Key, Is.EqualTo("A1C1G1V1"), "Expected balanced composition to rank first.");
                Assert.That(best.explainedCount, Is.EqualTo(4), "All four peaks should be explained by descendants.");
                Assert.That(best.totalError, Is.LessThanOrEqualTo(4 * 0.01 + 1e-9), "Total error should be consistent with tolerance.");
            });
        }
    }
}