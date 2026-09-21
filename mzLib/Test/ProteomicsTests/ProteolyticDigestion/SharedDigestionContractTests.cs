using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// One protein and one RNA walk through the same SearchModeType x FragmentationTerminus matrix. The point is to pin,
    /// in one shared fixture, that both analytes resolve their effective digestion agent, keep the named agent as
    /// <see cref="IDigestionParams.SpecificDigestionAgent"/> for reporting, and return the same kind of product (peptides
    /// vs seeds) for the same request. Per-peptide exhaustiveness is covered by SearchModeTypeDigestionTests and
    /// SemiSpecificDigestionTests for protein and TestDigestion for RNA.
    /// </summary>
    public class SharedDigestionContractTests
    {
        private static readonly Protein Protein =
            new("MAAKCCKDDKEEKFFKGGPEPTIDERSSTTAAK", "contract-protein");

        private static readonly RNA Rna =
            new("GUACUGGUACUG");

        private static readonly IDigestionParams ProteinParams =
            new DigestionParams("trypsin", 2, minPeptideLength: 5, maxPeptideLength: 30);

        private static readonly IDigestionParams RnaParams =
            new RnaDigestionParams("RNase T1", 2, 3);

        private static List<(int Start, int End, CleavageSpecificity Label)> Digest(IBioPolymer polymer, IDigestionParams dp) =>
            polymer.Digest(dp, new List<Modification>(), new List<Modification>())
                .Select(p => (p.OneBasedStartResidue, p.OneBasedEndResidue, p.CleavageSpecificityForFdrCategory))
                .OrderBy(p => p).ToList();

        private static IDigestionParams ParamsFor(IBioPolymer polymer, CleavageSpecificity mode, FragmentationTerminus terminus) =>
            polymer is Protein
                ? new DigestionParams("trypsin", 2, 5, 30, 1024, InitiatorMethionineBehavior.Variable, 2, mode, terminus)
                : new RnaDigestionParams("RNase T1", 2, 3, int.MaxValue, 1024, 2, terminus, mode);

        [TestCaseSource(nameof(MatrixCases))]
        public void Matrix_ReturnsTheRequestedKind(IBioPolymer polymer, IDigestionParams dp, string expectedKind)
        {
            if (expectedKind == "throws")
            {
                Assert.That(() => Digest(polymer, dp), Throws.ArgumentException,
                    "RNA Semi-specific digestion is deliberately limited to an anchored terminus (FivePrime or ThreePrime)");
                return;
            }

            var products = Digest(polymer, dp);
            Assert.That(products, Is.Not.Empty);

            Assert.That(dp.SpecificDigestionAgent, Is.Not.Null);
            Assert.That(dp.SpecificDigestionAgent.Name, Is.EqualTo(polymer is Protein ? "trypsin" : "RNase T1"),
                "the named agent must be preserved for every search mode");

            switch (expectedKind)
            {
                case "full":
                    Assert.That(products.Select(p => p.Label), Is.All.EqualTo(CleavageSpecificity.Full));
                    break;
                case "semi-full":
                    AllLabelsAreFullOrSemi(products);
                    Assert.That(products.Any(p => p.Label == CleavageSpecificity.Semi), Is.True);
                    break;
                case "singleN":
                    Assert.That(products.Select(p => p.Start).Distinct(), Has.Member(1));
                    Assert.That(dp.DigestionAgent.Name, Is.EqualTo("singleN"));
                    break;
                case "singleC":
                    Assert.That(dp.DigestionAgent.Name, Is.EqualTo("singleC"));
                    break;
                default:
                    Assert.Fail($"unexpected matrix case: {expectedKind}");
                    break;
            }
        }

        private static IEnumerable<TestCaseData> MatrixCases()
        {
            foreach (var polymer in new IBioPolymer[] { Protein, Rna })
            {
                yield return new TestCaseData(polymer, ParamsFor(polymer, CleavageSpecificity.Full, FragmentationTerminus.Both), "full");
                yield return new TestCaseData(polymer, ParamsFor(polymer, CleavageSpecificity.Full, FragmentationTerminus.N), "full");
            }
            yield return new TestCaseData(Protein, ParamsFor(Protein, CleavageSpecificity.Semi, FragmentationTerminus.Both), "semi-full");
            yield return new TestCaseData(Rna, ParamsFor(Rna, CleavageSpecificity.Semi, FragmentationTerminus.Both), "throws");
        }

        [Test]
        public void FullDigestion_PeptideCountIsStableAcrossTermini()
        {
            var byBoth = Digest(Protein, ParamsFor(Protein, CleavageSpecificity.Full, FragmentationTerminus.Both));
            var byN = Digest(Protein, ParamsFor(Protein, CleavageSpecificity.Full, FragmentationTerminus.N));
            var byC = Digest(Protein, ParamsFor(Protein, CleavageSpecificity.Full, FragmentationTerminus.C));
            Assert.That(byN.Count, Is.EqualTo(byBoth.Count));
            Assert.That(byC.Count, Is.EqualTo(byBoth.Count));
        }

        [Test]
        public void SemiSemi_BothTermini_IsASupersetOfFull()
        {
            var full = Digest(Protein, ParamsFor(Protein, CleavageSpecificity.Full, FragmentationTerminus.Both))
                .Select(p => (p.Start, p.End)).ToList();
            var semi = Digest(Protein, ParamsFor(Protein, CleavageSpecificity.Semi, FragmentationTerminus.Both))
                .Select(p => (p.Start, p.End)).ToList();
            Assert.That(semi, Is.SupersetOf(full));
        }

        [Test]
        public void Semi_Anchored_Rna_ProducesSeeds()
        {
            var fivePrime = Digest(Rna, ParamsFor(Rna, CleavageSpecificity.Semi, FragmentationTerminus.FivePrime));
            var threePrime = Digest(Rna, ParamsFor(Rna, CleavageSpecificity.Semi, FragmentationTerminus.ThreePrime));

            Assert.That(fivePrime, Is.Not.Empty);
            Assert.That(threePrime, Is.Not.Empty);
            Assert.That(fivePrime.Select(p => p.Start).Distinct(), Has.Member(1),
                "five-prime seeds are anchored at the 5' end");
            Assert.That(threePrime.Select(p => p.End).Distinct(), Has.Member(Rna.Length),
                "three-prime seeds are anchored at the 3' end");
        }

        [Test]
        public void Clone_Protein_EffectiveAgentFollowsTerminus()
        {
            var nParams = new DigestionParams("trypsin", 2, 5, 30, 1024, InitiatorMethionineBehavior.Variable, 2,
                CleavageSpecificity.None, FragmentationTerminus.N);
            var cClone = (DigestionParams)nParams.Clone(FragmentationTerminus.C);

            Assert.That(cClone.SpecificDigestionAgent, Is.EqualTo(nParams.SpecificDigestionAgent));
            Assert.That(cClone.SpecificDigestionAgent.Name, Is.EqualTo("trypsin"));
            Assert.That(cClone.DigestionAgent.Name, Is.EqualTo("singleC"));
            Assert.That(cClone.SearchModeType, Is.EqualTo(CleavageSpecificity.None));
            Assert.That(cClone.FragmentationTerminus, Is.EqualTo(FragmentationTerminus.C));
        }

        [Test]
        public void Clone_SpecificDigestionAgentIsPreserved_EffectiveAgentFollowsTerminus()
        {
            var rnaParams = new RnaDigestionParams("RNase T1", 2, 3, int.MaxValue, 1024, 2,
                FragmentationTerminus.FivePrime, CleavageSpecificity.None);
            var rnaClone = (RnaDigestionParams)rnaParams.Clone(FragmentationTerminus.ThreePrime);

            Assert.That(rnaClone.SpecificDigestionAgent, Is.EqualTo(rnaParams.SpecificDigestionAgent));
            Assert.That(rnaClone.SpecificDigestionAgent.Name, Is.EqualTo("RNase T1"));
            Assert.That(rnaClone.DigestionAgent.Name, Is.EqualTo("singleC"));
            Assert.That(rnaClone.SearchModeType, Is.EqualTo(CleavageSpecificity.None));
            Assert.That(rnaClone.FragmentationTerminus, Is.EqualTo(FragmentationTerminus.ThreePrime));
        }

        /// <summary>
        /// An implementation that does not perform the non-specific agent swap does not have to say so. The member is
        /// defaulted, not abstract, so a type that implements the rest of the interface compiles and reports its own
        /// agent.
        /// </summary>
        /// <remarks>
        /// This test DOES NOT COMPILE against the parent commit, which is the point: <see cref="AgentlessDigestionParams"/>
        /// deliberately omits <c>SpecificDigestionAgent</c>, so before the default existed it failed CS0535 -- the same
        /// error MetaMorpheus's two hand-written doubles fail with, which has reddened mzLib's integration job since
        /// #1318 merged.
        /// </remarks>
        [Test]
        public void AnImplementationThatSwapsNoAgent_NeedNotRestateSpecificDigestionAgent()
        {
            IDigestionParams minimal = new MinimalDigestionParams(ProteaseDictionary.Dictionary["trypsin"]);

            Assert.That(minimal.SpecificDigestionAgent, Is.SameAs(minimal.DigestionAgent),
                "with no swap performed, the named agent and the effective agent are one object");
            Assert.That(minimal.SpecificDigestionAgent.Name, Is.EqualTo("trypsin"));
        }

        /// <summary>
        /// An implementation that DOES swap still has to override, and both of mzLib's own do. Without this, the default
        /// above could be wrong for the case the member exists to describe and nothing here would say so.
        /// </summary>
        [Test]
        public void AnImplementationThatSwapsAgents_StillReportsTheNamedOne()
        {
            var nonSpecific = new DigestionParams("trypsin", searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.N);

            Assert.That(nonSpecific.DigestionAgent.Name, Is.EqualTo("singleN"));
            Assert.That(((IDigestionParams)nonSpecific).SpecificDigestionAgent.Name, Is.EqualTo("trypsin"),
                "the swap is exactly the case the member exists for, so the default must not reach it");
        }

        /// <summary>
        /// Stands in for the minimal hand-written doubles MetaMorpheus keeps -- ParameterTest.BadDigestionParams and
        /// DigestionAgentNameTests.AgentlessDigestionParams. It implements every member of IDigestionParams EXCEPT
        /// SpecificDigestionAgent, on purpose.
        /// </summary>
        private class MinimalDigestionParams : IDigestionParams
        {
            private readonly DigestionAgent _agent;

            internal MinimalDigestionParams(DigestionAgent agent) => _agent = agent;

            public int MaxMissedCleavages { get; set; }
            public int MinLength { get; set; }
            public int MaxLength { get; set; }
            public int MaxModificationIsoforms { get; set; }
            public int MaxMods { get; set; }
            public DigestionAgent DigestionAgent => _agent;
            public FragmentationTerminus FragmentationTerminus => FragmentationTerminus.Both;
            public CleavageSpecificity SearchModeType => CleavageSpecificity.Full;
            public IDigestionParams Clone(FragmentationTerminus? newTerminus = null) => this;
            public bool Equals(IDigestionParams other) => ReferenceEquals(this, other);
        }

        private static void AllLabelsAreFullOrSemi(List<(int Start, int End, CleavageSpecificity Label)> products)
        {
            foreach (var product in products)
            {
                Assert.That(product.Label, Is.AnyOf(CleavageSpecificity.Full, CleavageSpecificity.Semi));
            }
        }
    }
}