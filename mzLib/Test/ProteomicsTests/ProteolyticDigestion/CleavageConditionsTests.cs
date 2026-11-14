using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
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
    /// <see cref="DigestionProduct.GetCleavageConditions"/>: what a product's ends need, or must not carry, if the protease
    /// cut them there. Each condition is keyed by the parent's one-based residue, which may lie on the neighbouring product.
    /// </summary>
    [TestFixture]
    public static class CleavageConditionsTests
    {
        private static PeptideWithSetModifications Product(Protein protein, string protease, int start, int end) =>
            new(protein, new DigestionParams(protease), start, end, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

        // "5" is a glycan required on residue 5, "!4" one forbidden on residue 4.
        private static List<string> Conditions(DigestionProduct product, string protease) =>
            product.GetCleavageConditions(ProteaseDictionary.Dictionary[protease])
                .Select(c => (c.Requirement.IsForbidden ? "!" : "") + c.OneBasedResidue).OrderBy(s => s).ToList();

        /// <summary>RPPITQSSL: StcE cuts T5-Q6|S7 (TX|S), so its P2 is T5.</summary>
        [Test]
        public static void StcE_EndCut_RequiresP2InsideTheProduct() =>
            Assert.That(Conditions(Product(new Protein("RPPITQSSL", "P1"), "StcE", 1, 6), "StcE"), Is.EqualTo(new[] { "5" }));

        [Test]
        public static void StcE_StartCut_RequiresP2OnTheNeighbouringProduct() =>
            Assert.That(Conditions(Product(new Protein("RPPITQSSL", "P1"), "StcE", 7, 9), "StcE"), Is.EqualTo(new[] { "5" }));

        /// <summary>GGGGTGGGG: OpeRATOR cuts before T5 (|T), so its P1' is T5.</summary>
        [Test]
        public static void OpeRATOR_StartCut_RequiresTheFirstResidue_AndEndCut_TheNextOne()
        {
            var protein = new Protein("GGGGTGGGG", "P1");
            Assert.That(Conditions(Product(protein, "OpeRATOR", 5, 9), "OpeRATOR"), Is.EqualTo(new[] { "5" }));
            Assert.That(Conditions(Product(protein, "OpeRATOR", 1, 4), "OpeRATOR"), Is.EqualTo(new[] { "5" }));
        }

        /// <summary>GGGATGGG: IMPa ([D]|T) cuts before T5, needing a glycan on T5 and none on A4.</summary>
        [Test]
        public static void IMPa_RequiresP1Prime_AndForbidsP1_OnBothSidesOfItsBond()
        {
            var protein = new Protein("GGGATGGG", "P1");
            Assert.That(Conditions(Product(protein, "IMPa", 5, 8), "IMPa"), Is.EqualTo(new[] { "!4", "5" }));
            Assert.That(Conditions(Product(protein, "IMPa", 1, 4), "IMPa"), Is.EqualTo(new[] { "!4", "5" }));
        }

        /// <summary>GGGDTGGG: IMPa's [D]|T does not cut before T5 when P1 is aspartate, so nothing is imposed there.</summary>
        [Test]
        public static void AResidueThatPreventsTheCut_LeavesNoConditions() =>
            Assert.That(Conditions(Product(new Protein("GGGDTGGG", "P1"), "IMPa", 5, 8), "IMPa"), Is.Empty);

        [Test]
        public static void ProteinEnds_AnInitiatorMet_AndATruncationBoundary_AreNotCuts()
        {
            Assert.That(Conditions(Product(new Protein("TGGGGTGGG", "P1"), "OpeRATOR", 1, 9), "OpeRATOR"), Is.Empty, "the protein's ends");
            Assert.That(Conditions(Product(new Protein("MTGGGG", "P1"), "OpeRATOR", 2, 6), "OpeRATOR"), Is.Empty, "after the initiator Met");
            var truncated = new Protein("GGGGTGGGG", "P1", proteolysisProducts: new List<TruncationProduct> { new(5, 9, "chain") });
            Assert.That(Conditions(Product(truncated, "OpeRATOR", 5, 9), "OpeRATOR"), Is.Empty, "a truncation boundary");
            var unknownStart = new Protein("GGGGTGGGG", "P1", proteolysisProducts: new List<TruncationProduct> { new(null, 4, "propeptide") });
            Assert.That(Conditions(Product(unknownStart, "OpeRATOR", 5, 9), "OpeRATOR"), Is.Empty, "the end of a product whose start is unknown");
        }

        /// <summary>GTKSGG: T2-K3|S4 is a StcE cut (TX|S) and a trypsin cut (K|), and trypsin needs no glycan.</summary>
        [Test]
        public static void ACutARequirementFreeMotifAlsoMakes_ImposesNothing()
        {
            Assert.That(Conditions(Product(new Protein("GTKSGG", "P1"), "StcE-trypsin", 4, 6), "StcE-trypsin"), Is.Empty);
            Assert.That(Conditions(Product(new Protein("GGGKSGGG", "P1"), "OpeRATOR-trypsin|P", 5, 8), "OpeRATOR-trypsin|P"), Is.Empty);
        }

        /// <summary>
        /// Two motifs that both cut AAG|TAA, |T needing a glycan on T4 (P1'): a condition holds only if every motif cutting
        /// there imposes it, so G| needing P1 instead leaves none, and G| needing P1' as well leaves T4.
        /// </summary>
        [Test]
        public static void WhereTwoMotifsCut_OnlyWhatBothImposeIsACondition()
        {
            var p1Prime = new CleavageRequirement(1, isPrimeSide: true);
            var p1 = new CleavageRequirement(1, isPrimeSide: false);
            Protease Agent(params CleavageRequirement[] onG) => new("two-motif-test", CleavageSpecificity.Full, null, null,
                new List<DigestionMotif> { new("T", null, 0, null, new[] { p1Prime }), new("G", null, 1, null, onG) });
            var product = Product(new Protein("AAGTAA", "P1"), "trypsin", 4, 6); // the agent below decides, not these params

            Assert.That(product.GetCleavageConditions(Agent(p1)), Is.Empty, "the two motifs need different residues");
            Assert.That(product.GetCleavageConditions(Agent(p1Prime, p1)).Select(c => c.OneBasedResidue), Is.EqualTo(new[] { 4 }),
                "both need T4");
        }

        /// <summary>PITQ in RPPITQSSL: StcE cut its C-terminus but not its N-terminus, which a semi-specific search allows.</summary>
        [Test]
        public static void ASemiSpecificProduct_ReportsOnlyTheEndTheProteaseCut() =>
            Assert.That(Conditions(Product(new Protein("RPPITQSSL", "P1"), "StcE", 3, 6), "StcE"), Is.EqualTo(new[] { "5" }));

        [Test]
        public static void AProductWithoutAProtein_HasNoConditions() =>
            Assert.That(new PeptideWithSetModifications("TQSSL").GetCleavageConditions(ProteaseDictionary.Dictionary["StcE"]), Is.Empty);

        /// <summary>SQSSL from residue 2: StcE's TX|S needs two residues before the bond, so it cannot have cut after S1.</summary>
        [Test]
        public static void StcE_CannotHaveCutAfterTheFirstResidue() =>
            Assert.That(Conditions(Product(new Protein("SQSSL", "P1"), "StcE", 2, 5), "StcE"), Is.Empty);

        /// <summary>
        /// GTGGT, cut before each T, with conditions on P3 and P3': TGG's start cut puts P3 before the protein and its end cut
        /// puts P3' after it, so only the start cut's P3' (residue 4) and the end cut's P3 (residue 2) remain.
        /// </summary>
        [Test]
        public static void AConditionOnAResidueBeyondTheProtein_IsDropped()
        {
            var agent = new Protease("p3-test", CleavageSpecificity.Full, null, null, new List<DigestionMotif>
                { new("T", null, 0, null, new[] { new CleavageRequirement(3, isPrimeSide: false), new CleavageRequirement(3, isPrimeSide: true) }) });
            var product = Product(new Protein("GTGGT", "P1"), "trypsin", 2, 4); // the agent below decides, not these params
            Assert.That(product.GetCleavageConditions(agent).Select(c => c.OneBasedResidue).OrderBy(r => r), Is.EqualTo(new[] { 2, 4 }));
        }

        [Test]
        public static void AnRnaseProduct_HasNoConditions()
        {
            var oligos = new RNA("GUACUGCCUCUAGUGAAGCA").Digest(new RnaDigestionParams("RNase T1", maxMissedCleavages: 1, minLength: 1),
                new List<Modification>(), new List<Modification>()).ToList();
            Assert.That(oligos, Is.Not.Empty, "premise");
            Assert.That(oligos.SelectMany(o => o.GetCleavageConditions(RnaseDictionary.Dictionary["RNase T1"])), Is.Empty);
        }

        /// <summary>
        /// Over a semi-specific digest of leukosialin (P16150, residues 1-120), the required residue each product carries
        /// is the one the glyco-search experiment hard-coded (MetaMorpheus exp/pd-child-precursor,
        /// ExperimentalObligatedKey), whose benefit was measured on 26 files. That experiment used OpeRATOR with plain K|,R|;
        /// the shipped K[P]|,R[P]| differs only before a proline, which is never OpeRATOR's S/T.
        /// </summary>
        [Test]
        public static void OverAMucinDigest_TheRequiredResidueInsideEachProduct_IsTheMeasuredRule([Values("StcE-trypsin", "OpeRATOR-trypsin|P")] string protease)
        {
            const string leukosialin = "MATLLLLLGVLVVSPDALGSTTAVQTPTSGEPLVSTSEPLSSKMYTTSITSDPKADSTGD" +
                                       "QTSALPPSTSINEGSPLWTSIGASTGSPLPEPTTYQEVSIKMSSVPQETPHATSHPAVPI";
            var products = new Protein(leukosialin, "P16150").Digest(new DigestionParams(protease, maxMissedCleavages: 20, minPeptideLength: 5,
                    maxPeptideLength: 40, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.Both),
                new List<Modification>(), new List<Modification>()).ToList();

            int obligated = 0;
            foreach (var product in products)
            {
                int? expected = MeasuredRule(protease, leukosialin, product.OneBasedStartResidue, product.OneBasedEndResidue);
                obligated += expected.HasValue ? 1 : 0;
                var inside = product.GetCleavageConditions(ProteaseDictionary.Dictionary[protease])
                    .Where(c => !c.Requirement.IsForbidden && c.OneBasedResidue >= product.OneBasedStartResidue && c.OneBasedResidue <= product.OneBasedEndResidue)
                    .Select(c => c.OneBasedResidue).ToList();
                Assert.That(inside, Is.EqualTo(expected.HasValue ? new[] { expected.Value } : new int[0]),
                    $"{product.BaseSequence} [{product.OneBasedStartResidue}-{product.OneBasedEndResidue}]");
            }
            Assert.That(obligated, Is.GreaterThan(10), "premise: the digest has products the rule applies to");
        }

        // The experiment's rule, as parent residues (ExperimentalObligatedKey returned peptide keys).
        private static int? MeasuredRule(string protease, string q, int start, int end)
        {
            if (protease.StartsWith("StcE"))
            {
                if (end < 2 || end >= q.Length || !"ST".Contains(q[end - 2]) || !"ST".Contains(q[end]) || q[end - 1] is 'K' or 'R')
                    return null;
                return end - 1;
            }
            if (start < 2 || (start == 2 && q[0] == 'M') || !"ST".Contains(q[start - 1]) || q[start - 2] is 'K' or 'R')
                return null;
            return start;
        }
    }
}
