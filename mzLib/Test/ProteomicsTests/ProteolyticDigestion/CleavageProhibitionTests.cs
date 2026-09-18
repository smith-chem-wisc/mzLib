using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// The two rules that are about a subsite being the WRONG thing rather than the right one: a
    /// modification that forbids the cut, and a residue that forbids it.
    /// </summary>
    /// <remarks>
    /// Both come from IMPa, and together they are what separates it from SmE on identical substrate.
    /// IMPa cleaves N-terminal to a glycosylated Ser/Thr, refuses when P1 is aspartate, and refuses when
    /// P1 is itself glycosylated; SmE does none of the refusing. Truth-set IMPA-02, IMPA-10 and SME-02.
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CleavageProhibitionTests
    {
        private static Modification OGlycan(string residue)
        {
            Assert.IsTrue(ModificationMotif.TryGetMotif(residue, out ModificationMotif motif));
            return new Modification(_originalId: "H1N1", _modificationType: "O-linked glycosylation",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373);
        }

        private static List<string> Digest(string sequence, string protease, List<Modification> fixedMods)
        {
            var parameters = new DigestionParams(protease: protease, maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                respectCleavagePromotingModifications: true);

            return new Protein(sequence, "PROHIBIT")
                .Digest(parameters, fixedMods ?? new List<Modification>(), new List<Modification>())
                .Select(p => p.BaseSequence)
                .Distinct()
                .OrderBy(s => s, System.StringComparer.Ordinal)
                .ToList();
        }

        // ---------------------------------------------------------------------------------------
        // A modification that ABOLISHES the cut
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void AForbiddenModificationAtASubsiteAbolishesTheCleavage()
        {
            // PPDATSAAPLR with Thr5 AND Ser6 both glycosylated. IMPa cleaves N-terminal to a glycosylated
            // Ser/Thr, so on motif alone both bonds qualify -- but it "produced no reaction products
            // indicating cleavage between O-glycosylated Thr and Ser even after 24 hours". The bond before
            // Ser6 has a glycosylated P1, and that forbids it.
            List<string> products = Digest("PPDATSAAPLR", "IMPa",
                new List<Modification> { OGlycan("T"), OGlycan("S") });

            CollectionAssert.AreEqual(new[] { "PPDA", "TSAAPLR" }, products,
                "only the bond before Thr5 may be cut; a glycan at P1 forbids the bond before Ser6");
        }

        [Test]
        public static void TheSameSubstrateIsCutTwiceByAnEnzymeWithoutTheProhibition()
        {
            // SmE against IMPa on one substrate, which is the sharpest statement of what the prohibition
            // buys. SmE cleaves N-terminal to EVERY glycosylated Ser/Thr and has no P1 rule at all.
            var glycans = new List<Modification> { OGlycan("T") };

            List<string> bySmE = Digest("AVFTTA", "SmE", glycans);
            List<string> byIMPa = Digest("AVFTTA", "IMPa", glycans);

            CollectionAssert.AreEqual(new[] { "AVF", "T", "TA" }, bySmE,
                "SmE cuts before both glycosylated threonines, leaving Thr4 alone in the middle");
            CollectionAssert.AreEqual(new[] { "AVF", "TTA" }, byIMPa,
                "IMPa is forbidden the second bond, whose P1 is the glycosylated Thr");
            Assert.AreNotEqual(bySmE.Count, byIMPa.Count,
                "the two enzymes must not agree here, or the prohibition is doing nothing");
        }

        [Test]
        public static void AProhibitionIsNotAnObligationToLocalize()
        {
            // The obligation API reports where a glycan MUST be. A forbidden condition says the opposite,
            // and emitting it would force a glycan onto the one residue the enzyme says cannot carry one.
            var protein = new Protein("PPDATSAAPLR", "PROHIBIT");
            Protease impa = ProteaseDictionary.Dictionary["IMPa"];

            var peptide = new PeptideWithSetModifications(protein, null, 5, 11,
                CleavageSpecificity.Full, "test", 0, new Dictionary<int, Modification>(), 0);

            List<int> obligated = peptide.GetCleavageObligatedSites(impa);

            CollectionAssert.AreEqual(new[] { 2 }, obligated,
                "the required glycan at P1' is the peptide's first residue, and the forbidden one at P1 "
                + "must contribute nothing");
        }

        // ---------------------------------------------------------------------------------------
        // A residue that abolishes the cut, with syntax that already existed
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void AspartateAtP1AbolishesTheCleavage()
        {
            // The single flat refusal in the 20-member P1 panel. Expressed with the bracket syntax
            // proteases.tsv already had: for a cut-BEFORE motif the bracketed residues are read BEFORE
            // the recognition sequence, so "[D]|S" is exactly "not when P1 is aspartate".
            var glycan = new List<Modification> { OGlycan("S") };

            CollectionAssert.AreEqual(new[] { "PPDADSAAPLR" }, Digest("PPDADSAAPLR", "IMPa", glycan),
                "P1 is aspartate, so the bond is refused and the substrate stays intact");
            CollectionAssert.AreEqual(new[] { "PPDAA", "SAAPLR" }, Digest("PPDAASAAPLR", "IMPa", glycan),
                "the same substrate with alanine at P1 is cleaved");
        }

        // ---------------------------------------------------------------------------------------
        // The bounds bug the above uncovered, which has nothing to do with glycans
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void APreventingResidueIsHonouredAtTheStartOfTheSequence()
        {
            // SHIPPED trypsin|P, no glycoprotease involved. DigestionMotif.Fits bounds-checked BOTH the
            // forward and the backward index regardless of which one it was about to read, so a
            // cut-after motif at position 0 computed location - 1 = -1, abandoned the prevention, and cut
            // anyway. trypsin|P digested RPAAAAK into "R" + "PAAAAK" -- precisely the R|P cleavage the
            // entry exists to forbid -- while the same motif one residue later behaved correctly.
            CollectionAssert.AreEqual(new[] { "RPAAAAK" }, Digest("RPAAAAK", "trypsin|P", null),
                "proline follows the arginine, so trypsin|P must not cut, even at position 1");

            CollectionAssert.AreEqual(new[] { "AARPAAAK" }, Digest("AARPAAAK", "trypsin|P", null),
                "the same bond away from the terminus was always handled correctly");

            // And the rule must still only block where proline actually follows.
            CollectionAssert.AreEqual(new[] { "AAAAK", "R" }, Digest("RAAAAK", "trypsin|P", null),
                "no proline after the arginine, so the cut is made");
        }

        [Test]
        public static void TheProlineAwareCompositeDoesNotOverDigest()
        {
            // StcE-trypsin|P, added alongside the proline-blind StcE-trypsin exactly as mzLib already
            // ships trypsin beside trypsin|P. Truth-set COMBO-01: StcE cuts before Ser7 because Thr5
            // carries a glycan, trypsin is blocked by the proline after Arg1, and the C-terminal Arg
            // severs nothing.
            List<string> products = Digest("RPPITQSSLR", "StcE-trypsin|P",
                new List<Modification> { OGlycan("T") });

            CollectionAssert.AreEqual(new[] { "RPPITQ", "SSLR" }, products,
                "two products: the glycan-justified StcE cut and nothing else");

            // COMBO-02: no glycan, so StcE cannot fire either and nothing is cut at all.
            CollectionAssert.AreEqual(new[] { "RPPITQSSLR" }, Digest("RPPITQSSLR", "StcE-trypsin|P", null),
                "with no glycan the StcE rule is silent and the proline rule blocks trypsin");
        }
    }
}
