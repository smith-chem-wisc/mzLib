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
    /// A preventing-cleavage rule must be honoured everywhere in the sequence, including at the very
    /// ends.
    /// </summary>
    /// <remarks>
    /// <para><see cref="DigestionMotif.Fits"/> bounds-checked BOTH the forward index
    /// (<c>location + m + n</c>) and the backward one
    /// (<c>location - PreventingCleavage.Length + n</c>) before deciding which one it was about to read.
    /// Whenever either was out of range the prevention was silently abandoned and the bond was cut.</para>
    ///
    /// <para>For the shipped proline entries -- all of them cut-AFTER motifs written <c>X[P]|</c> -- that
    /// means the rule never applied at sequence position 0, because the backward index came out at -1.
    /// <c>trypsin|P</c> therefore digested <c>RPAAAAK</c> into <c>R</c> + <c>PAAAAK</c>, making exactly the
    /// R|P cleavage the entry exists to forbid, while behaving correctly one residue later. The same
    /// applies to <c>chymotrypsin|P</c>, <c>elastase|P</c>, <c>Lys-C|P</c> and <c>subtilisin|P</c>.</para>
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class PreventingCleavageBoundsTests
    {
        private static List<string> Digest(string sequence, string protease)
        {
            var parameters = new DigestionParams(protease: protease, maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            return new Protein(sequence, "BOUNDS")
                .Digest(parameters, new List<Modification>(), new List<Modification>())
                .Select(p => p.BaseSequence)
                .Distinct()
                .OrderBy(s => s, System.StringComparer.Ordinal)
                .ToList();
        }

        [Test]
        public static void ACutAfterMotifHonoursItsPreventingResidueAtPositionOne()
        {
            // The bug, in shipped configuration and with no test fixture involved.
            CollectionAssert.AreEqual(new[] { "RPAAAAK" }, Digest("RPAAAAK", "trypsin|P"),
                "proline follows the arginine, so trypsin|P must not cut -- position in the sequence "
                + "cannot change whether a rule applies");
        }

        [Test]
        public static void TheSameBondAwayFromTheTerminusWasAlwaysCorrect()
        {
            // The control that localises the fault: identical bond, two residues further along.
            CollectionAssert.AreEqual(new[] { "AARPAAAK" }, Digest("AARPAAAK", "trypsin|P"));
        }

        [Test]
        public static void ThePreventingRuleStillOnlyBlocksWhereItShould()
        {
            // The fix must not turn the prevention into a blanket refusal.
            CollectionAssert.AreEqual(new[] { "AAAAK", "R" }, Digest("RAAAAK", "trypsin|P"),
                "no proline after the arginine, so the cut is made even at position 1");

            // And the proline-blind entry is unaffected, since it has no preventing rule at all.
            CollectionAssert.AreEqual(new[] { "PAAAAK", "R" }, Digest("RPAAAAK", "trypsin"),
                "trypsin is proline-blind by design and must keep cutting here");
        }

        [Test]
        [TestCase("chymotrypsin|P", "FPAAAAK", TestName = "chymotrypsin_P honours its proline rule at position one")]
        // elastase|P recognises almost every residue, so the rest of the sequence uses ones it does not
        // (W, M) -- otherwise the test would be measuring its other sites rather than this rule.
        [TestCase("elastase|P", "IPWWWM", TestName = "elastase_P honours its proline rule at position one")]
        [TestCase("Lys-C|P", "KPAAAAR", TestName = "Lys-C_P honours its proline rule at position one")]
        public static void EveryShippedProlineEntryIsAffectedTheSameWay(string protease, string sequence)
        {
            CollectionAssert.AreEqual(new[] { sequence }, Digest(sequence, protease),
                protease + " must not cut before the proline at position 2");
        }

        [Test]
        public static void ACutBeforeMotifHonoursItsPreventingResidueAtTheLastResidue()
        {
            // The mirror fault. No shipped protease is a cut-BEFORE motif with a preventing rule -- every
            // "|P" entry is cut-after -- so this builds one. For CutIndex 0 the bracketed residues are
            // read BEFORE the recognition sequence, and it was the FORWARD index that went out of range
            // at the C-terminus and abandoned the rule.
            var motifs = DigestionMotif.ParseDigestionMotifsFromString("[D]|T");
            var protease = new Protease("bounds-cut-before", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary["bounds-cut-before"] = protease;

            CollectionAssert.AreEqual(new[] { "AAADT" }, Digest("AAADT", "bounds-cut-before"),
                "the residue before the final Thr is the forbidden Asp, so the bond must not be cut");

            CollectionAssert.AreEqual(new[] { "AAA", "TA" }, Digest("AAATA", "bounds-cut-before"),
                "with alanine before it the same bond is cut, so the rule is still doing its job");
        }
    }
}
