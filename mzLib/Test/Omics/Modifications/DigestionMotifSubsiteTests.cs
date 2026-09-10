using System.Collections.Generic;
using NUnit.Framework;
using Omics.Digestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.Omics.Modifications
{
    /// <summary>
    /// A digestion motif's sidedness lives in its CutIndex, and the codebase has historically read that
    /// by comparing it against a literal -- "CutIndex == 1 means the protease cuts C-terminal",
    /// "CutIndex == 0 means N-terminal". That reading has a hole: a motif may STRADDLE the bond, and
    /// four such motifs ship today (collagenase "GPX|GPX" at CutIndex 3, StcE-trypsin "TX|T" at 2,
    /// colicin_E5 "G|U" at 1, plus the "SX|T" family members). They are C-terminal cutters that the
    /// "== 1" test misses entirely.
    ///
    /// These tests pin CleavageSide's three-way classification against the real shipped motifs, pin the
    /// Schechter-Berger subsite accessors that replace index arithmetic at the call site, and pin that
    /// CleavesCTerminalTo still answers exactly what it answered before it was re-expressed as the P1
    /// question.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class DigestionMotifSubsiteTests
    {
        private const string AllResidues = "ACDEFGHIKLMNPQRSTVWY";

        private static DigestionMotif Motif(string motifString)
        {
            List<DigestionMotif> parsed = DigestionMotif.ParseDigestionMotifsFromString(motifString);
            Assert.AreEqual(1, parsed.Count, $"test setup: '{motifString}' should parse to one motif");
            return parsed[0];
        }

        // --- Side: the three regimes, against motifs that actually ship ---

        [Test]
        [TestCase("K|", CleavageSide.CTerminal)]          // trypsin, Lys-C
        [TestCase("E|", CleavageSide.CTerminal)]          // Glu-C
        [TestCase("|D", CleavageSide.NTerminal)]          // Asp-N
        [TestCase("|K", CleavageSide.NTerminal)]          // Lys-N
        [TestCase("|M", CleavageSide.NTerminal)]          // CNBr_N
        [TestCase("|U", CleavageSide.NTerminal)]          // RNase_MC1
        [TestCase("GPX|GPX", CleavageSide.Straddling)]    // collagenase
        [TestCase("TX|T", CleavageSide.Straddling)]       // StcE-trypsin
        [TestCase("SX|S", CleavageSide.Straddling)]       // StcE-trypsin
        [TestCase("G|U", CleavageSide.Straddling)]        // colicin_E5
        public static void Side_ClassifiesShippedMotifs(string motifString, CleavageSide expected)
        {
            Assert.AreEqual(expected, Motif(motifString).Side);
        }

        [Test]
        public static void Side_StraddlingMotifs_AreWhatTheMagicNumberTestsMiss()
        {
            // The regression guard. Sidedness-by-literal classifies with "CutIndex == 1" for C-terminal
            // and "CutIndex == 0" for N-terminal; these two motifs match NEITHER, which is precisely how
            // they fall through both arms and get classified as neither side.
            DigestionMotif collagenase = Motif("GPX|GPX");
            DigestionMotif stcETrypsin = Motif("TX|T");

            Assert.AreEqual(3, collagenase.CutIndex);
            Assert.AreEqual(2, stcETrypsin.CutIndex);

            foreach (DigestionMotif straddling in new[] { collagenase, stcETrypsin })
            {
                Assert.AreNotEqual(0, straddling.CutIndex);
                Assert.AreNotEqual(1, straddling.CutIndex);
                Assert.AreEqual(CleavageSide.Straddling, straddling.Side);

                // ...and yet they do sever a bond after a residue they name, so they are genuinely
                // C-terminal cutters. Side and CleavesCTerminalTo answer different questions.
                Assert.IsTrue(straddling.CleavesCTerminalTo('T') || straddling.CleavesCTerminalTo('X'));
            }
        }

        // --- Subsite accessors ---

        [Test]
        public static void Subsites_Trypsin_NamesP1Only()
        {
            DigestionMotif trypsin = Motif("K|");
            Assert.AreEqual('K', trypsin.NonPrimeSubsite(1));
            Assert.AreEqual('\0', trypsin.NonPrimeSubsite(2));
            Assert.AreEqual('\0', trypsin.PrimeSubsite(1), "trypsin's motif says nothing about P1'");
        }

        [Test]
        public static void Subsites_AspN_NamesP1PrimeOnly()
        {
            DigestionMotif aspN = Motif("|D");
            Assert.AreEqual('\0', aspN.NonPrimeSubsite(1));
            Assert.AreEqual('D', aspN.PrimeSubsite(1));
            Assert.AreEqual('\0', aspN.PrimeSubsite(2));
        }

        [Test]
        public static void Subsites_StcETrypsin_NamesBothSides_WithTheGlycosylatableResidueAtP2()
        {
            // "TX|T" -- the literature places StcE's glycan requirement at P2, and P2 is exactly where
            // this motif names a Thr. That correspondence is the reason subsite addressing is worth
            // having: the rule can be stated at the subsite the enzymology names.
            DigestionMotif stcE = Motif("TX|T");
            Assert.AreEqual('X', stcE.NonPrimeSubsite(1));
            Assert.AreEqual('T', stcE.NonPrimeSubsite(2));
            Assert.AreEqual('\0', stcE.NonPrimeSubsite(3));
            Assert.AreEqual('T', stcE.PrimeSubsite(1));
            Assert.AreEqual('\0', stcE.PrimeSubsite(2));
        }

        [Test]
        public static void Subsites_Collagenase_WalksBothDirectionsAcrossTheBond()
        {
            DigestionMotif collagenase = Motif("GPX|GPX");
            Assert.AreEqual('X', collagenase.NonPrimeSubsite(1));
            Assert.AreEqual('P', collagenase.NonPrimeSubsite(2));
            Assert.AreEqual('G', collagenase.NonPrimeSubsite(3));
            Assert.AreEqual('\0', collagenase.NonPrimeSubsite(4));
            Assert.AreEqual('G', collagenase.PrimeSubsite(1));
            Assert.AreEqual('P', collagenase.PrimeSubsite(2));
            Assert.AreEqual('X', collagenase.PrimeSubsite(3));
            Assert.AreEqual('\0', collagenase.PrimeSubsite(4));
        }

        [Test]
        public static void Subsites_ColicinE5_OneResidueEitherSide()
        {
            DigestionMotif colicin = Motif("G|U");
            Assert.AreEqual('G', colicin.NonPrimeSubsite(1));
            Assert.AreEqual('U', colicin.PrimeSubsite(1));
        }

        [Test]
        [TestCase(0)]
        [TestCase(-1)]
        public static void Subsites_PositionBelowOne_IsNotASubsite(int position)
        {
            // Guard against off-by-one misuse: P0 does not exist, and without this guard
            // NonPrimeSubsite(0) would silently return the P1' residue.
            DigestionMotif stcE = Motif("TX|T");
            Assert.AreEqual('\0', stcE.NonPrimeSubsite(position));
            Assert.AreEqual('\0', stcE.PrimeSubsite(position));
        }

        // --- Accepts: ambiguity codes honoured, unconstrained subsites never match ---

        [Test]
        public static void SubsiteAccepts_Wildcard_AcceptsEveryResidue()
        {
            DigestionMotif stcE = Motif("TX|T");
            foreach (char residue in AllResidues)
            {
                Assert.IsTrue(stcE.NonPrimeSubsiteAccepts(1, residue),
                    $"P1 of TX|T is the wildcard X, so it should accept {residue}");
            }
        }

        [Test]
        [TestCase('B', 'D', true)]
        [TestCase('B', 'N', true)]
        [TestCase('B', 'E', false)]
        [TestCase('J', 'I', true)]
        [TestCase('J', 'L', true)]
        [TestCase('Z', 'E', true)]
        [TestCase('Z', 'Q', true)]
        [TestCase('Z', 'D', false)]
        public static void SubsiteAccepts_HonoursAmbiguityCodes(char motifChar, char residue, bool expected)
        {
            // The subsite matcher must go through the same matcher digestion uses, or a rule keyed on a
            // subsite would disagree with the site list the protease actually produces.
            DigestionMotif ambiguous = Motif($"{motifChar}|");
            Assert.AreEqual(expected, ambiguous.NonPrimeSubsiteAccepts(1, residue));
        }

        [Test]
        public static void SubsiteAccepts_UnconstrainedSubsite_NeverMatches()
        {
            // Trypsin constrains nothing at P1', so asking whether P1' demands alanine is false -- not
            // vacuously true. A rule that requires something at an unconstrained subsite must fail
            // closed.
            DigestionMotif trypsin = Motif("K|");
            foreach (char residue in AllResidues)
            {
                Assert.IsFalse(trypsin.PrimeSubsiteAccepts(1, residue));
            }

            Assert.IsFalse(trypsin.NonPrimeSubsiteAccepts(2, 'K'));
            Assert.IsFalse(trypsin.NonPrimeSubsiteAccepts(1, '\0'),
                "the null character is the absence of a residue and must not match a real motif char");
        }

        // --- CleavesCTerminalTo: behaviour preserved across the re-expression ---

        [Test]
        public static void CleavesCTerminalTo_IsExactlyTheP1Question_ForEveryShippedMotif()
        {
            // The behaviour-preservation evidence for re-expressing CleavesCTerminalTo in terms of the
            // new accessors. Checked over every motif shape and every residue, not a sample.
            foreach (string motifString in new[] { "K|", "R|", "E|", "|D", "|K", "|M", "GPX|GPX", "TX|T", "G|U", "X|" })
            {
                DigestionMotif motif = Motif(motifString);
                foreach (char residue in AllResidues)
                {
                    Assert.AreEqual(motif.NonPrimeSubsiteAccepts(1, residue), motif.CleavesCTerminalTo(residue),
                        $"{motifString} disagreed at {residue}");
                }
            }
        }

        [Test]
        public static void CleavesCTerminalTo_KeepsItsDocumentedAnswers()
        {
            Assert.IsTrue(Motif("K|").CleavesCTerminalTo('K'));
            Assert.IsFalse(Motif("K|").CleavesCTerminalTo('R'));

            // An N-terminal cutter severs nothing after its recognition residue, so it reports false
            // for every residue including its own.
            foreach (char residue in AllResidues)
            {
                Assert.IsFalse(Motif("|D").CleavesCTerminalTo(residue));
            }

            // The wildcard reports true for every residue, as its own remarks state.
            foreach (char residue in AllResidues)
            {
                Assert.IsTrue(Motif("X|").CleavesCTerminalTo(residue));
            }
        }

        [Test]
        public static void Side_AgreesWithCleavesCTerminalTo_OnWhetherAnythingIsSeveredAfterAResidue()
        {
            // An N-terminal motif can never sever a bond after a residue it names; the other two sides
            // always can. This is the invariant that makes Side safe to branch on.
            foreach (string motifString in new[] { "K|", "R|", "E|", "|D", "|K", "|M", "GPX|GPX", "TX|T", "G|U", "X|" })
            {
                DigestionMotif motif = Motif(motifString);
                bool severedSomethingAfterAResidue = false;
                foreach (char residue in AllResidues)
                {
                    severedSomethingAfterAResidue |= motif.CleavesCTerminalTo(residue);
                }

                Assert.AreEqual(motif.Side != CleavageSide.NTerminal, severedSomethingAfterAResidue,
                    $"{motifString} ({motif.Side})");
            }
        }
    }
}
