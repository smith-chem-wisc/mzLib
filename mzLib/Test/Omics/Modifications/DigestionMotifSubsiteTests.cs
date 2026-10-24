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

                // ...and yet each names a P1 subsite, so it severs a bond after a residue and is
                // genuinely a C-terminal cutter. Both P1s are the wildcard, so the literal
                // CleavesCTerminalTo reports false for them; the subsite accessors are the question here.
                Assert.AreNotEqual('\0', straddling.NonPrimeSubsite(1));
                Assert.IsTrue(straddling.NonPrimeSubsiteAccepts(1, 'T'));
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
        [TestCase("|D", 1, 'D', true)]
        [TestCase("|D", 1, 'E', false)]
        [TestCase("TX|T", 1, 'T', true)]
        [TestCase("TX|T", 1, 'S', false)]
        [TestCase("GPX|GPX", 1, 'G', true)]
        [TestCase("GPX|GPX", 2, 'P', true)]
        [TestCase("GPX|GPX", 2, 'G', false)]
        [TestCase("GPX|GPX", 3, 'W', true)]   // the wildcard sits at P3'
        [TestCase("|X{P}", 1, 'P', false)]   // an excluded wildcard at P1'
        [TestCase("|X{P}", 1, 'A', true)]
        [TestCase("|B", 1, 'N', true)]
        [TestCase("|B", 1, 'E', false)]
        [TestCase("|Z", 1, 'Q', true)]
        [TestCase("|Z", 1, 'D', false)]
        public static void PrimeSubsiteAccepts_MatchesTheResidueTheMotifNames(string motifString, int position, char residue, bool expected)
        {
            Assert.AreEqual(expected, Motif(motifString).PrimeSubsiteAccepts(position, residue));
        }

        [Test]
        [TestCase("TX|T", 2, 'T', true)]      // StcE's glycosylatable P2
        [TestCase("TX|T", 2, 'S', false)]
        [TestCase("SX|S", 2, 'S', true)]
        [TestCase("GPX|GPX", 2, 'P', true)]
        [TestCase("GPX|GPX", 3, 'G', true)]
        [TestCase("GPX|GPX", 3, 'P', false)]
        [TestCase("RX{P}|", 2, 'R', true)]
        [TestCase("RX{P}|", 2, 'K', false)]
        public static void NonPrimeSubsiteAccepts_BeyondP1_MatchesTheResidueTheMotifNames(string motifString, int position, char residue, bool expected)
        {
            Assert.AreEqual(expected, Motif(motifString).NonPrimeSubsiteAccepts(position, residue));
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

        /// <summary>
        /// CleavesCTerminalTo as master wrote it before it was re-expressed through NonPrimeSubsite:
        /// bounds check, then compare InducingCleavage[CutIndex - 1] with the residue LITERALLY. The
        /// wildcard and the ambiguity codes deliberately do not match (mzLib #1119), so the shared
        /// matcher is not consulted. Restated here rather than reached through the motif, so that an
        /// error in the subsite arithmetic cannot move both sides of the comparison together.
        /// </summary>
        private static bool CleavesCTerminalToAsOriginallyWritten(DigestionMotif motif, char residue)
        {
            if (motif.CutIndex < 1 || motif.CutIndex > motif.InducingCleavage.Length)
            {
                return false;
            }

            return motif.InducingCleavage[motif.CutIndex - 1] == residue;
        }

        [Test]
        public static void CleavesCTerminalTo_AgreesWithTheOriginalFormula_ForEveryMotifShape()
        {
            // The behaviour-preservation evidence for re-expressing CleavesCTerminalTo as the P1
            // subsite, checked against the formula it replaced over every residue. The shapes cover
            // each side, each straddling offset, the wildcard with and without an exclusion, and the
            // ambiguity codes -- the last two are where reading P1 through the matcher instead of
            // literally would disagree.
            foreach (string motifString in new[]
                     {
                         "K|", "R|", "E|", "|D", "|K", "|M", "GPX|GPX", "|X{P}", "TX|T", "G|U",
                         "X|", "X|U", "U|A", "RX{P}|", "B|", "|Z",
                     })
            {
                DigestionMotif motif = Motif(motifString);
                foreach (char residue in AllResidues + "U")
                {
                    Assert.AreEqual(CleavesCTerminalToAsOriginallyWritten(motif, residue), motif.CleavesCTerminalTo(residue),
                        $"{motifString} disagreed with the original formula at {residue}");
                }
            }
        }

        [Test]
        public static void CleavesCTerminalTo_EmptyMotif_CleavesAfterNothing()
        {
            // The original bounds check rejected every CutIndex against an empty recognition sequence.
            DigestionMotif empty = new DigestionMotif(string.Empty, null, 0, null);
            foreach (char residue in AllResidues)
            {
                Assert.IsFalse(empty.CleavesCTerminalTo(residue));
                Assert.IsFalse(empty.NonPrimeSubsiteAccepts(1, residue));
                Assert.IsFalse(empty.PrimeSubsiteAccepts(1, residue));
            }
        }

        [Test]
        [TestCase("|D")]
        [TestCase("K|")]
        [TestCase("TX|T")]
        public static void CleavesCTerminalTo_TheNullCharacter_IsNoResidue(string motifString)
        {
            // NonPrimeSubsite answers '\0' wherever the motif names nothing at P1 -- every N-terminal
            // motif -- so a bare comparison would report "|D" as cleaving after the null character.
            Assert.IsFalse(Motif(motifString).CleavesCTerminalTo('\0'));
        }

        [Test]
        public static void NoRecognitionSequence_NamesNoSubsite_AndDoesNotThrow()
        {
            // The public constructor accepts a null recognition sequence. A motif that names no residue
            // cannot straddle a bond, so the cut index alone decides its side, and every subsite question
            // answers "nothing" rather than throwing.
            DigestionMotif motif = new DigestionMotif(null, null, 1, null);
            Assert.AreEqual(CleavageSide.CTerminal, motif.Side);
            Assert.AreEqual('\0', motif.NonPrimeSubsite(1));
            Assert.AreEqual('\0', motif.PrimeSubsite(1));
            foreach (char residue in AllResidues)
            {
                Assert.IsFalse(motif.CleavesCTerminalTo(residue));
                Assert.IsFalse(motif.NonPrimeSubsiteAccepts(1, residue));
                Assert.IsFalse(motif.PrimeSubsiteAccepts(1, residue));
            }
        }

        [Test]
        [TestCase("GPX|GPX", 'A', false)]  // P1 is the wildcard, which names no residue
        [TestCase("TX|T", 'T', false)]     // P1 is the wildcard; the T is at P2
        [TestCase("RX{P}|", 'P', false)]
        [TestCase("RX{P}|", 'R', false)]   // R is P2, not P1
        [TestCase("B|", 'D', false)]       // an ambiguity code at P1 is not a literal residue
        [TestCase("G|U", 'G', true)]
        [TestCase("G|U", 'U', false)]      // U is P1', not P1
        [TestCase("U|A", 'U', true)]
        public static void CleavesCTerminalTo_TabulatedP1Answers(string motifString, char residue, bool expected)
        {
            Assert.AreEqual(expected, Motif(motifString).CleavesCTerminalTo(residue));
        }

        [Test]
        [TestCase("GPX|GPX", 'A')]
        [TestCase("TX|T", 'T')]
        [TestCase("RX{P}|", 'R')]
        [TestCase("B|", 'D')]
        public static void CleavesCTerminalTo_IsLiteral_WhereTheP1SubsiteMatcherAccepts(string motifString, char residue)
        {
            // The two P1 questions differ exactly here, and the difference is deliberate: the subsite
            // matcher says the motif tolerates this residue at P1, but the residue's side chain does
            // not direct the cut, so a modification on it abolishes nothing (mzLib #1119).
            DigestionMotif motif = Motif(motifString);
            Assert.IsTrue(motif.NonPrimeSubsiteAccepts(1, residue));
            Assert.IsFalse(motif.CleavesCTerminalTo(residue));
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

            // The wildcard names no residue literally, so it reports false for every residue, as its
            // remarks state.
            foreach (char residue in AllResidues)
            {
                Assert.IsFalse(Motif("X|").CleavesCTerminalTo(residue));
            }
        }

        [Test]
        public static void Side_AgreesWithTheP1Subsite_OnWhetherAnythingIsSeveredAfterAResidue()
        {
            // An N-terminal motif can never sever a bond after a residue it names; the other two sides
            // always can. This is the invariant that makes Side safe to branch on. It is asked of the
            // P1 subsite matcher rather than of CleavesCTerminalTo, which is literal and so reports
            // nothing for a wildcard P1 ("GPX|GPX", "TX|T", "X|").
            foreach (string motifString in new[] { "K|", "R|", "E|", "|D", "|K", "|M", "GPX|GPX", "TX|T", "G|U", "X|" })
            {
                DigestionMotif motif = Motif(motifString);
                bool severedSomethingAfterAResidue = false;
                foreach (char residue in AllResidues)
                {
                    severedSomethingAfterAResidue |= motif.NonPrimeSubsiteAccepts(1, residue);
                    if (motif.Side == CleavageSide.NTerminal)
                    {
                        Assert.IsFalse(motif.CleavesCTerminalTo(residue), $"{motifString} at {residue}");
                    }
                }

                Assert.AreEqual(motif.Side != CleavageSide.NTerminal, severedSomethingAfterAResidue,
                    $"{motifString} ({motif.Side})");
            }
        }
    }
}
