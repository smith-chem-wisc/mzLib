using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Omics.Modifications.IO;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Telling one glycan from another by what it is built from, and being honest about where that stops
    /// working.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class MonosaccharideCompositionTests
    {
        private static Modification Glycan(string composition, string residue = "T")
        {
            Assert.IsTrue(ModificationMotif.TryGetMotif(residue, out ModificationMotif motif));
            return new Modification(_originalId: composition ?? "unknown",
                _modificationType: "O-linked glycosylation", _target: motif,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373,
                _monosaccharideComposition: composition is null ? null : MonosaccharideComposition.Parse(composition));
        }

        private static List<string> DigestWithOpeRATOR(string sequence, Modification glycan)
        {
            var parameters = new DigestionParams(protease: "OpeRATOR", maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                respectCleavagePromotingModifications: true);

            var fixedMods = glycan is null ? new List<Modification>() : new List<Modification> { glycan };
            return new Protein(sequence, "COMP")
                .Digest(parameters, fixedMods, new List<Modification>())
                .Select(p => p.BaseSequence)
                .Distinct()
                .OrderBy(s => s, System.StringComparer.Ordinal)
                .ToList();
        }

        // ---------------------------------------------------------------------------------------
        // Reading a composition
        // ---------------------------------------------------------------------------------------

        [Test]
        [TestCase("HexNAc1", TestName = "Parse a bare monosaccharide and count")]
        [TestCase("HexNAc(1)", TestName = "Parse the parenthesised form")]
        [TestCase("hexnac1", TestName = "Parse is case-insensitive")]
        public static void TheSameCompositionCanBeWrittenSeveralWays(string text)
        {
            MonosaccharideComposition composition = MonosaccharideComposition.Parse(text);
            Assert.AreEqual(1, composition["HexNAc"]);
            Assert.AreEqual(0, composition["Hex"]);
            Assert.AreEqual(1, composition.TotalUnits);
        }

        [Test]
        public static void AMultiSugarCompositionCountsEachSugar()
        {
            // The greedy name match must not read "Hex1HexNAc1" as Hex + NAc.
            MonosaccharideComposition coreTwo = MonosaccharideComposition.Parse("Hex1HexNAc2");
            Assert.AreEqual(1, coreTwo["Hex"]);
            Assert.AreEqual(2, coreTwo["HexNAc"]);
            Assert.AreEqual(3, coreTwo.TotalUnits);
        }

        [Test]
        public static void TheOneLetterDatabaseShorthandIsTheSameVocabulary()
        {
            // The glycan databases write a core 1 glycan as "H1N1", and MetaMorpheus's Glycan.Composition
            // returns exactly that. Reading both spellings is what lets a composition cross from a glycan
            // database into a cleavage rule without a translation step that could drift.
            MonosaccharideComposition shorthand = MonosaccharideComposition.Parse("H1N1");
            MonosaccharideComposition longhand = MonosaccharideComposition.Parse("Hex1HexNAc1");

            Assert.AreEqual(longhand.ToString(), shorthand.ToString(), "the two spellings are one composition");
            Assert.IsTrue(shorthand.IsSupersetOf(longhand) && longhand.IsSupersetOf(shorthand));

            // And the greedy name match must still prefer the long name where one exists: "N" is HexNAc,
            // but "NeuAc1" must not be read as N + euAc.
            Assert.AreEqual(1, MonosaccharideComposition.Parse("NeuAc1")["NeuAc"]);
            Assert.AreEqual(0, MonosaccharideComposition.Parse("NeuAc1")["HexNAc"]);
        }

        [Test]
        [TestCase("Glucose1", TestName = "Reject a monosaccharide that is not a known slot")]
        [TestCase("HexNAc", TestName = "Reject a monosaccharide with no count")]
        [TestCase("1HexNAc", TestName = "Reject a count with no monosaccharide")]
        [TestCase("HexNAc1junk", TestName = "Reject trailing text the parser cannot account for")]
        public static void AMalformedCompositionIsRefusedLoudly(string text)
        {
            // Loud, because a composition that failed to parse and was ignored would silently relax an
            // enzyme's rule -- the exact over-digestion this whole feature exists to stop.
            Assert.Throws<MzLibException>(() => MonosaccharideComposition.Parse(text));
            Assert.IsFalse(MonosaccharideComposition.TryParse(text, out _));
        }

        [Test]
        public static void IncludesIsAFloorOnEveryComponentAtOnce()
        {
            MonosaccharideComposition coreOne = MonosaccharideComposition.Parse("Hex1HexNAc1");

            Assert.IsTrue(MonosaccharideComposition.Parse("Hex1HexNAc1").IsSupersetOf(coreOne), "equal counts qualify");
            Assert.IsTrue(MonosaccharideComposition.Parse("Hex1HexNAc2").IsSupersetOf(coreOne), "core 2 is larger in one component");
            Assert.IsTrue(MonosaccharideComposition.Parse("NeuAc1Hex1HexNAc1").IsSupersetOf(coreOne), "extra sugars do not disqualify");
            Assert.IsFalse(MonosaccharideComposition.Parse("HexNAc1").IsSupersetOf(coreOne), "Tn is short of the Hex");
            Assert.IsFalse(MonosaccharideComposition.Parse("Hex2").IsSupersetOf(coreOne), "more of one sugar cannot pay for another");
            Assert.IsTrue(MonosaccharideComposition.Parse("HexNAc1").IsSupersetOf(null), "no floor is no constraint");
        }

        // ---------------------------------------------------------------------------------------
        // What it buys, end to end
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void OpeRATORTellsTnFromCoreOneFromCoreTwo()
        {
            // The published Table 1 result, in one test. OpeRATOR needs at least core 1 and is blocked
            // again by the second HexNAc that makes core 2.
            const string substrate = "GKPRPYSPRPTSH";

            CollectionAssert.AreEqual(new[] { "GKPRPYSPRPTSH" }, DigestWithOpeRATOR(substrate, Glycan("HexNAc1")),
                "the Tn antigen is a lone HexNAc and does not reach core 1, so there is no cleavage");

            CollectionAssert.AreEqual(new[] { "GKPRPYSPRP", "TSH" }, DigestWithOpeRATOR(substrate, Glycan("Hex1HexNAc1")),
                "core 1 is exactly the floor, so the bond is cut");

            CollectionAssert.AreEqual(new[] { "GKPRPYSPRPTSH" }, DigestWithOpeRATOR(substrate, Glycan("Hex1HexNAc2")),
                "core 2 adds a second HexNAc, which forbids the cleavage again");
        }

        // ---------------------------------------------------------------------------------------
        // Where it stops working, pinned so nobody assumes otherwise
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void CompositionCannotSeparateTwoGlycansThatDifferOnlyInLinkage()
        {
            // alpha-2,6-sialyl core 1 blocks OpeRATOR outright; alpha-2,3 merely slows it. They are the
            // same counts. This asserts the LIMIT, so that a future reader does not conclude composition
            // solved the structure problem -- truth-set OGPA-05 stays a known gap because of exactly this.
            MonosaccharideComposition alpha23 = MonosaccharideComposition.Parse("NeuAc1Hex1HexNAc1");
            MonosaccharideComposition alpha26 = MonosaccharideComposition.Parse("NeuAc1Hex1HexNAc1");

            Assert.IsTrue(alpha23.IsSupersetOf(alpha26) && alpha26.IsSupersetOf(alpha23),
                "the two are indistinguishable by composition, and no rule written against counts can part them");

            CollectionAssert.AreEqual(new[] { "GKPRPYSPRP", "TSH" },
                DigestWithOpeRATOR("GKPRPYSPRPTSH", Glycan("NeuAc1Hex1HexNAc1")),
                "so the sialylated form is cleaved, which is right for alpha-2,3 and wrong for alpha-2,6");
        }

        [Test]
        public static void TnAndOGlcNAcAreAlsoIndistinguishable()
        {
            // The other pair, and the reason STCE-02 stays a gap: alpha-O-GalNAc and beta-O-GlcNAc are
            // both one HexNAc, differing in sugar and anomer.
            Assert.AreEqual(MonosaccharideComposition.Parse("HexNAc1").ToString(),
                MonosaccharideComposition.Parse("hexnac(1)").ToString());
        }

        // ---------------------------------------------------------------------------------------
        // Backward compatibility
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void AGlycanWithNoCompositionDigestsExactlyAsItDidBefore()
        {
            // Nothing populates MonosaccharideComposition yet, so this is the case every existing glycan database
            // is in. Read uniformly as "no composition, so no match", OpeRATOR's own required floor would
            // never be met and the enzyme would silently cleave NOTHING. The unknown case has to fall back
            // to the answer the rule gave before it named a floor.
            CollectionAssert.AreEqual(new[] { "GKPRPYSPRP", "TSH" },
                DigestWithOpeRATOR("GKPRPYSPRPTSH", Glycan(null)),
                "an O-glycan of unknown composition must still satisfy the requirement, as it did before");
        }

        [Test]
        public static void ACompositionSurvivesADatabaseWriteAndRead()
        {
            // BINDING, and the reason this test asserts the field by hand: ProteinDbWriter stores a
            // modification as Modification.ToString() and ProteinDbLoader reads it back through
            // ModificationLoader, so a field missing from either side is silently gone after any
            // write-then-read. Modification.Equals compares only id, type and mass, so the existing
            // Equals-based round-trip test would still pass with the composition dropped -- it is not a
            // regression net for this.
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            var original = new Modification(_originalId: "core1", _modificationType: "O-linked glycosylation",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 365.1322,
                _monosaccharideComposition: MonosaccharideComposition.Parse("Hex1HexNAc1"));

            string record = original.ToString();
            Assert.IsTrue(record.Contains("GC   "), "the record must carry the composition: " + record);

            var reloaded = ModificationLoader
                .ReadModsFromString(record + System.Environment.NewLine + "//", out var errors)
                .First() as Modification;

            Assert.AreEqual(0, errors.Count, "the record this wrote must read back without error");
            Assert.IsNotNull(reloaded.MonosaccharideComposition, "the composition must survive the round trip");
            Assert.AreEqual("Hex(1)HexNAc(1)", reloaded.MonosaccharideComposition.ToString());
            Assert.IsTrue(reloaded.MonosaccharideComposition.IsSupersetOf(MonosaccharideComposition.Parse("Hex1HexNAc1")));
        }

        [Test]
        public static void AModificationWithNoCompositionRoundTripsWithoutTheKey()
        {
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            var plain = new Modification(_originalId: "Phospho", _modificationType: "Common Biological",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);

            string record = plain.ToString();
            Assert.IsFalse(record.Contains("GC   "),
                "a modification with no composition must not emit the key at all");

            var reloaded = ModificationLoader
                .ReadModsFromString(record + System.Environment.NewLine + "//", out var errors)
                .First() as Modification;
            Assert.AreEqual(0, errors.Count);
            Assert.IsNull(reloaded.MonosaccharideComposition);
        }

        [Test]
        public static void AMalformedCompositionInADatabaseIsRefusedRatherThanIgnored()
        {
            // The loader's switch ignores keys it does not know, but a key it DOES know with a value it
            // cannot read must not be skipped -- that would load a glycan whose rule silently does not
            // apply.
            string[] record =
            {
                "ID   bad", "TG   T", "PP   Anywhere.", "MT   test", "MM   100",
                "GC   Glucose1", "//",
            };

            // It THROWS rather than collecting into the warnings list, which is what the CF and BL cases
            // already do for a value they cannot read (ModificationLoader rethrows as MzLibException).
            // The warnings list is for modifications that parse but fail validation; a value that cannot
            // be read at all aborts the file, and that is the established split.
            var thrown = Assert.Throws<MzLibException>(() => ModificationLoader
                .ReadModsFromString(string.Join(System.Environment.NewLine, record), out _)
                .ToList());

            Assert.IsTrue(thrown.Message.Contains("Glucose1"),
                "the message must name the value that could not be read: " + thrown.Message);
        }

        [Test]
        public static void AnOrdinaryModificationCarriesNoCompositionAndIsUnaffected()
        {
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            var phospho = new Modification(_originalId: "Phospho", _modificationType: "Common Biological",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);

            Assert.IsNull(phospho.MonosaccharideComposition,
                "the property is opt-in and every existing modification must be unchanged by it");
        }
    }
}
