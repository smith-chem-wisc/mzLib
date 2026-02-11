using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestDecoyProteinGenerator
    {
        private static SequenceVariation ReverseOne(Protein protein, SequenceVariation sv)
        {
            // Call private static ReverseSequenceVariations via reflection
            var method = typeof(DecoyProteinGenerator).GetMethod(
                "ReverseSequenceVariations",
                BindingFlags.Static | BindingFlags.NonPublic);

            Assert.That(method, Is.Not.Null, "Could not reflect ReverseSequenceVariations");

            string reversedSequence = new string(protein.BaseSequence.Reverse().ToArray());
            var result = (List<SequenceVariation>)method.Invoke(
                null,
                new object[] { new List<SequenceVariation> { sv }, protein, reversedSequence, "TEST" });

            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            return result[0];
        }

        private static Dictionary<int, List<Modification>> Mods(params (int pos, string id)[] pairs)
        {
            var dict = new Dictionary<int, List<Modification>>();
            foreach (var (pos, id) in pairs)
                dict[pos] = new List<Modification> { new Modification(_originalId: id) };
            return dict;
        }

        /// <summary>
        /// CRITICAL: Tests the stop-gain modification mapping branch (DecoyProteinGenerator.cs lines 199-202).
        /// Stop-gain variants (ending with '*') require special handling to preserve decoy length consistency.
        /// Without this test, stop-gain modifications could be incorrectly mapped, causing decoy peptides
        /// to have different lengths than their target counterparts, breaking FDR calculations.
        /// Verifies: (1) original key is preserved, (2) default reverse mapping also applies.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_StopsGain_KeepsOriginalKey_AndAddsDefaultReverse()
        {
            // startsWithM = false; stopGain = true
            var protein = new Protein("APEPTIDE", "P1"); // length=8, no leading M
            // Point edit at pos 3; E -> E* (stop gained). Variant has a mod at variant-index 3.
            var sv = new SequenceVariation(
                oneBasedBeginPosition: 3,
                oneBasedEndPosition: 3,
                originalSequence: "E",
                variantSequence: "E*",
                description: "stop",
                oneBasedModifications: Mods((3, "m3")));

            var decoy = ReverseOne(protein, sv);

            // For stop-gain:
            // - First add keeps key (3)
            // - Then default branch adds (variantLen = 8 + 2 - 1 = 9) ? 9 - 3 + 1 = 7
            Assert.That(decoy.OneBasedModifications.ContainsKey(3), "stopGain should keep original key");
            Assert.That(decoy.OneBasedModifications.ContainsKey(7), "default reverse mapping for stopGain");
            Assert.That(decoy.OneBasedModifications.Count, Is.EqualTo(2));
        }

        /// <summary>
        /// CRITICAL: Tests the methionine-retention modification mapping branch (lines 204-207).
        /// When the protein starts with 'M', the first residue is retained during reversal to maintain
        /// biological plausibility (initiating methionine). Modifications at positions > 1 must use
        /// the formula: variantSeqLength - key + 2. Without this test, modifications on proteins with
        /// initiating methionine could be mapped to wrong positions, corrupting decoy modification sites.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_StartsWithM_KeyGreaterThanOne_UsesPlusTwoFormula()
        {
            // startsWithM = true; kvp.Key > 1
            var protein = new Protein("MPEPTIDE", "P2"); // length=8, leading M
            // Substitution at 4..4; mod at variant-index 3
            var sv = new SequenceVariation(
                4, 4, "P", "V", "sub", oneBasedModifications: Mods((3, "m3")));

            var decoy = ReverseOne(protein, sv);

            // variantSeqLength = 8 + 1 - 1 = 8 ? key = 8 - 3 + 2 = 7
            Assert.That(decoy.OneBasedModifications.Keys.Single(), Is.EqualTo(7));
        }

        /// <summary>
        /// CRITICAL: Tests the "modification on starting methionine" branch (lines 209-212).
        /// When a variant introduces a methionine at position 1, modifications at that position
        /// must remain at position 1 (since M is retained at the N-terminus during reversal).
        /// This ensures N-terminal modifications on initiating methionines are correctly preserved
        /// in decoy proteins, which is essential for matching N-terminal modification searches.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_ModOnStartingMethionine_VariantStartsWithM_MapsToOne()
        {
            // kvp.Key == 1 && variant starts with 'M'
            var protein = new Protein("APEPTIDE", "P3"); // no leading M required
            var sv = new SequenceVariation(
                1, 1, "A", "M", "gainM", oneBasedModifications: Mods((1, "mN")));

            var decoy = ReverseOne(protein, sv);

            Assert.That(decoy.OneBasedModifications.Keys.Single(), Is.EqualTo(1));
        }

        /// <summary>
        /// CRITICAL: Tests the "modification on starting non-methionine" branch (lines 214-217).
        /// When the variant does NOT start with 'M' but has a modification at position 1, the
        /// modification should map to the protein's C-terminus (protein.Length) after reversal.
        /// This is essential for correct placement of what was originally an N-terminal modification
        /// when the sequence is reversed without methionine retention.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_ModOnStart_NonMethionine_MapsToProteinLength()
        {
            // kvp.Key == 1 && variant does NOT start with 'M'
            var protein = new Protein("APEPTIDE", "P4"); // length = 8
            var sv = new SequenceVariation(
                1, 1, "A", "V", "nonMStart", oneBasedModifications: Mods((1, "mN")));

            var decoy = ReverseOne(protein, sv);

            Assert.That(decoy.OneBasedModifications.Keys.Single(), Is.EqualTo(protein.BaseSequence.Length));
        }

        /// <summary>
        /// CRITICAL: Tests the default modification mapping branch (lines 218-221).
        /// This is the most common case: protein doesn't start with M, modification key > 1.
        /// Uses formula: variantSeqLength - key + 1. This branch handles the majority of
        /// modification mappings and must be correct for standard decoy generation to work.
        /// An error here would affect most decoy proteins in a typical database search.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_DefaultElse_UsesPlusOneFormula()
        {
            // startsWithM = false; kvp.Key > 1; not stopGain
            var protein = new Protein("APEPTIDE", "P5"); // length=8, no leading M
            // Substitution; mod at variant-index 3
            var sv = new SequenceVariation(
                4, 4, "P", "V", "default", oneBasedModifications: Mods((3, "m3")));

            var decoy = ReverseOne(protein, sv);

            // variantSeqLength = 8 + 1 - 1 = 8 ? key = 8 - 3 + 1 = 6
            Assert.That(decoy.OneBasedModifications.Keys.Single(), Is.EqualTo(6));
        }

        /// <summary>
        /// CRITICAL: Tests the stop-gain sequence reversal branch (lines 237-243).
        /// Stop-gain variants (premature stop codons) require special handling: positions are
        /// preserved and the variant sequence is rotated (not simply reversed). This ensures
        /// decoy proteins with stop-gain variants maintain the same effective length as targets,
        /// which is crucial for accurate mass calculations and peptide identification.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_StopGain_BuildsSegmentAndRotatedVariant()
        {
            // Protein does NOT start with M (startsWithM == false)
            var protein = new Protein("APEPTIDE", "P_stop"); // len = 8
            // Point variant at pos=3; original E -> variant E* (stop-gain)
            var sv = new SequenceVariation(3, 3, "E", "E*", "stopgain", (System.Collections.Generic.Dictionary<int, System.Collections.Generic.List<Modification>>)null);

            var decoy = ReverseOne(protein, sv);

            // For stopGain branch:
            // - begin == original begin
            // - end   == original end
            // - variant is rotation of reversed array which yields original "E*" back for this 2-char example
            Assert.That(decoy.OneBasedBeginPosition, Is.EqualTo(3));
            Assert.That(decoy.OneBasedEndPosition, Is.EqualTo(3));
            Assert.That(decoy.VariantSequence, Is.EqualTo("E*"));

            // ReverseOne() passes decoyIdentifier = "TEST"
            StringAssert.StartsWith("TEST VARIANT:", decoy.Description);
        }

        /// <summary>
        /// CRITICAL: Tests the start-loss sequence reversal branch (lines 245-248).
        /// Start-loss occurs when the original sequence begins with 'M' but the variant does not
        /// (e.g., frameshift removing initiating methionine). The variant must be placed at the
        /// C-terminus of the decoy. This test ensures correct handling of variants that affect
        /// protein translation initiation, which are biologically significant mutations.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_StartLoss_VariantAtEndOfDecoy()
        {
            // startLoss: original starts with M at pos 1, variant does NOT start with M
            var protein = new Protein("APEPTIDE", "P_startloss"); // len=8
            // Original "MA" (1..2) -> Variant "A"
            var sv = new SequenceVariation(1, 2, "MA", "A", "startloss", (System.Collections.Generic.Dictionary<int, System.Collections.Generic.List<Modification>>)null);

            var decoy = ReverseOne(protein, sv);

            // begin = L - end + 2 = 8 - 2 + 2 = 8, end = L = 8
            // original = reverse("MA").Substring(0, len-1) = "AM" -> "A"
            // variant  = reverse("A") = "A"
            Assert.That(decoy.OneBasedBeginPosition, Is.EqualTo(8));
            Assert.That(decoy.OneBasedEndPosition, Is.EqualTo(8));
            Assert.That(decoy.OriginalSequence, Is.EqualTo("A"));
            Assert.That(decoy.VariantSequence, Is.EqualTo("A"));
        }

        /// <summary>
        /// CRITICAL: Tests the "both start with M, but longer" branch (lines 250-255).
        /// When both original and variant sequences start with 'M' and have length > 1, the
        /// reversed sequences must trim the last character (which was the 'M' before reversal)
        /// to avoid duplicating the retained N-terminal methionine. Without this test, decoys
        /// could have incorrect sequence lengths or duplicated methionines.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_BothStartWithM_ButLonger_TrimsLastChar()
        {
            // both start with M, but length > 1 (hits the 'both start with M, but there’s more' branch)
            var protein = new Protein("APEPTIDE", "P_bothM"); // len=8
            // Original "MA" (1..2) -> Variant "MI"
            var sv = new SequenceVariation(1, 2, "MA", "MI", "bothM", (System.Collections.Generic.Dictionary<int, System.Collections.Generic.List<Modification>>)null);

            var decoy = ReverseOne(protein, sv);

            // begin = 8 - 2 + 2 = 8, end = 8
            // original = reverse("MA") = "AM" -> trim last => "A"
            // variant  = reverse("MI") = "IM" -> trim last => "I"
            Assert.That(decoy.OneBasedBeginPosition, Is.EqualTo(8));
            Assert.That(decoy.OneBasedEndPosition, Is.EqualTo(8));
            Assert.That(decoy.OriginalSequence, Is.EqualTo("A"));
            Assert.That(decoy.VariantSequence, Is.EqualTo("I"));
        }

        /// <summary>
        /// CRITICAL: Tests the "gained initiating methionine" branch (lines 257-260).
        /// When a variant introduces a new initiating methionine at position 1 (single residue),
        /// it must stay at position 1 in the decoy. This handles the case of variants that add
        /// an alternative translation start site. Incorrect handling would misplace important
        /// N-terminal variants that affect protein expression and localization.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_GainedInitiatingMethionine_AtPos1()
        {
            // gained an initiating methionine (variant starts with M, begin==1) single-length case
            var protein = new Protein("APEPTIDE", "P_gainM");
            var sv = new SequenceVariation(1, 1, "A", "M", "gainM", (System.Collections.Generic.Dictionary<int, System.Collections.Generic.List<Modification>>)null);

            var decoy = ReverseOne(protein, sv);

            // begin=end=1; sequences reversed character-wise (single-char => unchanged)
            Assert.That(decoy.OneBasedBeginPosition, Is.EqualTo(1));
            Assert.That(decoy.OneBasedEndPosition, Is.EqualTo(1));
            Assert.That(decoy.OriginalSequence, Is.EqualTo("A"));
            Assert.That(decoy.VariantSequence, Is.EqualTo("M"));
        }

        /// <summary>
        /// CRITICAL: Tests the "protein starts with M, variation elsewhere" branch (lines 262-265).
        /// When the protein has an initiating methionine but the variant is at a different position,
        /// the reversal must account for methionine retention using the +2 offset formula.
        /// This is the common case for missense mutations in proteins with standard initiation,
        /// making it essential for accurate decoy generation in most proteomics experiments.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_ProteinStartsWithM_NoVariationOnStart()
        {
            // Protein starts with M => hits 'startsWithM' branch if earlier conditions don't apply
            var protein = new Protein("MPEPTIDE", "P_startsM"); // len=8, startsWithM=true
            // Simple substitution away from start: 4..4 P -> V
            var sv = new SequenceVariation(4, 4, "P", "V", "middle", (System.Collections.Generic.Dictionary<int, System.Collections.Generic.List<Modification>>)null);

            var decoy = ReverseOne(protein, sv);

            // begin = L - end + 2 = 8 - 4 + 2 = 6
            // end   = L - begin + 2 = 6
            // sequences are reversed char-wise (single-char => unchanged)
            Assert.That(decoy.OneBasedBeginPosition, Is.EqualTo(6));
            Assert.That(decoy.OneBasedEndPosition, Is.EqualTo(6));
            Assert.That(decoy.OriginalSequence, Is.EqualTo("P"));
            Assert.That(decoy.VariantSequence, Is.EqualTo("V"));
        }

        /// <summary>
        /// CRITICAL: Tests the default sequence reversal branch (lines 267-270).
        /// When the protein does NOT start with 'M' (e.g., after signal peptide cleavage or
        /// non-canonical start), standard reversal with +1 offset is used. This is the fallback
        /// case covering proteins without retained initiating methionine and must correctly
        /// reverse variant positions to maintain search accuracy.
        /// </summary>
        [Test]
        public void ReverseSequenceVariations_NoStartingMethionine_DefaultElseBranch()
        {
            // Protein does NOT start with M => final else branch (no special start handling)
            var protein = new Protein("APEPTIDE", "P_else"); // len=8
            // Simple substitution: 4..4 P -> V
            var sv = new SequenceVariation(4, 4, "P", "V", "middle", (System.Collections.Generic.Dictionary<int, System.Collections.Generic.List<Modification>>)null);

            var decoy = ReverseOne(protein, sv);

            // begin = L - end + 1 = 8 - 4 + 1 = 5
            // end   = L - begin + 1 = 5
            Assert.That(decoy.OneBasedBeginPosition, Is.EqualTo(5));
            Assert.That(decoy.OneBasedEndPosition, Is.EqualTo(5));
            Assert.That(decoy.OriginalSequence, Is.EqualTo("P"));
            Assert.That(decoy.VariantSequence, Is.EqualTo("V"));
        }
    }
}