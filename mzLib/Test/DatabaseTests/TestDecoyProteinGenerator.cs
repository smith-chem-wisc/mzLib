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