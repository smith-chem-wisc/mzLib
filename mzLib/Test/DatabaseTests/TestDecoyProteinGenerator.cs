using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
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
    }
}