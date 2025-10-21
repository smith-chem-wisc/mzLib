using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Test.DatabaseTests.VariantTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceVariationSplitPerGenotypeZygosityBranchTests
    {
        private static SequenceVariation Make(string vcf, Dictionary<int, List<Modification>> mods = null) =>
            new SequenceVariation(
                oneBasedPosition: 25,
                originalSequence: "M",
                variantSequence: "K",
                description: "ZygoVar",
                variantCallFormatDataString: vcf,
                oneBasedModifications: mods);

        // Helper: assert single variant with expected Mode substring
        private static void AssertSingleMode(List<SequenceVariation> list, string modeContains)
        {
            Assert.That(list.Count, Is.EqualTo(1), "Expected exactly one variant.");
            Assert.That(list[0].Description.Contains(modeContains), Is.True, $"Mode tag '{modeContains}' missing.");
        }

        [Test]
        public void ZygosityAlreadyPresent_KeyExists_NoRecalc()
        {
            // Heterozygous 0/1; storedAltIndex = -1 (ANN=.) so Mode=HeterozygousAlt
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=.\tGT:AD:DP\t0/1:5,6:11";
            var sv = Make(vcf);
            // Key "0" should already exist
            Assert.That(sv.VariantCallFormatData.ZygosityBySample.ContainsKey("0"), Is.True);
            var split = sv.SplitPerGenotype();
            AssertSingleMode(split, "HeterozygousAlt");
        }

        [Test]
        public void ZygosityFallback_Recomputed_AfterRemovingEntry()
        {
            // Remove zygosity entry to force fallback path
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=.\tGT:AD:DP\t0/1:4,5:9";
            var sv = Make(vcf);
            sv.VariantCallFormatData.ZygosityBySample.Remove("0");
            var split = sv.SplitPerGenotype();
            AssertSingleMode(split, "HeterozygousAlt");
        }

        [Test]
        public void ZygosityFallback_Unknown_NoCalledAlleles_Skipped()
        {
            // GT ./.
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=.\tGT:AD:DP\t./.:3,0:3";
            var sv = Make(vcf);
            // Remove key so fallback occurs, producing Unknown then numericAlleles empty => continue
            sv.VariantCallFormatData.ZygosityBySample.Remove("0");
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void ParseError_SkipsSample()
        {
            // Non-numeric allele token 'X' => parseError => continue
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/X:5,5:10";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void AllReference_AllRefTrue_NoVariantAdded()
        {
            // 0/0 homozygous reference; even with emitReferenceForHomozygousRef true, no-op variant invalid -> none returned
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/0:8,0:8";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(emitReferenceForHomozygousRef: true);
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void HomozygousAlt_allStoredAltPath()
        {
            // ANN allele = T => storedAltIndex = 1; genotype 1/1 => HomozygousAlt
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t1/1:0,9:9";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            AssertSingleMode(split, "HomozygousAlt");
        }

        [Test]
        public void MixedAltIndex_SkipDueToFlag()
        {
            // ALT T,G; storedAltIndex=1; genotype 0/2 containsDifferentAlt and skipIfAltIndexMismatch default true => skipped
            string vcf = "1\t200\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:4,0,5:9";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void MixedAltIndex_AddedWhenFlagFalse()
        {
            string vcf = "1\t200\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:3,0,4:7";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(skipIfAltIndexMismatch: false);
            AssertSingleMode(split, "MixedAltIndex(StoredAltOnly)");
        }

        [Test]
        public void Heterozygous_WithIncludeReference_AttemptsRefAndAddsAlt()
        {
            // includeReferenceForHeterozygous true requests HeterozygousRef (no-op dropped) + HeterozygousAlt (kept)
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,6:11";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(includeReferenceForHeterozygous: true);
            AssertSingleMode(split, "HeterozygousAlt");
            Assert.That(split[0].Description.Contains("HeterozygousRef"), Is.False);
        }

        [Test]
        public void CloneMods_CreatesIndependentDictionary()
        {
            var mods = new Dictionary<int, List<Modification>>
            {
                { 25, new List<Modification>{ new Modification(_originalId:"ModA", _modificationType:"TestType") } }
            };
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:4,5:9";
            var sv = Make(vcf, mods);
            var split = sv.SplitPerGenotype();
            Assert.That(split.Count, Is.EqualTo(1));
            Assert.That(split[0].OneBasedModifications, Is.Not.Null);
            Assert.That(split[0].OneBasedModifications.Count, Is.EqualTo(1));
            Assert.That(ReferenceEquals(split[0].OneBasedModifications, sv.OneBasedModifications), Is.False,
                "Expected cloned modification dictionary, not original reference.");
        }

        [Test]
        public void AlleleIndexZero_NoAllStoredAltBranch()
        {
            // ANN allele = REF (A) => storedAltIndex=0; genotype 1/1 but allStoredAlt false => heterozygous path yields HeterozygousAlt
            string vcf = "1\t200\t.\tA\tT\t.\tPASS\tANN=A|.\tGT:AD:DP\t1/1:0,10:10";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            AssertSingleMode(split, "HeterozygousAlt");
        }

        [Test]
        public void ContainingDifferentAlt_NoSkipWhenFlagFalse()
        {
            // ContainsDifferentAlt true (0/2 with storedAltIndex=1) skipIfAltIndexMismatch false => MixedAltIndex variant
            string vcf = "1\t200\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:2,0,6:8";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(skipIfAltIndexMismatch: false);
            AssertSingleMode(split, "MixedAltIndex(StoredAltOnly)");
        }
    }
}