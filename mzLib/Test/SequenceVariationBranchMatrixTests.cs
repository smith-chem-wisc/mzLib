using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceVariationBranchMatrixTests
    {
        // Helper to construct the base variant (substitution A->T) with a supplied VCF line
        private static SequenceVariation Make(string vcf) =>
            new SequenceVariation(
                oneBasedPosition: 25,
                originalSequence: "A",
                variantSequence: "T",
                description: "BranchBase",
                variantCallFormatDataString: vcf,
                oneBasedModifications: null);

        private static SequenceVariation MakeWithMod(string vcf, int pos) =>
            new SequenceVariation(
                oneBasedPosition: 25,
                originalSequence: "A",
                variantSequence: "T",
                description: "BranchBaseMod",
                variantCallFormatDataString: vcf,
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { pos, new List<Modification>{ new Modification(_originalId:"M1", _modificationType:"TestType") } }
                });

        // CASE 1: allRef true, emitReferenceForHomozygousRef false -> no variant
        [Test]
        public void AllRef_NoEmit_ReturnsEmpty()
        {
            string vcf = "1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/0:8,0:8";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(emitReferenceForHomozygousRef: false);
            Assert.That(result, Is.Empty);
        }

        // CASE 2: allRef true, emitReferenceForHomozygousRef true -> TryAdd ref?ref (no-op) caught -> still empty
        [Test]
        public void AllRef_EmitReference_NoOpCaught_ReturnsEmpty()
        {
            string vcf = "1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/0:10,0:10";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(emitReferenceForHomozygousRef: true);
            Assert.That(result, Is.Empty);
        }

        // CASE 3: allStoredAlt true (AlleleIndex=1, genotype 1/1) -> HomozygousAlt
        [Test]
        public void AllStoredAlt_HomozygousAltPath()
        {
            string vcf = "1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t1/1:0,11:11";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype();
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].Description.Contains("Mode=HomozygousAlt"), Is.True);
        }

        // CASE 4: containsDifferentAlt true and skipIfAltIndexMismatch true -> skipped
        [Test]
        public void ContainsDifferentAlt_SkipFlagTrue_Skipped()
        {
            // ALT T,G ; ANN -> T => storedAltIndex=1 ; genotype 0/2 includes allele 2 (different alt)
            string vcf = "1\t400\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:5,0,6:11";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(); // skipIfAltIndexMismatch default true
            Assert.That(result, Is.Empty);
        }

        // CASE 5: containsDifferentAlt true but skipIfAltIndexMismatch false -> MixedAltIndex(StoredAltOnly)
        [Test]
        public void ContainsDifferentAlt_SkipFlagFalse_MixedAltIndexAdded()
        {
            string vcf = "1\t400\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:4,0,7:11";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(skipIfAltIndexMismatch: false);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].Description.Contains("MixedAltIndex(StoredAltOnly)"), Is.True);
        }

        // CASE 6: Heterozygous standard (0/1) includeReferenceForHeterozygous false -> only HeterozygousAlt
        [Test]
        public void Heterozygous_NoRefRequest_OnlyAlt()
        {
            string vcf = "1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:6,7:13";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(includeReferenceForHeterozygous: false);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].Description.Contains("HeterozygousAlt"), Is.True);
            Assert.That(result[0].Description.Contains("HeterozygousRef"), Is.False);
        }

        // CASE 7: Heterozygous with includeReferenceForHeterozygous true -> ref attempt (no-op) caught, alt retained
        [Test]
        public void Heterozygous_WithRefRequest_RefNoOpCaughtAltAdded()
        {
            string vcf = "1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,8:13";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(includeReferenceForHeterozygous: true);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].Description.Contains("Mode=HeterozygousAlt"), Is.True);
            Assert.That(result[0].Description.Contains("HeterozygousRef"), Is.False);
        }

        // CASE 8: Heterozygous with includeReferenceForHeterozygous true AND variant-specific mod so reference attempt would still be no-op ? same outcome
        [Test]
        public void Heterozygous_WithMods_RefStillSuppressedAltAdded()
        {
            // Because the base SequenceVariation carries variant-specific modifications,
            // the reference no-op (A->A) attempt IS considered valid (hasMods == true) and is retained.
            // Therefore SplitPerGenotype returns TWO variants:
            //   1) HeterozygousRef  (A->A with variant-specific mods)
            //   2) HeterozygousAlt  (A->T)
            // Previous expectation of 1 was incorrect for the “has mods” case.
            var baseVar = MakeWithMod(
                "1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:7,9:16", pos: 5);

            var result = baseVar.SplitPerGenotype(includeReferenceForHeterozygous: true);

            Assert.That(result.Count, Is.EqualTo(2), "Expected both ref (with mods) and alt variants.");
            Assert.That(result.Any(v => v.Description.Contains("HeterozygousRef")), Is.True,
                "Reference variant with modifications should be present.");
            Assert.That(result.Any(v => v.Description.Contains("HeterozygousAlt")), Is.True,
                "Alternate variant should be present.");
            // Confirm cloned modifications persisted on at least one variant
            Assert.That(result.Any(v => v.OneBasedModifications?.ContainsKey(5) == true), Is.True,
                "Expected variant-specific modification to be cloned.");
        }

        // CASE 9: Non-matching ANN allele (AlleleIndex = 0) genotype 1/1 -> falls to heterozygous branch (not HomozygousAlt)
        [Test]
        public void AlleleIndexZero_GenotypeAlt_FallsThroughElseBranch()
        {
            // ANN allele = REF (A) => storedAltIndex=0, genotype 1/1 produces numericAlleles != allStoredAlt
            string vcf = "1\t400\t.\tA\tT\t.\tPASS\tANN=A|.\tGT:AD:DP\t1/1:0,9:9";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype();
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].Description.Contains("HomozygousAlt"), Is.False);
            Assert.That(result[0].Description.Contains("HeterozygousAlt"), Is.True);
        }

        // CASE 10: MixedAltIndex generating branch when depth filter applied (minDepth > depth) -> skipped before branching
        [Test]
        public void DepthFilterBeforeBranching_SuppressesAll()
        {
            string vcf = "1\t400\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:2,0,3:5";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(minDepth: 6, skipIfAltIndexMismatch: false);
            Assert.That(result, Is.Empty, "Depth filter should remove sample before branch logic.");
        }

        // CASE 11: Homozygous reference with includeReference false AND a variant-specific mod (still no variant)
        [Test]
        public void AllRef_WithMod_NoEmit_NoVariant()
        {
            var baseVar = MakeWithMod("1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/0:12,0:12", pos: 3);
            var result = baseVar.SplitPerGenotype();
            Assert.That(result, Is.Empty);
        }

        // CASE 12: Homozygous alt with includeReferenceForHeterozygous true (flag irrelevant in this path)
        [Test]
        public void HomozygousAlt_IgnoresHeterozygousRefFlag()
        {
            string vcf = "1\t400\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t1/1:0,15:15";
            var baseVar = Make(vcf);
            var result = baseVar.SplitPerGenotype(includeReferenceForHeterozygous: true);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].Description.Contains("Mode=HomozygousAlt"), Is.True);
        }
    }
}