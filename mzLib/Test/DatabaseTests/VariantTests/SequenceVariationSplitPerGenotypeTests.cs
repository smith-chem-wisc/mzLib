using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;

namespace Test.DatabaseTests.VariantTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceVariationSplitPerGenotypeTests
    {
        /*
         VCF (4 samples):
           0: 0/0 depth 10  (pure reference) ? would create a no?op; constructor rejects it (AreValid=false) ? excluded
           1: 0/1 depth 11  (heterozygous)   ? yields one alt variant (Mode=HeterozygousAlt)
           2: 1/1 depth 12  (homozygous alt) BUT storedAltIndex = -1 (ANN=.) so logic routes through heterozygous branch ? Mode=HeterozygousAlt
           3: 0/2 depth 9   (mixed alleles)  storedAltIndex = -1 so treated same (hetero path) ? Mode=HeterozygousAlt if depth passes filter

         Depth thresholds:
           minDepth 0 or 1  ? samples 1,2,3 pass (depths 11,12,9) ? 3 variants
           minDepth 10      ? samples 1,2 pass (11,12)            ? 2 variants
           
         Flags includeReferenceForHeterozygous / emitReferenceForHomozygousRef attempt to add ref “variants”
         but those are no?ops and SequenceVariation constructor rejects them ? no effect on output.
         skipIfAltIndexMismatch also no effect because storedAltIndex = -1 (guard requires >0).
        */
        private const string MultiSampleVcf =
            "1\t1000\trsX\tA\tT,G\t.\tPASS\tANN=.\tGT:AD:DP\t0/0:10,0,0:10\t0/1:5,6,0:11\t1/1:0,12,0:12\t0/2:4,0,5:9";

        private static SequenceVariation MakeBaseVariant() =>
            new SequenceVariation(
                oneBasedPosition: 10,
                originalSequence: "A",
                variantSequence: "T",
                description: "BaseVariant",
                variantCallFormatDataString: MultiSampleVcf,
                oneBasedModifications: null);

        private static IEnumerable<TestCaseData> Matrix()
        {
            int[] depths = { 0, 1, 10 };
            bool[] bools = { false, true };
            foreach (var minDepth in depths)
            {
                // Expected variant count based solely on depth (see comment above)
                int expected = (11 >= minDepth ? 1 : 0) + (12 >= minDepth ? 1 : 0) + (9 >= minDepth ? 1 : 0);
                foreach (var includeRefHet in bools)
                foreach (var emitRefHomRef in bools)
                foreach (var skipAltMismatch in bools)
                {
                    yield return new TestCaseData(minDepth, includeRefHet, emitRefHomRef, skipAltMismatch, expected)
                        .SetName($"MinDepth={minDepth},IncludeHetRef={includeRefHet},EmitHomRef={emitRefHomRef},SkipAltMismatch={skipAltMismatch},Expected={expected}");
                }
            }
        }

        [TestCaseSource(nameof(Matrix))]
        public void SplitPerGenotype_AdjustedExpectations(
            int minDepth,
            bool includeReferenceForHeterozygous,
            bool emitReferenceForHomozygousRef,
            bool skipIfAltIndexMismatch,
            int expectedCount)
        {
            var baseVar = MakeBaseVariant();

            var split = baseVar.SplitPerGenotype(
                minDepth: minDepth,
                includeReferenceForHeterozygous: includeReferenceForHeterozygous,
                emitReferenceForHomozygousRef: emitReferenceForHomozygousRef,
                skipIfAltIndexMismatch: skipIfAltIndexMismatch);

            // Count check
            Assert.That(split.Count, Is.EqualTo(expectedCount), "Variant count mismatch.");

            // All variants must represent a sequence change (no no-ops)
            Assert.That(split.Any(v => v.OriginalSequence == v.VariantSequence), Is.False, "Found unexpected no-op variant.");

            // Because AlleleIndex == -1 (ANN=.), every alt follows heterozygous branch ? Mode=HeterozygousAlt
            Assert.That(split.All(v => v.Description.Contains("Mode=HeterozygousAlt")),
                Is.True, "Expected only Mode=HeterozygousAlt due to AlleleIndex=-1 routing.");

            // Ensure no HomozygousAlt or MixedAltIndex modes appear
            Assert.That(split.Any(v => v.Description.Contains("HomozygousAlt")), Is.False);
            Assert.That(split.Any(v => v.Description.Contains("MixedAltIndex")), Is.False);
        }
    }
}