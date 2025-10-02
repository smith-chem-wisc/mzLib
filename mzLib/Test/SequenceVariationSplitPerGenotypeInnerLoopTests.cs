using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.BioPolymer;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceVariationSplitPerGenotypeInnerLoopTests
    {
        private static SequenceVariation Make(string vcf) =>
            new SequenceVariation(
                oneBasedPosition: 25,
                originalSequence: "M",
                variantSequence: "K",
                description: "InnerLoopVariant",
                variantCallFormatDataString: vcf,
                oneBasedModifications: null);

        [Test]
        public void MissingGenotypeKey_Continues()
        {
            // Single-sample VCF, then remove genotype key -> loop sees missing -> continue -> no variants
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,7:12";
            var sv = Make(vcf);
            sv.VariantCallFormatData.Genotypes.Remove("0");
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void EmptyGenotypeTokens_Continues()
        {
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,6:11";
            var sv = Make(vcf);
            sv.VariantCallFormatData.Genotypes["0"] = Array.Empty<string>();
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void AlleleDepthSummation_SkipsDots_AndWhitespace()
        {
            // AD tokens include '.', whitespace, and valid ints. Depth = 4 + 3 + 2 = 9
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:4,.,  ,3,2:20";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 0);
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Depth=9", split[0].Description);
        }

        [Test]
        public void AlleleDepthAllDots_DepthZero_PassesWhenMinDepthZero()
        {
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:.,.,.:15";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 0);
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Depth=0", split[0].Description);
        }

        [Test]
        public void AlleleDepthNegativeValues_FallbacksToDP()
        {
            // Negative AD token makes entire AD invalid per ADvaluesAreValid (all tokens must be '.' or non-negative ints).
            // Implementation discards AD and falls back to DP=30.
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:6,-3,2:30";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Depth=30", split[0].Description, "Expected DP fallback when AD contains a negative value.");
            // Ensure no AD-based partial accumulation occurred
            Assert.That(split[0].Description.Contains("Depth=8"), Is.False, "AD summation should NOT occur when AD is invalid.");
        }

        [Test]
        public void DpFallbackUsed_WhenNoADFieldInFormat()
        {
            // Format excludes AD; dpIndex resolves; depth = 14
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:DP\t0/1:14";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Depth=14", split[0].Description);
        }

        [Test]
        public void DpFallback_NotApplied_WhenTokenCountMismatch()
        {
            // The VariantCallFormat parser enforces that FORMAT token count matches sample column token count.
            // This VCF line has FORMAT GT:AD:DP (3 fields) but the sample column only has 2 (0/1:5,6) ? constructor throws.
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,6";

            Assert.Throws<ArgumentException>(() => Make(vcf),
                "Expected an ArgumentException due to genotype / FORMAT token count mismatch.");
        }

        [Test]
        public void DepthBelowMinDepth_Continues()
        {
            // Depth from AD = 5 + 2 =7; minDepth=8 => variant skipped
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,2:20";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 8);
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void DepthExactlyMinDepth_Passes()
        {
            // Depth = 6; minDepth=6 -> included
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:1,5:20";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 6);
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Depth=6", split[0].Description);
        }
            
        [Test]
        public void MultipleSamples_MixedPaths_ADAndDP()
        {
            // Sample0: GT:AD:DP -> AD valid => depth = 3+4=7 (meets minDepth 5)
            // Sample1: GT:AD:DP -> AD token "." (length>0) => AD branch runs, all skipped => depth=0 (no DP fallback) -> excluded
            // Sample2: GT:AD:DP -> AD contains invalid token 'X' -> AD invalid ? stored as empty array ? AD branch skipped ? DP fallback depth=25
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:3,4:15\t0/1:.:9\t0/1:.,X,8:25";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 5);

            Assert.That(split.Count, Is.EqualTo(2), "Expected only samples 0 and 2 to pass depth filter.");

            // Sample 0
            Assert.That(split.Any(v => v.Description.Contains("Sample=0") && v.Description.Contains("Depth=7")),
                Is.True, "Sample 0 (depth 7) should be present.");

            // Sample 2 (DP fallback = 25, not partial AD sum)
            Assert.That(split.Any(v => v.Description.Contains("Sample=2") && v.Description.Contains("Depth=25")),
                Is.True, "Sample 2 should use DP fallback (25) after invalid AD.");

            // Ensure sample 1 excluded
            Assert.That(split.Any(v => v.Description.Contains("Sample=1")), Is.False, "Sample 1 depth=0 should be excluded.");
        }

        [Test]
        public void GenotypeParseError_SkipsSample()
        {
            // Introduce an invalid token in GT (non-numeric letter 'X') so numericAlleles remains maybe partial but parseError triggers continue
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/X:5,5:10";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 0);
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void NoCalledAlleles_SkipsSample()
        {
            // GT is './.' -> gtTokens are ['.','.'] -> numericAlleles empty -> continue
            string vcf = "1\t101\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t./.:5,5:10";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }
    }
}