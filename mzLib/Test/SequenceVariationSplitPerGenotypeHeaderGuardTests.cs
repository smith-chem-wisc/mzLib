using System;
using System.Linq;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.BioPolymer;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceVariationSplitPerGenotypeHeaderGuardTests
    {
        private static SequenceVariation Make(string vcf) =>
            new SequenceVariation(
                oneBasedPosition: 10,
                originalSequence: "A",
                variantSequence: "T",
                description: "Var",
                variantCallFormatDataString: vcf,
                oneBasedModifications: null);

        [Test]
        public void SplitPerGenotype_ReturnsEmpty_WhenNoVcfData()
        {
            // Variant created without a VCF line
            var sv = new SequenceVariation(10, "A", "T", "NoVcf");
            var list = sv.SplitPerGenotype();
            Assert.That(list, Is.Empty);
        }

        [Test]
        public void SplitPerGenotype_ReturnsEmpty_WhenGenotypesMissing()
        {
            // <10 columns (only 9) ? parsing aborts; Genotypes null/empty triggers first early return
            string vcfNoSamples = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=.\tGT:AD";
            var sv = Make(vcfNoSamples);
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void SplitPerGenotype_ReturnsEmpty_WhenFieldsBelowThresholdWithGenotypeCheck()
        {
            // Same as above; documents unreachable second guard (vcfFields.Length < 10) because initial genotype guard fires first.
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=.\tGT";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void SplitPerGenotype_NoDPToken_DepthFromAD()
        {
            // FORMAT excludes DP ? dpIndex = -1; depth calculated from AD sum (5+4=9)
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=T|missense_variant\tGT:AD\t0/1:5,4";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 0);
            Assert.That(split.Count, Is.EqualTo(1));
            var d = split[0].Description;
            StringAssert.Contains("Depth=9", d);
            StringAssert.Contains("Mode=HeterozygousAlt", d);
        }

        [Test]
        public void SplitPerGenotype_WithDPToken_NoAD_UsesDP()
        {
            // FORMAT has GT:DP, no AD. dpIndex valid. Depth=14.
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=T|.\tGT:DP\t0/1:14";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Depth=14", split[0].Description);
        }

        [Test]
        public void SplitPerGenotype_HomozygousAlt_StoredAltIndexPositive()
        {
            // ANN allele = T (ALT1) => AlleleIndex=1; genotype 1/1 => HomozygousAlt path
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t1/1:0,8:8";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Mode=HomozygousAlt", split[0].Description);
        }

        [Test]
        public void SplitPerGenotype_HomozygousAlt_ButAlleleIndexZero_TreatedAsHeterozygousAltPath()
        {
            // ANN allele = REF (A) => storedAltIndex=0 ? allStoredAlt false even for 1/1 => falls through heterozygous branch
            // genotype 1/1 still includes only alt allele index 1, but code uses storedAltIndex (0) so "HomozygousAlt" not used.
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=A|.\tGT:AD:DP\t1/1:0,9:9";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Mode=HeterozygousAlt", split[0].Description);
            Assert.False(split[0].Description.Contains("HomozygousAlt"));
        }

        [Test]
        public void SplitPerGenotype_AlleleIndexUnknown_NegativeOne()
        {
            // ANN=.; AlleleIndex = -1; heterozygous 0/1
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=.\tGT:AD:DP\t0/1:4,7:11";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype();
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Mode=HeterozygousAlt", split[0].Description);
        }

        [Test]
        public void SplitPerGenotype_MixedAltIndex_SkippedWhenFlagTrue()
        {
            // ALT = T,G ; ANN allele = T -> storedAltIndex=1; sample genotype 0/2 (containsDifferentAlt).
            // skipIfAltIndexMismatch = true (default) => no variant yielded for sample 0/2
            string vcf = "1\t1000\trsX\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:3,0,5:8";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 0); // depth 8 passes
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void SplitPerGenotype_MixedAltIndex_YieldsWhenFlagFalse()
        {
            string vcf = "1\t1000\trsX\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:3,0,5:8";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 0, skipIfAltIndexMismatch: false);
            Assert.That(split.Count, Is.EqualTo(1));
            StringAssert.Contains("Mode=MixedAltIndex(StoredAltOnly)", split[0].Description);
        }

        [Test]
        public void SplitPerGenotype_IncludeReferenceForHeterozygous_NoOpFiltered()
        {
            // includeReferenceForHeterozygous tries to add a ref variant (no-op) which will fail validation; only alt remains.
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:6,7:13";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(includeReferenceForHeterozygous: true);
            Assert.That(split.Count, Is.EqualTo(1));
            Assert.False(split.Any(v => v.Description.Contains("HeterozygousRef")));
            StringAssert.Contains("HeterozygousAlt", split[0].Description);
        }

        [Test]
        public void SplitPerGenotype_EmitReferenceHomozygousRef_NoOpFiltered()
        {
            // Homozygous reference sample only: attempt to emit ref variant but it's a no-op; result empty.
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/0:8,0:8";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(emitReferenceForHomozygousRef: true);
            Assert.That(split, Is.Empty);
        }

        [Test]
        public void SplitPerGenotype_DepthFilterApplied()
        {
            // depth = 9; minDepth = 10 => excluded
            string vcf = "1\t1000\trsX\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:4,5:9";
            var sv = Make(vcf);
            var split = sv.SplitPerGenotype(minDepth: 10);
            Assert.That(split, Is.Empty);
        }
    }
}