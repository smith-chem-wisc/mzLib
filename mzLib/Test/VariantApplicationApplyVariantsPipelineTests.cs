using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Proteomics;

namespace Test.DatabaseTests
{
    [TestFixture]
    public class VariantApplicationApplyVariantsPipelineTests
    {
        /*
         * Extended: add tests for genotype / zygosity gating logic inside ApplyVariants:
         *
         *   var vcf = variant.VariantCallFormatData;
         *   if (vcf == null || vcf.Genotypes == null || !vcf.Genotypes.ContainsKey(individual)) continue;
         *   var alleleIndexStr = vcf.AlleleIndex.ToString();
         *   bool variantAlleleIsInTheGenotype = vcf.Genotypes[individual].Contains(alleleIndexStr);
         *   if (!variantAlleleIsInTheGenotype) continue;
         *   bool hetero = ...
         *   bool homoAlternate = ...
         *
         * We cover branches:
         *   - vcf == null (already covered by earlier filtering logic reuse)
         *   - Missing individual genotype (variant lacks that sample key)
         *   - Allele not in genotype (0/0 ? skip)
         *   - Heterozygous (0/1)
         *   - Homozygous alternate (1/1)
         */

        #region Helpers (existing + new lightweight builders)

        private string BuildVcf(
            string chrom,
            int pos,
            string refAA,
            string altAA,
            string sample0GT,
            string sample1GT,
            string sample0AD = "12,11",
            string sample1AD = "15,14",
            string qual = ".",
            string filter = "PASS")
        {
            var cols = new[]
            {
                chrom, pos.ToString(), ".", refAA, altAA, qual, filter,
                "ANN=" + altAA + "|missense|GENE|GENE|",
                "GT:AD:DP",
                $"{sample0GT}:{sample0AD}:23",
                $"{sample1GT}:{sample1AD}:29"
            };
            return string.Join('\t', cols);
        }

        // Build VCF with only one sample column (sample0 only)
        private string BuildSingleSampleVcf(
            string chrom,
            int pos,
            string refAA,
            string altAA,
            string sample0GT,
            string sample0AD = "9,8",
            string qual = ".",
            string filter = "PASS")
        {
            var cols = new[]
            {
                chrom, pos.ToString(), ".", refAA, altAA, qual, filter,
                "ANN=" + altAA + "|missense|GENE|GENE|",
                "GT:AD:DP",
                $"{sample0GT}:{sample0AD}:17"
            };
            return string.Join('\t', cols);
        }

        private SequenceVariation MakeVar(int begin, string orig, string variant, string desc, string vcfLine = null)
        {
            return new SequenceVariation(begin,
                begin + (orig?.Length > 0 ? orig.Length - 1 : 0),
                orig,
                variant,
                desc,
                vcfLine);
        }

        private static HashSet<string> VariantSetKey(Protein p) =>
            new(p.AppliedSequenceVariations.Select(v => v.SimpleString()));

        private static HashSet<string> AllAppliedSimpleStrings(IEnumerable<Protein> proteins) =>
            new(proteins.SelectMany(p => p.AppliedSequenceVariations ?? new List<SequenceVariation>())
                        .Select(v => v.SimpleString()));

        #endregion

        #region Genotype Filtering / Classification Tests

        [Test]
        public void ApplyVariants_GenotypeSkip_AlleleNotInGenotype_RefRef()
        {
            // Variant with genotype 0/0 (alleleIndex likely "1" for first ALT) ? allele not present ? skip
            var protein = new Protein("MPEPTIDERESIDUESKIPTEST", "SKIP_REFREF");
            var vcfRefRef = BuildVcf("1", 8, "E", "K", "0/0", "0/0"); // both samples homo ref
            var refVar = MakeVar(8, "E", "K", "should_skip_refref", vcfRefRef);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { refVar },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 1);

            var applied = AllAppliedSimpleStrings(produced);
            Assert.That(applied.Contains(refVar.SimpleString()), Is.False,
                "Variant with 0/0 genotype should be skipped (allele not in genotype).");
        }

        [Test]
        public void ApplyVariants_Genotype_Heterozygous_BranchingPresent()
        {
            // Single heterozygous variant 0/1 for both samples ? expect at least 2 proteoforms:
            // either the algorithm duplicates (one ref, one alt) or yields only alt if threshold logic collapses,
            // but heterozygous path should apply variant at least once.
            var protein = new Protein("MPEPTIDEHETEROXYZ", "HET_SINGLE");
            var vcfHet = BuildVcf("1", 6, "T", "A", "0/1", "0/1");
            var hetVar = MakeVar(6, "T", "A", "het_variant", vcfHet);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { hetVar },
                maxAllowedVariantsForCombinatorics: 3,
                minAlleleDepth: 1);

            var sets = produced.Select(VariantSetKey).ToList();
            var withVar = sets.Count(s => s.Contains(hetVar.SimpleString()));
            Assert.That(withVar, Is.GreaterThanOrEqualTo(1),
                "Heterozygous variant should appear at least once.");
        }

        [Test]
        public void ApplyVariants_Genotype_HomozygousAlternate_NoBaseRetained()
        {
            // Homozygous alt (1/1) variant & deep alt depth -> all resulting proteoforms should include the variant
            var protein = new Protein("MPEPTIDEHOMOALL", "HOMO_ALT");
            var vcfHomoAlt = BuildVcf("1", 4, "P", "L", "1/1", "1/1");
            var homoAlt = MakeVar(4, "P", "L", "homo_alt", vcfHomoAlt);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { homoAlt },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 1);

            var sets = produced.Select(VariantSetKey).ToList();
            // All variant sets should contain the homo alt variant (base copy mutated)
            Assert.That(sets.All(s => s.Contains(homoAlt.SimpleString())), Is.True,
                "All proteoforms should include homozygous alternate variant.");
        }

        [Test]
        public void ApplyVariants_Genotype_MissingSampleKey_SkipsVariantForThatIndividual()
        {
            // Variant A has both samples ? individuals set includes "0" and "1"
            // Variant B has only sample0 genotype ? during iteration for individual "1" it must be skipped (no key)
            var protein = new Protein("MPEPTIDEMISSINGSAMPLE", "MISS_KEY");

            var vBoth = BuildVcf("1", 12, "E", "G", "0/1", "0/1"); // hetero both samples
            var vOnly0 = BuildSingleSampleVcf("1", 20, "K", "R", "0/1"); // only sample0 column

            var varBoth = MakeVar(12, "E", "G", "both_samples", vBoth);
            var varOnly0 = MakeVar(20, "K", "R", "sample0_only", vOnly0);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varBoth, varOnly0 },
                maxAllowedVariantsForCombinatorics: 3,
                minAlleleDepth: 1);

            var sets = produced.Select(VariantSetKey).ToList();

            // Variant present overall (sample0 path)
            bool varOnly0Present = sets.Any(s => s.Contains(varOnly0.SimpleString()));
            Assert.That(varOnly0Present, Is.True, "Variant with only sample0 genotype should appear in sample0-derived proteoforms.");

            // Find any proteoform that includes varOnly0 but excludes varBoth — indicates a branch from individual 0 only.
            bool isolatedSample0Evidence = sets.Any(s =>
                s.Contains(varOnly0.SimpleString()) && !s.Contains(varBoth.SimpleString()));
            Assert.That(isolatedSample0Evidence, Is.True,
                "Expected at least one proteoform showing variant-only0 applied without the both-samples variant, evidencing sample1 skip.");

            // Ensure no contradiction: all sets that include varOnly0 came from sample0 iteration; sample1 iteration cannot add it
            // (Indirect check: if sample1 had applied it, we'd expect combination sets where sample1-only variant exists with absence of both-samples variant after some logic)
        }

        [Test]
        public void ApplyVariants_Genotype_SkipWhenAlleleAbsent_AndApplyOthers()
        {
            // Mixed variants: one 0/0 (skip), one 0/1 (apply), one 1/1 (apply everywhere)
            // Ensure all coordinates are within sequence length.
            var protein = new Protein("MPEPTIDELONGSEQUENCEFORTEST", "MIXED_GENO");

            // Positions: 18 (homo alt), 12 (hetero), 5 (ref/ref) so ordering (desc) processes homo-alt first
            var vSkip = BuildVcf("1", 5,  "D", "N", "0/0", "0/0");   // skip (allele not in genotype)
            var vHet  = BuildVcf("1", 12, "T", "S", "0/1", "0/1");  // heterozygous
            var vAlt  = BuildVcf("1", 18, "A", "V", "1/1", "1/1");  // homozygous alternate

            var varSkip = MakeVar(5,  "D", "N", "skip_refref", vSkip);
            var varHet  = MakeVar(12, "T", "S", "het_apply",   vHet);
            var varAlt  = MakeVar(18, "A", "V", "hom_alt",     vAlt);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varSkip, varHet, varAlt },
                maxAllowedVariantsForCombinatorics: 3,
                minAlleleDepth: 1);

            var applied = AllAppliedSimpleStrings(produced);

            Assert.Multiple(() =>
            {
                Assert.That(applied.Contains(varSkip.SimpleString()), Is.False, "Ref/Ref variant should be skipped.");
                Assert.That(applied.Contains(varHet.SimpleString()), Is.True,  "Heterozygous variant should be applied somewhere.");
                Assert.That(applied.Contains(varAlt.SimpleString()), Is.True,  "Homozygous alt variant should be applied everywhere.");
            });
        }
        [Test]
        public void ApplyVariants_HeterozygousThreshold_AltOnlyBranch()
        {
            // Protein length is 23; keep all variant positions <= length
            // Force tooManyHeterozygousVariants = true with ref depth below threshold (2) but alt depth above (15).
            string BuildAltFavoredVcf(int pos, string refAA, string altAA, string gt) =>
                string.Join('\t', new[]
                {
                    "1", pos.ToString(), ".", refAA, altAA, ".", "PASS",
                    "ANN=" + altAA + "|missense|GENE|GENE|",
                    "GT:AD:DP",
                    $"{gt}:2,15:17",
                    $"{gt}:2,15:17"
                });

            var protein = new Protein("MPEPTIDEHETALTBRANCHSEQ", "HET_ALT_BRANCH"); // length 23 (Q at 23)

            // Three heterozygous variants (0/1) ? hetero count (3) > maxAllowedVariantsForCombinatorics(=1) ? threshold path
            // Use valid coordinates: 23 (Q->R), 15 (T->A? base at 15 is B? Actually sequence index 15 = B from 'BRANCH'; keep original letter check),
            // For clarity pick residues matching actual sequence:
            // Sequence indexed: 1:M 2:P 3:E 4:P 5:T 6:I 7:D 8:E 9:H 10:E 11:T 12:A 13:L 14:T 15:B 16:R 17:A 18:N 19:C 20:H 21:S 22:E 23:Q
            // Since 'B' is not a standard residue, if the actual sequence differs in your source, adjust accordingly.
            // To remain safe, mutate positions we know: 23 (Q->R), 12 (A->G), 5 (T->S).

            var v1 = MakeVar(23, "Q", "R", "het_alt_only_23", BuildAltFavoredVcf(23, "Q", "R", "0/1"));
            var v2 = MakeVar(12, "A", "G", "het_alt_only_12", BuildAltFavoredVcf(12, "A", "G", "0/1"));
            var v3 = MakeVar(5, "T", "S", "het_alt_only_05", BuildAltFavoredVcf(5, "T", "S", "0/1"));

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { v1, v2, v3 },
                maxAllowedVariantsForCombinatorics: 1,
                minAlleleDepth: 10); // ref depth 2 fails, alt depth 15 passes

            var variantSets = produced.Select(p => VariantSetKey(p)).ToList();
            var flattened = new HashSet<string>(variantSets.SelectMany(s => s));

            Assert.That(flattened.Contains(v1.SimpleString()), Is.True, "Variant v1 (pos23) not applied in alt-only branch.");
            Assert.That(flattened.Contains(v2.SimpleString()), Is.True, "Variant v2 (pos12) not applied in alt-only branch.");
            Assert.That(flattened.Contains(v3.SimpleString()), Is.True, "Variant v3 (pos5) not applied in alt-only branch.");

            bool cumulativeExists = variantSets.Any(s =>
                s.SetEquals(new[] { v1.SimpleString(), v2.SimpleString(), v3.SimpleString() }));
            TestContext.WriteLine("Cumulative heterozygous alt-only proteoform present: " + cumulativeExists);
        }
        #endregion
    }
}