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
        #region Allele Depth (isDeepReferenceAllele / isDeepAlternateAllele) Tests

        private string BuildDepthVcf(string chrom, int pos, string refAA, string altAA,
            string sample0GT, string sample1GT,
            string sample0AD, string sample1AD,
            string format = "GT:AD:DP", string extraInfo = null)
        {
            // INFO with ANN so AlleleIndex resolves (single ALT ? index 1)
            var info = extraInfo ?? $"ANN={altAA}|missense|GENE|GENE|";
            return string.Join('\t', new[]
            {
                chrom, pos.ToString(), ".", refAA, altAA, ".", "PASS",
                info,
                format,
                $"{sample0GT}:{sample0AD}:20",
                $"{sample1GT}:{sample1AD}:22"
            });
        }

        private string BuildDepthVcfSingleSample(string chrom, int pos, string refAA, string altAA,
            string sample0GT, string sample0AD, string format = "GT:AD:DP", string extraInfo = null)
        {
            var info = extraInfo ?? $"ANN={altAA}|missense|GENE|GENE|";
            return string.Join('\t', new[]
            {
                chrom, pos.ToString(), ".", refAA, altAA, ".", "PASS",
                info,
                format,
                $"{sample0GT}:{sample0AD}:20"
            });
        }

        private Protein MakeBaseDepthProtein() => new Protein("MPEPTIDESEQVARTEST", "DEPTH_BASE"); // length 17

        private static bool VariantApplied(IEnumerable<Protein> proteins, SequenceVariation v) =>
            proteins.SelectMany(p => p.AppliedSequenceVariations ?? new List<SequenceVariation>())
                    .Any(ap => ap.SimpleString() == v.SimpleString());

        [Test]
        public void ApplyVariants_Depth_HomoAlt_AltPasses_RefPasses()
        {
            // Both ref & alt depths >= minAlleleDepth (10) ? homozygous alt variant applied
            var baseProt = MakeBaseDepthProtein();
            var vcf = BuildDepthVcf("1", 8, "E", "K", "1/1", "1/1", "12,15", "14,18"); // ref=12/14 alt=15/18
            var varAlt = MakeVar(8, "E", "K", "homoAltBothDeep", vcf);

            var produced = VariantApplication.ApplyVariants(
                baseProt,
                new[] { varAlt },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 10);

            Assert.That(VariantApplied(produced, varAlt), Is.True,
                "Homozygous alt variant should be applied when alt depth passes.");
        }

        [Test]
        public void ApplyVariants_Depth_HomoAlt_AltBelowThreshold_NotApplied()
        {
            // Alt depth < threshold (alt=5 < 10) ? not applied
            var baseProt = MakeBaseDepthProtein();
            var vcf = BuildDepthVcf("1", 8, "E", "K", "1/1", "1/1", "12,5", "14,7"); // alt depths below
            var varAlt = MakeVar(8, "E", "K", "homoAltAltTooShallow", vcf);

            var produced = VariantApplication.ApplyVariants(
                baseProt,
                new[] { varAlt },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 10);

            Assert.That(VariantApplied(produced, varAlt), Is.False,
                "Variant should not be applied when alt depth is below threshold.");
        }

        [Test]
        public void ApplyVariants_Depth_HomoAlt_AltDepthNonNumeric_NotApplied()
        {
            // Non-numeric alt depth token ? int.TryParse fails ? not applied
            var baseProt = MakeBaseDepthProtein();
            var vcf = BuildDepthVcf("1", 8, "E", "K", "1/1", "1/1", "12,XYZ", "11,QQ"); // alt tokens invalid
            var varAlt = MakeVar(8, "E", "K", "homoAltAltNonNumeric", vcf);

            var produced = VariantApplication.ApplyVariants(
                baseProt,
                new[] { varAlt },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 10);

            Assert.That(VariantApplied(produced, varAlt), Is.False,
                "Variant should not be applied when alt depth token is non-numeric.");
        }

        [Test]
        public void ApplyVariants_Depth_HomoAlt_AlleleDepthsMissing_NotApplied()
        {
            // Remove AD field entirely (FORMAT GT:DP only) ? AlleleDepths empty ? not applied
            var baseProt = MakeBaseDepthProtein();
            string vcf = string.Join('\t', new[]
            {
                "1","8",".","E","K",".","PASS",
                "ANN=K|missense|GENE|GENE|",
                "GT:DP",
                "1/1:20",
                "1/1:22"
            });
            var varAlt = MakeVar(8, "E", "K", "homoAltNoAD", vcf);

            var produced = VariantApplication.ApplyVariants(
                baseProt,
                new[] { varAlt },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 5);

            Assert.That(VariantApplied(produced, varAlt), Is.False,
                "Variant should not be applied when AD field absent.");
        }

        [Test]
        public void ApplyVariants_Depth_HomoAlt_AlleleIndexOutOfRange_NotApplied()
        {
            // AD has only one value (ref) so alt index (1) is out of range ? alt depth check fails
            var baseProt = MakeBaseDepthProtein();
            var vcf = BuildDepthVcf("1", 8, "E", "K", "1/1", "1/1", "12", "11"); // AD arrays length 1
            var varAlt = MakeVar(8, "E", "K", "homoAltAltIndexOutOfRange", vcf);

            var produced = VariantApplication.ApplyVariants(
                baseProt,
                new[] { varAlt },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 5);

            Assert.That(VariantApplied(produced, varAlt), Is.False,
                "Variant should not be applied when alt index is out of AD range.");
        }

        [Test]
        public void ApplyVariants_Depth_HomoAlt_MissingSampleKey_OnlyPresentSampleConsidered()
        {
            // Only sample0 present. Verify variant applied for sample0 path (alt deep), no error for missing sample1.
            var baseProt = MakeBaseDepthProtein();
            var vcf = BuildDepthVcfSingleSample("1", 8, "E", "K", "1/1", "12,15"); // single sample
            var varAlt = MakeVar(8, "E", "K", "homoAltSingleSample", vcf);

            var produced = VariantApplication.ApplyVariants(
                baseProt,
                new[] { varAlt },
                maxAllowedVariantsForCombinatorics: 2,
                minAlleleDepth: 10);

            Assert.That(VariantApplied(produced, varAlt), Is.True,
                "Variant should be applied for existing sample; missing other sample key should not block application.");
        }

        [Test]
        public void ApplyVariants_Depth_Hetero_AltDeep_RefShallow_AltPathOnly()
        {
            // Heterozygous 0/1, alt deep (15), ref shallow (2). Should still allow application via alt path.
            var baseProt = MakeBaseDepthProtein();
            string vcf = BuildDepthVcf("1", 8, "E", "K", "0/1", "0/1", "2,15", "3,14");
            var hetVar = MakeVar(8, "E", "K", "heteroAltOnlyDepth", vcf);

            var produced = VariantApplication.ApplyVariants(
                baseProt,
                new[] { hetVar },
                maxAllowedVariantsForCombinatorics: 3,
                minAlleleDepth: 10);

            Assert.That(VariantApplied(produced, hetVar), Is.True,
                "Heterozygous variant with alt deep / ref shallow should still be applied (alt path).");
        }

        #endregion
        #region Heterozygous Threshold Internal Branch Tests

        private string BuildThresholdVcf(int pos, string refAA, string altAA,
            string sample0GT, string sample1GT,
            string sample0AD, string sample1AD)
        {
            // GT:AD:DP with ANN annotation (single ALT)
            return string.Join('\t', new[]
            {
                "1", pos.ToString(), ".", refAA, altAA, ".", "PASS",
                $"ANN={altAA}|missense|GENE|GENE|",
                "GT:AD:DP",
                $"{sample0GT}:{sample0AD}:25",
                $"{sample1GT}:{sample1AD}:27"
            });
        }

        private Protein MakeBaseThresholdProtein() => new Protein("MPEPTIDEVARIANTBRANCHSEQ", "HET_THRESH_BASE");

        private static HashSet<string> VariantSimpleSets(IEnumerable<Protein> proteins) =>
            new(proteins.Select(p =>
                string.Join("|", (p.AppliedSequenceVariations ?? new List<SequenceVariation>())
                    .Select(v => v.SimpleString())
                    .OrderBy(s => s))));

        [Test]
        public void ApplyVariants_HeteroThreshold_AddsSecondProtein_ThenUpdatesSecond()
        {
            // Two heterozygous variants; maxAllowedVariantsForCombinatorics=1 ? threshold triggers (2 > 1)
            // Both ref & alt depths >= minDepth ? isDeepReferenceAllele && isDeepAlternateAllele
            // First variant (count==1) adds second protein; second variant (count>1) updates second protein only.
            var protein = MakeBaseThresholdProtein();

            var vcfHighA = BuildThresholdVcf(18, "E", "K", "0/1", "0/1", "12,14", "11,13");
            var vcfHighB = BuildThresholdVcf(10, "T", "A", "0/1", "0/1", "15,16", "14,15");

            var varA = MakeVar(18, "E", "K", "hetA_bothDeep", vcfHighA);
            var varB = MakeVar(10, "T", "A", "hetB_bothDeep", vcfHighB);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varA, varB },
                maxAllowedVariantsForCombinatorics: 1,
                minAlleleDepth: 10);

            var setStrings = VariantSimpleSets(produced);

            // Expect only:
            //  "" (base) and "T10A|E18K"  (second protein accumulates both after update)
            Assert.That(setStrings.Contains(""), Is.True, "Base branch (unmodified) missing.");
            Assert.That(setStrings.Contains($"{varB.SimpleString()}|{varA.SimpleString()}") ||
                        setStrings.Contains($"{varA.SimpleString()}|{varB.SimpleString()}"),
                Is.True, "Combined variant branch (both variants) missing.");

            // No intermediate single-variant proteoform should remain after second variant updates slot
            bool singleVariantPresent = setStrings.Any(s =>
                !string.IsNullOrEmpty(s) &&
                (s.Split('|').Length == 1));
            Assert.That(singleVariantPresent, Is.False,
                "Found a single-variant proteoform; expected replacement of second branch.");
        }

        [Test]
        public void ApplyVariants_HeteroThreshold_AltDeepRefShallow_AppliesToAllExistingProteins()
        {
            // Two heterozygous variants; each alt deep (>=10), ref shallow (<10).
            // threshold path; internal alt-only branch (isDeepAlternateAllele && !isDeepReferenceAllele)
            // Each variant maps across all current newVariantProteins (size stays 1, mutated sequentially).
            var protein = MakeBaseThresholdProtein();

            var vcfAltOnly1 = BuildThresholdVcf(16, "P", "L", "0/1", "0/1", "3,15", "2,14");
            var vcfAltOnly2 = BuildThresholdVcf(7, "D", "N", "0/1", "0/1", "4,13", "3,12");

            var var1 = MakeVar(16, "P", "L", "het_altOnly16", vcfAltOnly1);
            var var2 = MakeVar(7, "D", "N", "het_altOnly07", vcfAltOnly2);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { var1, var2 },
                maxAllowedVariantsForCombinatorics: 1,
                minAlleleDepth: 10);

            var sets = VariantSimpleSets(produced);

            // Expect exactly one proteoform with both variants applied; base eliminated.
            Assert.That(sets.Contains(""), Is.False,
                "Base proteoform should be absent (alt-only mapping replaced it).");

            string combinedKey1 = $"{var2.SimpleString()}|{var1.SimpleString()}";
            string combinedKey2 = $"{var1.SimpleString()}|{var2.SimpleString()}";
            Assert.That(sets.Contains(combinedKey1) || sets.Contains(combinedKey2), Is.True,
                "Combined alt-only heterozygous proteoform missing.");
            Assert.That(sets.Count, Is.EqualTo(1),
                "Unexpected additional proteoforms present for alt-only threshold scenario.");
        }

        [Test]
        public void ApplyVariants_HeteroThreshold_AltDeepRefDeep_FirstAddsSecond_SecondAltOnly_RewritesBoth()
        {
            // Mixed case: first variant both deep (adds second branch),
            // second variant alt-only (ref shallow) => alt-only branch applies to ALL existing branches,
            // producing two proteoforms each now carrying the second variant; first branch remains base-only + second variant.
            var protein = MakeBaseThresholdProtein();

            var vcfBoth = BuildThresholdVcf(14, "A", "V", "0/1", "0/1", "11,12", "10,11");
            var vcfAltOnly = BuildThresholdVcf(6, "T", "S", "0/1", "0/1", "3,14", "2,15"); // ref shallow

            var varBoth = MakeVar(14, "A", "V", "het_bothDeep14", vcfBoth);
            var varAltOnly = MakeVar(6, "T", "S", "het_altOnly06", vcfAltOnly);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varBoth, varAltOnly },
                maxAllowedVariantsForCombinatorics: 1,
                minAlleleDepth: 10);

            var variantSets = produced
                .Select(p => p.AppliedSequenceVariations.Select(v => v.SimpleString()).OrderBy(s => s))
                .Select(s => string.Join("|", s))
                .ToHashSet();

            // After first: {"" , "A14V"}
            // After second (alt-only applies to ALL): {"T6S", "A14V|T6S"}
            Assert.That(variantSets.Contains("T6S"), Is.True,
                "Expected modified base branch with only alt-only second variant.");
            Assert.That(variantSets.Contains($"{varAltOnly.SimpleString()}|{varBoth.SimpleString()}") ||
                        variantSets.Contains($"{varBoth.SimpleString()}|{varAltOnly.SimpleString()}"),
                Is.True,
                "Expected cumulative branch (both variants) missing.");
            Assert.That(variantSets.Contains(""), Is.False,
                "Base (unmodified) branch should have been replaced by alt-only mapping.");
            Assert.That(variantSets.Contains(varBoth.SimpleString()), Is.False,
                "Intermediate single first variant branch should have been overwritten.");
        }

        [Test]
        public void ApplyVariants_HeteroThreshold_LimitZero_NoApplication()
        {
            // With maxAllowedVariantsForCombinatorics=0 internal blocks are guarded; no variant application.
            var protein = MakeBaseThresholdProtein();
            var vcfDeep = BuildThresholdVcf(12, "E", "G", "0/1", "0/1", "11,14", "10,13");
            var varDeep = MakeVar(12, "E", "G", "het_deep_limit0", vcfDeep);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varDeep },
                maxAllowedVariantsForCombinatorics: 0,
                minAlleleDepth: 5);

            // Expect only base proteoform (sequence identical, no variants applied)
            Assert.That(produced.Count, Is.EqualTo(1), "Unexpected additional proteoforms created with limit zero.");
            Assert.That(produced[0].AppliedSequenceVariations.Count, Is.EqualTo(0),
                "No variants should be applied when maxAllowedVariantsForCombinatorics=0.");
        }
        #endregion
        #region Heterozygous Combinatorics Branch Tests

        private string BuildSingleSampleCombinatoricsVcf(
            int pos,
            string refAA,
            string altAA,
            string genotype,
            int refDepth,
            int altDepth)
        {
            // Single-sample, GT:AD:DP format. ANN ensures AlleleIndex resolves (single ALT -> index 1).
            return string.Join('\t', new[]
            {
                "1", pos.ToString(), ".", refAA, altAA, ".", "PASS",
                $"ANN={altAA}|missense|GENE|GENE|",
                "GT:AD:DP",
                $"{genotype}:{refDepth},{altDepth}:{refDepth + altDepth + 5}"
            });
        }

        private Protein MakeCombinatoricsProtein() => new Protein("MPEPTIDEVARIANTCOMBINATORICSEQ", "HET_COMB_BASE"); // length >= positions used

        private HashSet<string> ProteoformVariantSetStrings(IEnumerable<Protein> proteins) =>
            proteins.Select(p =>
                    string.Join("|",
                        (p.AppliedSequenceVariations ?? new List<SequenceVariation>())
                        .Select(v => v.SimpleString())
                        .OrderBy(s => s)))
                    .ToHashSet();

        [Test]
        public void ApplyVariants_Combinatorics_Hetero_BothDeep_TwoVariants_AllSubsets()
        {
            // Two heterozygous variants, both ref & alt depths pass ? should produce 2^2 = 4 subsets
            var protein = MakeCombinatoricsProtein();
            int minDepth = 10;
            // Provide positions so ordering (desc) = varHigh then varLow
            var vcfHigh = BuildSingleSampleCombinatoricsVcf(18, "E", "K", "0/1", 12, 15);
            var vcfLow = BuildSingleSampleCombinatoricsVcf(7, "D", "N", "0/1", 11, 13);

            var varHigh = MakeVar(18, "E", "K", "bothDeep_high", vcfHigh);
            var varLow = MakeVar(7, "D", "N", "bothDeep_low", vcfLow);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varLow, varHigh }, // input order irrelevant; pipeline sorts descending
                maxAllowedVariantsForCombinatorics: 5,
                minAlleleDepth: minDepth);

            var sets = ProteoformVariantSetStrings(produced);

            // Expected subsets: "", "E18K", "D7N", "D7N|E18K"
            Assert.That(sets.Contains(""), Is.True);
            Assert.That(sets.Contains(varHigh.SimpleString()), Is.True);
            Assert.That(sets.Contains(varLow.SimpleString()), Is.True);
            Assert.That(sets.Contains($"{varLow.SimpleString()}|{varHigh.SimpleString()}") ||
                        sets.Contains($"{varHigh.SimpleString()}|{varLow.SimpleString()}"),
                Is.True, "Combined variant subset missing.");
            Assert.That(sets.Count, Is.EqualTo(4), "Unexpected number of combinatoric subsets for two variants.");
        }

        [Test]
        public void ApplyVariants_Combinatorics_Hetero_AltOnlyThenBothDeep()
        {
            // First variant: alt deep / ref shallow ? only alt path (replaces base with 1 proteoform)
            // Second variant: both deep ? combinatorics on existing proteoform (gives two subsets: with first only, with first+second)
            var protein = MakeCombinatoricsProtein();
            int minDepth = 10;

            var vcfAltOnly = BuildSingleSampleCombinatoricsVcf(16, "P", "L", "0/1", 3, 18);   // ref < minDepth, alt >= minDepth
            var vcfBoth = BuildSingleSampleCombinatoricsVcf(8, "T", "A", "0/1", 12, 14);  // both deep

            var varAltOnly = MakeVar(16, "P", "L", "altOnly_first", vcfAltOnly);
            var varBoth = MakeVar(8, "T", "A", "bothDeep_second", vcfBoth);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varBoth, varAltOnly }, // order doesn't matter; sorted descending => varAltOnly applied first
                maxAllowedVariantsForCombinatorics: 5,
                minAlleleDepth: minDepth);

            var sets = ProteoformVariantSetStrings(produced);

            // After first (alt-only): only one proteoform: "P16L"
            // After second (both deep): two proteoforms: "P16L" and "P16L|T8A"
            Assert.That(sets.Contains(""), Is.False, "Base should be replaced by alt-only first variant.");
            Assert.That(sets.Contains(varAltOnly.SimpleString()), Is.True, "First (alt-only) variant subset missing.");
            string combined = $"{varAltOnly.SimpleString()}|{varBoth.SimpleString()}";
            string combinedAlt = $"{varBoth.SimpleString()}|{varAltOnly.SimpleString()}";
            Assert.That(sets.Contains(combined) || sets.Contains(combinedAlt), Is.True,
                "Combined alt-only + both-deep variant subset missing.");
            Assert.That(sets.Count, Is.EqualTo(2),
                "Unexpected number of proteoforms after alt-only then both-deep application.");
        }

        [Test]
        public void ApplyVariants_Combinatorics_Hetero_AltShallow_SkipsVariant()
        {
            // First variant alt shallow / ref deep ? third internal branch (add only reference ppp) ? effectively skip
            // Second variant both deep ? classic combinatorics on original base
            var protein = MakeCombinatoricsProtein();
            int minDepth = 10;

            var vcfSkip = BuildSingleSampleCombinatoricsVcf(14, "A", "V", "0/1", 14, 5); // alt < minDepth -> isDeepAlternate=false
            var vcfBoth = BuildSingleSampleCombinatoricsVcf(6, "K", "R", "0/1", 11, 12);

            var varSkip = MakeVar(14, "A", "V", "altShallow_skip", vcfSkip);
            var varBoth = MakeVar(6, "K", "R", "bothDeep_apply", vcfBoth);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varSkip, varBoth }, // sorted desc -> varSkip first
                maxAllowedVariantsForCombinatorics: 5,
                minAlleleDepth: minDepth);

            var sets = ProteoformVariantSetStrings(produced);

            // varSkip never applied; only combinatorics of varBoth => "", "K6R"
            Assert.That(sets.Contains(varSkip.SimpleString()), Is.False,
                "Alt-shallow heterozygous variant should not appear in any proteoform.");
            Assert.That(sets.Contains(""), Is.True);
            Assert.That(sets.Contains(varBoth.SimpleString()), Is.True);
            Assert.That(sets.Count, Is.EqualTo(2));
        }

        [Test]
        public void ApplyVariants_Combinatorics_Hetero_MixedThreePaths()
        {
            // Three variants descending positions:
            // 1) Both deep (duplicating base) -> subsets: "" , A
            // 2) Alt-only (ref shallow) applies to all existing proteoforms -> subsets: A, A|B (base replaced by B alone)
            // 3) Alt-shallow (skip branch) - should not modify sets
            var protein = MakeCombinatoricsProtein();
            int minDepth = 10;

            var vcfBoth = BuildSingleSampleCombinatoricsVcf(20, "E", "K", "0/1", 11, 13); // both deep
            var vcfAltOnly = BuildSingleSampleCombinatoricsVcf(12, "P", "L", "0/1", 4, 16);  // alt-only
            var vcfSkip = BuildSingleSampleCombinatoricsVcf(5, "T", "S", "0/1", 12, 4);  // alt shallow

            var varBoth = MakeVar(20, "E", "K", "bothDeep20", vcfBoth);
            var varAltOnly = MakeVar(12, "P", "L", "altOnly12", vcfAltOnly);
            var varSkip = MakeVar(5, "T", "S", "skip5", vcfSkip);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varSkip, varBoth, varAltOnly }, // sorted desc => varBoth, varAltOnly, varSkip
                maxAllowedVariantsForCombinatorics: 5,
                minAlleleDepth: minDepth);

            var sets = ProteoformVariantSetStrings(produced);

            // Expect only: "P12L" (alt-only applied to base path) and "E20K|P12L"
            string keyAlt = varAltOnly.SimpleString();
            string keyBoth = varBoth.SimpleString();
            Assert.That(sets.Contains(keyAlt), Is.True,
                "Alt-only variant subset missing.");
            Assert.That(sets.Contains($"{keyBoth}|{keyAlt}") || sets.Contains($"{keyAlt}|{keyBoth}"),
                Is.True, "Combined bothDeep + altOnly subset missing.");
            Assert.That(sets.Contains(varSkip.SimpleString()), Is.False,
                "Alt-shallow variant (skip) should not appear.");
            Assert.That(sets.Contains(""), Is.False,
                "Base subset should have been replaced by alt-only mapping.");
            Assert.That(sets.Count, Is.EqualTo(2),
                "Unexpected number of proteoforms after mixed three-path scenario.");
        }

        [Test]
        public void ApplyVariants_Combinatorics_Hetero_RefOnlyBranch_AllRefsRetained()
        {
            // All three variants alt shallow (isDeepAlternate=false, ref deep)
            // Each should pass through without creating variant-applied proteoforms
            var protein = MakeCombinatoricsProtein();
            int minDepth = 10;

            var vcfSkipHigh = BuildSingleSampleCombinatoricsVcf(19, "M", "V", "0/1", 15, 5);
            var vcfSkipMid = BuildSingleSampleCombinatoricsVcf(11, "E", "D", "0/1", 14, 3);
            var vcfSkipLow = BuildSingleSampleCombinatoricsVcf(4, "A", "G", "0/1", 13, 4);

            var varHigh = MakeVar(19, "M", "V", "skipHigh", vcfSkipHigh);
            var varMid = MakeVar(11, "E", "D", "skipMid", vcfSkipMid);
            var varLow = MakeVar(4, "A", "G", "skipLow", vcfSkipLow);

            var produced = VariantApplication.ApplyVariants(
                protein,
                new[] { varLow, varMid, varHigh },
                maxAllowedVariantsForCombinatorics: 5,
                minAlleleDepth: minDepth);

            var sets = ProteoformVariantSetStrings(produced);

            // Only base proteoform expected
            Assert.That(sets.SetEquals(new[] { "" }), Is.True,
                "No variant should have been applied when alt depths are shallow for all heterozygous variants.");
        }

        #endregion
    }
}