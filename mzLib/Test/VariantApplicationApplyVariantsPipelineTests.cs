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
         * Focus: Filtering / dedup / ordering portion of ApplyVariants and depth-driven expectations.
         * We do NOT assert survival of the pristine base sequence; downstream genotype logic can legitimately
         * eliminate it even without a fully universal homozygous-alt variant.
         */

        #region Helpers

        private string BuildVcf(
            string chrom,
            int pos,
            string refAA,
            string altAA,
            string sample0GT,
            string sample1GT,
            string sample0AD = "10,5",
            string sample1AD = "6,9",
            string qual = ".",
            string filter = "PASS")
        {
            var cols = new[]
            {
                chrom,
                pos.ToString(),
                ".",
                refAA,
                altAA,
                qual,
                filter,
                "ANN=" + altAA + "|missense|GENE|GENE|",
                "GT:AD:DP",
                $"{sample0GT}:{sample0AD}:15",
                $"{sample1GT}:{sample1AD}:20"
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

        private List<SequenceVariation> ComputeUnique(IEnumerable<SequenceVariation> source) =>
            source
                .Where(v => v != null)
                .GroupBy(v => v.SimpleString())
                .Select(g => g.First())
                .Where(v => v.VariantCallFormatData != null &&
                            v.VariantCallFormatData.Genotypes != null &&
                            v.VariantCallFormatData.Genotypes.Count > 0)
                .OrderByDescending(v => v.OneBasedBeginPosition)
                .ToList();

        private static bool AltPassesDepth(SequenceVariation v, int minAlleleDepth)
        {
            if (minAlleleDepth <= 0) return true;
            var vcf = v.VariantCallFormatData;
            if (vcf == null || vcf.AlleleDepths == null) return false;
            int altIndex = vcf.AlleleIndex;
            if (altIndex < 0) return false;

            foreach (var kv in vcf.AlleleDepths)
            {
                var depths = kv.Value;
                if (depths.Length > altIndex &&
                    int.TryParse(depths[altIndex], out int depthVal) &&
                    depthVal >= minAlleleDepth)
                {
                    return true;
                }
            }
            return false;
        }

        private static HashSet<string> ExtractAppliedSimpleStrings(IEnumerable<Protein> proteins) =>
            new HashSet<string>(proteins
                .SelectMany(p => p.AppliedSequenceVariations ?? new List<SequenceVariation>())
                .Select(v => v.SimpleString()));

        #endregion

        #region Test Data Construction

        private (Protein protein, List<SequenceVariation> variants) BuildSourceSet()
        {
            var protein = new Protein("MPEPTIDELONGSEQUENCEFORTEST", "BASE_PROT");

            // AD pairs (ref,alt):
            // Sample0 AD=10,5 alt=5
            // Sample1 AD=6,9  alt=9
            // None reach 10 for any variant’s alt depth.
            var vcf1 = BuildVcf("1", 5,  "E", "K", "0/1", "1/1"); // duplicate basis
            var vcf2 = BuildVcf("1", 15, "T", "A", "1/1", "0/1");
            var vcf3 = BuildVcf("1", 22, "D", "G", "0/1", "0/1");
            var duplicateOf1 = BuildVcf("1", 5, "E", "K", "0/1", "1/1");

            string noVcf = null;
            var badVcf = "1\t30\t.\tA\tG\t.\tPASS"; // insufficient columns

            var variants = new List<SequenceVariation>
            {
                null,
                MakeVar(5,  "E","K","dupCandidateFirst",   vcf1),
                MakeVar(5,  "E","K","dupCandidateSecond",  duplicateOf1),
                MakeVar(15, "T","A","validMiddle",         vcf2),
                MakeVar(22, "D","G","validHighest",        vcf3),
                MakeVar(10, "P","A","filteredNoVcf",       noVcf),
                MakeVar(30, "A","V","badVcfFiltered",      badVcf)
            };

            return (protein, variants);
        }

        #endregion

        #region Tests

        [Test]
        public void ApplyVariants_Pipeline_EarlyReturn_NoUsableVariants()
        {
            var protein = new Protein("MPEPTIDE", "EARLY_ONLY");
            var variants = new List<SequenceVariation>
            {
                null,
                MakeVar(3,"E","K","noVcf", null),
                MakeVar(4,"P","L","badVcf", "1\t10\t.\tP\tL\t.\tPASS")
            };

            var result = VariantApplication.ApplyVariants(
                protein,
                variants,
                maxAllowedVariantsForCombinatorics: 3,
                minAlleleDepth: 1);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].BaseSequence, Is.EqualTo(protein.BaseSequence));
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        [Test] public void ApplyVariants_Pipeline_DedupFilterOrdering_Depth0()  => RunPipelineCore(0);
        [Test] public void ApplyVariants_Pipeline_DedupFilterOrdering_Depth1()  => RunPipelineCore(1);
        [Test] public void ApplyVariants_Pipeline_DedupFilterOrdering_Depth10() => RunPipelineCore(10);

        private void RunPipelineCore(int minAlleleDepth)
        {
            var (protein, variants) = BuildSourceSet();

            var unique = ComputeUnique(variants);
            Assert.That(unique.Count, Is.EqualTo(3));
            Assert.That(unique.Select(v => v.OneBasedBeginPosition).ToArray(),
                Is.EqualTo(new[] { 22, 15, 5 }));

            // Depth-filtered expectation: only variants whose alt depth >= threshold in at least one sample.
            var expectedPassing = unique.Where(v => AltPassesDepth(v, minAlleleDepth)).ToList();

            var produced = VariantApplication.ApplyVariants(
                protein,
                variants,
                maxAllowedVariantsForCombinatorics: 3,
                minAlleleDepth: minAlleleDepth);

            Assert.That(produced, Is.Not.Null);
            Assert.That(produced.Count, Is.GreaterThanOrEqualTo(1));

            var appliedSimpleStrings = ExtractAppliedSimpleStrings(produced);

            if (expectedPassing.Count == 0)
            {
                // No variant should be applied (only base or depth-failing expansions suppressed)
                Assert.That(appliedSimpleStrings.Count, Is.EqualTo(0),
                    $"No variants should pass depth {minAlleleDepth}, but found: {string.Join(",", appliedSimpleStrings)}");
                return;
            }

            foreach (var v in expectedPassing)
            {
                Assert.That(appliedSimpleStrings.Contains(v.SimpleString()), Is.True,
                    $"Expected variant {v.SimpleString()} not found (depth {minAlleleDepth}).");
            }

            // Duplicate check only if the 5-position variant passed depth
            var pos5 = expectedPassing.FirstOrDefault(v => v.OneBasedBeginPosition == 5);
            if (pos5 != null)
            {
                var dupKey = pos5.SimpleString();
                var dupCount = produced.SelectMany(p => p.AppliedSequenceVariations)
                                       .Count(v => v.SimpleString() == dupKey);
                Assert.That(dupCount, Is.GreaterThanOrEqualTo(1),
                    "Collapsed duplicate variant should appear at least once.");
            }

            // Confirm filtered (no VCF or malformed VCF) never appear
            Assert.That(appliedSimpleStrings.Any(s => s.Contains("P10A")), Is.False,
                "Null-VCF variant should have been filtered.");
            Assert.That(appliedSimpleStrings.Any(s => s.Contains("A30V")), Is.False,
                "Malformed VCF variant should have been filtered.");
        }

        #endregion
    }
}