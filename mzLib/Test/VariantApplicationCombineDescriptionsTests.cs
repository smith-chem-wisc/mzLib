using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.BioPolymer;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class VariantApplicationCombineDescriptionsTests
    {
        private static SequenceVariation MakeVar(int pos, string orig, string variant, string desc, string vcf = null)
            => new SequenceVariation(pos,
                pos + (orig?.Length > 0 ? orig.Length - 1 : 0),
                orig,
                variant,
                desc,
                vcf);

        private static string DeriveToken(SequenceVariation v)
        {
            if (v == null) return null;
            // VCF precedence
            if (v.VariantCallFormatData?.Description is string d) return d;
            // Fallback: Description if non-whitespace, else SimpleString
            return string.IsNullOrWhiteSpace(v.Description) ? v.SimpleString() : v.Description;
        }

        private static List<string> ExpectedTokens(IEnumerable<SequenceVariation> vars) =>
            vars?
                .Where(v => v != null)
                .Select(DeriveToken)
                .Where(s => !string.IsNullOrWhiteSpace(s))
                .Distinct()
                .Take(10)
                .ToList()
            ?? new List<string>();

        [Test]
        public void CombineDescriptions_Comprehensive()
        {
            // Shared VCF token (duplicate across v1 & v3)
            string tokenA_Vcf =
                "1\t100\t.\tA\tG\t.\tPASS\tANN=A|missense|X|GENE|\tGT:AD:DP\t0/1:10,12:22";
            // Second distinct VCF token (v6)
            string tokenE_Vcf =
                "1\t200\t.\tC\tT\t.\tPASS\tANN=C|synonymous|Y|GENE2|\tGT:AD:DP\t0/1:5,9:14";

            // 12 variants:
            // v1: VCF token A (preempts description)
            var v1 = MakeVar(10, "M", "V", "DescIgnoredByVCF", tokenA_Vcf);
            // v2: Plain description (B)
            var v2 = MakeVar(20, "P", "A", "B_desc");
            // v3: Duplicate VCF token A (must deduplicate)
            var v3 = MakeVar(30, "K", "R", "AnotherIgnored", tokenA_Vcf);
            // v4: Whitespace description but real change (insertion) -> fallback to SimpleString
            var v4 = MakeVar(40, "L", "LL", "   ");
            // v5: Plain description (D)
            var v5 = MakeVar(50, "S", "T", "D_desc");
            // v6: Second VCF token (E)
            var v6 = MakeVar(60, "Q", "E", "IgnoredVCF2", tokenE_Vcf);
            // v7: NEW unique description (X13) to push unique count above 10
            var v7 = MakeVar(70, "A", "G", "X13");
            var v8 = MakeVar(80, "R", "K", "X8");
            var v9 = MakeVar(90, "H", "Y", "X9");
            var v10 = MakeVar(100, "N", "D", "X10");
            var v11 = MakeVar(110, "F", "S", "X11");
            var v12 = MakeVar(120, "C", "W", "X12"); // 11th unique token (should be truncated out)

            var all = new List<SequenceVariation>
            {
                v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12
            };

            // Subsets
            List<SequenceVariation> subsetNull = null;
            var subsetEmpty = new List<SequenceVariation>();
            var subset1 = all.Take(1).ToList();   // 1
            var subset5 = all.Take(5).ToList();   // up to v5
            var subset10 = all.Take(10).ToList();  // up to v10
            var subset12 = all.ToList();           // full set

            // 0 (null)
            Assert.That(VariantApplication.CombineDescriptions(subsetNull), Is.EqualTo(string.Empty));
            // 0 (empty)
            Assert.That(VariantApplication.CombineDescriptions(subsetEmpty), Is.EqualTo(string.Empty));

            // 1
            var expected1 = ExpectedTokens(subset1);
            var got1 = VariantApplication.CombineDescriptions(subset1);
            Assert.That(got1, Is.EqualTo(expected1.Single()));
            Assert.That(got1.Contains(", variant:"), Is.False);

            // 5
            var expected5 = ExpectedTokens(subset5);
            var got5 = VariantApplication.CombineDescriptions(subset5);
            var tokens5 = got5.Split(new[] { ", variant:" }, StringSplitOptions.None);
            Assert.That(tokens5.Length, Is.EqualTo(expected5.Count));
            CollectionAssert.AreEqual(expected5, tokens5);

            // 10
            var expected10 = ExpectedTokens(subset10);
            var got10 = VariantApplication.CombineDescriptions(subset10);
            var tokens10 = got10.Split(new[] { ", variant:" }, StringSplitOptions.None);
            Assert.That(tokens10.Length, Is.EqualTo(expected10.Count));
            Assert.That(tokens10.Length, Is.LessThanOrEqualTo(10));
            CollectionAssert.AreEqual(expected10, tokens10);

            // 12 (trigger truncation: 11 distinct -> keep first 10)
            var expected12 = ExpectedTokens(subset12); // already applies Distinct().Take(10)
            var got12 = VariantApplication.CombineDescriptions(subset12);
            var tokens12 = got12.Split(new[] { ", variant:" }, StringSplitOptions.None);
            Assert.That(tokens12.Length, Is.EqualTo(expected12.Count));
            Assert.That(tokens12.Length, Is.EqualTo(10), "Should truncate to 10 tokens when >10 unique encountered.");
            CollectionAssert.AreEqual(expected12, tokens12, "Truncated token ordering/content mismatch.");

            // Branch / behavior verifications:

            // VCF precedence: Description ignored when VCF present
            Assert.That(DeriveToken(v1), Is.EqualTo(tokenA_Vcf));
            // Duplicate VCF token only once after distinct
            Assert.That(expected12.Count(t => t == tokenA_Vcf), Is.EqualTo(1));

            // Whitespace description fallback (v4)
            Assert.That(string.IsNullOrWhiteSpace(v4.Description), Is.True);
            Assert.That(expected12.Contains(v4.SimpleString()), Is.True, "Whitespace fallback token missing.");

            // Truncation: ensure last distinct (X12) excluded (since it would be the 11th)
            var fullDistinct = all.Select(DeriveToken)
                                  .Where(s => !string.IsNullOrWhiteSpace(s))
                                  .Distinct()
                                  .ToList();
            if (fullDistinct.Count > 10)
            {
                var eleventh = fullDistinct[10];
                Assert.That(tokens12.Contains(eleventh), Is.False, "11th token should be truncated.");
            }
        }
    }
}