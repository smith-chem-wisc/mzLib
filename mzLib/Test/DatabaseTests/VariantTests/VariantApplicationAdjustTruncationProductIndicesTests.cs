using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.BioPolymer;
using Proteomics;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.DatabaseTests.VariantTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class VariantApplicationAdjustTruncationProductIndicesTests
    {
        private static readonly MethodInfo AdjustMethod =
            typeof(VariantApplication).GetMethod("AdjustTruncationProductIndices",
                BindingFlags.NonPublic | BindingFlags.Static)
            ?? throw new InvalidOperationException("Could not locate AdjustTruncationProductIndices via reflection.");

        private static List<TruncationProduct> InvokeAdjust(
            SequenceVariation variant,
            string variantAppliedSequence,
            Protein protein,
            IEnumerable<TruncationProduct> products)
        {
            return (List<TruncationProduct>)AdjustMethod.Invoke(
                null,
                new object[] { variant, variantAppliedSequence, protein, products })!;
        }

        private Protein MakeProtein(string accession = "BASE") => new Protein("MPEPTIDESEQX", accession); // length 12

        private static SequenceVariation MakeVar(int begin, string original, string variant, string desc) =>
            new SequenceVariation(begin,
                                  begin + original.Length - 1,
                                  original,
                                  variant,
                                  desc,
                                  variantCallFormatDataString: null,
                                  oneBasedModifications: null);

        // Light coverage test already added previously (left for context)
        [Test]
        public void AdjustTruncationProducts_LightCoverage_InsertionAndStopGain()
        {
            var baseProducts = new List<TruncationProduct>
            {
                new TruncationProduct(1, 3, "before"),
                new TruncationProduct(2, 10, "spanning"),
                new TruncationProduct(8, 12, "after"),
                new TruncationProduct(1, 12, "full")
            };

            // Insertion (+2)
            var prot = MakeProtein("INS");
            var insVar = MakeVar(5, "TI", "TAAI", "Insertion");
            string appliedIns = "MPEPTAAIDESEQX"; // length 14
            var adjustedIns = InvokeAdjust(insVar, appliedIns, prot, baseProducts);
            Assert.Contains(new TruncationProduct(1, 3, "before"), adjustedIns);
            Assert.Contains(new TruncationProduct(2, 12, "spanning"), adjustedIns);
            Assert.Contains(new TruncationProduct(10, 14, "after"), adjustedIns);
            Assert.Contains(new TruncationProduct(1, 14, "full"), adjustedIns);

            // Stop gain
            var protStop = MakeProtein("STOP");
            var stopVar = MakeVar(5, "TIDES", "T*", "StopGain");
            string appliedStop = "MPEPT"; // truncated at stop (len 5)
            var adjustedStop = InvokeAdjust(stopVar, appliedStop, protStop, baseProducts);
            NUnit.Framework.Assert.That(adjustedStop.Count, Is.EqualTo(3));
            Assert.Contains(new TruncationProduct(1, 3, "before"), adjustedStop);
            Assert.Contains(new TruncationProduct(2, 5, "spanning"), adjustedStop);
            Assert.Contains(new TruncationProduct(1, 5, "full"), adjustedStop);
        }

        // ========= Targeted branch tests for the specified if / else-if block ==========
        [Test]
        public void TruncationProducts_Branch_EntirelyBeforeVariant_Unchanged()
        {
            var prot = MakeProtein("BEFORE");
            // Variant starts at position 8 (ESEQ -> KSEQ) – valid substitution (not a no-op)
            var variant = MakeVar(8, "ESEQ", "KSEQ", "Substitution");

            // Apply change (replace residue at 8 with K, keep rest)
            string applied = prot.BaseSequence.Substring(0, 7) + "K" + prot.BaseSequence.Substring(8);

            var products = new List<TruncationProduct>
            {
                new TruncationProduct(1,5,"before"),   // entirely before variant region (positions 8–11)
                new TruncationProduct(2,11,"spanning"),
                new TruncationProduct(9,12,"after")
            };

            var adjusted = InvokeAdjust(variant, applied, prot, products);

            Assert.Contains(new TruncationProduct(1, 5, "before"), adjusted,
                "Product before variant should be retained unchanged.");
        }

        [Test]
        public void TruncationProducts_Branch_Spanning_StopGain_AdjustsToNewLength()
        {
            var prot = MakeProtein("SPAN_STOP");
            // Replace positions 5-7 (TID) with "A*" ? new sequence truncated to prefix + 'A' (positions 1..5)
            var variant = MakeVar(5, "TID", "A*", "Stop");
            string applied = prot.BaseSequence.Substring(0, 4) + "A"; // length = 5

            var spanning = new TruncationProduct(2, 11, "span");
            var products = new List<TruncationProduct> { spanning };

            var adjusted = InvokeAdjust(variant, applied, prot, products);

            // Expect new product from original begin to new truncated protein length
            NUnit.Framework.Assert.That(adjusted.Count, Is.EqualTo(1));
            Assert.Contains(new TruncationProduct(2, applied.Length, "span"), adjusted);
        }

        [Test]
        public void TruncationProducts_Branch_Spanning_Insertion_ShiftsEnd()
        {
            var prot = MakeProtein("SPAN_INS");
            // Insertion: positions 5-6 (TI) -> TAAI (+2)
            var variant = MakeVar(5, "TI", "TAAI", "Insertion");
            string applied = prot.BaseSequence.Substring(0, 4) + "TAAI" + prot.BaseSequence.Substring(6); // length 14

            var spanning = new TruncationProduct(2, 10, "span");
            var products = new List<TruncationProduct> { spanning };

            var adjusted = InvokeAdjust(variant, applied, prot, products);
            // End should shift +2: 10 -> 12
            Assert.Contains(new TruncationProduct(2, 12, "span"), adjusted);
        }

        [Test]
        public void TruncationProducts_Branch_FullLengthProduct_LeftClause_BeginEquals1()
        {
            var prot = MakeProtein("FULL_BEGIN1");
            // Simple substitution mid-protein (positions 6-6)
            var variant = MakeVar(6, "I", "K", "Sub");
            string applied = prot.BaseSequence.Substring(0, 5) + "K" + prot.BaseSequence.Substring(6);

            var full = new TruncationProduct(1, prot.BaseSequence.Length, "full");
            var products = new List<TruncationProduct> { full };

            var adjusted = InvokeAdjust(variant, applied, prot, products);
            // No length change, so end unchanged
            Assert.Contains(new TruncationProduct(1, prot.BaseSequence.Length, "full"), adjusted);
        }

        [Test]
        public void TruncationProducts_Branch_FullLengthProduct_LeftClause_BeginEquals2()
        {
            var prot = MakeProtein("FULL_BEGIN2");
            // Variant internal substitution
            var variant = MakeVar(7, "D", "N", "Sub");
            string applied = prot.BaseSequence.Substring(0, 6) + "N" + prot.BaseSequence.Substring(7);

            var fullFrom2 = new TruncationProduct(2, prot.BaseSequence.Length, "full2");
            var products = new List<TruncationProduct> { fullFrom2 };

            var adjusted = InvokeAdjust(variant, applied, prot, products);
            Assert.Contains(new TruncationProduct(2, prot.BaseSequence.Length, "full2"), adjusted);
        }

        [Test]
        public void TruncationProducts_Branch_LeftSideViaEquality_BeginEquals2_VariantStartsAt2()
        {
            var prot = MakeProtein("BEGIN_EQ2");
            // Variant starts at position 2 (P->L, single AA)
            var variant = MakeVar(2, "P", "L", "Sub");
            string applied = "ML" + prot.BaseSequence.Substring(2); // length unchanged

            var product = new TruncationProduct(2, prot.BaseSequence.Length, "edge");
            var products = new List<TruncationProduct> { product };

            var adjusted = InvokeAdjust(variant, applied, prot, products);
            Assert.Contains(new TruncationProduct(2, prot.BaseSequence.Length, "edge"), adjusted);
        }

        [Test]
        public void TruncationProducts_Branch_Spanning_NoStop_FullEndCondition()
        {
            var prot = MakeProtein("SPAN_FULLEND");
            // Variant internal substitution positions 5-7 "TID" -> "KID" (length unchanged)
            var variant = MakeVar(5, "TID", "KID", "Sub");
            string applied = prot.BaseSequence.Substring(0, 4) + "KID" + prot.BaseSequence.Substring(7);

            // Product begins before variant (2) and extends to full length (end == base length) satisfying right side via equality.
            var product = new TruncationProduct(2, prot.BaseSequence.Length, "span_full");
            var products = new List<TruncationProduct> { product };

            var adjusted = InvokeAdjust(variant, applied, prot, products);
            Assert.Contains(new TruncationProduct(2, prot.BaseSequence.Length, "span_full"), adjusted);
        }
        // (append to existing file)

        #region After-Variant Branch Tests (final else-if)

        // Helpers local to this region
        private Protein MakeProteinCustom(string seq, string acc) => new Protein(seq, acc);

        [Test]
        public void AfterVariant_Substitution_NoLengthChange_ShiftZero()
        {
            // Base length 12
            var prot = MakeProtein("AFTER_SUB_ZERO");
            // Variant: single AA substitution at position 5 (T->A), length change = 0
            var variant = MakeVar(5, "T", "A", "Sub");
            string applied = prot.BaseSequence.Substring(0, 4) + "A" + prot.BaseSequence.Substring(5); // length unchanged

            // Product entirely after variant (variant spans 5..5; product starts 7)
            var productAfter = new TruncationProduct(7, 12, "after");
            var adjusted = InvokeAdjust(variant, applied, prot, new[] { productAfter });

            // lengthChange = 0 ? coordinates unchanged
            NUnit.Framework.Assert.That(adjusted, Has.Count.EqualTo(1));
            Assert.Contains(new TruncationProduct(7, 12, "after"), adjusted);
        }

        [Test]
        public void AfterVariant_Insertion_PositiveShift()
        {
            var prot = MakeProtein("AFTER_INS");
            // Insertion at 5-6: "TI" -> "TAAAI" (+3 length; original TI len=2, inserted len=5; +3)
            var variant = MakeVar(5, "TI", "TAAAI", "Insertion");
            string applied = prot.BaseSequence.Substring(0, 4) + "TAAAI" + prot.BaseSequence.Substring(6); // new length 12+3=15
            int lengthChange = 3;

            // Product after variant (variant end = 6). Pick original coordinates 8-12.
            var productAfter = new TruncationProduct(8, 12, "after");
            var adjusted = InvokeAdjust(variant, applied, prot, new[] { productAfter });

            // Expect begin/end shifted forward by +3
            NUnit.Framework.Assert.That(adjusted, Has.Count.EqualTo(1));
            Assert.Contains(new TruncationProduct(8 + lengthChange, 12 + lengthChange, "after"), adjusted);
        }

        [Test]
        public void AfterVariant_Deletion_NegativeShift()
        {
            var prot = MakeProtein("AFTER_DEL");
            // Deletion at 5-6: "TI" -> ""  (length change = -2)
            var variant = MakeVar(5, "TI", "", "Deletion");
            string applied = prot.BaseSequence.Remove(4, 2); // remove indices 4..5 (0-based) => new length 10
            int lengthChange = -2;

            // Product after variant (variant end=6). Original product 8-12.
            var productAfter = new TruncationProduct(8, 12, "after");
            var adjusted = InvokeAdjust(variant, applied, prot, new[] { productAfter });

            // Shift backward by 2: 8->6, 12->10
            NUnit.Framework.Assert.That(adjusted, Has.Count.EqualTo(1));
            Assert.Contains(new TruncationProduct(6, 10, "after"), adjusted);
        }

        [Test]
        public void AfterVariant_StopGain_NotAdded()
        {
            var prot = MakeProtein("AFTER_STOP");
            // Stop gain at 5-8: replace "TIDE" with "T*" (truncation).
            var variant = MakeVar(5, "TIDE", "T*", "Stop");
            string applied = prot.BaseSequence.Substring(0, 4) + "T"; // truncated length = 5

            // Product originally after variant (variant end = 8) -> choose 9-12
            var productAfter = new TruncationProduct(9, 12, "after");

            var adjusted = InvokeAdjust(variant, applied, prot, new[] { productAfter });

            // Since variant introduces stop (*), after-variant products are NOT added.
            NUnit.Framework.Assert.That(adjusted, Is.Empty);
        }

        [Test]
        public void AfterVariant_NotStrictlyAfter_FirstConditionFails_NotAdded()
        {
            var prot = MakeProtein("AFTER_FAIL");
            // Substitution at 8-9: "ES" -> "KS"
            var variant = MakeVar(8, "ES", "KS", "Sub");
            string applied = prot.BaseSequence.Substring(0, 7) + "KS" + prot.BaseSequence.Substring(9);

            // Product begins at 9 (variant spans 8..9); condition requires begin > variant end (9 > 9 false)
            var productAdjacent = new TruncationProduct(9, 12, "adjacent");

            var adjusted = InvokeAdjust(variant, applied, prot, new[] { productAdjacent });

            NUnit.Framework.Assert.That(adjusted, Is.Empty, "Product starting at variant end should not be treated as strictly after variant.");
        }

        [Test]
        public void AfterVariant_MultipleProducts_Mixed_AddsOnlyAfterOnes()
        {
            var prot = MakeProtein("AFTER_MIX");
            // Insertion at 5-6: "TI" -> "TIQQ" (+2 length change)
            var variant = MakeVar(5, "TI", "TIQQ", "Insertion");
            string applied = prot.BaseSequence.Substring(0, 4) + "TIQQ" + prot.BaseSequence.Substring(6); // new length +2
            int lengthChange = 2;

            var products = new List<TruncationProduct>
            {
                // This product straddles the variant (begin < variantBegin AND end > variantEnd) so it qualifies for the
                // second (straddling) branch and will have only its end extended by +lengthChange.
                new TruncationProduct(3,7,"straddling"),
                // These two are strictly after the variant (variant end = 6) and will be shifted by +lengthChange.
                new TruncationProduct(8,10,"after1"),
                new TruncationProduct(9,12,"after2")
            };

            var adjusted = InvokeAdjust(variant, applied, prot, products);

            // Expect:
            // straddling: (3, 7+2) = (3,9)
            // after1: (8+2, 10+2) = (10,12)
            // after2: (9+2, 12+2) = (11,14)
            NUnit.Framework.Assert.That(adjusted.Count, Is.EqualTo(3), "Straddling product is also retained and adjusted.");
            Assert.Contains(new TruncationProduct(3, 9, "straddling"), adjusted);
            Assert.Contains(new TruncationProduct(10, 12, "after1"), adjusted);
            Assert.Contains(new TruncationProduct(11, 14, "after2"), adjusted);

            // Sanity: none of the original (unadjusted) coordinates should appear
            Assert.False(adjusted.Any(p => p.OneBasedBeginPosition == 3 && p.OneBasedEndPosition == 7 && p.Type == "straddling"));
            Assert.False(adjusted.Any(p => p.OneBasedBeginPosition == 8 && p.OneBasedEndPosition == 10));
            Assert.False(adjusted.Any(p => p.OneBasedBeginPosition == 9 && p.OneBasedEndPosition == 12));
        }

        #endregion
    }
}