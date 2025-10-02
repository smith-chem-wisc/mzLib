using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.BioPolymer;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public partial class VariantApplicationAdjustSequenceVariationIndicesTests
    {
        // (Reuse existing reflection + helpers if this is appended to previous file.
        // If file is standalone, duplicate helper definitions.)

        private static readonly MethodInfo AdjustMethod2 =
            typeof(VariantApplication).GetMethod("AdjustSequenceVariationIndices",
                BindingFlags.NonPublic | BindingFlags.Static)
            ?? throw new InvalidOperationException("Could not locate AdjustSequenceVariationIndices via reflection.");

        private static List<SequenceVariation> InvokeAdjust2(
            SequenceVariation variantGettingApplied,
            string variantAppliedProteinSequence,
            IEnumerable<SequenceVariation> alreadyApplied)
        {
            return (List<SequenceVariation>)AdjustMethod2.Invoke(
                null,
                new object[] { variantGettingApplied, variantAppliedProteinSequence, alreadyApplied })!;
        }

        private static SequenceVariation MkVar(int begin, string original, string variant, string desc) =>
            new SequenceVariation(begin,
                begin + (original?.Length ?? 0) - 1,
                original,
                variant,
                desc,
                variantCallFormatDataString: null,
                oneBasedModifications: null);

        [Test]
        public void AdjustSequenceVariationIndices_NullCollection_ReturnsEmpty()
        {
            // Applied variant (simple substitution)
            var applied = MkVar(4, "P", "K", "Applied_P4K");
            string newSeq = "MPEK TIDESEQX".Replace(" ", ""); // base mutated (original base assumed MPEPTIDESEQX)

            var result = InvokeAdjust2(applied, newSeq, null);

            Assert.That(result, Is.Empty, "Expected empty list when alreadyAppliedVariations is null.");
        }

        [Test]
        public void AdjustSequenceVariationIndices_ContainsNulls_SkipsThem()
        {
            // Base sequence
            const string baseSeq = "MPEPTIDESEQX";

            // Applied variant: insertion (I -> IL)
            var applied = MkVar(6, "I", "IL", "Applied_Insertion_I6IL");
            string mutated = baseSeq.Substring(0, 5) + "IL" + baseSeq.Substring(6); // length +1

            // Another existing variation (after region) to ensure normal processing path (not reference equal)
            var other = MkVar(10, "E", "Q", "Other_E10Q");

            var list = new List<SequenceVariation>
            {
                null,
                applied,   // reference equal -> should be added directly and continue
                null,
                other
            };

            var result = InvokeAdjust2(applied, mutated, list);

            // Nulls skipped
            Assert.That(result.Count, Is.EqualTo(2));

            var appliedOut = result.Single(v => v.Description == "Applied_Insertion_I6IL");
            Assert.That(ReferenceEquals(appliedOut, applied), Is.True, "Applied variant should be added unchanged by reference.");

            var otherOut = result.Single(v => v.Description == "Other_E10Q");
            // After insertion (+1), original 10 shifts to 11 (no overlap, no subtraction)
            Assert.That(otherOut.OneBasedBeginPosition, Is.EqualTo(11));
            Assert.That(otherOut.OneBasedEndPosition, Is.EqualTo(11));
        }
        [Test]
        public void AdjustSequenceVariationIndices_VariantNotInList_NoReferenceEquality()
        {
            const string baseSeq = "MPEPTIDESEQX";

            // Variant getting applied (substitution) NOT present in alreadyApplied list
            var applied = MkVar(3, "E", "K", "Applied_E3K");
            string mutated = "MPK" + baseSeq.Substring(3); // position 3 altered

            // Only unrelated variants
            var v1 = MkVar(1, "M", "A", "Other_M1A");     // before applied
            var v2 = MkVar(8, "E", "Q", "Other_E8Q");     // after applied
            var v3 = MkVar(3, "E", "L", "Overlap_Alt");   // overlaps applied coordinates but is a different object

            var list = new List<SequenceVariation> { v1, v2, v3 };

            var result = InvokeAdjust2(applied, mutated, list);

            Assert.That(result.Count, Is.EqualTo(3));
            Assert.That(result.Any(v => v.Description == "Applied_E3K"), Is.False);

            var r1 = result.Single(v => v.Description == "Other_M1A");
            Assert.That(r1.OneBasedBeginPosition, Is.EqualTo(1));
            Assert.That(r1.OneBasedEndPosition, Is.EqualTo(1));

            var r2 = result.Single(v => v.Description == "Other_E8Q");
            Assert.That(r2.OneBasedBeginPosition, Is.EqualTo(8));
            Assert.That(r2.OneBasedEndPosition, Is.EqualTo(8));

            var r3 = result.Single(v => v.Description == "Overlap_Alt");
            // Because the overlapping region (position 3) is shared with the applied variant and overlap=1,
            // the algorithm shifts begin/end: new = old + seqLenChange (0) - overlap (1) => 2.
            Assert.That(r3.OneBasedBeginPosition, Is.EqualTo(2));
            Assert.That(r3.OneBasedEndPosition, Is.EqualTo(2));
        }

        [Test]
        public void AdjustSequenceVariationIndices_AllNullExceptApplied()
        {
            const string baseSeq = "MPEPTIDESEQX";
            var applied = MkVar(7, "D", "N", "Applied_D7N");
            string mutated = baseSeq.Substring(0, 6) + "N" + baseSeq.Substring(7);

            var list = new List<SequenceVariation> { null, applied, null };

            var result = InvokeAdjust2(applied, mutated, list);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(ReferenceEquals(result[0], applied), Is.True);
            Assert.That(result[0].OneBasedBeginPosition, Is.EqualTo(7));
            Assert.That(result[0].OneBasedEndPosition, Is.EqualTo(7));
        }
        #region Branch tests: sameVcfRecord / effective-before (addedIdx) early-continue logic

        private static SequenceVariation MkVarVcf(int begin, string orig, string varSeq, string desc, string vcfLine) =>
            new SequenceVariation(begin,
                begin + (orig?.Length ?? 0) - 1,
                orig,
                varSeq,
                desc,
                variantCallFormatDataString: vcfLine,
                oneBasedModifications: null);

        [Test]
        public void AdjustSequenceVariationIndices_SameVcfRecord_AddedUnmodified()
        {
            // Applied variant with VCF
            string vcf = "1\t1000\trs1\tA\tT\t.\tPASS\tANN=.\tGT:AD\t0/1:10,8";
            var applied = MkVarVcf(6, "I", "K", "Applied_I6K", vcf);
            string mutated = "MPEPTKDESEQX"; // base after substitution

            // Another variant sharing identical VCF record but overlapping AND after (so second condition would be false if evaluated)
            var sameVcfDifferentCoords = MkVarVcf(8, "E", "Q", "SameVCF_E8Q", vcf);

            var list = new List<SequenceVariation> { applied, sameVcfDifferentCoords };

            var result = InvokeAdjust2(applied, mutated, list);

            // Expect both variants present, the second added via sameVcfRecord early path (no coordinate shift)
            var outVar = result.Single(v => v.Description == "SameVCF_E8Q");
            Assert.That(outVar.OneBasedBeginPosition, Is.EqualTo(sameVcfDifferentCoords.OneBasedBeginPosition));
            Assert.That(outVar.OneBasedEndPosition, Is.EqualTo(sameVcfDifferentCoords.OneBasedEndPosition));
        }
        [Test]
        public void AdjustSequenceVariationIndices_EntirelyBefore_AfterPositiveAddedIdx_AddedUnmodified()
        {
            // Applied variant later in sequence (position 10)
            var applied = MkVar(10, "E", "Q", "Applied_E10Q");

            // Earlier insertion that contributes positive length change (+2) fully before 'beforeVariant'.
            // IMPORTANT: For an insertion, use the single?position constructor with originalSequence = null
            // so that (variant length - original length) contributes correctly and coordinates are valid.
            var earlyInsertion = new SequenceVariation(2, null, "AA", "Ins_Pos2");

            // Variant before applied (ends at 5). addedIdx from earlyInsertion = +2.
            // Effective end for comparison logic: 5 - 2 = 3 which is < applied begin (10) ? early-continue path.
            var beforeVariant = MkVar(5, "T", "A", "Before_T5A");

            // Mutated sequence reflecting the insertion (length base 12 + 2 = 14) and later substitution at pos 10.
            // Base: M P E P T I D E S E Q X
            // After insertion at pos2: M A A P E P T I D E S E Q X
            // After substitution at pos10 (E->Q): M A A P E P T I D Q S E Q X
            string mutated = "MAAPEPTIDQSEQX";

            var list = new List<SequenceVariation> { earlyInsertion, beforeVariant, applied };

            var result = InvokeAdjust2(applied, mutated, list);

            var outBefore = result.Single(v => v.Description == "Before_T5A");
            Assert.That(outBefore.OneBasedBeginPosition, Is.EqualTo(beforeVariant.OneBasedBeginPosition));
            Assert.That(outBefore.OneBasedEndPosition, Is.EqualTo(beforeVariant.OneBasedEndPosition));
        }

        [Test]
        public void AdjustSequenceVariationIndices_DeletionEarlier_NegativeAddedIdx_ForcesAdjustPath()
        {
            // Applied variant starts at 8
            var applied = MkVar(8, "E", "K", "Applied_E8K");

            // Earlier deletion spanning positions 2-4 (orig 'PEP' -> '') length change -3
            var earlyDeletion = MkVar(2, "PEP", "", "Del_2_4");

            // Overlapping candidate variant whose end is not strictly before applied when adjusted (should NOT early-continue)
            // Coordinates 5..6; after deletion addedIdx is -3, so effective end = 6 - (-3) = 9 which is NOT < 8
            var overlapping = MkVar(5, "TI", "TA", "Overlap_TI5_6");

            string baseSeq = "MPEPTIDESEQX";
            // Apply deletion (remove positions 2-4) => M + TIDESEQX
            string afterDeletion = "M" + baseSeq.Substring(4);
            // Apply substitution at (original) 8; due to deletion shift, adjust manually (not strictly needed for test)
            string mutated = afterDeletion.Substring(0, 6) + "K" + afterDeletion.Substring(7);

            var list = new List<SequenceVariation> { earlyDeletion, overlapping, applied };

            var result = InvokeAdjust2(applied, mutated, list);

            // overlapping should have passed through adjustment path (coordinates changed)
            var outOverlap = result.Single(v => v.Description == "Overlap_TI5_6");
            // Expect begin shifted: seqLenChange for applied (K vs E is 0), overlap with applied variant? They don't overlap (5..6 vs applied 8)
            // But addedIdx = (-3) from deletion; condition failed so it enters adjust block:
            // overlap = 0 (no direct intersection with applied variant range 8..8)
            // begin = 5 + 0 - 0? + seqLenChange(applied)=0 - overlap=0 -> 5
            // Because addedIdx only influences early-continue decision; coordinates remain same here.
            // Validate it was NOT added via early path by verifying object reference (it is a new instance, not original)
            Assert.That(ReferenceEquals(outOverlap, overlapping), Is.False);
            Assert.That(outOverlap.OneBasedBeginPosition, Is.EqualTo(5));
            Assert.That(outOverlap.OneBasedEndPosition, Is.EqualTo(6));
        }

        [Test]
        public void AdjustSequenceVariationIndices_SameVcfRecord_TakesPrecedenceOverBeforeLogic()
        {
            // Applied variant at 7 with VCF
            string vcf = "1\t2000\trs2\tA\tG\t.\tPASS\tANN=.\tGT:AD\t0/1:15,5";
            var applied = MkVarVcf(7, "D", "N", "Applied_D7N", vcf);

            // Another variant ending after applied (would not satisfy before condition) but same VCF ensures early add
            var followerSameVcf = MkVarVcf(9, "S", "T", "FollowerSameVCF_S9T", vcf);

            string mutated = "MPEPTINSEQX"; // approximate after substitution

            var list = new List<SequenceVariation> { applied, followerSameVcf };

            var result = InvokeAdjust2(applied, mutated, list);

            var outVar = result.Single(v => v.Description == "FollowerSameVCF_S9T");
            Assert.That(ReferenceEquals(outVar, followerSameVcf), Is.True, "Must be added via sameVcfRecord early path without cloning.");
        }

        #endregion
        #region Branch tests: overlap / shifting / begin-skip / end-clamp logic

        [Test]
        public void AdjustSequenceVariationIndices_NoOverlap_PositiveSeqLenChange_ShiftsForward()
        {
            // Applied insertion at position 6 (I -> ILM) delta +2
            var applied = MkVar(6, "I", "ILM", "Applied_I6ILM");
            string mutated = "MPEPTILMDESEQX"; // length 14 (base 12 +2)

            // Variant v entirely after applied (no overlap) original coords 10..11
            var after = MkVar(10, "ES", "QT", "After_ES10_11QT"); // span 10..11

            var list = new List<SequenceVariation> { after, applied };
            var result = InvokeAdjust2(applied, mutated, list);

            var adj = result.Single(v => v.Description == "After_ES10_11QT");
            // Shifted by +2 (delta) because overlap=0
            Assert.That(adj.OneBasedBeginPosition, Is.EqualTo(12));
            Assert.That(adj.OneBasedEndPosition, Is.EqualTo(13));
        }

        [Test]
        public void AdjustSequenceVariationIndices_NoOverlap_NegativeSeqLenChange_ShiftsBackward()
        {
            // Applied deletion 6..8 (IDE -> '') delta -3
            var applied = MkVar(6, "IDE", "", "Applied_Del_6_8");
            string mutated = "MPEPTSEQX"; // original 12 -> new 9

            // Variant after region (positions 9..10 originally; note original end 10 inside base)
            var after = MkVar(9, "SE", "QT", "After_SE9_10QT");

            var list = new List<SequenceVariation> { after, applied };
            var result = InvokeAdjust2(applied, mutated, list);

            var adj = result.Single(v => v.Description == "After_SE9_10QT");
            // Shift by -3 (delta), overlap=0 -> 9-3=6, 10-3=7
            Assert.That(adj.OneBasedBeginPosition, Is.EqualTo(6));
            Assert.That(adj.OneBasedEndPosition, Is.EqualTo(7));
        }

        [Test]
        public void AdjustSequenceVariationIndices_PartialOverlap_PositiveDelta()
        {
            // Applied insertion at single residue 6 (I -> IL) delta +1
            var applied = MkVar(6, "I", "IL", "Applied_I6IL");
            string mutated = "MPEPTILDESEQX"; // base 12 +1

            // Variant spanning 5..7 overlaps applied at position 6 (overlap=1)
            var span = MkVar(5, "TID", "TMD", "Span_5_7");

            var list = new List<SequenceVariation> { span, applied };
            var result = InvokeAdjust2(applied, mutated, list);

            var adj = result.Single(v => v.Description == "Span_5_7");
            // begin = 5 +1 -1 =5; end=7 +1 -1 =7 (net unchanged because overlap absorbed delta)
            Assert.That(adj.OneBasedBeginPosition, Is.EqualTo(5));
            Assert.That(adj.OneBasedEndPosition, Is.EqualTo(7));
        }

        [Test]
        public void AdjustSequenceVariationIndices_FullContainment_PositiveDelta_ShiftsBackWithinApplied()
        {
            // Applied insertion enlarging region 6 (I -> ILL) delta +2
            var applied = MkVar(6, "I", "ILL", "Applied_I6ILL");
            string mutated = "MPEPTILLDESEQX"; // len 14

            // Variant fully inside applied original 6..6 (point change same site) but distinct object
            var inside = MkVar(6, "I", "K", "Inside_I6K");

            var list = new List<SequenceVariation> { inside, applied };
            var result = InvokeAdjust2(applied, mutated, list);

            // Because same coordinates but not same reference & not sameVcfRecord ? overlap = 1; begin=6+2-1=7; end=6+2-1=7
            // This shows containment adjustment (shifts forward by delta-overlap)
            var adj = result.Single(v => v.Description == "Inside_I6K");
            Assert.That(adj.OneBasedBeginPosition, Is.EqualTo(7));
            Assert.That(adj.OneBasedEndPosition, Is.EqualTo(7));
        }

        [Test]
        public void AdjustSequenceVariationIndices_BeginBeyondLength_SkippedByStopTruncation()
        {
            // Applied variant introduces early stop: replace 6..10 (length 5) with "K*" (length 2) delta -3
            // New truncated sequence length = 6 (positions 1..6 kept)
            var applied = MkVar(6, "IDESE", "K*", "Applied_Stop_6_10");
            string mutated = "MPEPTK"; // truncated at stop

            // Variant after original region 11..11 (E->Q) original coordinate now beyond truncated length
            var after = MkVar(11, "Q", "R", "After_Q11R");

            var list = new List<SequenceVariation> { after, applied };
            var result = InvokeAdjust2(applied, mutated, list);

            // 'after' should be skipped (not present) because begin > truncated length
            Assert.That(result.Any(v => v.Description == "After_Q11R"), Is.False);
        }
        [Test]
        public void AdjustSequenceVariationIndices_EndClamped_ByStopTruncation_NoException()
        {
            // Applied stop variant: 5..7 (len 3) -> "K*" (len 2) delta -1; resulting truncated sequence length = 5
            var applied = MkVar(5, "TID", "K*", "Applied_Stop_5_7");
            string mutated = "MPEPK"; // truncated sequence

            // Long span variant starting at 5 extending beyond new sequence end
            // Original span 5..10 (len 6)
            // Overlap with applied = 3 (5..7)
            // seqLenChange (applied) = -1
            // begin = 5 -1 -3 = 1
            // end = 10 -1 -3 = 6 > truncated len (5) ? clamped to 5
            // Constructor does NOT throw; it produces a SequenceVariation whose (end - begin + 1) < original sequence length.
            var longSpan = MkVar(5, "TIDESE", "XTIDESE", "LongSpan_5_10");

            var list = new List<SequenceVariation> { longSpan, applied };

            var result = InvokeAdjust2(applied, mutated, list);

            // Verify adjusted variant exists (no exception was thrown)
            var adj = result.Single(v => v.Description == "LongSpan_5_10");
            Assert.That(adj.OneBasedBeginPosition, Is.EqualTo(1));
            Assert.That(adj.OneBasedEndPosition, Is.EqualTo(5));
            Assert.That(adj.OriginalSequence, Is.EqualTo("TIDESE"));
            Assert.That(adj.VariantSequence, Is.EqualTo("XTIDESE"));

            // Document current behavior: coordinate span (5) shorter than original sequence length (6)
            Assert.That(adj.OneBasedEndPosition - adj.OneBasedBeginPosition + 1, Is.LessThan(adj.OriginalSequence.Length),
                "Current implementation allows truncation producing a shorter coordinate span than OriginalSequence length.");
        }

        #endregion
    }
}