using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class VariantApplicationApplySingleVariantTests
    {
        private static MethodInfo _applySingleVariantGeneric;

        [OneTimeSetUp]
        public void LocateMethod()
        {
            _applySingleVariantGeneric = typeof(VariantApplication)
                .GetMethods(BindingFlags.NonPublic | BindingFlags.Static)
                .FirstOrDefault(m =>
                    m.Name == "ApplySingleVariant"
                    && m.IsGenericMethodDefinition
                    && m.GetParameters().Length == 3
                    && m.GetParameters()[0].ParameterType == typeof(SequenceVariation))
                ?? throw new InvalidOperationException("Unable to locate ApplySingleVariant<> by reflection.");
        }

        private static Protein InvokeApplySingleVariant(SequenceVariation variant, Protein protein, string individual = "")
        {
            var mi = _applySingleVariantGeneric.MakeGenericMethod(typeof(Protein));
            return (Protein)mi.Invoke(null, new object[] { variant, protein, individual })!;
        }

        private static SequenceVariation Var(int begin, string original, string variant, string desc) =>
            new SequenceVariation(begin,
                begin + (original?.Length ?? 0) - 1,
                original,
                variant,
                desc,
                variantCallFormatDataString: null,
                oneBasedModifications: null);

        private Protein MakeBaseProtein(string accession = "BASE_APPLY", string seq = "MPEPTIDESEQX")
        {
            var p = new Protein(seq, accession);
            p.TruncationProducts.AddRange(new[]
            {
                new TruncationProduct(1,3,"before"),
                new TruncationProduct(4,10,"span"),
                new TruncationProduct(8,12,"after")
            });
            return p;
        }

        private Modification DummyMod(string id = "Mod1") =>
            new Modification(_originalId: id, _accession: "ACC", _modificationType: "TestType");

        [Test]
        public void ApplySingleVariant_Insertion_AdjustsSequence_Variants_TruncationProducts()
        {
            var baseProtein = MakeBaseProtein();
            var insertion = Var(6, "I", "ILM", "Insertion_I6_to_ILM");
            var variantProtein = InvokeApplySingleVariant(insertion, baseProtein);

            Assert.That(variantProtein.BaseSequence, Is.EqualTo("MPEPTILMDESEQX"));
            var applied = variantProtein.AppliedSequenceVariations.Single(v => v.Description == "Insertion_I6_to_ILM");
            Assert.Multiple(() =>
            {
                Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(6));
                Assert.That(applied.OneBasedEndPosition, Is.EqualTo(6));
                Assert.That(applied.OriginalSequence, Is.EqualTo("I"));
                Assert.That(applied.VariantSequence, Is.EqualTo("ILM"));
            });

            var tps = variantProtein.TruncationProducts;
            Assert.That(tps.Any(tp => tp.OneBasedBeginPosition == 1 && tp.OneBasedEndPosition == 3), Is.True);
            Assert.That(tps.Any(tp => tp.OneBasedBeginPosition == 4 && tp.OneBasedEndPosition == 12), Is.True);
            Assert.That(tps.Any(tp => tp.OneBasedBeginPosition == 10 && tp.OneBasedEndPosition == 14), Is.True);

            Assert.That(baseProtein.BaseSequence, Is.EqualTo("MPEPTIDESEQX")); // unchanged
        }

        [Test]
        public void ApplySingleVariant_NullVariant_ReturnsOriginal()
        {
            var baseProtein = MakeBaseProtein();
            var result = InvokeApplySingleVariant(null, baseProtein);
            Assert.That(ReferenceEquals(result, baseProtein), Is.True);
            Assert.That(result.AppliedSequenceVariations, Is.Empty);
        }

        [Test]
        public void ApplySingleVariant_NullProtein_ReturnsNull()
        {
            var variant = Var(3, "E", "K", "Sub_E3K");
            var mi = _applySingleVariantGeneric.MakeGenericMethod(typeof(Protein));
            var result = mi.Invoke(null, new object[] { variant, null, "" });
            Assert.That(result, Is.Null);
        }

        [Test]
        public void ApplySingleVariant_InvalidBeginPastLengthPlusOne_ReturnsOriginal()
        {
            var baseProtein = MakeBaseProtein();
            int invalidBegin = baseProtein.BaseSequence.Length + 2; // length+2 triggers guard
            var variant = new SequenceVariation(invalidBegin, null, "AA", "OutOfRangeInsertion"); // insertion form
            var result = InvokeApplySingleVariant(variant, baseProtein);
            Assert.That(ReferenceEquals(result, baseProtein), Is.True);
            Assert.That(result.AppliedSequenceVariations, Is.Empty);
            Assert.That(result.BaseSequence, Is.EqualTo(baseProtein.BaseSequence));
        }

        [Test]
        public void ApplySingleVariant_InsertionAtLengthPlusOne_AppendsSequence()
        {
            var baseProtein = MakeBaseProtein();
            int appendPos = baseProtein.BaseSequence.Length + 1; // legal insertion site
            var variant = new SequenceVariation(appendPos, null, "AA", "TailInsertion");
            var result = InvokeApplySingleVariant(variant, baseProtein);

            Assert.That(result.BaseSequence, Is.EqualTo(baseProtein.BaseSequence + "AA"));
            var applied = result.AppliedSequenceVariations.Single(v => v.Description == "TailInsertion");
            Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(appendPos));
            Assert.That(applied.OneBasedEndPosition, Is.EqualTo(appendPos));
        }

        [Test]
        public void ApplySingleVariant_OverrunOriginalSequence_AdjustsReplacedLength()
        {
            var baseProtein = MakeBaseProtein(); // length 12
            // Begin at 11 with original length 5 (runs past end)
            var overrun = Var(11, "EQXZZ", "K", "OverrunNearEnd");
            // Manually craft variant that tries to replace beyond end:
            // originalSeq length 5 -> afterIdx=11+5-1=15>12 triggers adjust path
            var result = InvokeApplySingleVariant(overrun, baseProtein);

            // New sequence: first 10 chars + 'K' (since afterIdx clipped) => positions 11..12 replaced by original substring clipped to remaining (2 residues)
            // Original tail starting at 11 (index 10 zero-based) is 'QX'
            // Replacement logic: seqBefore = first 10 chars, seqAfter becomes empty (afterIdx >= length)
            // seqBefore = MPEPTIDESE (first 10)
            Assert.That(result.BaseSequence, Is.EqualTo("MPEPTIDESEK"));

            var applied = result.AppliedSequenceVariations.Single(v => v.Description == "OverrunNearEnd");
            Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(11));
            Assert.That(applied.OneBasedEndPosition, Is.EqualTo(11 + ("EQXZZ".Length - 1))); // end based on original span (even if overrun)
        }

        [Test]
        public void ApplySingleVariant_Deletion_RemovesSegment()
        {
            var baseProtein = MakeBaseProtein(); // MPEPTIDESEQX
            // Delete 'TID' at positions 5..7 -> variantSeq empty
            var deletion = Var(5, "TID", "", "Del_5_7");
            var result = InvokeApplySingleVariant(deletion, baseProtein);

            // Expected sequence: positions 1-4 + positions 8-12 => MPEP ESEQX
            Assert.That(result.BaseSequence, Is.EqualTo("MPEPESEQX"));

            var applied = result.AppliedSequenceVariations.Single(v => v.Description == "Del_5_7");
            Assert.Multiple(() =>
            {
                Assert.That(applied.OriginalSequence, Is.EqualTo("TID"));
                Assert.That(applied.VariantSequence, Is.Empty);
                Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(5));
                Assert.That(applied.OneBasedEndPosition, Is.EqualTo(7));
            });
        }

        [Test]
        public void ApplySingleVariant_VariantSpecificModifications_Copied()
        {
            var baseProtein = MakeBaseProtein();
            // Substitution with variant-specific modification at position 6 (global)
            var mods = new Dictionary<int, List<Modification>>
            {
                { 6, new List<Modification>{ DummyMod("VarMod1") } }
            };
            var variantWithMods = new SequenceVariation(6, 6, "I", "K", "Sub_I6K_WithMod", variantCallFormatDataString: null, oneBasedModifications: mods);

            var result = InvokeApplySingleVariant(variantWithMods, baseProtein);
            var applied = result.AppliedSequenceVariations.Single(v => v.Description == "Sub_I6K_WithMod");

            Assert.That(applied.OneBasedModifications, Is.Not.Null);
            Assert.That(applied.OneBasedModifications.ContainsKey(6), Is.True);
            Assert.That(applied.OneBasedModifications[6].Count, Is.EqualTo(1));
            Assert.That(result.BaseSequence[5], Is.EqualTo('K'));
        }

        [Test]
        public void ApplySingleVariant_PointSubstitution_NoLengthChange()
        {
            var baseProtein = MakeBaseProtein();
            var sub = Var(3, "E", "K", "Sub_E3K");
            var result = InvokeApplySingleVariant(sub, baseProtein);

            Assert.That(result.BaseSequence, Is.EqualTo("MPKPTIDESEQX"));
            var applied = result.AppliedSequenceVariations.Single();
            Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(3));
            Assert.That(applied.OneBasedEndPosition, Is.EqualTo(3));
        }

        [Test]
        public void ApplySingleVariant_InsertionCreatesStop_TruncatesSequence()
        {
            var baseProtein = MakeBaseProtein(); // length 12
            // Insert "AA*" at position 6 replacing "I" (stop terminates after concatenation and split('*')[0])
            var stopIns = Var(6, "I", "AA*", "InsertionWithStop");
            var result = InvokeApplySingleVariant(stopIns, baseProtein);

            // New sequence should truncate before '*' : prefix (positions 1..5) + "AA" => MPEPTAA
            Assert.That(result.BaseSequence, Is.EqualTo("MPEPTAA"));

            var applied = result.AppliedSequenceVariations.Single(v => v.Description == "InsertionWithStop");
            Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(6));
        }
        #region Branch tests: intersectsAppliedRegionIncompletely path coverage

        [Test]
        public void ApplySingleVariant_IncompleteIntersection_DropsPreviousAndUsesConsensusSeqAfter()
        {
            // Base protein
            var baseProt = MakeBaseProtein();

            // First variant (variant1) spanning 5..9: "TIDES" -> "QQQQQ" (same length substitution)
            var variant1 = Var(5, "TIDES", "QQQQQ", "Span_5_9_Qs");
            var protAfterV1 = InvokeApplySingleVariant(variant1, baseProt);

            Assert.That(protAfterV1.BaseSequence, Is.EqualTo("MPEPQQQQQEQX"), "Precondition altered sequence unexpected.");

            // Second variant (variant2) fully INSIDE variant1 span but NOT including variant1:
            // Replace positions 7..8 (currently 'QQ') with 'KK'.
            var variant2 = Var(7, "QQ", "KK", "Inner_7_8_KK");

            // Because variant2 is strictly inside variant1 (variant2 does NOT include variant1;
            // variant2 span 7..8, variant1 span 5..9) AND they intersect, the condition:
            // Intersects && !Includes == true ? intersectsAppliedRegionIncompletely = true
            var protAfterV2 = InvokeApplySingleVariant(variant2, protAfterV1);

            // Expected sequence logic:
            // seqBefore (1..6) from protAfterV1: MPEPQQ
            // replaced segment (7..8) => KK
            // seqAfter (override uses consensus ORIGINAL base, not protAfterV1) from position 9 onward of consensus: S E Q X
            // Final: M P E P Q Q K K S E Q X
            Assert.That(protAfterV2.BaseSequence, Is.EqualTo("MPEPQQKKSEQX"), "Sequence did not reflect consensus-based seqAfter override.");

            // Applied variations: ONLY the second variant (previous one not merged due to incomplete intersection)
            var appliedDescs = protAfterV2.AppliedSequenceVariations.Select(v => v.Description).ToList();
            Assert.That(appliedDescs, Is.EquivalentTo(new[] { "Inner_7_8_KK" }),
                "Previous intersecting variant should not be merged when intersection is incomplete.");
        }

        [Test]
        public void ApplySingleVariant_PreviousVariantFullyIncluded_RemovedFromMergedAppliedList()
        {
            var baseProt = MakeBaseProtein();

            // Small prior variant (variant1) inside the region of the next larger variant
            var variant1 = Var(6, "I", "L", "Point_I6L");
            var protAfterV1 = InvokeApplySingleVariant(variant1, baseProt);
            Assert.That(protAfterV1.AppliedSequenceVariations.Count, Is.EqualTo(1));

            // Larger variant (variant2) completely includes variant1 span: 5..9
            var variant2 = Var(5, "TIDES", "AAAAA", "Block_5_9_AAAAA");
            var protAfterV2 = InvokeApplySingleVariant(variant2, protAfterV1);

            // Since variant2 includes variant1, intersectsAppliedRegionIncompletely == false
            // Merge path excludes included variant (filter !Includes)
            Assert.That(protAfterV2.AppliedSequenceVariations.Count, Is.EqualTo(1), "Included prior variant should have been excluded.");
            Assert.That(protAfterV2.AppliedSequenceVariations.Single().Description, Is.EqualTo("Block_5_9_AAAAA"));

            // Sequence positions 5..9 replaced with AAAAA
            Assert.That(protAfterV2.BaseSequence, Is.EqualTo("MPEPAAAAAEQX"));
        }

        [Test]
        public void ApplySingleVariant_NoIncompleteIntersection_MergesNonOverlappingPriorVariants()
        {
            var baseProt = MakeBaseProtein();

            // Manually seed two non-overlapping prior applied variations (simulate earlier applications)
            var prior1 = Var(2, "P", "A", "Prior_P2A");   // span 2..2
            var prior2 = Var(11, "Q", "R", "Prior_Q11R"); // span 11..11
            baseProt.AppliedSequenceVariations.Add(prior1);
            baseProt.AppliedSequenceVariations.Add(prior2);

            // New variant does not intersect either (substitution at position 6)
            var newVar = Var(6, "I", "K", "Central_I6K");
            var result = InvokeApplySingleVariant(newVar, baseProt);

            // Merge path, keep those not included
            var descs = result.AppliedSequenceVariations.Select(v => v.Description).OrderBy(s => s).ToList();
            Assert.That(descs, Is.EqualTo(new[] { "Central_I6K", "Prior_P2A", "Prior_Q11R" }.OrderBy(s => s)));

            // Base sequence updated only at position 6
            Assert.That(result.BaseSequence, Is.EqualTo("MPEPTKDESEQX".Replace("P2A", "")), "Sequence mismatch (only position 6 substitution expected).");
            Assert.That(result.BaseSequence[5], Is.EqualTo('K'));
        }

        [Test]
        public void ApplySingleVariant_IncompleteIntersection_PriorVariantExtendsRight()
        {
            var baseProt = MakeBaseProtein();
            // Prior variant extends beyond the new variant to the right: prior 6..10, new 6..7
            var prior = Var(6, "IDESE", "AAAAA", "Prior_6_10_AAAAA");
            var protAfterPrior = InvokeApplySingleVariant(prior, baseProt);

            var newVar = Var(6, "ID", "KK", "New_6_7_KK"); // inside left portion of prior, not including its full span
            var protAfterNew = InvokeApplySingleVariant(newVar, protAfterPrior);

            // Incomplete overlap ? previous not merged
            Assert.That(protAfterNew.AppliedSequenceVariations.Count, Is.EqualTo(1));
            Assert.That(protAfterNew.AppliedSequenceVariations.Single().Description, Is.EqualTo("New_6_7_KK"));

            // Sequence rebuild uses consensus tail (original base) after position 7
            // Consensus (original) after position 7 => positions 8..12: E S E Q X
            // Sequence prefix (positions 1..5) from prior variant base: M P E P T (prior replaced 6..10 so position 5 remains T)
            // Insert new variant at 6..7 'KK'
            // Final: M P E P T K K E S E Q X
            Assert.That(protAfterNew.BaseSequence, Is.EqualTo("MPEPTKKESEQX"));
        }
        [Test]
        public void ApplySingleVariant_IncompleteIntersection_PriorVariantExtendsLeft()
        {
            var baseProt = MakeBaseProtein();
            // Prior variant spanning 4..8 replaces 'PTIDE' (positions 4..8) with AAAAA
            var prior = Var(4, "PTIDE", "AAAAA", "Prior_4_8_AAAAA");
            var protAfterPrior = InvokeApplySingleVariant(prior, baseProt);

            // New variant fully inside prior span (5..6) and does not include full prior region -> incomplete intersection
            var newVar = Var(5, "TI", "KK", "New_5_6_KK");
            var protAfterNew = InvokeApplySingleVariant(newVar, protAfterPrior);

            // Only the new inner variant should remain (prior is discarded due to incomplete intersection rule)
            Assert.That(protAfterNew.AppliedSequenceVariations.Select(v => v.Description),
                Is.EquivalentTo(new[] { "New_5_6_KK" }));

            // Explanation:
            // After prior: M P E A A A A A S E Q X  (position 4 changed to 'A')
            // New variant (5..6) ? seqBefore = first 4 residues of modified sequence = M P E A
            // Variant seq = KK
            // seqAfter sourced from consensus (original) starting at afterIdx (6) ? original positions 7..12 = D E S E Q X
            // Final: M P E A K K D E S E Q X
            Assert.That(protAfterNew.BaseSequence, Is.EqualTo("MPEAKKDESEQX"));
        }
        #endregion

    }
}