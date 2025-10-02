using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;

namespace Test.DatabaseTests
{
    [TestFixture]
    public class VariantApplicationGetVariantBioPolymersExitTests
    {
        /// <summary>
        /// Dummy IHasSequenceVariants; can optionally simulate a null SequenceVariations list.
        /// </summary>
        private sealed class NullVariantsProtein : IHasSequenceVariants
        {
            private readonly Protein _consensus;
            private readonly bool _returnNullSequenceVariations;
            private readonly List<SequenceVariation>? _seqVars;

            public NullVariantsProtein(string sequence,
                                       string accession,
                                       bool returnNullSequenceVariations = true)
            {
                BaseSequence = sequence;
                _consensus = new Protein(sequence, accession + "_CONS");
                AppliedSequenceVariations = new List<SequenceVariation>();
                OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
                OriginalNonVariantModifications = new Dictionary<int, List<Modification>>();
                TruncationProducts = new List<TruncationProduct>();
                _returnNullSequenceVariations = returnNullSequenceVariations;
                if (!returnNullSequenceVariations)
                {
                    _seqVars = new List<SequenceVariation>();
                }
            }

            public string BaseSequence { get; }
            public string SampleNameForVariants => string.Empty;
            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
            public IBioPolymer ConsensusVariant => _consensus;
            public List<SequenceVariation> AppliedSequenceVariations { get; }
            public List<TruncationProduct> TruncationProducts { get; }

#pragma warning disable CS8603
            public List<SequenceVariation> SequenceVariations =>
                _returnNullSequenceVariations ? null : _seqVars;
#pragma warning restore CS8603

            public TBioPolymerType CreateVariant<TBioPolymerType>(
                string variantBaseSequence,
                TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants,
                IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications,
                string sampleNameForVariants) where TBioPolymerType : IHasSequenceVariants
            {
                return original;
            }
        }

        private Protein CreateProteinWithVariants(string accession, params SequenceVariation[] vars)
        {
            var p = new Protein("MPEPTIDESEQ", accession);
            if (vars != null && vars.Length > 0)
            {
                p.SequenceVariations.AddRange(vars);
            }
            return p;
        }

        private SequenceVariation Sub(int pos, char from, char to, string desc = null)
            => new SequenceVariation(pos, from.ToString(), to.ToString(), desc ?? $"{from}{pos}{to}");

        private Modification MakeMod(string id) =>
            new Modification(_originalId: id, _accession: id, _modificationType: "unit-test", _featureType: "ft", _target: null);

        #region Guard: (maxSequenceVariantsPerIsoform == 0 || maxSequenceVariantIsoforms == 1)

        [TestCase(0, 0)]
        [TestCase(0, 1)]
        [TestCase(0, 2)]
        [TestCase(0, 10)]
        [TestCase(1, 1)]
        [TestCase(4, 1)]
        public void GetVariantBioPolymers_Exit_CombinatoricsDisabled(int maxVariantsPerIsoform, int maxIsoforms)
        {
            var v1 = Sub(3, 'E', 'K');
            var v2 = Sub(7, 'D', 'N');
            var protein = CreateProteinWithVariants($"P_{maxVariantsPerIsoform}_{maxIsoforms}", v1, v2);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: maxVariantsPerIsoform,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: maxIsoforms);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(ReferenceEquals(result[0], protein), Is.True);
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        #endregion

        #region Guard: (all.Count == 0) with non-guard combinatorics settings

        [TestCase(1, 0)]
        [TestCase(1, 2)]
        [TestCase(1, 10)]
        [TestCase(4, 0)]
        [TestCase(4, 2)]
        [TestCase(4, 10)]
        public void GetVariantBioPolymers_NoVariants_ListEmpty(int maxVariantsPerIsoform, int maxIsoforms)
        {
            var protein = CreateProteinWithVariants($"EMPTY_{maxVariantsPerIsoform}_{maxIsoforms}");

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: maxVariantsPerIsoform,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: maxIsoforms);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(ReferenceEquals(result[0], protein), Is.True);
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        [TestCase(1, 0)]
        [TestCase(4, 0)]
        public void GetVariantBioPolymers_NoVariants_IsoformsZero(int maxVariantsPerIsoform, int maxIsoforms)
        {
            var protein = CreateProteinWithVariants($"EMPTY_ISO0_{maxVariantsPerIsoform}", Array.Empty<SequenceVariation>());

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: maxVariantsPerIsoform,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: maxIsoforms);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(ReferenceEquals(result[0], protein), Is.True);
        }

        [TestCase(1, 2)]
        [TestCase(4, 10)]
        public void GetVariantBioPolymers_NullSequenceVariations(int maxVariantsPerIsoform, int maxIsoforms)
        {
            var nullProt = new NullVariantsProtein("MPEPTIDESEQ",
                                                   $"NULLSEQ_{maxVariantsPerIsoform}_{maxIsoforms}",
                                                   returnNullSequenceVariations: true);

            var result = nullProt.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: maxVariantsPerIsoform,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: maxIsoforms);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(ReferenceEquals(result[0], nullProt), Is.True);
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        #endregion

        #region Non-guard path sanity

        [Test]
        public void GetVariantBioPolymers_VariantsApplied()
        {
            var v1 = Sub(3, 'E', 'K');
            var v2 = Sub(7, 'D', 'N');
            var protein = CreateProteinWithVariants("APPLY_OK", v1, v2);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 10);

            Assert.That(result.Count, Is.GreaterThanOrEqualTo(3));
            Assert.That(result.Any(p => p.AppliedSequenceVariations.Count > 0), Is.True);
            Assert.That(result.First().AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        [Test]
        public void GetVariantBioPolymers_IsoformLimitRestricts()
        {
            var v1 = Sub(3, 'E', 'K');
            var v2 = Sub(7, 'D', 'N');
            var protein = CreateProteinWithVariants("LIMITED", v1, v2);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 2);

            Assert.That(result.Count, Is.EqualTo(2));
            Assert.That(result.Count(p => p.AppliedSequenceVariations.Count > 0), Is.EqualTo(1));
        }

        #endregion

        #region Validation Loop Branch Tests

        /// <summary>
        /// Hits: v == null (failed++ & continue) AND final early base return (valid.Count == 0 after fallback).
        /// </summary>
        [Test]
        public void ValidationLoop_NullOnlyVariant_ListContainsNull_ReturnsBase()
        {
            var protein = new Protein("MPEPTIDESEQ", "NULL_ONLY_CASE");
            protein.SequenceVariations.Add(null);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 2,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(ReferenceEquals(result[0], protein), Is.True);
        }

        /// <summary>
        /// Hits: AreValid() true path ? variant added to 'valid'.
        /// </summary>
        [Test]
        public void ValidationLoop_ValidVariant_AddedToValidList()
        {
            var v1 = Sub(4, 'P', 'L', "valid");
            var protein = CreateProteinWithVariants("VALID_ONLY", v1);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 2,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Any(p => p.BaseSequence != protein.BaseSequence), Is.True,
                "Expected at least one variant protein derived from a valid variant.");
        }

        /// <summary>
        /// Hits: AreValid() returns false (after mutation) ? failed++ branch (no exception).
        /// </summary>
        [Test]
        public void ValidationLoop_InvalidAfterMutation_FailedBranch()
        {
            // Create a variant that is only valid because it has a modification while Original==Variant.
            int pos = 5;
            var modVariant = new SequenceVariation(
                pos,
                pos,
                "T",
                "T",
                "noop_with_mod",
                variantCallFormatDataString: null,
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { pos, new List<Modification>{ MakeMod("TempMod") } }
                });
            // Now mutate to invalid (no-op without modifications).
            modVariant.OneBasedModifications.Clear();

            var protein = CreateProteinWithVariants("INVALID_MUTATED", modVariant);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 3,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            // Only base expected (variant filtered later as pure no-op).
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        /// <summary>
        /// Hits: catch path (AreValid throws) ? ok forced true, threw++, variant retained.
        /// Reflection corrupts backing field for OneBasedModifications to force InvalidCastException inside AreValid.
        /// </summary>
        [Test]
        public void ValidationLoop_AreValidThrows_CatchBranchTreatsAsValid()
        {
            var throwingVar = Sub(6, 'E', 'K', "throw_test");

            // Corrupt backing field to force runtime cast failure on property access inside AreValid.
            var fld = typeof(SequenceVariation).GetField("<OneBasedModifications>k__BackingField",
                BindingFlags.Instance | BindingFlags.NonPublic);
            Assert.That(fld, Is.Not.Null, "Could not locate backing field for OneBasedModifications via reflection.");
            fld!.SetValue(throwingVar, new object()); // incompatible type

            var protein = CreateProteinWithVariants("THROWING", throwingVar);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 2,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 10);

            // Expect at least one variant produced (the substitution) despite exception.
            Assert.That(result.Any(p => p.BaseSequence != protein.BaseSequence), Is.True,
                "Exception during AreValid() should not prevent variant inclusion (catch branch).");
        }
        // ... (unchanged using directives and class header)

        [Test]
        public void ValidationLoop_Mixed_AllBranchesCoveredSimultaneously()
        {
            var protein = new Protein("MPEPTIDESEQ", "MIXED_BRANCHES");

            // 1. Null
            protein.SequenceVariations.Add(null);

            // 2. Mutated invalid (initially valid via mod)
            int pos = 3;
            var modVar = new SequenceVariation(
                pos,
                pos,
                "E",
                "E",
                "noop_mod",
                variantCallFormatDataString: null,
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { pos, new List<Modification>{ MakeMod("Keep") } }
                });
            modVar.OneBasedModifications.Clear(); // now invalid (AreValid false)
            protein.SequenceVariations.Add(modVar);

            // 3. Throwing variant
            var throwVar = Sub(5, 'T', 'A', "thrower");
            var fld = typeof(SequenceVariation).GetField("<OneBasedModifications>k__BackingField",
                BindingFlags.Instance | BindingFlags.NonPublic);
            fld!.SetValue(throwVar, new object());
            protein.SequenceVariations.Add(throwVar);

            // 4. Normal valid variant
            var goodVar = Sub(8, 'D', 'N', "good");
            protein.SequenceVariations.Add(goodVar);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 10);

            var variantSeqs = result.Where(p => p.BaseSequence != protein.BaseSequence).ToList();
            Assert.That(variantSeqs.Count, Is.GreaterThanOrEqualTo(1));

            bool containsMutatedOnly = variantSeqs.Any(p =>
                p.AppliedSequenceVariations.Count == 1 &&
                p.AppliedSequenceVariations[0].SimpleString().Contains("E3E"));
            Assert.That(containsMutatedOnly, Is.False);
        }

        [Test]
        public void ValidationLoop_FallbackAfterEmptyValidList_NoUsableVariants()
        {
            var protein = new Protein("MPEPTIDESEQ", "FALLBACK_CASE");

            // Null variant
            protein.SequenceVariations.Add(null);

            // Mutated invalid variant (initially valid with mod)
            int pos = 4;
            var temp = new SequenceVariation(
                pos,
                pos,
                "P",
                "P",
                "noop_temp",
                variantCallFormatDataString: null,
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { pos, new List<Modification>{ MakeMod("TempMod") } }
                });
            temp.OneBasedModifications.Clear(); // now invalid
            protein.SequenceVariations.Add(temp);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 3,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].BaseSequence, Is.EqualTo(protein.BaseSequence));
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        // ... (rest of file unchanged)
        #endregion
    }
}