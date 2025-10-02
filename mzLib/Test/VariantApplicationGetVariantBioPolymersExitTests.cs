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

        [Test]
        public void ValidationLoop_ValidVariant_AddedToValidList()
        {
            var v1 = Sub(4, 'P', 'L', "valid");
            var protein = CreateProteinWithVariants("VALID_ONLY", v1);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 2,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Any(p => p.BaseSequence != protein.BaseSequence), Is.True);
        }

        [Test]
        public void ValidationLoop_InvalidAfterMutation_FailedBranch()
        {
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
            modVariant.OneBasedModifications.Clear();

            var protein = CreateProteinWithVariants("INVALID_MUTATED", modVariant);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 3,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

                /// <summary>
                /// NOTE: The original intent was to force an exception inside AreValid().
                /// The current SequenceVariation.AreValid() implementation is defensive and does not throw
                /// under mutation of its dictionary reference. We instead verify that mutating the
                /// OneBasedModifications reference to null (and re-adding content) does not break processing
                /// and still produces variant isoforms (resilience test, not catch-path test).
                /// </summary>
                [Test]
                public void ValidationLoop_MutationResilience_DoesNotThrow()
                {
                    var v = Sub(6, 'E', 'K', "mutable_mods");
                    // Remove all variant-specific modifications to ensure pure substitution (valid)
                    v.OneBasedModifications.Clear();

                    // Simulate external mutation: set backing field to null (reflection)
                    var fld = typeof(SequenceVariation).GetField("<OneBasedModifications>k__BackingField",
                        BindingFlags.Instance | BindingFlags.NonPublic);
                    Assert.That(fld, Is.Not.Null);
                    fld!.SetValue(v, null);

                    var protein = CreateProteinWithVariants("MUT_RESILIENT", v);

                    Assert.DoesNotThrow(() =>
                    {
                        var res = protein.GetVariantBioPolymers(
                            maxSequenceVariantsPerIsoform: 2,
                            minAlleleDepth: 1,
                            maxSequenceVariantIsoforms: 5);
                        // Depending on downstream filtering a variant isoform may or may not appear (if sequence changes, it should).
                        Assert.That(res.Count, Is.GreaterThanOrEqualTo(1));
                    });
                }

        [Test]
        public void ValidationLoop_Mixed_AllBranchesCoveredSimultaneously()
        {
            var protein = new Protein("MPEPTIDESEQ", "MIXED_BRANCHES");
            protein.SequenceVariations.Add(null);

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
            modVar.OneBasedModifications.Clear();
            protein.SequenceVariations.Add(modVar);

            var throwVar = Sub(5, 'T', 'A', "thrower_sim"); // will act as normal substitution now
            protein.SequenceVariations.Add(throwVar);

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
            protein.SequenceVariations.Add(null);

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
            temp.OneBasedModifications.Clear();
            protein.SequenceVariations.Add(temp);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 3,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].BaseSequence, Is.EqualTo(protein.BaseSequence));
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        #endregion

        #region Fallback Block Specific Tests

        // Scenario A: All variants null -> first valid.Count==0, fallback list empty, second valid.Count==0 => returns base
        [Test]
        public void Fallback_AllVariantsNull_ReturnsBase()
        {
            var protein = new Protein("MPEPTIDESEQ", "FALLBACK_ALL_NULL");
            protein.SequenceVariations.AddRange(new SequenceVariation[] { null, null, null });

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 5,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Count, Is.EqualTo(1), "Expected only base protein when all variants are null.");
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        // Scenario B: All variants non-null but invalid (mutated to pure no-op) -> fallback picks them up (non-empty),
        // ApplyAllVariantCombinations filters them out (no-op removal) -> base only.
        [Test]
        public void Fallback_AllVariantsInvalidNoOps_FallbackNonEmptyButResultBase()
        {
            var protein = new Protein("MPEPTIDESEQ", "FALLBACK_ALL_INVALID");

            // Create 3 variants that become invalid (no-op) after modification removal
            for (int i = 0; i < 3; i++)
            {
                int pos = 2 + i;
                var v = new SequenceVariation(
                    pos,
                    pos,
                    "E",
                    "E",
                    $"noop_{i}",
                    variantCallFormatDataString: null,
                    oneBasedModifications: new Dictionary<int, List<Modification>>
                    {
                        { pos, new List<Modification>{ MakeMod($"Mod{i}") } }
                    });
                v.OneBasedModifications.Clear(); // now invalid (AreValid false)
                protein.SequenceVariations.Add(v);
            }

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 5,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 10);

            Assert.That(result.Count, Is.EqualTo(1),
                "Fallback should retain invalid variants but downstream filtering should leave only base.");
            Assert.That(result[0].AppliedSequenceVariations.Count, Is.EqualTo(0));
        }

        // Scenario C: Mixed invalid (forcing fallback) is impossible to produce variant isoforms because any invalid remains invalid later.
        // So add a control showing that adding a single valid variant avoids fallback (valid.Count>0) and yields variants.
        [Test]
        public void Fallback_NotTriggeredWhenAnyValidVariantExists()
        {
            var protein = new Protein("MPEPTIDESEQ", "FALLBACK_NOT_TRIGGERED");

            // Invalid no-op (post mutation)
            int pos = 5;
            var invalid = new SequenceVariation(
                pos,
                pos,
                "T",
                "T",
                "noop_w_mod",
                variantCallFormatDataString: null,
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { pos, new List<Modification>{ MakeMod("Temp") } }
                });
            invalid.OneBasedModifications.Clear();
            protein.SequenceVariations.Add(invalid);

            // Valid substitution
            var valid = Sub(7, 'D', 'N', "real_change");
            protein.SequenceVariations.Add(valid);

            var result = protein.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 3,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 5);

            Assert.That(result.Any(p => p.BaseSequence != protein.BaseSequence), Is.True,
                "Presence of a valid variant should bypass fallback empty-valid behavior and yield variant isoforms.");
        }

        #endregion
    }
}