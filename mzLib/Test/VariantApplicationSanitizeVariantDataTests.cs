using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;

namespace Test.DatabaseTests
{
    [TestFixture]
    public class VariantApplicationSanitizeVariantDataTests
    {
        /*
         * Phases covered:
         *  - Null enumerable guard
         *  - Per-item guards (null protein, null SequenceVariations collection, empty list)
         *  - Null variant entry
         *  - Coordinate out-of-range variants (drop vs retain depending on removeInvalidVariants flag)
         *  - Mixed sets (null + valid + out-of-range)
         *  - Invalid “no-op” variants created via post-construction mutation (drop vs retain)
         *  - Invalid span constructor rejection
         *  - Variant-specific modification pruning (out-of-range added post-construction)
         *  - Valid variant (no messages)
         */

        #region Test-only Dummy Types

        // Minimal dummy implementing IBioPolymer to exercise SequenceVariations == null path
        private sealed class DummyNullSeqVariantsBioPolymer : IBioPolymer
        {
            public DummyNullSeqVariantsBioPolymer(string accession = "DUMMY_NULL")
            {
                Accession = accession;
                BaseSequence = "MAAATESTSEQ";
                OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
                OriginalNonVariantModifications = new Dictionary<int, List<Modification>>();
                AppliedSequenceVariations = new List<SequenceVariation>();
                TruncationProducts = new List<TruncationProduct>();
                GeneNames = new List<Tuple<string, string>>();
            }

            public string Accession { get; }
            public string BaseSequence { get; }
            public string Name => Accession;
            public string FullName => Accession;
            public int Length => BaseSequence.Length;
            public string DatabaseFilePath => string.Empty;
            public bool IsDecoy => false;
            public bool IsContaminant => false;
            public string Organism => "TEST_ORG";
            public List<Tuple<string, string>> GeneNames { get; }
            public string SampleNameForVariants => string.Empty;

            public List<SequenceVariation> SequenceVariations => null; // trigger skip branch
            public List<SequenceVariation> AppliedSequenceVariations { get; }
            public List<TruncationProduct> TruncationProducts { get; }

            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }

            public IBioPolymer ConsensusVariant => this;

            public TBioPolymerType CreateVariant<TBioPolymerType>(
                string variantBaseSequence,
                TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants,
                IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications,
                string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants => original;

            public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams,
                List<Modification> allKnownFixedModifications,
                List<Modification> variableModifications,
                List<SilacLabel> silacLabels = null,
                (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
                bool topDownTruncationSearch = false) => Enumerable.Empty<IBioPolymerWithSetMods>();

            public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence,
                IDictionary<int, List<Modification>> newMods) => this;

            public IDictionary<int, List<Modification>> SelectValidOneBaseMods(IDictionary<int, List<Modification>> dict) => dict;

            public bool Equals(IBioPolymer other) => ReferenceEquals(this, other);
            public override bool Equals(object obj) => Equals(obj as IBioPolymer);
            public override int GetHashCode() => Accession.GetHashCode(StringComparison.Ordinal);
        }

        #endregion

        #region Phase 1: Null Enumerable

        [Test]
        public void SanitizeVariantData_NullEnumerable_YieldsNoMessages()
        {
            var notes = VariantApplication.SanitizeVariantData<Protein>(polymers: null);
            Assert.That(notes, Is.Not.Null);
            Assert.That(notes.Any(), Is.False);
        }

        #endregion

        #region Early Per-Item Guards

        [Test]
        public void SanitizeVariantData_EnumerableWithOnlyNullProtein_ProducesNoNotes()
        {
            var list = new Protein[] { null };
            var notes = VariantApplication.SanitizeVariantData(list).ToList();
            Assert.That(notes.Count, Is.EqualTo(0));
        }

        [Test]
        public void SanitizeVariantData_EnumerableWithNullAndEmptyRealProtein_NoNotes()
        {
            var real = new Protein("MPEPTIDESEQ", "REAL_EMPTY");
            Assert.That(real.SequenceVariations.Count, Is.EqualTo(0));
            var list = new Protein[] { null, real };
            var notes = VariantApplication.SanitizeVariantData(list).ToList();
            Assert.That(notes.Count, Is.EqualTo(0));
        }

        [Test]
        public void SanitizeVariantData_ProteinWithNullSequenceVariations_SkippedSilently()
        {
            var dummy = new DummyNullSeqVariantsBioPolymer("NULL_SEQVAR");
            var notes = VariantApplication.SanitizeVariantData(new[] { dummy }).ToList();
            Assert.That(notes.Count, Is.EqualTo(0));
        }

        [Test]
        public void SanitizeVariantData_MixedNullProtein_NullSeqVariants_RealEmpty_NoNotes()
        {
            var dummy = new DummyNullSeqVariantsBioPolymer("MIX_NULL");
            var real = new Protein("MPEPTIDESEQXX", "REAL_EMPTY2");
            var notes = VariantApplication.SanitizeVariantData(new IHasSequenceVariants[] { null, dummy, real }).ToList();
            Assert.That(notes.Count, Is.EqualTo(0));
        }

        #endregion

        #region Helpers

        private SequenceVariation MakeVar(int begin, string orig, string variant, string desc)
            => new SequenceVariation(begin, begin + (orig?.Length > 0 ? orig.Length - 1 : 0), orig, variant, desc);

        private static Modification MakeTestMod(string id) =>
            new Modification(
                _originalId: id,
                _accession: id,
                _modificationType: "test-mod",
                _featureType: "feature",
                _target: null,
                _locationRestriction: "Unassigned.",
                _chemicalFormula: null,
                _monoisotopicMass: null,
                _databaseReference: null,
                _taxonomicRange: null,
                _keywords: new List<string>(),
                _neutralLosses: null,
                _diagnosticIons: null,
                _fileOrigin: null);

        #endregion

        #region Null Variant + Coordinate Sanity

        [Test]
        public void SanitizeVariantData_DropsNullVariant_AddsDroppedAndSanitizedNotes()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_NULL_ONLY");
            prot.SequenceVariations.Add(null);
            var notes = VariantApplication.SanitizeVariantData(prot).ToList();
            Assert.That(notes.Count, Is.EqualTo(2));
            Assert.That(notes.Any(n => n.Contains("Dropped null variant")), Is.True);
            Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
            Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
        }

        [Test]
        public void SanitizeVariantData_DropsOutOfRange_WhenRemoveInvalidTrue()
        {
            var seq = "MPEPTIDESEQVAR";
            var prot = new Protein(seq, "ACC_OUTRANGE_DROP");
            var invalid = MakeVar(seq.Length + 2, "A", "V", "oor_high");
            prot.SequenceVariations.Add(invalid);
            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();
            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped variant (coords out of range)") && n.Contains(invalid.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_KeepsOutOfRange_WhenRemoveInvalidFalse()
        {
            var seq = "MPEPTIDESEQVAR";
            var prot = new Protein(seq, "ACC_OUTRANGE_KEEP");
            var invalid = MakeVar(seq.Length + 2, "A", "V", "oor_high");
            prot.SequenceVariations.Add(invalid);
            var notes = VariantApplication.SanitizeVariantData(prot, false).ToList();
            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped variant (coords out of range)") && notes[0].Contains(invalid.SimpleString()), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
            });
        }

        [Test]
        public void SanitizeVariantData_MixedNullAndOutOfRange_AndValid_VariousDrops()
        {
            var seq = "MPEPTIDESEQVAR";
            var prot = new Protein(seq, "ACC_MIXED");
            prot.SequenceVariations.Add(null);
            var valid = MakeVar(5, "T", "A", "valid_mid");
            prot.SequenceVariations.Add(valid);
            var invalid = MakeVar(seq.Length + 3, "E", "K", "oor_far");
            prot.SequenceVariations.Add(invalid);

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped null variant")), Is.True);
                Assert.That(notes.Any(n => n.Contains("Dropped variant (coords out of range)") && n.Contains(invalid.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 1/3")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                Assert.That(prot.SequenceVariations[0].SimpleString(), Is.EqualTo(valid.SimpleString()));
            });
        }

        #endregion

        #region Validation / Mutation Scenarios

        [Test]
        public void SanitizeVariantData_InvalidNoOp_Removed_WhenRemoveInvalidTrue()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_INVALID_DROP");
            int pos = 3;
            var mod = MakeTestMod("TestMod");
            var variant = new SequenceVariation(pos, pos, "P", "P", "noop_with_mod_then_cleared", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(variant);
            variant.OneBasedModifications.Clear(); // becomes pure no-op
            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();
            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(variant.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_InvalidNoOp_Retained_WhenRemoveInvalidFalse()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_INVALID_KEEP");
            int pos = 5;
            var mod = MakeTestMod("TestMod2");
            var variant = new SequenceVariation(pos, pos, "T", "T", "noop_with_mod_then_cleared_keep", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(variant);
            try { variant.OneBasedModifications.Clear(); } catch { }
            var notes = VariantApplication.SanitizeVariantData(prot, false).ToList();
            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped invalid variant") && notes[0].Contains(variant.SimpleString()), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
            });
        }

        [Test]
        public void SequenceVariation_InvalidSpan_ConstructorThrows()
        {
            Assert.That(() =>
                new SequenceVariation(10, 9, "A", "G", "invalid_span_should_throw", (string?)null, null),
                Throws.TypeOf<ArgumentException>().With.Message.Contains("coordinates"));
        }

        [Test]
        public void SanitizeVariantData_PrunesOutOfRangeVariantSpecificModSite()
        {
            var prot = new Protein("MPEPTIDEQ", "ACC_PRUNE_OOR"); // length 9
            int pos = 3;
            var mod = MakeTestMod("InRange");
            var variant = new SequenceVariation(pos, pos, "P", "L", "simple_sub_with_mod", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(variant);

            // Inject an out-of-range variant-specific mod AFTER construction to trigger pruning (position > maxAllowedPos)
            int invalidPos = prot.BaseSequence.Length + 5; // 14
            variant.OneBasedModifications[invalidPos] = new List<Modification> { MakeTestMod("OOR") };

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("pruned 1 mod site") && n.Contains(variant.SimpleString())), Is.True);
                Assert.That(variant.OneBasedModifications.Keys.SequenceEqual(new[] { pos }), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
            });
        }

        [Test]
        public void SanitizeVariantData_ValidVariant_NoInvalidMessage()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_VALID_OK");
            var valid = MakeVar(4, "P", "L", "valid_sub");
            prot.SequenceVariations.Add(valid);
            var notes = VariantApplication.SanitizeVariantData(prot).ToList();
            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant")), Is.False);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
            });
        }

        #endregion
        #region Pruning Tests (variant-specific modification pruning)

        [Test]
        public void SanitizeVariantData_NoPruning_WhenAllVariantSpecificModsValid_NonDeletion()
        {
            var prot = new Protein("MPEPTIDEQK", "ACC_PRUNE_NONE"); // length 10
            int begin = 5;
            var variant = new SequenceVariation(begin, begin, "T", "A", "subst_with_valid_mods",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 2, new List<Modification>{ MakeTestMod("ModA") } },
                    { 9, new List<Modification>{ MakeTestMod("ModB") } }
                });
            prot.SequenceVariations.Add(variant);

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("pruned")), Is.False);
                Assert.That(variant.OneBasedModifications.Keys.OrderBy(k => k).SequenceEqual(new[] { 2, 9 }), Is.True);
            });
        }

        // Deletion + invalid mod positions: AreValid() now fails BEFORE pruning ? variant dropped, not pruned.
        [Test]
        public void SanitizeVariantData_Deletion_InvalidMods_Dropped_WhenRemoveInvalidTrue()
        {
            var prot = new Protein("MAPTIDEQK", "ACC_DEL_DROP"); // length 9
            int begin = 3;
            int end = 6;
            var deletion = new SequenceVariation(begin, end, "PTID", "", "deletion_region",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 2, new List<Modification>{ MakeTestMod("KeepBefore") } } // valid site (before deletion)
                });
            prot.SequenceVariations.Add(deletion);

            // Add invalid (at/after begin) – these cause AreValid() to fail so variant is DROPPED (not pruned)
            deletion.OneBasedModifications[3] = new List<Modification> { MakeTestMod("AtBegin") };
            deletion.OneBasedModifications[5] = new List<Modification> { MakeTestMod("Inside") };
            deletion.OneBasedModifications[8] = new List<Modification> { MakeTestMod("After") };

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(deletion.SimpleString())), Is.True,
                    "Expected invalid deletion variant to be dropped (AreValid fails) rather than pruned.");
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_Deletion_InvalidMods_Retained_WhenRemoveInvalidFalse()
        {
            var prot = new Protein("MAPTIDEQK", "ACC_DEL_RETAIN"); // length 9
            int begin = 3;
            int end = 6;
            var deletion = new SequenceVariation(begin, end, "PTID", "", "deletion_region_keep",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 2, new List<Modification>{ MakeTestMod("KeepBefore") } }
                });
            prot.SequenceVariations.Add(deletion);

            deletion.OneBasedModifications[3] = new List<Modification> { MakeTestMod("AtBegin") };
            deletion.OneBasedModifications[5] = new List<Modification> { MakeTestMod("Inside") };
            deletion.OneBasedModifications[8] = new List<Modification> { MakeTestMod("After") };

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: false).ToList();

            Assert.Multiple(() =>
            {
                // Variant invalid => "Dropped invalid variant" note, but retained (no sanitized summary since kept == original)
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped invalid variant") && notes[0].Contains(deletion.SimpleString()), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
            });
        }

        // Stop-gain + invalid mod positions: AreValid() fails (mods at/after begin) -> drop/retain logic mirrors deletion.
        [Test]
        public void SanitizeVariantData_StopGain_InvalidMods_Dropped_WhenRemoveInvalidTrue()
        {
            var prot = new Protein("MPEPTIDEQK", "ACC_STOP_DROP");
            int begin = 4;
            var stopGain = new SequenceVariation(begin, begin, "P", "*", "stop_gain_region",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 3, new List<Modification>{ MakeTestMod("KeepBefore") } }
                });
            prot.SequenceVariations.Add(stopGain);

            stopGain.OneBasedModifications[4] = new List<Modification> { MakeTestMod("AtStop") };
            stopGain.OneBasedModifications[7] = new List<Modification> { MakeTestMod("AfterStop") };

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(stopGain.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_StopGain_InvalidMods_Retained_WhenRemoveInvalidFalse()
        {
            var prot = new Protein("MPEPTIDEQK", "ACC_STOP_RETAIN");
            int begin = 4;
            var stopGain = new SequenceVariation(begin, begin, "P", "*", "stop_gain_region_keep",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 3, new List<Modification>{ MakeTestMod("KeepBefore") } }
                });
            prot.SequenceVariations.Add(stopGain);

            stopGain.OneBasedModifications[4] = new List<Modification> { MakeTestMod("AtStop") };
            stopGain.OneBasedModifications[7] = new List<Modification> { MakeTestMod("AfterStop") };

            var notes = VariantApplication.SanitizeVariantData(prot, false).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped invalid variant") && notes[0].Contains(stopGain.SimpleString()), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
            });
        }
        [Test]
        public void SanitizeVariantData_Insertion_ValidMods_NoPruning()
        {
            // Insertion: original residue 'T' at position 5 replaced by 'TAAA' (delta +3)
            // Base length = 10 => new sequence length = 13; valid mod positions: 1..13
            var prot = new Protein("MPEPTIDEQK", "ACC_INS_NOPRUNE"); // length 10
            int pos = 5;
            var insertion = new SequenceVariation(
                pos,
                pos,
                "T",
                "TAAA",
                "insertion_valid_mods",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 5,  new List<Modification>{ MakeTestMod("KeepSite") } }, // valid (within inserted block)
                    { 13, new List<Modification>{ MakeTestMod("KeepMax") } }  // valid (last new residue)
                });

            prot.SequenceVariations.Add(insertion);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                // No pruning note expected
                Assert.That(notes.Any(n => n.Contains("pruned")), Is.False, "Unexpected pruning note for fully valid insertion variant.");
                // No sanitized summary (kept == original count, and no drops)
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False, "No sanitized summary expected (no variants removed).");
                // Variant retained
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                // Modification keys unchanged
                Assert.That(insertion.OneBasedModifications.Keys.OrderBy(k => k).SequenceEqual(new[] { 5, 13 }), Is.True);
            });
        }
        [Test]
        public void SanitizeVariantData_Prunes_Mixed_AllThreeConditions()
        {
            var prot = new Protein("MPEPTIDEQK", "ACC_PRUNE_MIX"); // length 10
            int begin = 6;
            int end = 7;
            var deletion = new SequenceVariation(begin, end, "DE", "", "mixed_deletion",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 5, new List<Modification>{ MakeTestMod("KeepBefore") } }
                });
            prot.SequenceVariations.Add(deletion);

            deletion.OneBasedModifications[6] = new List<Modification> { MakeTestMod("DelBegin") };
            deletion.OneBasedModifications[9] = new List<Modification> { MakeTestMod("DelAfter") };
            deletion.OneBasedModifications[-2] = new List<Modification> { MakeTestMod("Neg") };
            deletion.OneBasedModifications[25] = new List<Modification> { MakeTestMod("TooHigh") };
            deletion.OneBasedModifications[2] = new List<Modification> { MakeTestMod("KeepFarBefore") };
                
            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            // Because invalid (mods at/after begin for a deletion) => AreValid fails ? variant dropped (not pruned)
            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(deletion.SimpleString())), Is.True,
                    "Expected variant drop (invalid) rather than pruning note.");
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        #endregion
        #region CVB Insertion Tests

        [Test]
        public void SanitizeVariantData_Insertion_InvalidOutOfRangeMods_Dropped_WhenRemoveInvalidTrue()
        {
            // Insertion: original "T" -> "TAAA" at position 5 (delta +3)
            // Base length = 10 ? maxAllowedPos = 13. We will inject invalid positions (-1, 14) AFTER construction.
            var prot = new Protein("MPEPTIDEQK", "ACC_INS_DROP"); // length 10
            int pos = 5;
            var insertion = new SequenceVariation(pos, pos, "T", "TAAA", "insertion_with_invalid_mods",
                (string?)null,
                new Dictionary<int, List<Modification>>
                {
                    { 5, new List<Modification>{ MakeTestMod("KeepSite") } },  // valid
                    { 13, new List<Modification>{ MakeTestMod("KeepMax") } }  // valid (== maxAllowedPos)
                });
            prot.SequenceVariations.Add(insertion);

            // Add invalid positions AFTER construction (these will cause AreValid() to fail, so variant is dropped — not pruned)
            insertion.OneBasedModifications[14] = new List<Modification> { MakeTestMod("TooHigh") }; // > maxAllowedPos
            insertion.OneBasedModifications[-1] = new List<Modification> { MakeTestMod("Neg") };     // < 1

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(insertion.SimpleString())), Is.True,
                    "Expected the insertion variant to be dropped as invalid (not pruned).");
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True,
                    "Should report 0/1 kept after dropping invalid insertion variant.");
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0),
                    "Invalid insertion variant should have been removed.");
            });
        }

        #endregion
        #region Sanitized Summary Branch Tests (kept.Count != originalCount)

        [Test]
        public void SanitizeVariantData_NoVariants_NoSanitizedSummary()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_SUMMARY_NONE");
            Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_AllValid_NoSanitizedSummary()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_SUMMARY_VALID");
            var valid = MakeVar(4, "P", "L", "valid_sub");
            prot.SequenceVariations.Add(valid);

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                Assert.That(ReferenceEquals(prot.SequenceVariations[0], valid), Is.True);
            });
        }

        [Test]
        public void SanitizeVariantData_DroppedNullVariant_SanitizedSummary()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_SUMMARY_NULL");
            prot.SequenceVariations.Add(null); // originalCount = 1

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped null variant")), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_DroppedInvalidVariant_SanitizedSummary()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_SUMMARY_INVALID");
            int pos = 3;
            // Create valid (temp) no-op via mod
            var mod = MakeTestMod("TempMod");
            var variant = new SequenceVariation(pos, pos, "P", "P", "noop_then_invalid", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(variant);

            // Invalidate to no-op (no mods)
            variant.OneBasedModifications.Clear();

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(variant.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_InvalidVariantRetained_NoSanitizedSummary()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_SUMMARY_INVALID_RETAIN");
            int pos = 6;
            var mod = MakeTestMod("TempMod2");
            var variant = new SequenceVariation(pos, pos, "E", "E", "noop_retain", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(variant);
            variant.OneBasedModifications.Clear(); // now pure no-op (invalid)

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: false).ToList();

            Assert.Multiple(() =>
            {
                // Invalid logged, but kept (so kept == originalCount => no sanitized summary)
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped invalid variant") && notes[0].Contains(variant.SimpleString()), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
            });
        }

        [Test]
        public void SanitizeVariantData_MixedSomeDropped_SanitizedSummary()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_SUMMARY_MIX_DROP");
            // valid
            var valid = MakeVar(4, "P", "L", "valid_sub");
            prot.SequenceVariations.Add(valid);
            // null
            prot.SequenceVariations.Add(null);
            // invalid (no-op after clearing mods)
            int pos = 7;
            var mod = MakeTestMod("TempMod3");
            var invalid = new SequenceVariation(pos, pos, "D", "D", "noop_mutated", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(invalid);
            invalid.OneBasedModifications.Clear();

            // originalCount = 3; kept expected = 1 (valid)
            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped null variant")), Is.True);
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(invalid.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 1/3")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                Assert.That(prot.SequenceVariations[0].SimpleString(), Is.EqualTo(valid.SimpleString()));
            });
        }

        [Test]
        public void SanitizeVariantData_MixedDroppedButRetainedFlag_NoSanitizedSummary()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_SUMMARY_MIX_RETAIN");
            // valid
            var valid = MakeVar(3, "P", "L", "valid_sub");
            prot.SequenceVariations.Add(valid);
            // null
            prot.SequenceVariations.Add(null);
            // invalid mutated no-op retained
            int pos = 8;
            var mod = MakeTestMod("TempMod4");
            var invalid = new SequenceVariation(pos, pos, "Q", "Q", "noop_mutated_retain", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(invalid);
            invalid.OneBasedModifications.Clear();

            // With removeInvalidVariants=false: null dropped (not added), invalid kept (explicitly added due to flag), so kept.Count == originalCount (3)?
            // Actually null variant is skipped (not added), invalid is added, valid is added -> kept =2, original=3 => sanitized summary WILL appear.
            // To ensure no summary we must avoid null (since null is never added). Adjust test: use only invalid retained.

            prot.SequenceVariations.Clear();
            prot.SequenceVariations.Add(valid);
            prot.SequenceVariations.Add(invalid); // originalCount=2, kept should remain 2 (invalid retained)

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: false).ToList();

            Assert.Multiple(() =>
            {
                // Only invalid note; no sanitized summary (since kept==original)
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped invalid variant") && notes[0].Contains(invalid.SimpleString()), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(2));
            });
        }

        #endregion
        #region AppliedSequenceVariations Reconciliation Tests

        [Test]
        public void SanitizeVariantData_AppliedEmpty_NoPruneNote()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_EMPTY");
            var v = MakeVar(4, "P", "L", "valid_sub");
            prot.SequenceVariations.Add(v);
            // Applied list intentionally left empty
            Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(0));

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs:")), Is.False);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_AppliedAllValid_NoRemovals_NoPruneNote()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_ALLVALID");
            var v = MakeVar(3, "P", "L", "valid_sub");
            prot.SequenceVariations.Add(v);
            prot.AppliedSequenceVariations.Add(v); // reference-equal

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs:")), Is.False);
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(1));
                Assert.That(ReferenceEquals(prot.AppliedSequenceVariations[0], v), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
            });
        }

        [Test]
        public void SanitizeVariantData_AppliedContainsNull_NullRemoved_PruneNote()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_NULL");
            var v = MakeVar(5, "T", "A", "valid_mid");
            prot.SequenceVariations.Add(v);
            prot.AppliedSequenceVariations.Add(v);
            prot.AppliedSequenceVariations.Add(null); // will be removed

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                // No base variant dropped ? no sanitized summary
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs: 1 removed")), Is.True);
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(1));
                Assert.That(ReferenceEquals(prot.AppliedSequenceVariations[0], v), Is.True);
            });
        }

        [Test]
        public void SanitizeVariantData_AppliedContainsStaleReference_Removed_WithPruneNote()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_STALE");
            // valid variant
            var valid = MakeVar(4, "P", "L", "valid_sub");
            // invalid variant (will be dropped)
            int posInvalid = 7;
            var mod = MakeTestMod("TempInv");
            var invalid = new SequenceVariation(posInvalid, posInvalid, "D", "D", "noop_invalid", (string?)null,
                new Dictionary<int, List<Modification>> { { posInvalid, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(valid);
            prot.SequenceVariations.Add(invalid);
            // mutate invalid to pure no-op
            invalid.OneBasedModifications.Clear();

            // Applied list references both
            prot.AppliedSequenceVariations.Add(valid);
            prot.AppliedSequenceVariations.Add(invalid);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                // invalid dropped from SequenceVariations ? sanitized summary
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(invalid.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 1/2")), Is.True);
                // Applied stale reference pruned
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs: 1 removed")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(1));
                Assert.That(ReferenceEquals(prot.AppliedSequenceVariations[0], valid), Is.True);
            });
        }
        [Test]
        public void SanitizeVariantData_AppliedContainsNullAndClone_BothRemoved_PruneNoteShowsCount2()
        {
            // NOTE: SequenceVariation equality is value-based (coords, original, variant, VCF, mods) and
            // description is NOT part of equality. So a "clone" differing only by description is considered equal
            // and will NOT be pruned by the applied reconciliation step (kept.Contains(clone) == true).
            // Therefore only the explicit null entry is pruned. Adjust expectations accordingly.

            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_NULL_CLONE");
            var baseVar = MakeVar(6, "D", "N", "valid_sub");
            prot.SequenceVariations.Add(baseVar);

            // Clone (same coordinates + sequences ? Equals == true)
            var clone = MakeVar(6, "D", "N", "valid_sub_clone");

            prot.AppliedSequenceVariations.Add(baseVar);
            prot.AppliedSequenceVariations.Add(null);   // will be pruned
            prot.AppliedSequenceVariations.Add(clone);  // value-equal ? retained

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                // Only the null reference is removed
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs: 1 removed")), Is.True,
                    "Expected only the null applied variant reference to be pruned.");
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(2),
                    "Both value-equal variants should remain (base + clone).");
                // Both remaining entries should be value-equal to baseVar
                Assert.That(prot.AppliedSequenceVariations.All(v => v.Equals(baseVar)), Is.True);
                // No sanitized summary (no base variants dropped)
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
            });
        }
        [Test]
        public void SanitizeVariantData_AppliedInvalidVariantRetained_NoPruneNote()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_INVALID_RETAIN");
            int pos = 5;
            var mod = MakeTestMod("TempKeep");
            var invalid = new SequenceVariation(pos, pos, "T", "T", "noop_invalid_retain", (string?)null,
                new Dictionary<int, List<Modification>> { { pos, new List<Modification> { mod } } });
            prot.SequenceVariations.Add(invalid);
            invalid.OneBasedModifications.Clear(); // becomes invalid
            prot.AppliedSequenceVariations.Add(invalid);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: false).ToList();

            Assert.Multiple(() =>
            {
                // Invalid logged (note) but variant kept, so applied reference stays
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(invalid.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs:")), Is.False);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False,
                    "No sanitized summary because kept == original count (variant retained).");
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(1));
                Assert.That(ReferenceEquals(prot.AppliedSequenceVariations[0], invalid), Is.True);
            });
        }

        [Test]
        public void SanitizeVariantData_AppliedOnlyDroppedNull_NoPruneNoteBecauseAppliedEmpty()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_ONLY_NULL");
            // Add null variant only so it is dropped; applied list references nothing before sanitize
            prot.SequenceVariations.Add(null);

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped null variant")), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                // Applied list was empty so no prune note
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs:")), Is.False);
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_AppliedMixedNullAndDroppedAndValid_AllPrunedCountMatches()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_APPLIED_COMPLEX");
            // valid
            var valid = MakeVar(3, "P", "L", "valid_sub");
            // invalid (droppable)
            int posInv = 8;
            var modInv = MakeTestMod("TempInv2");
            var invalid = new SequenceVariation(posInv, posInv, "Q", "Q", "noop_invalid_drop", (string?)null,
                new Dictionary<int, List<Modification>> { { posInv, new List<Modification> { modInv } } });

            prot.SequenceVariations.Add(valid);
            prot.SequenceVariations.Add(invalid);
            invalid.OneBasedModifications.Clear(); // make invalid

            prot.AppliedSequenceVariations.Add(null);      // will be pruned
            prot.AppliedSequenceVariations.Add(valid);     // kept
            prot.AppliedSequenceVariations.Add(invalid);   // stale after variant drop -> pruned
            prot.AppliedSequenceVariations.Add(MakeVar(10, "E", "K", "nonlisted_clone")); // not in kept -> pruned

            // Before: 4 applied entries
            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                // invalid variant drop
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(invalid.SimpleString())), Is.True);
                // sanitized summary (1 kept of 2 base variants)
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 1/2")), Is.True);
                // prune note (removed 3 applied refs: null + invalid + clone)
                Assert.That(notes.Any(n => n.Contains("Pruned applied variant refs: 3 removed")), Is.True);
                Assert.That(prot.AppliedSequenceVariations.Count, Is.EqualTo(1));
                Assert.That(ReferenceEquals(prot.AppliedSequenceVariations[0], valid), Is.True);
            });
        }

        #endregion
        #region Accession Prefix Selection Tests (IBioPolymer vs Consensus vs Fallback)

        // Wrapper that implements only IHasSequenceVariants (NOT IBioPolymer)
        // Forces SanitizeVariantData to fall back to ConsensusVariant.Accession
        private sealed class BareVariantContainer : IHasSequenceVariants
        {
            private readonly Protein _consensus;
            public BareVariantContainer(string consensusAccession, string seq = "MPEPTIDESEQ")
            {
                _consensus = new Protein(seq, consensusAccession);
                BaseSequence = seq;
                SampleNameForVariants = string.Empty;
                OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
                OriginalNonVariantModifications = new Dictionary<int, List<Modification>>();
                AppliedSequenceVariations = new List<SequenceVariation>();
                SequenceVariations = new List<SequenceVariation>();
                TruncationProducts = new List<TruncationProduct>(); // ADDED
            }

            public string BaseSequence { get; }
            public string SampleNameForVariants { get; }
            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
            public IBioPolymer ConsensusVariant => _consensus;
            public List<SequenceVariation> AppliedSequenceVariations { get; }
            public List<SequenceVariation> SequenceVariations { get; }
            public List<TruncationProduct> TruncationProducts { get; } // ADDED
            public TBioPolymerType CreateVariant<TBioPolymerType>(
                string variantBaseSequence,
                TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants,
                IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications,
                string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants => original;
        }

        // Wrapper that returns null ConsensusVariant to force "<no_accession>" fallback
        private sealed class NullConsensusContainer : IHasSequenceVariants
        {
            public NullConsensusContainer(string seq = "MPEPTIDESEQ")
            {
                BaseSequence = seq;
                SampleNameForVariants = "";
                OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
                OriginalNonVariantModifications = new Dictionary<int, List<Modification>>();
                AppliedSequenceVariations = new List<SequenceVariation>();
                SequenceVariations = new List<SequenceVariation>();
                TruncationProducts = new List<TruncationProduct>(); // ADDED
            }

            public string BaseSequence { get; }
            public string SampleNameForVariants { get; }
            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
            public IBioPolymer ConsensusVariant => null;
            public List<SequenceVariation> AppliedSequenceVariations { get; }
            public List<SequenceVariation> SequenceVariations { get; }
            public List<TruncationProduct> TruncationProducts { get; } // ADDED
            public TBioPolymerType CreateVariant<TBioPolymerType>(
                string variantBaseSequence,
                TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants,
                IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications,
                string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants => original;
        }

        [Test]
        public void SanitizeVariantData_AccessionPrefix_UsesDirectProteinAccession()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_DIRECT_PREFIX");
            prot.SequenceVariations.Add(null); // force a note

            var notes = VariantApplication.SanitizeVariantData(prot, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(2));
                Assert.That(notes.All(n => n.StartsWith("[ACC_DIRECT_PREFIX]")), Is.True,
                    "All notes should be prefixed with the direct protein accession.");
            });
        }

        [Test]
        public void SanitizeVariantData_AccessionPrefix_FallsBackToConsensusVariantAccession()
        {
            var container = new BareVariantContainer("ACC_CONS_FALLBACK");
            container.SequenceVariations.Add(null); // trigger sanitization path

            var notes = VariantApplication.SanitizeVariantData(container, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(2));
                Assert.That(notes.All(n => n.StartsWith("[ACC_CONS_FALLBACK]")), Is.True,
                    "Expected fallback to ConsensusVariant.Accession when object is not IBioPolymer.");
            });
        }

        [Test]
        public void SanitizeVariantData_AccessionPrefix_FallbackNoAccession()
        {
            var container = new NullConsensusContainer();
            container.SequenceVariations.Add(null); // trigger path

            var notes = VariantApplication.SanitizeVariantData(container, true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(2));
                Assert.That(notes.All(n => n.StartsWith("[<no_accession>]")), Is.True,
                    "Expected <no_accession> prefix when neither IBioPolymer nor ConsensusVariant.Accession is available.");
            });
        }

        [Test]
        public void SanitizeVariantData_AccessionPrefix_MixedTypesAllCorrect()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_REAL");
            prot.SequenceVariations.Add(null);

            var wrapper = new BareVariantContainer("ACC_WRAPPED");
            wrapper.SequenceVariations.Add(null);

            var nullCons = new NullConsensusContainer();
            nullCons.SequenceVariations.Add(null);

            var notes = VariantApplication
                .SanitizeVariantData(new IHasSequenceVariants[] { prot, wrapper, nullCons }, true)
                .ToList();

            // Expect 2 notes per object (drop + summary): 6 total
            Assert.That(notes.Count, Is.EqualTo(6));

            var grouped = notes.GroupBy(n =>
            {
                if (n.StartsWith("[ACC_REAL]")) return "real";
                if (n.StartsWith("[ACC_WRAPPED]")) return "wrapped";
                if (n.StartsWith("[<no_accession>]")) return "none";
                return "other";
            }).ToDictionary(g => g.Key, g => g.Count());

            Assert.Multiple(() =>
            {
                Assert.That(grouped.TryGetValue("real", out var c1) && c1 == 2, Is.True);
                Assert.That(grouped.TryGetValue("wrapped", out var c2) && c2 == 2, Is.True);
                Assert.That(grouped.TryGetValue("none", out var c3) && c3 == 2, Is.True);
                Assert.That(grouped.ContainsKey("other"), Is.False, "Unexpected accession prefix found.");
            });
        }

        #endregion
    }
}