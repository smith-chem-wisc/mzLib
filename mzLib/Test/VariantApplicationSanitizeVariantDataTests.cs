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
         * Phase 1 already covered: null enumerable ? yield break (no messages).
         * Phase 2 here: cover early per-item guards inside the foreach:
         *
         *   if (prot == null) continue;
         *   if (prot.SequenceVariations == null) continue;
         *
         * Also cover "real" Protein with an empty (but non-null) SequenceVariations list
         * which should produce no notes, exercising originalCount == 0 with no mutations.
         *
         * Because real Proteomics.Protein never exposes a null SequenceVariations collection,
         * we introduce a minimal test-only dummy biopolymer that implements IHasSequenceVariants
         * (via IBioPolymer) returning null for SequenceVariations to exercise that path.
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

            // IBioPolymer / IHasSequenceVariants core
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

            // Variant-related collections (simulate null SequenceVariations path)
            public List<SequenceVariation> SequenceVariations => null; // intentionally null to trigger skip
            public List<SequenceVariation> AppliedSequenceVariations { get; }
            public List<TruncationProduct> TruncationProducts { get; }

            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }

            public IBioPolymer ConsensusVariant => this;

            // Generic variant factory (return original unchanged for sanitation tests)
            public TBioPolymerType CreateVariant<TBioPolymerType>(
                string variantBaseSequence,
                TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants,
                IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications,
                string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants
            {
                return original;
            }

            public IEnumerable<IBioPolymerWithSetMods> Digest(
                IDigestionParams digestionParams,
                List<Modification> allKnownFixedModifications,
                List<Modification> variableModifications,
                List<SilacLabel> silacLabels = null,
                (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
                bool topDownTruncationSearch = false) =>
                Enumerable.Empty<IBioPolymerWithSetMods>();

            public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence,
                IDictionary<int, List<Modification>> newMods) => this;

            public IDictionary<int, List<Modification>> SelectValidOneBaseMods(IDictionary<int, List<Modification>> dict) => dict;

            public bool Equals(IBioPolymer other) => ReferenceEquals(this, other);
            public override bool Equals(object obj) => Equals(obj as IBioPolymer);
            public override int GetHashCode() => Accession.GetHashCode(StringComparison.Ordinal);
        }

        #endregion

        #region Existing Phase-1 Test

        [Test]
        public void SanitizeVariantData_NullEnumerable_YieldsNoMessages()
        {
            var notes = VariantApplication.SanitizeVariantData<Protein>(polymers: null);

            Assert.That(notes, Is.Not.Null);
            Assert.That(notes.Any(), Is.False);
        }

        #endregion

        #region New Early-Loop Guard Tests

        [Test]
        public void SanitizeVariantData_EnumerableWithOnlyNullProtein_ProducesNoNotes()
        {
            var list = new Protein[] { null };
            var notes = VariantApplication.SanitizeVariantData(list).ToList();

            Assert.That(notes.Count, Is.EqualTo(0),
                "Null protein entries should be skipped silently.");
        }

        [Test]
        public void SanitizeVariantData_EnumerableWithNullAndEmptyRealProtein_NoNotes()
        {
            // Real Protein initializes SequenceVariations to an empty non-null list
            var real = new Protein("MPEPTIDESEQ", "REAL_EMPTY");
            // Ensure it has no sequence variations
            Assert.That(real.SequenceVariations, Is.Not.Null);
            Assert.That(real.SequenceVariations.Count, Is.EqualTo(0));

            var list = new Protein[] { null, real };
            var notes = VariantApplication.SanitizeVariantData(list).ToList();

            Assert.That(notes.Count, Is.EqualTo(0),
                "Empty SequenceVariations (non-null) should produce no sanitation notes.");
        }

        [Test]
        public void SanitizeVariantData_ProteinWithNullSequenceVariations_SkippedSilently()
        {
            var dummy = new DummyNullSeqVariantsBioPolymer("NULL_SEQVAR");
            var notes = VariantApplication.SanitizeVariantData(new[] { dummy }).ToList();

            Assert.That(notes.Count, Is.EqualTo(0),
                "Protein with null SequenceVariations should be skipped with no messages.");
        }

        [Test]
        public void SanitizeVariantData_MixedNullProtein_NullSeqVariants_RealEmpty_NoNotes()
        {
            var dummy = new DummyNullSeqVariantsBioPolymer("MIX_NULL");
            var real = new Protein("MPEPTIDESEQXX", "REAL_EMPTY2");
            var list = new IHasSequenceVariants[] { null, dummy, real };

            var notes = VariantApplication.SanitizeVariantData(list).ToList();

            Assert.That(notes.Count, Is.EqualTo(0),
                "All early-guarded entries (null, null SequenceVariations, empty list) should yield no notes.");
        }

        #endregion
        #region Variant Loop: Null Variant + Coordinate Sanity Tests

        private SequenceVariation MakeVar(int begin, string orig, string variant, string desc) =>
            new SequenceVariation(begin, begin + (orig?.Length > 0 ? orig.Length - 1 : 0), orig, variant, desc);

        [Test]
        public void SanitizeVariantData_DropsNullVariant_AddsDroppedAndSanitizedNotes()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_NULL_ONLY");
            prot.SequenceVariations.Add(null); // single null entry

            var notes = VariantApplication.SanitizeVariantData(prot).ToList();

            Assert.That(notes.Count, Is.EqualTo(2), "Expected a drop note and a sanitized summary note.");
            Assert.That(notes.Any(n => n.Contains("Dropped null variant")), Is.True);
            Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
            Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0), "Null variant should have been removed.");
        }

        [Test]
        public void SanitizeVariantData_DropsOutOfRange_WhenRemoveInvalidTrue()
        {
            var baseSeq = "MPEPTIDESEQVAR"; // length = 14
            var prot = new Protein(baseSeq, "ACC_OUTRANGE_DROP");
            // Out-of-range: begin > length + 1 (length+1 allowed for insertion; need >)
            int invalidBegin = baseSeq.Length + 2; // 16
            var outOfRange = MakeVar(invalidBegin, "A", "V", "oor_high");
            prot.SequenceVariations.Add(outOfRange);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped variant (coords out of range)") && n.Contains(outOfRange.SimpleString())),
                    Is.True, "Missing coordinate out-of-range drop message.");
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True,
                    "Expected sanitized summary after dropping only variant.");
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0),
                    "Out-of-range variant should have been removed when removeInvalidVariants=true.");
            });
        }

        [Test]
        public void SanitizeVariantData_KeepsOutOfRange_WhenRemoveInvalidFalse()
        {
            var baseSeq = "MPEPTIDESEQVAR"; // length = 14
            var prot = new Protein(baseSeq, "ACC_OUTRANGE_KEEP");
            int invalidBegin = baseSeq.Length + 2;
            var outOfRange = MakeVar(invalidBegin, "A", "V", "oor_high");
            prot.SequenceVariations.Add(outOfRange);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: false).ToList();

            // Expect ONLY the drop message (variant kept, so no sanitized variants note)
            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(1), "Only the drop message should be emitted.");
                Assert.That(notes[0].Contains("Dropped variant (coords out of range)") &&
                            notes[0].Contains(outOfRange.SimpleString()), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1),
                    "Variant should be retained when removeInvalidVariants=false.");
                Assert.That(prot.SequenceVariations[0], Is.SameAs(outOfRange));
            });
        }

        [Test]
        public void SanitizeVariantData_MixedNullAndOutOfRange_AndValid_VariousDrops()
        {
            var baseSeq = "MPEPTIDESEQVAR"; // length = 14
            var prot = new Protein(baseSeq, "ACC_MIXED");
            // Add: null, valid in-range, out-of-range
            prot.SequenceVariations.Add(null);

            var valid = MakeVar(5, "T", "A", "valid_mid"); // in-range substitution
            prot.SequenceVariations.Add(valid);

            int invalidBegin = baseSeq.Length + 3; // further out-of-range
            var outOfRange = MakeVar(invalidBegin, "E", "K", "oor_far");
            prot.SequenceVariations.Add(outOfRange);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped null variant")), Is.True);
                Assert.That(notes.Any(n => n.Contains("Dropped variant (coords out of range)") &&
                                           n.Contains(outOfRange.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 1/3")), Is.True,
                    "Expected sanitized summary with 1 kept of 3 original.");
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                Assert.That(prot.SequenceVariations[0].SimpleString(), Is.EqualTo(valid.SimpleString()),
                    "Only the valid in-range variant should remain.");
            });
        }

        #endregion
        #region Variant Loop: Validation (AreValid true/false + exception path)

        [Test]
        public void SanitizeVariantData_InvalidNoOp_Removed_WhenRemoveInvalidTrue()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_INVALID_DROP");
            // No-op (Original == Variant, no mods) => AreValid() should return false
            var noOp = MakeVar(3, "P", "P", "noop_same");
            prot.SequenceVariations.Add(noOp);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(noOp.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_InvalidNoOp_Retained_WhenRemoveInvalidFalse()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_INVALID_KEEP");
            var noOp = MakeVar(5, "T", "T", "noop_same2");
            prot.SequenceVariations.Add(noOp);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: false).ToList();

            Assert.Multiple(() =>
            {
                // Drop message still logged
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped invalid variant") && notes[0].Contains(noOp.SimpleString()), Is.True);
                // Variant retained
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                Assert.That(prot.SequenceVariations[0], Is.SameAs(noOp));
                // No sanitized summary because kept == original
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
            });
        }

        [Test]
        public void SanitizeVariantData_AreValidThrows_Removed_WhenRemoveInvalidTrue()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_THROW_DROP");
            // Craft variant likely to throw in AreValid by supplying a null modification list entry
            var badMods = new Dictionary<int, List<Modification>> { { 2, null } };
            var throwing = new SequenceVariation(2, 2, "P", "A", "throw_mod_null", variantCallFormatDataString: null, oneBasedModifications: badMods);
            prot.SequenceVariations.Add(throwing);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: true).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Any(n => n.Contains("Dropped invalid variant") && n.Contains(throwing.SimpleString())), Is.True);
                Assert.That(notes.Any(n => n.Contains("Sanitized variants: kept 0/1")), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(0));
            });
        }

        [Test]
        public void SanitizeVariantData_AreValidThrows_Retained_WhenRemoveInvalidFalse()
        {
            var prot = new Protein("MPEPTIDESEQ", "ACC_THROW_KEEP");
            var badMods = new Dictionary<int, List<Modification>> { { 3, null } };
            var throwing = new SequenceVariation(3, 3, "E", "K", "throw_mod_null_keep", variantCallFormatDataString: null, oneBasedModifications: badMods);
            prot.SequenceVariations.Add(throwing);

            var notes = VariantApplication.SanitizeVariantData(prot, removeInvalidVariants: false).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(notes.Count, Is.EqualTo(1));
                Assert.That(notes[0].Contains("Dropped invalid variant") && notes[0].Contains(throwing.SimpleString()), Is.True);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
                Assert.That(prot.SequenceVariations[0], Is.SameAs(throwing));
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
                // No changes ? no sanitized summary
                Assert.That(notes.Any(n => n.Contains("Sanitized variants:")), Is.False);
                Assert.That(prot.SequenceVariations.Count, Is.EqualTo(1));
            });
        }

        #endregion
    }
}