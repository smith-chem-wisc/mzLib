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
    }
}