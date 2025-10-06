using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Test.DatabaseTests
{
    [TestFixture]
    internal class SequenceVariationRandomTests
    {
        // ---------------- Existing Tests ----------------

        [Test]
        public void Constructor_InvalidCoordinates_ThrowsArgumentException()
        {
            // Minimal valid VCF line (10 columns) so VariantCallFormat parses without truncation.
            string vcf =
                "1\t100\t.\tA\tT\t.\tPASS\t" +
                "ANN=T|missense_variant|MODERATE|GENE1|GENE1|transcript|TX1|protein_coding|1/1|c.100A>T|p.K34N|100/1000|34/300|34/100|0|\t" +
                "GT:AD:DP\t0/1:5,4:9";

            var parsedVcf = new VariantCallFormat(vcf);

            // Intentionally invalid: end < begin (5,4) triggers AreValid() == false
            Assert.That(
                () => new SequenceVariation(
                    oneBasedBeginPosition: 5,
                    oneBasedEndPosition: 4,
                    originalSequence: "A",
                    variantSequence: "V",
                    description: "invalid-coords",
                    vcf: parsedVcf),
                Throws.TypeOf<ArgumentException>()
                      .With.Message.EqualTo("SequenceVariation coordinates are invalid."));
        }

        [Test]
        public void Equals_ReturnsFalse_ForNonSequenceVariationObjects()
        {
            // Valid point substitution so Equals reaches the type check cleanly
            var sv = new SequenceVariation(
                oneBasedBeginPosition: 5,
                oneBasedEndPosition: 5,
                originalSequence: "A",
                variantSequence: "V",
                description: "point");

            Assert.Multiple(() =>
            {
                // Different runtime type
                Assert.That(sv.Equals("not a variation"), Is.False);
                // Null
                Assert.That(sv.Equals(null), Is.False);
                // Different type but structurally similar data holder
                var anonymous = new { OneBasedBeginPosition = 5, OneBasedEndPosition = 5, OriginalSequence = "A", VariantSequence = "V" };
                Assert.That(sv.Equals(anonymous), Is.False);
            });
        }

        // ---------------- New Tests For ModificationDictionariesEqual ----------------

        private MethodInfo _modDictEqualMethod;
        private ModificationMotif _motifA;
        private ModificationMotif _motifC;

        [OneTimeSetUp]
        public void OneTimeSetUp()
        {
            _modDictEqualMethod = typeof(SequenceVariation)
                .GetMethod("ModificationDictionariesEqual", BindingFlags.NonPublic | BindingFlags.Static)
                ?? throw new InvalidOperationException("Could not reflect ModificationDictionariesEqual.");

            Assert.That(ModificationMotif.TryGetMotif("A", out _motifA), Is.True);
            Assert.That(ModificationMotif.TryGetMotif("C", out _motifC), Is.True);
        }

        private bool InvokeCompare(Dictionary<int, List<Modification>> a, Dictionary<int, List<Modification>> b)
            => (bool)_modDictEqualMethod.Invoke(null, new object[] { a, b });

        private static Modification MakeMod(string id, ModificationMotif motif) =>
            new Modification(_originalId: id, _modificationType: "TestType", _target: motif, _locationRestriction: "Anywhere.");

        [Test]
        public void ModDictEqual_ReturnsFalse_WhenOneDictionaryIsNull()
        {
            var b = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA) } }
            };

            Assert.That(InvokeCompare(null, b), Is.False);
            Assert.That(InvokeCompare(b, null), Is.False);
        }

        [Test]
        public void ModDictEqual_ReturnsFalse_WhenCountDiffers()
        {
            var a = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA) } }
            };
            var b = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA) }},
                {2, new List<Modification>{ MakeMod("M2", _motifA) }}
            };

            Assert.That(InvokeCompare(a, b), Is.False);
        }

        [Test]
        public void ModDictEqual_ReturnsFalse_WhenKeySetsDiffer()
        {
            var a = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA) }},
                {2, new List<Modification>{ MakeMod("M2", _motifA) }}
            };
            var b = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA) }},
                {3, new List<Modification>{ MakeMod("M3", _motifA) }}
            };

            Assert.That(InvokeCompare(a, b), Is.False);
        }

        [Test]
        public void ModDictEqual_ReturnsFalse_WhenOneListIsNull()
        {
            var a = new Dictionary<int, List<Modification>>
            {
                {1, null}
            };
            var b = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA) }}
            };

            Assert.That(InvokeCompare(a, b), Is.False);
        }

        [Test]
        public void ModDictEqual_ReturnsFalse_WhenListCountsDiffer()
        {
            var a = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA), MakeMod("M2", _motifA) }}
            };
            var b = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA) }}
            };

            Assert.That(InvokeCompare(a, b), Is.False);
        }

        [Test]
        public void ModDictEqual_ReturnsFalse_WhenDistinctKeyCountsDiffer()
        {
            var a = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA), MakeMod("M2", _motifA) }}
            };
            var b = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA), MakeMod("M1", _motifA) }}
            };

            Assert.That(InvokeCompare(a, b), Is.False);
        }

        [Test]
        public void ModDictEqual_ReturnsFalse_WhenFrequencyMismatchForSameDistinctCount()
        {
            var a = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("AX", _motifA), MakeMod("AY", _motifA) }}
            };
            var b = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("BX", _motifA), MakeMod("BY", _motifA) }}
            };

            Assert.That(InvokeCompare(a, b), Is.False);
        }

        [Test]
        public void ModDictEqual_Control_ReturnsTrue_ForEquivalentDictionaries()
        {
            var a = new Dictionary<int, List<Modification>>
            {
                {1, new List<Modification>{ MakeMod("M1", _motifA), MakeMod("M2", _motifA) }},
                {3, new List<Modification>{ MakeMod("M3", _motifC) }}
            };
            var b = new Dictionary<int, List<Modification>>
            {
                {3, new List<Modification>{ MakeMod("M3", _motifC) }},
                {1, new List<Modification>{ MakeMod("M2", _motifA), MakeMod("M1", _motifA) }}
            };

            Assert.That(InvokeCompare(a, b), Is.True);
        }
        private static SequenceVariation MakeSpanVar(int begin, int end)
        {
            // length = end - begin + 1
            int len = end - begin + 1;
            string original = new string('A', len);
            string variant = new string('V', len); // ensure sequence actually changes so AreValid passes
            return new SequenceVariation(begin, end, original, variant, "span-var");
        }

        [Test]
        public void Intersects_TruncationProduct_TrueAndFalse()
        {
            var sv = MakeSpanVar(10, 20);

            // Build truncation products
            var overlapMiddle = new TruncationProduct(15, 25, "overlap"); // overlaps (15..20)
            var entirelyBefore = new TruncationProduct(1, 9, "before");   // ends just before
            var entirelyAfter = new TruncationProduct(21, 30, "after");   // starts just after
            var touchingLeftEdge = new TruncationProduct(1, 10, "touch-left"); // end == begin of sv => intersects
            var touchingRightEdge = new TruncationProduct(20, 40, "touch-right"); // begin == end of sv => intersects

            // Reflect internal Intersects(TruncationProduct)
            var intersectsTpMethod = typeof(SequenceVariation).GetMethod(
                "Intersects",
                BindingFlags.Instance | BindingFlags.NonPublic,
                binder: null,
                types: new[] { typeof(TruncationProduct) },
                modifiers: null);

            Assert.That(intersectsTpMethod, Is.Not.Null, "Could not reflect Intersects(TruncationProduct).");

            bool Invoke(TruncationProduct tp) => (bool)intersectsTpMethod.Invoke(sv, new object[] { tp });

            Assert.Multiple(() =>
            {
                Assert.That(Invoke(overlapMiddle), Is.True, "Expected overlap in middle");
                Assert.That(Invoke(entirelyBefore), Is.False, "Expected no overlap (before)");
                Assert.That(Invoke(entirelyAfter), Is.False, "Expected no overlap (after)");
                Assert.That(Invoke(touchingLeftEdge), Is.True, "Expected intersection at left boundary");
                Assert.That(Invoke(touchingRightEdge), Is.True, "Expected intersection at right boundary");
            });
        }

        [Test]
        public void Intersects_Position_TrueAndFalse()
        {
            var sv = MakeSpanVar(100, 110); // inclusive span 100-110

            // Reflect internal Intersects(int)
            var intersectsPosMethod = typeof(SequenceVariation).GetMethod(
                "Intersects",
                BindingFlags.Instance | BindingFlags.NonPublic,
                binder: null,
                types: new[] { typeof(int) },
                modifiers: null);

            Assert.That(intersectsPosMethod, Is.Not.Null, "Could not reflect Intersects(int).");

            bool Invoke(int pos) => (bool)intersectsPosMethod.Invoke(sv, new object[] { pos });

            Assert.Multiple(() =>
            {
                // Inside
                Assert.That(Invoke(100), Is.True, "Begin boundary");
                Assert.That(Invoke(105), Is.True, "Middle position");
                Assert.That(Invoke(110), Is.True, "End boundary");

                // Outside
                Assert.That(Invoke(99), Is.False, "Just before");
                Assert.That(Invoke(111), Is.False, "Just after");
            });
        }
        [Test]
        public void SplitPerGenotype_EarlyReturn_WhenVcfHasFewerThanTenColumns()
        {
            // Truncated VCF (8 columns: CHROM POS ID REF ALT QUAL FILTER INFO) – NO FORMAT/SAMPLES
            // This causes VariantCallFormat to mark IsTruncated and leave Genotypes empty.
            string truncatedVcf =
                "1\t100\t.\tA\tT\t.\tPASS\tANN=T|missense_variant";

            // Construct a valid SequenceVariation (sequence actually changes so AreValid passes)
            var sv = new SequenceVariation(
                oneBasedBeginPosition: 34,
                oneBasedEndPosition: 34,
                originalSequence: "K",
                variantSequence: "N",
                description: "truncated-vcf",
                truncatedVcf);

            Assert.That(sv.VariantCallFormatData, Is.Not.Null);
            // Normally this would trigger the FIRST early return (Genotypes empty).
            // To specifically cover the vcfFields.Length < 10 branch, we artificially add a fake genotype.
            sv.VariantCallFormatData.Genotypes.Add("0", new[] { "0", "1" });

            // Act
            var perSample = sv.SplitPerGenotype();

            // Because the underlying raw line still has <10 tab-delimited fields,
            // the method hits:
            // if (vcfFields.Length < 10) { return result; }
            // producing an empty list.
            Assert.That(perSample, Is.Empty);
        }
        [Test]
        public void SplitPerGenotype_TryAdd_Success_AddsVariant()
        {
            // Valid minimal VCF line with exactly 10 tab-delimited columns (single sample)
            // Columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
            string vcf =
                "1\t100\t.\tA\tT\t.\tPASS\tANN=T|missense_variant|MODERATE|GENE1|GENE1|transcript|TX1|protein_coding|1/1|" +
                "c.100A>T|p.K34N|100/1000|34/300|34/100|0|\tGT:AD:DP\t1/1:0,10:10";

            // Create base variation (single residue substitution K->N)
            var sv = new SequenceVariation(
                oneBasedBeginPosition: 34,
                oneBasedEndPosition: 34,
                originalSequence: "K",
                variantSequence: "N",
                description: "homozygous-alt",
                vcf);

            // Act
            var perSample = sv.SplitPerGenotype();

            Assert.Multiple(() =>
            {
                Assert.That(perSample, Has.Count.EqualTo(1), "Exactly one per-sample variant expected");
                Assert.That(perSample[0].Description, Does.Contain("HomozygousAlt"), "Expected HomozygousAlt mode");
                Assert.That(perSample[0].VariantSequence, Is.EqualTo("N"));
                Assert.That(perSample[0].OriginalSequence, Is.EqualTo("K"));
            });
        }

        [Test]
        public void SplitPerGenotype_TryAdd_Failure_NoOpReferenceNotAdded()
        {
            // Heterozygous sample (0/1). includeReferenceForHeterozygous=true will attempt:
            // 1) A ref->ref "no-op" variant (invalid; SequenceVariation constructor throws; caught and skipped)
            // 2) A ref->alt valid variant (added)
            string vcf =
                "1\t200\t.\tG\tC\t.\tPASS\tANN=C|missense_variant|MODERATE|GENE2|GENE2|transcript|TX2|protein_coding|1/1|" +
                "c.200G>C|p.R67P|200/1200|67/400|67/150|0|\tGT:AD:DP\t0/1:7,6:13";

            var sv = new SequenceVariation(
                oneBasedBeginPosition: 67,
                oneBasedEndPosition: 67,
                originalSequence: "R",
                variantSequence: "P",
                description: "heterozygous",
                vcf);

            var perSample = sv.SplitPerGenotype(
                minDepth: 0,
                includeReferenceForHeterozygous: true,
                emitReferenceForHomozygousRef: false);

            Assert.Multiple(() =>
            {
                // Only the alt variant should be present (reference no-op filtered by failed TryAdd)
                Assert.That(perSample, Has.Count.EqualTo(1));
                Assert.That(perSample[0].Description, Does.Contain("HeterozygousAlt"));
                Assert.That(perSample[0].Description, Does.Not.Contain("HeterozygousRef"));
            });
        }
        [Test]
        public void CombineEquivalent_NullInput_ReturnsEmptyList()
        {
            var combined = SequenceVariation.CombineEquivalent(null);
            Assert.That(combined, Is.Empty);
        }
        private static Modification CreateValidModification(string id = "TestMod")
        {
            Assert.That(ModificationMotif.TryGetMotif("A", out var motif), Is.True, "Failed to create motif 'A'");
            // Provide minimal valid fields: OriginalId, Type, Target motif, valid location, monoisotopic mass
            return new Modification(
                _originalId: id,
                _modificationType: "TestType",
                _target: motif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 42.010565); // arbitrary positive mass
        }

        private static SequenceVariation CreateSimpleVariation()
        {
            // Valid substitution (positions equal, sequence changes) so AreValid() passes
            return new SequenceVariation(
                oneBasedBeginPosition: 10,
                oneBasedEndPosition: 10,
                originalSequence: "K",
                variantSequence: "N",
                description: "simple-sub");
        }

        [Test]
        public void TryAddModification_ReturnsFalse_WhenModificationIsNull()
        {
            var sv = CreateSimpleVariation();

            var ok = sv.TryAddModification(oneBasedPosition: 5, modification: null, out string error);

            Assert.Multiple(() =>
            {
                Assert.That(ok, Is.False);
                Assert.That(error, Is.EqualTo("Modification is null."));
                Assert.That(sv.OneBasedModifications, Is.Empty, "No modification entries should be added");
            });
        }

        [Test]
        public void TryAddModification_ReturnsFalse_WhenPositionIsNonPositive()
        {
            var sv = CreateSimpleVariation();
            var mod = CreateValidModification();

            var okZero = sv.TryAddModification(0, mod, out string errorZero);
            var okNegative = sv.TryAddModification(-3, mod, out string errorNeg);

            Assert.Multiple(() =>
            {
                Assert.That(okZero, Is.False);
                Assert.That(errorZero, Is.EqualTo("Position must be > 0."));
                Assert.That(okNegative, Is.False);
                Assert.That(errorNeg, Is.EqualTo("Position must be > 0."));
                Assert.That(sv.OneBasedModifications, Is.Empty, "No modification entries should be added");
            });
        }
        [Test]
        public void AddModifications_NullEnumerable_ReturnsZeroAndNoChanges()
        {
            // Arrange: valid substitution so SequenceVariation is valid
            var sv = new SequenceVariation(
                oneBasedBeginPosition: 12,
                oneBasedEndPosition: 12,
                originalSequence: "A",
                variantSequence: "V",
                description: "valid-sub");

            // Act
            var added = sv.AddModifications(
                modifications: null,
                throwOnFirstInvalid: false,
                out var skipped);

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(added, Is.EqualTo(0), "Expected zero affected positions for null input");
                Assert.That(skipped, Is.Null, "Skipped list should remain null when nothing processed");
                Assert.That(sv.OneBasedModifications, Is.Empty, "No modifications should have been added");
            });
        }
        private static Modification MakeMod(string id, string motif = "A", double mass = 42.010565)
        {
            Assert.That(ModificationMotif.TryGetMotif(motif, out var m), Is.True, "Failed to get motif");
            return new Modification(
                _originalId: id,
                _modificationType: "TestType",
                _target: m,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: mass);
        }

        private static SequenceVariation MakeSubstitutionVar(int begin, int end)
        {
            int len = end - begin + 1;
            string orig = new string('K', len);
            string variant = new string('N', len);
            return new SequenceVariation(begin, end, orig, variant, "sub");
        }

        private static SequenceVariation MakeDeletionVar(int begin, int end)
        {
            string orig = new string('A', end - begin + 1);
            // Deletion: variant sequence empty
            return new SequenceVariation(begin, end, orig, string.Empty, "del");
        }

        [Test]
        public void AddModifications_ThrowOnFirstInvalid_Throws()
        {
            var sv = MakeSubstitutionVar(10, 15);
            var goodMod = MakeMod("Good1");

            // First tuple invalid because position <= 0; second would be valid but never reached
            var tuples = new List<(int position, Modification modification)>
            {
                (0, goodMod),
                (12, goodMod)
            };

            var ex = Assert.Throws<ArgumentException>(() =>
                sv.AddModifications(tuples, throwOnFirstInvalid: true, out var _));

            Assert.That(ex!.Message, Does.Contain("Invalid modification at position 0: Position must be > 0."));
            Assert.That(sv.OneBasedModifications, Is.Empty);
        }

        [Test]
        public void AddModifications_SkipInvalids_CollectsSkipped()
        {
            // Deletion variant: any position >= begin (10) invalid when variantSequence == "" (termination semantics)
            var sv = MakeDeletionVar(10, 12);

            var modA = MakeMod("mA");
            var modB = MakeMod("mB");
            var modC = MakeMod("mC");

            var batch = new List<(int position, Modification modification)>
            {
                // Invalid: deletion / termination prevents mod at or after begin
                (11, modA),
                // Invalid: position <= 0
                (0, modB),
                // Invalid: null modification
                (8, null),
                // Valid: position before begin on deletion variant
                (5, modC)
            };

            int added = sv.AddModifications(batch, throwOnFirstInvalid: false, out var skipped);

            Assert.Multiple(() =>
            {
                Assert.That(added, Is.EqualTo(1), "Only one valid position should have been added");
                Assert.That(skipped, Is.Not.Null);
                Assert.That(skipped, Has.Count.EqualTo(3));

                // Extract reasons
                var reasons = skipped!.Select(s => s.reason).ToList();

                Assert.That(reasons.Any(r => r == "Position invalid for a termination or deletion at/after the begin coordinate."), Is.True);
                Assert.That(reasons.Any(r => r == "Position must be > 0."), Is.True);
                Assert.That(reasons.Any(r => r == "Modification is null."), Is.True);

                // Current implementation always supplies a concrete reason; "Unknown reason" would only appear
                // if TryAddModification returned false with a null error (not possible at present).
                Assert.That(reasons.Any(r => r == "Unknown reason"), Is.False, "Fallback 'Unknown reason' path is unreachable with current logic");
            });

            // Confirm the valid modification stored under position 5
            Assert.That(sv.OneBasedModifications.ContainsKey(5), Is.True);
            Assert.That(sv.OneBasedModifications[5], Has.Count.EqualTo(1));
            Assert.That(sv.OneBasedModifications[5][0].OriginalId, Is.EqualTo("mC"));
        }
        [Test]
        public void GetInvalidModificationPositions_YieldsAndContinues_OnNonPositivePosition()
        {
            // Create a valid variation first (original length 3, variant length 2 ? frameshift but valid)
            // Begin=10 End=12; variant length=2 => newSpanEnd = 10 + 2 - 1 = 11
            var sv = new SequenceVariation(
                oneBasedBeginPosition: 10,
                oneBasedEndPosition: 12,
                originalSequence: "AAA",
                variantSequence: "VV",
                description: "frameshift");

            // Prepare real modification instances
            Assert.That(ModificationMotif.TryGetMotif("A", out var motif), Is.True);
            var mod1 = new Modification(_originalId: "ModX", _modificationType: "TestType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.9949);
            var mod2 = new Modification(_originalId: "ModY", _modificationType: "TestType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 42.0106);

            // Directly inject invalid modification positions (bypass TryAddModification which rejects them):
            // -1 (<=0) triggers: yield return pos; continue;
            // 12 is inside edited region (10–12) but > newSpanEnd (11) ? also invalid
            sv.OneBasedModifications[-1] = new List<Modification> { mod1 };
            sv.OneBasedModifications[12] = new List<Modification> { mod2 };

            // Reflect the private iterator method
            var method = typeof(SequenceVariation)
                .GetMethod("GetInvalidModificationPositions", BindingFlags.Instance | BindingFlags.NonPublic);
            Assert.That(method, Is.Not.Null);

            var enumerable = (IEnumerable<int>)method.Invoke(sv, Array.Empty<object>());
            var invalidList = enumerable.ToList();

            Assert.Multiple(() =>
            {
                Assert.That(invalidList, Has.Count.EqualTo(2), "Expected two invalid positions");
                Assert.That(invalidList, Does.Contain(-1), "Non-positive position should be reported");
                Assert.That(invalidList, Does.Contain(12), "Position beyond new variant span should be reported");
            });
        }
    }
}