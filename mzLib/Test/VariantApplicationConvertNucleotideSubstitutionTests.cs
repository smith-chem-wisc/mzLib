using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;

namespace Test.DatabaseTests
{
    [TestFixture]
    public class VariantApplicationConvertNucleotideSubstitutionTests
    {
        // Helper to create a minimal substitution modification matching the required detection pattern
        private static Modification Substitution(string idArrow)
            => new Modification(
                idArrow,                  // OriginalId must contain "X->Y"
                null,                     // accession
                "1 nucleotide substitution", // ModificationType must contain this substring
                null,                     // secondary accession / source
                null,                     // motif (irrelevant here)
                "Anywhere.",              // position restriction
                null,                     // feature type
                0,                        // mass delta (not relevant for test)
                null, null, null, null, null, null);

        // Non-substitution (should be ignored)
        private static Modification Other(string id, double mass = 15.9949)
            => new Modification(
                id,
                null,
                "oxidation",
                null,
                null,
                "Anywhere.",
                null,
                mass,
                null, null, null, null, null, null);

        // Malformed substitution (no "->" pattern) must be ignored
        private static Modification Malformed()
            => new Modification(
                "E>A",
                null,
                "1 nucleotide substitution",
                null,
                null,
                "Anywhere.",
                null,
                0,
                null, null, null, null, null, null);

        [Test]
        public void ConvertNucleotideSubstitutionModificationsToSequenceVariants_Comprehensive()
        {
            // Sequence indices (1-based):
            // 1 M, 2 A, 3 E, 4 W, 5 P, 6 Q, 7 K
            var protein = new Protein("MAEWPQK", "TEST_PROT");

            // Seed: ensure dictionaries exist (Protein constructor normally does this, but be defensive)
            Assert.That(protein.OneBasedPossibleLocalizedModifications, Is.Not.Null);
            Assert.That(protein.OriginalNonVariantModifications, Is.Not.Null);
            Assert.That(protein.ConsensusVariant, Is.Not.Null);

            // Substitution modifications to be converted
            var modEtoA = Substitution("E->A"); // position 3
            var modWtoK = Substitution("W->K"); // position 4

            // Non-substitution modification (should remain)
            var modOxidP = Other("Oxidation_P"); // position 5

            // Malformed substitution (contains correct modification type but no "->" pattern in OriginalId)
            var malformed = Malformed(); // position 6

            // Populate modification dictionaries (both possible localized & original non-variant)
            AddMod(protein, 3, modEtoA);
            AddMod(protein, 4, modWtoK);
            AddMod(protein, 5, modOxidP);
            AddMod(protein, 6, malformed);

            // Pre-existing variant matching W->K (should prevent duplicate)
            var preExistingWtoK = new SequenceVariation(4, 4, "W", "K", "Existing substitution");
            protein.SequenceVariations.Add(preExistingWtoK);
            Assert.That(protein.SequenceVariations.Count, Is.EqualTo(1), "Precondition failed: pre-existing variant not added.");

            // Capture snapshot counts
            int initialModKeyCount = protein.OneBasedPossibleLocalizedModifications.Count;
            Assert.That(initialModKeyCount, Is.EqualTo(4));

            // Invoke conversion
            protein.ConvertNucleotideSubstitutionModificationsToSequenceVariants();

            // EXPECTATIONS:
            // 1. A new variant for E3->A (position 3) added.
            // 2. No duplicate variant for W4->K (still exactly one at position 4).
            // 3. Modifications at positions 3 & 4 removed from:
            //    - OneBasedPossibleLocalizedModifications
            //    - OriginalNonVariantModifications
            //    - ConsensusVariant mirrored dictionaries
            // 4. Unrelated oxidation mod (position 5) retained.
            // 5. Malformed substitution (position 6) retained (not converted).
            // 6. Description of newly created SequenceVariation is "Putative GPTMD Substitution".

            // Variants present
            var variants = protein.SequenceVariations;
            Assert.That(variants.Count, Is.EqualTo(2), "Exactly two variants expected (pre-existing W->K + new E->A).");

            var eToAVariant = variants.SingleOrDefault(v => v.OneBasedBeginPosition == 3
                                                            && v.OneBasedEndPosition == 3
                                                            && v.OriginalSequence == "E"
                                                            && v.VariantSequence == "A");
            Assert.That(eToAVariant, Is.Not.Null, "E->A variant missing.");
            Assert.That(eToAVariant.Description, Is.EqualTo("Putative GPTMD Substitution"),
                "E->A variant should use standardized description.");

            var wToKVariantMatches = variants.Where(v => v.OneBasedBeginPosition == 4
                                                         && v.OneBasedEndPosition == 4
                                                         && v.OriginalSequence == "W"
                                                         && v.VariantSequence == "K")
                                             .ToList();
            Assert.That(wToKVariantMatches.Count, Is.EqualTo(1),
                "Pre-existing W->K variant should not be duplicated.");

            // Modifications removed at positions 3 and 4
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(3), Is.False,
                "Converted mod (E->A) should be removed from OneBasedPossibleLocalizedModifications.");
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(4), Is.False,
                "Converted mod (W->K) should be removed from OneBasedPossibleLocalizedModifications.");

            Assert.That(protein.OriginalNonVariantModifications.ContainsKey(3), Is.False,
                "Converted mod (E->A) should be removed from OriginalNonVariantModifications.");
            Assert.That(protein.OriginalNonVariantModifications.ContainsKey(4), Is.False,
                "Converted mod (W->K) should be removed from OriginalNonVariantModifications.");

            // Consensus variant dictionaries mirror removal
            Assert.That(protein.ConsensusVariant.OneBasedPossibleLocalizedModifications.ContainsKey(3), Is.False);
            Assert.That(protein.ConsensusVariant.OneBasedPossibleLocalizedModifications.ContainsKey(4), Is.False);
            Assert.That(protein.ConsensusVariant.OriginalNonVariantModifications.ContainsKey(3), Is.False);
            Assert.That(protein.ConsensusVariant.OriginalNonVariantModifications.ContainsKey(4), Is.False);

            // Unaffected modifications remain (position 5 & 6)
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(5), Is.True,
                "Non-substitution modification at position 5 should remain.");
            Assert.That(protein.OneBasedPossibleLocalizedModifications[5]
                .Any(m => m.OriginalId == "Oxidation_P"), Is.True);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(6), Is.True,
                "Malformed substitution at position 6 should remain.");
            Assert.That(protein.OneBasedPossibleLocalizedModifications[6]
                .Any(m => m.OriginalId == "E>A"), Is.True);

            // Ensure removal did not accidentally clear unrelated keys
            Assert.That(protein.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2),
                "Unexpected modification key removals (expected only positions 3 & 4 removed).");
        }

        private static void AddMod(Protein protein, int position, Modification mod)
        {
            if (!protein.OneBasedPossibleLocalizedModifications.TryGetValue(position, out var list1))
            {
                list1 = new List<Modification>();
                protein.OneBasedPossibleLocalizedModifications[position] = list1;
            }
            list1.Add(mod);

            if (!protein.OriginalNonVariantModifications.TryGetValue(position, out var list2))
            {
                list2 = new List<Modification>();
                protein.OriginalNonVariantModifications[position] = list2;
            }
            list2.Add(mod);

            // Mirror expected initial state in consensus variant as constructor usually does
            if (!protein.ConsensusVariant.OneBasedPossibleLocalizedModifications.TryGetValue(position, out var list3))
            {
                list3 = new List<Modification>();
                protein.ConsensusVariant.OneBasedPossibleLocalizedModifications[position] = list3;
            }
            list3.Add(mod);

            if (!protein.ConsensusVariant.OriginalNonVariantModifications.TryGetValue(position, out var list4))
            {
                list4 = new List<Modification>();
                protein.ConsensusVariant.OriginalNonVariantModifications[position] = list4;
            }
            list4.Add(mod);
        }
    }
}