using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Omics.BioPolymer;
using UsefulProteomicsDatabases;
using Proteomics;

namespace Test
{
    [TestFixture]
    public class TestDecoyProteinGenerator
    {
        //[Test]
        //public void TestReverseDecoySequenceVariationCoordinates()
        //{
        //    // Target sequence: M A C D E F G H I K (10 aa)
        //    string targetSequence = "MACDEFGHIK";

        //    // Sequence variations:
        //    // 1. On M at position 1: M -> MM
        //    var variationOnM = new SequenceVariation(
        //        1, 1, "M", "MM", "N-terminal extension"
        //    );

        //    // 2. On E at position 5: E -> Q
        //    var variationOnE = new SequenceVariation(
        //        5, 5, "E", "Q", "Middle substitution"
        //    );

        //    // 3. On K at position 10: K -> R
        //    var variationOnK = new SequenceVariation(
        //        10, 10, "K", "R", "C-terminal substitution"
        //    );

        //    // Create the target protein with the sequence variations
        //    var targetProtein = new Protein(
        //        targetSequence,
        //        "TestProtein",
        //        sequenceVariations: new List<SequenceVariation> { variationOnM, variationOnE, variationOnK }
        //    );

        //    // Generate the reverse decoy
        //    var decoys = DecoyProteinGenerator.GenerateDecoys(
        //        new List<Protein> { targetProtein },
        //        DecoyType.Reverse
        //    );

        //    // Validate the decoy
        //    Assert.That(decoys.Count, Is.EqualTo(1));
        //    var decoy = decoys[0];

        //    // Decoy sequence: M K I H G F E D C A (M at position 1, K at position 2, etc.)
        //    Assert.That(decoy.BaseSequence, Is.EqualTo("MKIHGFEDCA"));

        //    // Validate sequence variations in the decoy
        //    var decoyVariations = decoy.SequenceVariations;
        //    Assert.That(decoyVariations.Count, Is.EqualTo(3));

        //    // Use Assert.Multiple to evaluate all assertions
        //    Assert.Multiple(() =>
        //    {
        //        // 1. Variant on M at position 1 in target should remain at position 1 in decoy
        //        var decoyVariationOnM = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 1);
        //        Assert.That(decoyVariationOnM, Is.Not.Null, "Decoy variant on M at position 1 should exist.");
        //        Assert.That(decoyVariationOnM.OriginalSequence, Is.EqualTo("M"));
        //        Assert.That(decoyVariationOnM.VariantSequence, Is.EqualTo("MM"));

        //        // 2. Variant on E at position 5 in target should map to position 6 in decoy
        //        var decoyVariationOnE = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 6);
        //        Assert.That(decoyVariationOnE, Is.Not.Null, "Decoy variant on E at position 5 should map to position 6.");
        //        Assert.That(decoyVariationOnE.OriginalSequence, Is.EqualTo("E"));
        //        Assert.That(decoyVariationOnE.VariantSequence, Is.EqualTo("Q"));

        //        // 3. Variant on K at position 10 in target should map to position 2 in decoy
        //        var decoyVariationOnK = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 2);
        //        Assert.That(decoyVariationOnK, Is.Not.Null, "Decoy variant on K at position 10 should map to position 2.");
        //        Assert.That(decoyVariationOnK.OriginalSequence, Is.EqualTo("K"));
        //        Assert.That(decoyVariationOnK.VariantSequence, Is.EqualTo("R"));
        //    });
        //}
        [Test]
        public void TestReverseDecoySingleSequenceVariation()
        {
            // Target sequence: M A C D E F G H I K (10 aa)
            string targetSequence = "MACDEFGHIK";

            // Sequence variation: C -> Z at position 3
            var variationOnC = new SequenceVariation(
                3, 3, "C", "Z", "Single substitution"
            );

            // Create the target protein with the sequence variation
            var targetProtein = new Protein(
                targetSequence,
                "TestProtein",
                sequenceVariations: new List<SequenceVariation> { variationOnC }
            );

            // Generate the reverse decoy
            var decoys = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { targetProtein },
                DecoyType.Reverse
            );

            // Validate the decoy
            Assert.That(decoys.Count, Is.EqualTo(1));
            var decoy = decoys[0];

            // Expected reverse decoy sequence: M K I H G F E D C A
            Assert.That(decoy.BaseSequence, Is.EqualTo("MKIHGFEDCA"));

            // Validate sequence variations in the decoy
            var decoyVariations = decoy.SequenceVariations;
            Assert.That(decoyVariations.Count, Is.EqualTo(1));

            // Use Assert.Multiple to evaluate all assertions
            Assert.Multiple(() =>
            {
                // Variant on C at position 3 in target should map to position 8 in decoy
                var decoyVariationOnC = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 9);
                Assert.That(decoyVariationOnC, Is.Not.Null, "Decoy variant on C at position 3 should map to position 9.");
                Assert.That(decoyVariationOnC.OriginalSequence, Is.EqualTo("C"));
                Assert.That(decoyVariationOnC.VariantSequence, Is.EqualTo("Z"));
            });
        }
        [Test]
        public void TestReverseDecoySingleSequenceVariationWithInsertion()
        {
            // Target sequence: M A C D E F G H I K (10 aa)
            string targetSequence = "MACDEFGHIK";

            // Sequence variation: C -> ZZ at position 3
            var variationOnC = new SequenceVariation(
                3, 3, "C", "ZZ", "Single substitution with insertion"
            );

            // Create the target protein with the sequence variation
            var targetProtein = new Protein(
                targetSequence,
                "TestProtein",
                sequenceVariations: new List<SequenceVariation> { variationOnC }
            );

            // Generate the reverse decoy
            var decoys = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { targetProtein },
                DecoyType.Reverse
            );

            // Validate the decoy
            Assert.That(decoys.Count, Is.EqualTo(1));
            var decoy = decoys[0];

            // Expected reverse decoy sequence: M K I H G F E D C A
            Assert.That(decoy.BaseSequence, Is.EqualTo("MKIHGFEDCA"));

            // Validate sequence variations in the decoy
            var decoyVariations = decoy.SequenceVariations;
            Assert.That(decoyVariations.Count, Is.EqualTo(1));

            // Use Assert.Multiple to evaluate all assertions
            Assert.Multiple(() =>
            {
                // Variant on C at position 3 in target should map to position 9 in decoy
                var decoyVariationOnC = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 9);
                Assert.That(decoyVariationOnC, Is.Not.Null, "Decoy variant on C at position 3 should map to position 9.");
                Assert.That(decoyVariationOnC.OriginalSequence, Is.EqualTo("C"));
                Assert.That(decoyVariationOnC.VariantSequence, Is.EqualTo("ZZ"));
            });
        }
        [Test]
        public void TestReverseDecoySingleSequenceVariationWithAcetylation()
        {
            // Target sequence: M A C D E F G H I K (10 aa)
            string targetSequence = "MACDEFGHIK";

            // Sequence variation: C -> Z at position 3
            var variationOnC = new SequenceVariation(
                3, 3, "C", "Z", "Single substitution"
            );
            // Create a ModificationMotif for lysine (K)
            ModificationMotif.TryGetMotif("K", out var lysineMotif);
            // Add acetylation modification on K at position 10
            var acetylation = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Acetyl",
                _target: lysineMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 42.0106
            );
            var modifications = new Dictionary<int, List<Modification>>
            {
                { 10, new List<Modification> { acetylation } } // Lysine at position 10
            };

            // Create the target protein with the modification and sequence variation
            var targetProtein = new Protein(
                targetSequence,
                "TestProtein",
                oneBasedModifications: modifications, // Apply the modification to the target protein
                sequenceVariations: new List<SequenceVariation> { variationOnC }
            );

            // Generate the reverse decoy
            var decoys = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { targetProtein },
                DecoyType.Reverse
            );

            // Validate the decoy
            Assert.That(decoys.Count, Is.EqualTo(1));
            var decoy = decoys[0];

            // Expected reverse decoy sequence: M K I H G F E D C A
            Assert.That(decoy.BaseSequence, Is.EqualTo("MKIHGFEDCA"));

            // Validate sequence variations in the decoy
            var decoyVariations = decoy.SequenceVariations;
            Assert.That(decoyVariations.Count, Is.EqualTo(1));

            // Validate modifications in the decoy
            var decoyModifications = decoy.OneBasedPossibleLocalizedModifications;
            Assert.That(decoyModifications.Count, Is.EqualTo(1));

            // Use Assert.Multiple to evaluate all assertions
            Assert.Multiple(() =>
            {
                // Variant on C at position 3 in target should map to position 9 in decoy
                var decoyVariationOnC = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 9);
                Assert.That(decoyVariationOnC, Is.Not.Null, "Decoy variant on C at position 3 should map to position 9.");
                Assert.That(decoyVariationOnC.OriginalSequence, Is.EqualTo("C"));
                Assert.That(decoyVariationOnC.VariantSequence, Is.EqualTo("Z"));

                // Acetylation on K at position 10 in target should map to position 2 in decoy
                Assert.That(decoyModifications.ContainsKey(2), Is.True, "Acetylation on K at position 10 in target should map to position 2 in decoy.");
                Assert.That(decoyModifications[2].Any(mod => mod.ToString() == acetylation.ToString()), Is.True, "Decoy modification at position 2 should be acetylation.");
            });
        }
        [Test]
        public void TestReverseDecoySequenceVariationWithModificationOnVariant()
        {
            // Target sequence: M A C D E F G H I K (10 aa)
            string targetSequence = "MACDEFGHIK";

            // Sequence variation: C -> K at position 3
            // Create a ModificationMotif for lysine (K)
            ModificationMotif.TryGetMotif("K", out var lysineMotif);

            // Add acetylation modification on K at position 3
            var acetylation = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Acetyl",
                _target: lysineMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 42.0106
            );
            var variantModifications = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { acetylation } } // Acetylation on K at position 3
            };

            var variationOnC = new SequenceVariation(
                3, 3, "C", "K", "Single substitution with modification",
                oneBasedModifications: variantModifications
            );

            // Create the target protein with the sequence variation
            var targetProtein = new Protein(
                targetSequence,
                "TestProtein",
                sequenceVariations: new List<SequenceVariation> { variationOnC }
            );

            // Generate the reverse decoy
            var decoys = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { targetProtein },
                DecoyType.Reverse
            );

            // Validate the decoy
            Assert.That(decoys.Count, Is.EqualTo(1));
            var decoy = decoys[0];

            // Expected reverse decoy sequence: M K I H G F E D C A
            Assert.That(decoy.BaseSequence, Is.EqualTo("MKIHGFEDCA"));

            // Validate sequence variations in the decoy
            var decoyVariations = decoy.SequenceVariations;
            Assert.That(decoyVariations.Count, Is.EqualTo(1));

            // Use Assert.Multiple to evaluate all assertions
            Assert.Multiple(() =>
            {
                // Variant on C at position 3 in target should map to position 9 in decoy
                var decoyVariationOnC = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 9);
                Assert.That(decoyVariationOnC, Is.Not.Null, "Decoy variant on C at position 3 should map to position 9.");
                Assert.That(decoyVariationOnC.OriginalSequence, Is.EqualTo("C"));
                Assert.That(decoyVariationOnC.VariantSequence, Is.EqualTo("K"));

                // Validate the modification on the variant
                Assert.That(decoyVariationOnC.OneBasedModifications.ContainsKey(9), Is.True, "Acetylation on K at position 3 in target should map to position 9 in decoy.");
                Assert.That(decoyVariationOnC.OneBasedModifications[9].Any(mod => mod.ToString() == acetylation.ToString()), Is.True, "Decoy modification at position 9 should be acetylation.");
            });
        }
        [Test]
        public void TestReverseDecoySequenceVariationWithModificationOnVariantAndInsertion()
        {
            // Target sequence: M A C D E F G H I K (10 aa)
            string targetSequence = "MACDEFGHIK";

            // Sequence variation: C -> KR at position 3
            // Create a ModificationMotif for lysine (K)
            ModificationMotif.TryGetMotif("K", out var lysineMotif);

            // Add acetylation modification on K at position 3
            var acetylation = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Acetyl",
                _target: lysineMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 42.0106
            );
            var variantModifications = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { acetylation } } // Acetylation on K at position 3
            };

            var variationOnC = new SequenceVariation(
                3, 3, "C", "KR", "Single substitution with insertion and modification",
                oneBasedModifications: variantModifications
            );

            // Create the target protein with the sequence variation
            var targetProtein = new Protein(
                targetSequence,
                "TestProtein",
                sequenceVariations: new List<SequenceVariation> { variationOnC }
            );

            // Generate the reverse decoy
            var decoys = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { targetProtein },
                DecoyType.Reverse
            );

            // Validate the decoy
            Assert.That(decoys.Count, Is.EqualTo(1));
            var decoy = decoys[0];

            // Expected reverse decoy sequence: M K I H G F E D C A
            Assert.That(decoy.BaseSequence, Is.EqualTo("MKIHGFEDCA"));

            // Validate sequence variations in the decoy
            var decoyVariations = decoy.SequenceVariations;
            Assert.That(decoyVariations.Count, Is.EqualTo(1));

            // Use Assert.Multiple to evaluate all assertions
            Assert.Multiple(() =>
            {
                // Variant on C at position 3 in target should map to position 9 in decoy
                var decoyVariationOnC = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 9);
                Assert.That(decoyVariationOnC, Is.Not.Null, "Decoy variant on C at position 3 should map to position 9.");
                Assert.That(decoyVariationOnC.OriginalSequence, Is.EqualTo("C"));
                Assert.That(decoyVariationOnC.VariantSequence, Is.EqualTo("KR"));

                // Validate the modification on the variant
                Assert.That(decoyVariationOnC.OneBasedModifications.ContainsKey(9), Is.True, "Acetylation on K at position 3 in target should map to position 9 in decoy.");
                Assert.That(decoyVariationOnC.OneBasedModifications[9].Any(mod => mod.ToString() == acetylation.ToString()), Is.True, "Decoy modification at position 9 should be acetylation.");
            });
        }
        [Test]
        public void TestReverseDecoySequenceVariationWithModificationOnVariantAndProtein()
        {
            // Target sequence: M A C D E F G H I K (10 aa)
            string targetSequence = "MACDEFGHIK";

            // Sequence variation: C -> KR at position 3
            // Create a ModificationMotif for lysine (K)
            ModificationMotif.TryGetMotif("K", out var lysineMotif);

            // Add acetylation modification on K at position 3 (in the sequence variant)
            var acetylationOnVariant = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Acetyl",
                _target: lysineMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 42.0106
            );
            var variantModifications = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { acetylationOnVariant } } // Acetylation on K at position 3
            };

            var variationOnC = new SequenceVariation(
                3, 3, "C", "KR", "Single substitution with insertion and modification",
                oneBasedModifications: variantModifications
            );

            // Add acetylation modification on K at position 10 (in the protein)
            var acetylationOnProtein = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Acetyl",
                _target: lysineMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 42.0106
            );
            var proteinModifications = new Dictionary<int, List<Modification>>
            {
                { 10, new List<Modification> { acetylationOnProtein } } // Acetylation on K at position 10
            };

            // Create the target protein with the sequence variation and protein modification
            var targetProtein = new Protein(
                targetSequence,
                "TestProtein",
                oneBasedModifications: proteinModifications, // Apply the modification to the protein
                sequenceVariations: new List<SequenceVariation> { variationOnC }
            );

            // Generate the reverse decoy
            var decoys = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { targetProtein },
                DecoyType.Reverse
            );

            // Validate the decoy
            Assert.That(decoys.Count, Is.EqualTo(1));
            var decoy = decoys[0];

            // Expected reverse decoy sequence: M K I H G F E D C A
            Assert.That(decoy.BaseSequence, Is.EqualTo("MKIHGFEDCA"));

            // Validate sequence variations in the decoy
            var decoyVariations = decoy.SequenceVariations;
            Assert.That(decoyVariations.Count, Is.EqualTo(1));

            // Validate modifications in the decoy
            var decoyModifications = decoy.OneBasedPossibleLocalizedModifications;
            Assert.That(decoyModifications.Count, Is.EqualTo(1)); // one from the protein. THERE IS ALSO ONE ON THE VARIANT BUT IT HAS NOT BEEN APPLIED TO THE PROTEIN

            var sequenceVariantModifications = decoyVariations.SelectMany(v => v.OneBasedModifications).SelectMany(kvp => kvp.Value).Count();
            Assert.That(sequenceVariantModifications, Is.EqualTo(1));

            // Use Assert.Multiple to evaluate all assertions
            Assert.Multiple(() =>
            {
                // Variant on C at position 3 in target should map to position 9 in decoy
                var decoyVariationOnC = decoyVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 9);
                Assert.That(decoyVariationOnC, Is.Not.Null, "Decoy variant on C at position 3 should map to position 9.");
                Assert.That(decoyVariationOnC.OriginalSequence, Is.EqualTo("C"));
                Assert.That(decoyVariationOnC.VariantSequence, Is.EqualTo("KR"));

                // Validate the modification on the variant
                Assert.That(decoyVariationOnC.OneBasedModifications.ContainsKey(9), Is.True, "Acetylation on K at position 3 in target should map to position 9 in decoy.");
                Assert.That(decoyVariationOnC.OneBasedModifications[9].Any(mod => mod.ToString() == acetylationOnVariant.ToString()), Is.True, "Decoy modification at position 9 should be acetylation.");

                // Validate the modification on the protein
                Assert.That(decoyModifications.ContainsKey(2), Is.True, "Acetylation on K at position 10 in target should map to position 2 in decoy.");
                Assert.That(decoyModifications[2].Any(mod => mod.ToString() == acetylationOnProtein.ToString()), Is.True, "Decoy modification at position 2 should be acetylation.");
            });
        }
    }
}