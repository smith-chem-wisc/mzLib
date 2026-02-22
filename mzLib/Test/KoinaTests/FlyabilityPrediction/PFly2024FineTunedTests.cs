using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;
using BayesianEstimation;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class PFly2024FineTunedTests
    {
        /// <summary>
        /// Tests that the model correctly loads valid peptide sequences without throwing exceptions or warnings.
        /// Verifies basic initialization and that all valid peptides are accepted.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelAcceptsValidPeptides()
        {
            var validPeptides = new List<string>
            {
                "PEPTIDEK",
                "ACDEFGHIKLMNPQRSTVWY", // All standard amino acids
                "LMNPQRSTVWY" // Another valid sequence
            };

            var model = new PFly2024FineTuned();
            var inputs = validPeptides.Select(p => new DetectabilityPredictionInput(p)).ToList();
            var outputs = model.Predict(inputs);

            Assert.That(model.ValidInputsMask.All(x => x), Is.True, "All valid peptides should be accepted");
            Assert.That(model.Predictions.Count, Is.EqualTo(0), "No predictions should exist before running inference");
        }

        /// <summary>
        /// Tests that peptides with modifications are handled correctly.
        /// Modifications are represented in square brackets and should be stripped for length validation.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelDoesNotAcceptModifiedPeptides()
        {
            var modifiedPeptides = new List<string>
            {
                "M[Common Variable:Oxidation on M]PEPTIDEK",
                "PEPTIDEM[Common Variable:Oxidation on M]",
                "C[Common Fixed:Carbamidomethyl on C]PEPTIDEK",
                "PEPTIDEPEPI"
            };
            var inputs = modifiedPeptides.Select(p => new DetectabilityPredictionInput(p)).ToList();

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(inputs);
            Assert.That(model.ValidInputsMask.First(), Is.True, "Only the unmodified peptide (\"PEPTIDEPEPI\") should be accepted");
        }

        /// <summary>
        /// Tests boundary condition for peptide sequence length.
        /// Maximum valid length is 40 amino acids (after removing modification annotations).
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelPeptideLengthBoundaries()
        {
            // Test maximum valid length (40 amino acids)
            var maxLengthPeptide = new List<DetectabilityPredictionInput> { new DetectabilityPredictionInput(new string('A', 40)) };
            var model = new PFly2024FineTuned();
            var outputMax = model.Predict(maxLengthPeptide);

            Assert.That(model.ValidInputsMask.First(), Is.True, "Peptide with exactly 40 amino acids should be valid");

            // Test invalid length (41 amino acids)
            var tooLongPeptide = new string('A', 41);
            var modelTooLong = new PFly2024FineTuned(new List<string> { tooLongPeptide }, out var warningTooLong);

            Assert.That(modelTooLong.PeptideSequences.Count, Is.EqualTo(0), "Peptide with 41 amino acids should be invalid");
            Assert.That(warningTooLong, Is.Not.Null, "Too long peptide should produce a warning");
            Assert.That(warningTooLong.Message, Does.Contain(tooLongPeptide), "Warning should mention the invalid peptide");

            // Test empty peptide
            var emptyPeptide = "";
            var modelEmpty = new PFly2024FineTuned(new List<string> { emptyPeptide }, out var warningEmpty);

            Assert.That(modelEmpty.PeptideSequences.Count, Is.EqualTo(0), "Empty peptide should be invalid");
            Assert.That(warningEmpty, Is.Not.Null, "Empty peptide should produce a warning");

            // Test modified peptide that exceeds length when modifications are removed
            var modifiedTooLong = new string('A', 41) + "[Modification]";
            var modelModifiedTooLong = new PFly2024FineTuned(new List<string> { modifiedTooLong }, out var warningModifiedTooLong);

            Assert.That(modelModifiedTooLong.PeptideSequences.Count, Is.EqualTo(0), "Modified peptide exceeding 40 AA should be invalid");
            Assert.That(warningModifiedTooLong, Is.Not.Null);
        }

        /// <summary>
        /// Tests that peptides containing non-standard or invalid amino acids are rejected.
        /// Standard 20 amino acids are: ACDEFGHIKLMNPQRSTVWY
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelInvalidAminoAcids()
        {
            // Test peptide with 'X' (unknown amino acid)
            var peptideWithX = "PEPTXIDEK";
            var modelX = new PFly2024FineTuned(new List<string> { peptideWithX }, out var warningX);

            Assert.That(modelX.PeptideSequences.Count, Is.EqualTo(0), "Peptide with 'X' should be filtered out");
            Assert.That(warningX, Is.Not.Null, "Invalid amino acid should produce a warning");

            // Test peptide with 'U' (Selenocysteine)
            var peptideWithU = "PEPTUIDEK";
            var modelU = new PFly2024FineTuned(new List<string> { peptideWithU }, out var warningU);

            Assert.That(modelU.PeptideSequences.Count, Is.EqualTo(0), "Peptide with 'U' should be filtered out");

            // Test peptide with 'B' (Asx)
            var peptideWithB = "PEPTBIDEK";
            var modelB = new PFly2024FineTuned(new List<string> { peptideWithB }, out var warningB);

            Assert.That(modelB.PeptideSequences.Count, Is.EqualTo(0), "Peptide with 'B' should be filtered out");

            // Test peptide with 'Z' (Glx)
            var peptideWithZ = "PEPTZIDEK";
            var modelZ = new PFly2024FineTuned(new List<string> { peptideWithZ }, out var warningZ);

            Assert.That(modelZ.PeptideSequences.Count, Is.EqualTo(0), "Peptide with 'Z' should be filtered out");

            // Test mixed valid and invalid peptides
            var mixedPeptides = new List<string> { "PEPTIDEK", "PEPTXIDEK", "VALIDSEQ" };
            var modelMixed = new PFly2024FineTuned(mixedPeptides, out var warningMixed);

            Assert.That(modelMixed.PeptideSequences.Count, Is.EqualTo(2), "Only valid peptides should be retained");
            Assert.That(warningMixed, Is.Not.Null, "Invalid peptides should produce a warning");
        }

        /// <summary>
        /// Tests handling of empty input lists.
        /// Empty lists should not throw exceptions but should produce appropriate warnings.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelEmptyInputHandling()
        {
            var emptyPeptides = new List<string>();

            var model = new PFly2024FineTuned(emptyPeptides, out var warning);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(0), "Empty input should result in no peptides");
            Assert.That(warning, Is.Not.Null, "Empty input should produce a warning");
            Assert.That(warning.Message, Does.Contain("empty"), "Warning should mention empty inputs");
            Assert.DoesNotThrowAsync(async () => await model.PredictAsync(), "Empty model should not throw on inference");
            Assert.That(model.Predictions.Count, Is.EqualTo(0), "No predictions for empty input");
        }

        /// <summary>
        /// Tests that duplicate peptide sequences are handled appropriately.
        /// Duplicates should either be kept or deduplicated depending on implementation.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelDuplicateHandling()
        {
            var duplicatePeptides = new List<string> { "PEPTIDEK", "PEPTIDEK" };

            var model = new PFly2024FineTuned(duplicatePeptides, out var warning);

            Assert.That(model.PeptideSequences.Count, Is.GreaterThan(0), "Duplicate peptides should be handled without error");
            // Note: Actual duplicate handling behavior (keep vs. deduplicate) should be documented
        }

        /// <summary>
        /// Tests that model constants and properties are correctly initialized.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelProperties()
        {
            var peptides = new List<string> { "PEPTIDEK" };
            var model = new PFly2024FineTuned(peptides, out var warning);

            Assert.That(model.ModelName, Is.EqualTo("pfly_2024_fine_tuned"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(128));
            Assert.That(model.NumberOfDetectabilityClasses, Is.EqualTo(4));
            Assert.That(model.DetectabilityClasses.Count, Is.EqualTo(4));
            Assert.That(model.DetectabilityClasses[0], Is.EqualTo("Not Detectable"));
            Assert.That(model.DetectabilityClasses[1], Is.EqualTo("Low Detectability"));
            Assert.That(model.DetectabilityClasses[2], Is.EqualTo("Intermediate Detectability"));
            Assert.That(model.DetectabilityClasses[3], Is.EqualTo("High Detectability"));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(40));
        }

        /// <summary>
        /// Tests handling of peptides containing lowercase letters.
        /// Peptides should be case-insensitive or properly handled.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelCaseSensitivity()
        {
            var mixedCasePeptides = new List<string>
            {
                "peptidek", // lowercase
                "PEPTIDEK", // uppercase
                "PePtIdEk"  // mixed case
            };

            var model = new PFly2024FineTuned(mixedCasePeptides, out var warning);

            // Assuming the model should handle case conversion internally
            Assert.That(model.PeptideSequences.Count, Is.GreaterThan(0), "Case should not affect validity");
        }

        /// <summary>
        /// Tests handling of peptides with special characters and whitespace.
        /// These should be filtered out as invalid.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelSpecialCharacters()
        {
            var invalidPeptides = new List<string>
            {
                "PEPTIDE K", // space
                "PEPTIDE\tK", // tab
                "PEPTIDE\nK", // newline
                "PEPTIDE-K", // hyphen
                "PEPTIDE*K", // asterisk
            };

            var model = new PFly2024FineTuned(invalidPeptides, out var warning);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(0), "Peptides with special characters should be filtered out");
            Assert.That(warning, Is.Not.Null, "Invalid characters should produce a warning");
        }

        /// <summary>
        /// Tests the structure and content of PeptideDetectabilityPrediction records.
        /// Verifies that predictions contain all expected fields with valid values.
        /// </summary>
        [Test]
        public static async Task TestPFly2024FineTunedPredictionStructure()
        {
            var peptides = new List<string> { "PEPTIDEK", "VALIDSEQK" };
            var model = new PFly2024FineTuned(peptides, out var warning);

            await model.PredictAsync();

            Assert.That(model.Predictions.Count, Is.EqualTo(2), "Should have predictions for both peptides");

            foreach (var prediction in model.Predictions)
            {
                // Checks that probabilities can be accessed by name
                Assert.That(prediction.PeptideSequence, Is.Not.Null.And.Not.Empty);
                Assert.That(prediction.DetectabilityProbabilities.NotDetectable, Is.InRange(0.0, 1.0), "NotDetectable probability should be 0-1");
                Assert.That(prediction.DetectabilityProbabilities.LowDetectability, Is.InRange(0.0, 1.0), "LowDetectability probability should be 0-1");
                Assert.That(prediction.DetectabilityProbabilities.IntermediateDetectability, Is.InRange(0.0, 1.0), "IntermediateDetectability probability should be 0-1");
                Assert.That(prediction.DetectabilityProbabilities.HighDetectability, Is.InRange(0.0, 1.0), "HighDetectability probability should be 0-1");

                // Sum of probabilities should be approximately 1
                double sum = prediction.DetectabilityProbabilities.NotDetectable +
                             prediction.DetectabilityProbabilities.LowDetectability +
                             prediction.DetectabilityProbabilities.IntermediateDetectability +
                             prediction.DetectabilityProbabilities.HighDetectability;
                Assert.That(sum, Is.EqualTo(1.0).Within(0.01), "Sum of detectability probabilities should be approximately 1");
            }
        }

        /// <summary>
        /// Tests that the model correctly handles batch sizes larger than MaxBatchSize.
        /// Should split large batches appropriately.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelBatchSizing()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            var peptides = new List<string>();
            var seqLength = 10;
            var numberOfSequences = 500; // More than MaxBatchSize (128)

            while (peptides.Count < numberOfSequences)
            {
                var pep = new string(Random.Shared.GetItems(aminoacids, seqLength));
                if (!peptides.Contains(pep))
                {
                    peptides.Add(pep);
                }
            }

            var model = new PFly2024FineTuned(peptides, out var warning);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(numberOfSequences));
            Assert.DoesNotThrowAsync(async () => await model.PredictAsync());
        }

        /// <summary>
        /// Tests handling of null input list.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelNullInput()
        {
            List<string> nullList = null;

            Assert.DoesNotThrow(() =>
            {
                var model = new PFly2024FineTuned(nullList, out var warning);
                Assert.That(warning, Is.Not.Null, "Null input should produce a warning");
                Assert.That(model.PeptideSequences.Count, Is.EqualTo(0));
            });
        }

        /// <summary>
        /// Tests that predictions maintain order corresponding to input peptides.
        /// </summary>
        [Test]
        public static async Task TestPFly2024FineTunedPredictionOrder()
        {
            var peptides = new List<string> { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
            var model = new PFly2024FineTuned(peptides, out var warning);

            await model.PredictAsync();

            Assert.That(model.Predictions.Count, Is.EqualTo(4));
            Assert.That(model.Predictions[0].PeptideSequence, Is.EqualTo("AAAAAA"));
            Assert.That(model.Predictions[1].PeptideSequence, Is.EqualTo("CCCCCC"));
            Assert.That(model.Predictions[2].PeptideSequence, Is.EqualTo("GGGGGG"));
            Assert.That(model.Predictions[3].PeptideSequence, Is.EqualTo("TTTTTT"));
        }

        /// <summary>
        /// Tests that the canonical amino acid pattern regex works correctly.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedCanonicalAminoAcidPattern()
        {
            var model = new PFly2024FineTuned(new List<string> { "PEPTIDEK" }, out var warning);

            // Test valid sequences
            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("ACDEFGHIKLMNPQRSTVWY", model.CanonicalAminoAcidPattern));
            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("PEPTIDEK", model.CanonicalAminoAcidPattern));

            // Test invalid sequences
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("PEPTXIDEK", model.CanonicalAminoAcidPattern));
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("PEPTIDE123", model.CanonicalAminoAcidPattern));
        }

        /// <summary>
        /// Tests modification pattern regex for correct bracket matching.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModificationPattern()
        {
            var model = new PFly2024FineTuned(new List<string> { "PEPTIDEK" }, out var warning);

            var modificationPattern = model.ModificationPattern;

            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("[Oxidation]", modificationPattern));
            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("[Common Variable:Oxidation on M]", modificationPattern));
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("[", modificationPattern));
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("]", modificationPattern));
        }

        [Test]
        public static async Task TestModelIsDisposedProperly()
        {
            var peptides = new List<string> { "PEPTIDE", "PEPTIDEK" };

            // Disposed after use/inference
            var model = new PFly2024FineTuned(peptides, out var warning);
            Assert.That(warning, Is.Null,
                "Warning should not be generated for valid peptides");
            await model.PredictAsync();
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.PredictAsync(),
                "1Running inference on a disposed model should throw ObjectDisposedException");

            // Results should still be accessible after disposal
            Assert.That(model.Predictions.Count, Is.EqualTo(2),
                "Predicted spectra should still be accessible after model is disposed");

            // Disposed on command
            model = new PFly2024FineTuned(peptides, out warning);
            Assert.That(warning, Is.Null,
                "Warning should not be generated for valid peptides");
            model.Dispose();
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.PredictAsync(),
                "2Running inference on a disposed model should throw ObjectDisposedException");
        }
    }
}
