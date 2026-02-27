using NUnit.Framework;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;
using PredictionClients.Koina.Util;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class PFly2024FineTunedTests
    {
        /// <summary>
        /// Tests that the model correctly processes valid peptide sequences without throwing exceptions.
        /// Verifies basic prediction functionality and that all valid peptides return predictions.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelAcceptsValidPeptides()
        {
            var modelInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTIDEK"),
                new DetectabilityPredictionInput("ACDEFGHIKLMNPQRSTVWY"), // All standard amino acids
                new DetectabilityPredictionInput("LMNPQRSTVWY") // Another valid sequence
            };

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.DetectabilityProbabilities != null), Is.True, "All valid peptides should have predictions");
            Assert.That(predictions.All(p => p.Warning == null), Is.True, "Valid peptides should not produce warnings");
        }

        /// <summary>
        /// Tests that peptides with modifications are handled using UsePrimarySequence mode.
        /// PFly2024FineTuned does not support modified peptides for detectability prediction.
        /// ModHandlingMode.UsePrimarySequence strips all modifications and predicts on the base sequence.
        /// Modified peptides should still receive valid predictions (on their base sequence) with warnings.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelStripsModificationsAndPredictsBaseSequence()
        {
            var modelInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("M[Common Variable:Oxidation on M]PEPTIDEK"),
                new DetectabilityPredictionInput("PEPTIDEM[Common Variable:Oxidation on M]"),
                new DetectabilityPredictionInput("C[Common Fixed:Carbamidomethyl on C]PEPTIDEK"),
                new DetectabilityPredictionInput("PEPTIDEPEPI")
            };
   
            var model = new PFly2024FineTuned();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(4), "Should return predictions for all inputs");
            Assert.That(model.ValidInputsMask.Count, Is.EqualTo(4), "All inputs should be considered valid after stripping modifications");

            // Modified peptides should have valid predictions (on base sequence) with warnings
            Assert.That(predictions[0].DetectabilityProbabilities, Is.Not.Null, "Modified peptide should have predictions on base sequence");
            Assert.That(predictions[0].Warning, Is.Not.Null, "Modified peptide should have a warning");
            Assert.That(predictions[0].Warning.Message, Does.Contain("modifications"), "Warning should mention modifications");

            Assert.That(predictions[1].DetectabilityProbabilities, Is.Not.Null, "Modified peptide should have predictions on base sequence");
            Assert.That(predictions[1].Warning, Is.Not.Null, "Modified peptide should have a warning");

            Assert.That(predictions[2].DetectabilityProbabilities, Is.Not.Null, "Modified peptide should have predictions on base sequence");
            Assert.That(predictions[2].Warning, Is.Not.Null, "Modified peptide should have a warning");

            // Unmodified peptide should have no warnings
            Assert.That(predictions[3].DetectabilityProbabilities, Is.Not.Null, "Unmodified peptide should have valid predictions");
            Assert.That(predictions[3].Warning, Is.Null, "Unmodified peptide should not have warnings");
        }

        /// <summary>
        /// Tests boundary conditions for peptide sequence length.
        /// Valid lengths are 1-40 amino acids (after removing modification annotations).
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelPeptideLengthBoundaries()
        {
            var model = new PFly2024FineTuned();

            // Test minimum valid length (1 amino acid)
            var minLengthInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("K")
            };

            var predictionsMin = model.Predict(minLengthInputs);

            Assert.That(predictionsMin.Count, Is.EqualTo(1), "Should return a prediction entry");
            Assert.That(model.ValidInputsMask[0], Is.True, "Single amino acid should be considered valid");
            Assert.That(predictionsMin[0].DetectabilityProbabilities, Is.Not.Null, "Single amino acid should be valid");
            Assert.That(predictionsMin[0].Warning, Is.Null, "Minimum length peptide should produce a benign warning");

            // Test maximum valid length (40 amino acids)
            var maxLengthInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput(new string('A', 40))
            };

            var predictionsMax = model.Predict(maxLengthInputs);

            Assert.That(predictionsMax.Count, Is.EqualTo(1), "Should return a prediction entry");
            Assert.That(predictionsMax[0].DetectabilityProbabilities, Is.Not.Null, "Peptide with exactly 40 amino acids should be valid");
            Assert.That(predictionsMax[0].Warning, Is.Null, "Maximum length peptide should not produce a warning");

            // Test invalid length (0 amino acids)
            var tooShortInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("")
            };

            var predictionsTooShort = model.Predict(tooShortInputs);

            Assert.That(predictionsTooShort.Count, Is.EqualTo(1), "Should return a prediction entry for empty sequence");
            Assert.That(model.ValidInputsMask[0], Is.False, "Empty peptide should be considered invalid");
            Assert.That(predictionsTooShort[0].DetectabilityProbabilities, Is.Null, "Empty peptide should be invalid");
            Assert.That(predictionsTooShort[0].Warning, Is.Not.Null, "Empty peptide should produce a warning");

            // Test invalid length (41 amino acids)
            var tooLongInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput(new string('A', 41))
            };

            var predictionsTooLong = model.Predict(tooLongInputs);

            Assert.That(predictionsTooLong.Count, Is.EqualTo(1), "Should return a prediction entry for too long sequence");
            Assert.That(model.ValidInputsMask[0], Is.False, "Peptide with 41 amino acids should be considered invalid");
            Assert.That(predictionsTooLong[0].DetectabilityProbabilities, Is.Null, "Peptide with 41 amino acids should be invalid");
            Assert.That(predictionsTooLong[0].Warning, Is.Not.Null, "Too long peptide should produce a warning");
        }

        /// <summary>
        /// Tests that peptides containing non-standard or invalid amino acids are rejected.
        /// Standard 20 amino acids are: ACDEFGHIKLMNPQRSTVWY
        /// Characters like X (unknown), B (Asx), Z (Glx), U (Selenocysteine) are not supported.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelInvalidAminoAcids()
        {
            var model = new PFly2024FineTuned();

            // Test peptide with 'X' (unknown amino acid)
            var inputsWithX = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTXIDEK")
            };

            var predictionsX = model.Predict(inputsWithX);

            Assert.That(predictionsX.Count, Is.EqualTo(1), "Should return a prediction entry");
            Assert.That(predictionsX[0].DetectabilityProbabilities, Is.Null, "Peptide with 'X' should be filtered out");
            Assert.That(predictionsX[0].Warning, Is.Not.Null, "Invalid amino acid should produce a warning");

            // Test peptide with 'U' (Selenocysteine)
            var inputsWithU = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTUIDEK")
            };

            var predictionsU = model.Predict(inputsWithU);

            Assert.That(predictionsU.Count, Is.EqualTo(1), "Should return a prediction entry");
            Assert.That(predictionsU[0].DetectabilityProbabilities, Is.Null, "Peptide with 'U' should be filtered out");

            // Test peptide with 'B' (Asx)
            var inputsWithB = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTBIDEK")
            };

            var predictionsB = model.Predict(inputsWithB);

            Assert.That(predictionsB.Count, Is.EqualTo(1), "Should return a prediction entry");
            Assert.That(predictionsB[0].DetectabilityProbabilities, Is.Null, "Peptide with 'B' should be filtered out");

            // Test peptide with 'Z' (Glx)
            var inputsWithZ = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTZIDEK")
            };

            var predictionsZ = model.Predict(inputsWithZ);

            Assert.That(predictionsZ.Count, Is.EqualTo(1), "Should return a prediction entry");
            Assert.That(predictionsZ[0].DetectabilityProbabilities, Is.Null, "Peptide with 'Z' should be filtered out");

            // Test mixed valid and invalid peptides
            var mixedInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTIDEK"),
                new DetectabilityPredictionInput("PEPTXIDEK"),
                new DetectabilityPredictionInput("VALIDSEQ")
            };

            var predictionsMixed = model.Predict(mixedInputs);

            Assert.That(predictionsMixed.Count, Is.EqualTo(3), "Should return entries for all inputs");
            Assert.That(predictionsMixed.Count(p => p.DetectabilityProbabilities != null), Is.EqualTo(2), "Only valid peptides should have predictions");
        }

        /// <summary>
        /// Tests handling of empty input lists.
        /// Empty lists should not throw exceptions and should return empty predictions.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelEmptyInputHandling()
        {
            var emptyInputs = new List<DetectabilityPredictionInput>();

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(emptyInputs);

            Assert.That(predictions.Count, Is.EqualTo(0), "Empty input should result in no predictions");
            Assert.DoesNotThrow(() => model.Predict(emptyInputs), "Empty input should not throw exception");
        }

        /// <summary>
        /// Tests that duplicate peptide sequences are handled appropriately.
        /// Duplicates should both receive predictions without errors.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelDuplicateHandling()
        {
            var duplicateInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTIDEK"),
                new DetectabilityPredictionInput("PEPTIDEK")
            };

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(duplicateInputs);

            Assert.That(predictions.Count, Is.EqualTo(2), "Duplicate peptides should be handled without error");
            Assert.That(predictions.All(p => p.DetectabilityProbabilities != null), Is.True, "Both duplicates should have predictions");
        }

        /// <summary>
        /// Tests that model constants and properties are correctly initialized.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelProperties()
        {
            var model = new PFly2024FineTuned();

            Assert.That(model.ModelName, Is.EqualTo("pfly_2024_fine_tuned"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(128));
            Assert.That(model.NumberOfDetectabilityClasses, Is.EqualTo(4));
            Assert.That(model.DetectabilityClasses.Count, Is.EqualTo(4));
            Assert.That(model.DetectabilityClasses[0], Is.EqualTo("Not Detectable"));
            Assert.That(model.DetectabilityClasses[1], Is.EqualTo("Low Detectability"));
            Assert.That(model.DetectabilityClasses[2], Is.EqualTo("Intermediate Detectability"));
            Assert.That(model.DetectabilityClasses[3], Is.EqualTo("High Detectability"));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(40));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.ModHandlingMode, Is.EqualTo(IncompatibleModHandlingMode.UsePrimarySequence));
            Assert.That(model.ValidModificationUnimodMapping.Count, Is.EqualTo(0), "PFly does not support any modifications");
        }

        /// <summary>
        /// Tests handling of peptides containing lowercase letters.
        /// Tests whether the model/validation requires uppercase amino acids.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelCaseSensitivity()
        {
            var mixedCaseInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("peptidek"), // lowercase
                new DetectabilityPredictionInput("PEPTIDEK"), // uppercase
                new DetectabilityPredictionInput("PePtIdEk")  // mixed case
            };

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(mixedCaseInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), "Should return predictions for all inputs");
            // The AllowedAminoAcidPattern uses uppercase [ACDEFGHIKLMNPQRSTVWY], so lowercase may be invalid
            // If lowercase letters are rejected, predictions[0] and predictions[2] should have null probabilities
        }

        /// <summary>
        /// Tests handling of peptides with special characters and whitespace.
        /// These should be filtered out as invalid.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelSpecialCharacters()
        {
            var invalidInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTIDE K"), // space
                new DetectabilityPredictionInput("PEPTIDE\tK"), // tab
                new DetectabilityPredictionInput("PEPTIDE\nK"), // newline
                new DetectabilityPredictionInput("PEPTIDE-K"), // hyphen
                new DetectabilityPredictionInput("PEPTIDE*K"), // asterisk
            };

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(invalidInputs);

            Assert.That(predictions.Count, Is.EqualTo(5), "Should return entries for all inputs");
            Assert.That(predictions.All(p => p.DetectabilityProbabilities == null), Is.True, "Peptides with special characters should be filtered out");
            Assert.That(predictions.All(p => p.Warning != null), Is.True, "Invalid characters should produce warnings");
        }

        /// <summary>
        /// Tests the structure and content of PeptideDetectabilityPrediction records.
        /// Verifies that predictions contain all expected fields with valid values.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedPredictionStructure()
        {
            var modelInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTIDEK"),
                new DetectabilityPredictionInput("VALIDSEQK")
            };

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2), "Should have predictions for both peptides");

            foreach (var prediction in predictions.Where(p => p.DetectabilityProbabilities != null))
            {
                // Checks that probabilities can be accessed by name
                Assert.That(prediction.FullSequence, Is.Not.Null.And.Not.Empty);
                Assert.That(prediction.DetectabilityProbabilities.Value.NotDetectable, Is.InRange(0.0, 1.0), "NotDetectable probability should be 0-1");
                Assert.That(prediction.DetectabilityProbabilities.Value.LowDetectability, Is.InRange(0.0, 1.0), "LowDetectability probability should be 0-1");
                Assert.That(prediction.DetectabilityProbabilities.Value.IntermediateDetectability, Is.InRange(0.0, 1.0), "IntermediateDetectability probability should be 0-1");
                Assert.That(prediction.DetectabilityProbabilities.Value.HighDetectability, Is.InRange(0.0, 1.0), "HighDetectability probability should be 0-1");

                // Sum of probabilities should be approximately 1
                double sum = prediction.DetectabilityProbabilities.Value.NotDetectable +
                             prediction.DetectabilityProbabilities.Value.LowDetectability +
                             prediction.DetectabilityProbabilities.Value.IntermediateDetectability +
                             prediction.DetectabilityProbabilities.Value.HighDetectability;
                Assert.That(sum, Is.EqualTo(1.0).Within(0.01), "Sum of detectability probabilities should be approximately 1");
            }
        }

        /// <summary>
        /// Tests that the model correctly handles batch sizes larger than MaxBatchSize.
        /// Should split large batches appropriately and process all inputs.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelRequestBatching()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            var modelInputs = new List<DetectabilityPredictionInput>();
            var seqLength = 10;
            var numberOfSequences = 500; // More than MaxBatchSize (128)

            var peptides = new HashSet<string>();
            while (peptides.Count < numberOfSequences)
            {
                var pep = new string(Random.Shared.GetItems(aminoacids, seqLength));
                peptides.Add(pep);
            }

            foreach (var peptide in peptides)
            {
                modelInputs.Add(new DetectabilityPredictionInput(peptide));
            }

            var model = new PFly2024FineTuned();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));

            var predictions = model.Predict(modelInputs);
            Assert.That(predictions.Count, Is.EqualTo(numberOfSequences), "Should return predictions for all inputs");
        }

        /// <summary>
        /// Tests handling of null input list.
        /// Null list should be treated similar to empty list.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModelNullInput()
        {
            List<DetectabilityPredictionInput> nullList = null;

            var model = new PFly2024FineTuned();

            // This may throw or handle null gracefully depending on implementation
            // Test documents the expected behavior
            Assert.DoesNotThrow(() => model.Predict(nullList));
            Assert.That(model.Predictions.Count, Is.EqualTo(0), "Null input should result in no predictions");
            Assert.That(model.ValidInputsMask.Count, Is.EqualTo(0), "Null input should result in empty valid inputs mask");
        }

        /// <summary>
        /// Tests that predictions maintain order corresponding to input peptides.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedPredictionOrder()
        {
            var modelInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("AAAAAA"),
                new DetectabilityPredictionInput("CCCCCC"),
                new DetectabilityPredictionInput("GGGGGG"),
                new DetectabilityPredictionInput("TTTTTT")
            };

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(4));
            Assert.That(predictions[0].FullSequence, Is.EqualTo("AAAAAA"));
            Assert.That(predictions[1].FullSequence, Is.EqualTo("CCCCCC"));
            Assert.That(predictions[2].FullSequence, Is.EqualTo("GGGGGG"));
            Assert.That(predictions[3].FullSequence, Is.EqualTo("TTTTTT"));
        }

        /// <summary>
        /// Tests that the allowed amino acid pattern regex works correctly.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedAllowedAminoAcidPattern()
        {
            var model = new PFly2024FineTuned();

            // Test valid sequences
            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("ACDEFGHIKLMNPQRSTVWY", model.AllowedAminoAcidPattern));
            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("PEPTIDEK", model.AllowedAminoAcidPattern));

            // Test invalid sequences
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("PEPTXIDEK", model.AllowedAminoAcidPattern));
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("PEPTIDE123", model.AllowedAminoAcidPattern));
        }

        /// <summary>
        /// Tests modification pattern regex for correct bracket matching.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedModificationPattern()
        {
            var model = new PFly2024FineTuned();

            var modificationPattern = model.ModificationPattern;

            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("[Oxidation]", modificationPattern));
            Assert.That(System.Text.RegularExpressions.Regex.IsMatch("[Common Variable:Oxidation on M]", modificationPattern));
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("[", modificationPattern));
            Assert.That(!System.Text.RegularExpressions.Regex.IsMatch("]", modificationPattern));
        }

        /// <summary>
        /// Tests that predicted detectability values are chemically reasonable.
        /// Probability scores should sum to approximately 1.0 and be within valid ranges.
        /// </summary>
        [Test]
        public static void TestPFly2024FineTunedPredictionQuality()
        {
            var modelInputs = new List<DetectabilityPredictionInput>
            {
                new DetectabilityPredictionInput("PEPTIDEK"),
                new DetectabilityPredictionInput("ELVISLIVESK"),
                new DetectabilityPredictionInput("AAAAAAAAAA")
            };

            var model = new PFly2024FineTuned();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), "Should have predictions for all peptides");

            foreach (var prediction in predictions.Where(p => p.DetectabilityProbabilities != null))
            {
                var probs = prediction.DetectabilityProbabilities.Value;

                // Each probability should be in valid range
                Assert.That(probs.NotDetectable, Is.InRange(0.0, 1.0));
                Assert.That(probs.LowDetectability, Is.InRange(0.0, 1.0));
                Assert.That(probs.IntermediateDetectability, Is.InRange(0.0, 1.0));
                Assert.That(probs.HighDetectability, Is.InRange(0.0, 1.0));

                // Sum should be approximately 1.0
                double sum = probs.NotDetectable + probs.LowDetectability + 
                            probs.IntermediateDetectability + probs.HighDetectability;
                Assert.That(sum, Is.EqualTo(1.0).Within(0.01));
            }
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark test for the PFly2024FineTuned model with a large number of unique peptide sequences.
        /// This test is meant to evaluate the model's ability to handle large batch sizes and the efficiency of request batching logic
        /// 
        /// Notes: 
        ///  - 1 MaxBatchSize is 128 peptides and PFly2024FineTuned.Predict() takes about 0.75 seconds to run.
        ///  - With 4 million unique peptides of length 40, the test takes approximately 5 min minutes to run. 
        ///  - Timing is linearly proportional to peptide number due to the throttled approach.
        ///  - The default MaxNumberOfBatchesPerRequest is currently set to 500, which was tested up to 4 million peptides without hitting client issues. 
        /// </summary>
        public static void BenchmarkModelInputCountPerformance()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            var modelInputs = new List<DetectabilityPredictionInput>();
            var seqLength = 40; // max length increases combinatorial space to ensure we always get unique sequences and benchmark response length handling
            var numberOfSequences = 4000000;
            var peptides = new HashSet<string>();
            while (peptides.Count < numberOfSequences)
            {
                var pep = new string(Random.Shared.GetItems(aminoacids, seqLength));
                if (!peptides.Contains(pep))
                {
                    peptides.Add(pep);
                }
            }
            foreach (var peptide in peptides)
            {
                modelInputs.Add(new DetectabilityPredictionInput(peptide));
            }
            var model = new PFly2024FineTuned();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            var predictions = model.Predict(modelInputs);
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(numberOfSequences));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict {numberOfSequences:N0} peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }
    }
}
