using NUnit.Framework;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class Prosit2019iRTTests
    {
        /// <summary>
        /// Tests that the model correctly processes valid peptide sequences without throwing exceptions.
        /// Verifies basic prediction functionality and model properties.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelAcceptsValidPeptides()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("PEPTIDE"),
                new RetentionTimePredictionInput("TESTING")
            };

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All valid peptides should have predictions");
            Assert.That(predictions.All(p => p.Warning == null), Is.True, "Valid peptides should not produce warnings");

            // Verify model properties
            Assert.That(model.ModelName, Is.EqualTo("Prosit_2019_irt"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.IsIndexedRetentionTimeModel, Is.True);
        }

        /// <summary>
        /// Tests handling of empty input lists.
        /// Empty lists should not throw exceptions and should return empty predictions.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelEmptyInputHandling()
        {
            var emptyInputs = new List<RetentionTimePredictionInput>();

            var model = new Prosit2019iRT();
            var predictions = model.Predict(emptyInputs);

            Assert.That(predictions.Count, Is.EqualTo(0), "Empty input should result in no predictions");
            Assert.DoesNotThrow(() => model.Predict(emptyInputs), "Empty input should not throw exception");
        }

        /// <summary>
        /// Tests handling of null input list.
        /// Null list should throw ArgumentNullException.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelNullInput()
        {
            List<RetentionTimePredictionInput> nullList = null;

            var model = new Prosit2019iRT();

            Assert.DoesNotThrow(() => model.Predict(nullList));
            Assert.That(model.Predictions.Count, Is.EqualTo(0), "Null input should be treated as empty and return no predictions");
            Assert.That(model.ValidInputsMask.Count, Is.EqualTo(0), "Null input should not produce warnings");
        }

        /// <summary>
        /// Tests handling of peptides with invalid sequences.
        /// Invalid sequences include too long peptides, invalid characters, and non-canonical amino acids.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelInvalidSequences()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("PEPTIDE"),
                new RetentionTimePredictionInput("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), // Too long
                new RetentionTimePredictionInput("PEP*TIDE"), // Invalid character
                new RetentionTimePredictionInput("SEQUENS") // Noncanonical amino acid 'U'
            };

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(4), "Should return entries for all inputs");
            Assert.That(predictions[0].PredictedRetentionTime, Is.Not.Null, "Valid peptide should have predictions");
            Assert.That(predictions[0].Warning, Is.Null, "Valid peptide should not have warnings");

            Assert.That(predictions[1].PredictedRetentionTime, Is.Null, "Too long peptide should have null predictions");
            Assert.That(predictions[1].Warning, Is.Not.Null, "Too long peptide should have warning");

            Assert.That(predictions[2].PredictedRetentionTime, Is.Null, "Invalid character peptide should have null predictions");
            Assert.That(predictions[2].Warning, Is.Not.Null, "Invalid character should have warning");

            Assert.That(predictions[3].PredictedRetentionTime, Is.Null, "Non-canonical AA peptide should have null predictions");
            Assert.That(predictions[3].Warning, Is.Not.Null, "Non-canonical AA should have warning");
        }

        /// <summary>
        /// Tests the handling of valid post-translational modifications.
        /// Prosit2019iRT supports oxidation on M and carbamidomethyl on C.
        /// ModHandlingMode is RemoveIncompatibleMods by default.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelValidModifications()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("M[Common Variable:Oxidation on M]SEQENS"),
                new RetentionTimePredictionInput("TESTC[Common Fixed:Carbamidomethyl on C]ING")
            };

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2), "Should return predictions for all inputs");
            Assert.That(predictions[0].PredictedRetentionTime, Is.Not.Null, "Peptide with valid oxidation should have predictions");
            Assert.That(predictions[1].PredictedRetentionTime, Is.Not.Null, "Peptide with valid carbamidomethyl should have predictions");
        }

        /// <summary>
        /// Tests handling of unsupported modifications.
        /// ModHandlingMode.RemoveIncompatibleMods strips unsupported modifications and predicts with remaining valid mods/sequence.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelIncompatibleModificationsOrAminoAcids()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("PEPTIDE"),
                new RetentionTimePredictionInput("SEQENC[InvalidMod]E"),
                new RetentionTimePredictionInput("TESTING"),
                new RetentionTimePredictionInput("SEQUENS") // Non-canonical amino acid 'U'
            };

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(4), "Should return entries for all inputs");
            Assert.That(predictions[0].PredictedRetentionTime, Is.Not.Null, "Valid unmodified peptide should have predictions");
            Assert.That(predictions[0].Warning, Is.Null, "Valid unmodified peptide should not have warnings");

            // With RemoveIncompatibleMods, the invalid mod is stripped and prediction is made on the base sequence
            Assert.That(predictions[1].PredictedRetentionTime, Is.Not.Null, "Peptide with invalid mod should still have predictions after stripping");
            Assert.That(predictions[1].Warning, Is.Not.Null, "Peptide with stripped mod should have warning");

            Assert.That(predictions[2].PredictedRetentionTime, Is.Not.Null, "Valid unmodified peptide should have predictions");

            Assert.That(predictions[3].PredictedRetentionTime, Is.Null, "Peptide with non-canonical amino acid should have null predictions");
        }

        /// <summary>
        /// Tests batching behavior with a small input list (below max batch size).
        /// All sequences should be processed successfully.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelSmallBatchProcessing()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("PEPTIDE"),
                new RetentionTimePredictionInput("SEQENS"),
                new RetentionTimePredictionInput("TESTING")
            };

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All valid peptides should have predictions");
        }

        /// <summary>
        /// Tests batching behavior with a large input list (above max batch size).
        /// Should split into multiple batches and process all sequences successfully.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelRequestBatching()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            int numberOfSequences = 2500;
            int seqLength = 20;
            var peptides = new HashSet<string>();

            while (peptides.Count < numberOfSequences)
            {
                var pep = new string(Random.Shared.GetItems(aminoacids, seqLength));
                peptides.Add(pep);
            }

            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2019iRT();

            Assert.DoesNotThrow(() => model.Predict(modelInputs));

            var predictions = model.Predict(modelInputs);
            Assert.That(predictions.Count, Is.EqualTo(numberOfSequences), "Should return predictions for all inputs");
        }

        /// <summary>
        /// Tests batching behavior when input size equals max batch size.
        /// Should process all sequences in a single batch.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelExactlyAtMaxBatchSize()
        {
            var modelInputs = new List<RetentionTimePredictionInput>();
            for (int i = 0; i < 1000; i++)
            {
                modelInputs.Add(new RetentionTimePredictionInput("PEPTIDE"));
            }

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1000), "Should return predictions for all inputs");
        }

        /// <summary>
        /// Tests the properties of the Prosit2019iRT model.
        /// Validates model metadata and configuration.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelProperties()
        {
            var model = new Prosit2019iRT();

            Assert.That(model.ModelName, Is.EqualTo("Prosit_2019_irt"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.IsIndexedRetentionTimeModel, Is.True);
            Assert.That(model.ModHandlingMode, Is.EqualTo(IncompatibleModHandlingMode.RemoveIncompatibleMods));
            Assert.That(model.ValidModificationUnimodMapping, Is.Not.Null);
            Assert.That(model.ValidModificationUnimodMapping.Count, Is.EqualTo(2), "Should support Oxidation and Carbamidomethyl");
            Assert.That(model.ValidModificationUnimodMapping.ContainsKey("[Common Variable:Oxidation on M]"), Is.True);
            Assert.That(model.ValidModificationUnimodMapping.ContainsKey("[Common Fixed:Carbamidomethyl on C]"), Is.True);
        }

        /// <summary>
        /// Tests handling of mixed valid and invalid sequences.
        /// Invalid sequences should have null predictions with warnings.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelMixedValidAndInvalidSequences()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("PEPTIDE"), // Valid
                new RetentionTimePredictionInput("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), // Invalid - too long
                new RetentionTimePredictionInput("SEQENS"), // Valid
                new RetentionTimePredictionInput("INVALID*"), // Invalid - bad character
                new RetentionTimePredictionInput("TESTINGTWICE") // Valid
            };

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(5), "Should return entries for all inputs");
            Assert.That(predictions.Count(p => p.PredictedRetentionTime != null), Is.EqualTo(3), "Only valid sequences should have predictions");
            Assert.That(predictions.Count(p => p.Warning != null), Is.EqualTo(2), "Invalid sequences should have warnings");
        }

        /// <summary>
        /// Tests prediction order is maintained from input to output.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelPredictionOrder()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("AAAAAA"),
                new RetentionTimePredictionInput("CCCCCC"),
                new RetentionTimePredictionInput("GGGGGG")
            };

            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3));
            Assert.That(predictions[0].FullSequence, Is.EqualTo("AAAAAA"));
            Assert.That(predictions[1].FullSequence, Is.EqualTo("CCCCCC"));
            Assert.That(predictions[2].FullSequence, Is.EqualTo("GGGGGG"));
        }

        /// <summary>
        /// Tests custom constructor parameters for batching and throttling.
        /// </summary>
        [Test]
        public void TestProsit2019iRTModelCustomConstructorParameters()
        {
            var model = new Prosit2019iRT(
                modHandlingMode: IncompatibleModHandlingMode.RemoveIncompatibleMods,
                maxNumberOfBatchesPerRequest: 100,
                throttlingDelayInMilliseconds: 50
            );

            Assert.That(model.MaxNumberOfBatchesPerRequest, Is.EqualTo(100));
            Assert.That(model.ThrottlingDelayInMilliseconds, Is.EqualTo(50));
            Assert.That(model.ModHandlingMode, Is.EqualTo(IncompatibleModHandlingMode.RemoveIncompatibleMods));
        }
    }
}
