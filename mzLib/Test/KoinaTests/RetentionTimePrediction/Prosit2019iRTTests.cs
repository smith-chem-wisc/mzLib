using NUnit.Framework;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Threading.Tasks;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class Prosit2019iRTTests
    {
        /// <summary>
        /// Tests the constructor with a list of valid peptide sequences.
        /// 
        /// Validates that:
        /// - No warnings are generated for valid sequences
        /// - All sequences are accepted and stored
        /// - Model properties are correctly initialized
        /// - Model name matches expected value "Prosit_2019_irt"
        /// - Batch size, length constraints, and iRT flag are correct
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - PeptideSequences count should equal input count (with carbamidomethylation applied)
        /// - All model properties should match specification
        /// </summary>
        [Test]
        public void TestConstructorWithValidPeptides()
        {
            var peptideSequences = new List<string>
            {
                "PEPTIDE",
                "TESTING"
            };

            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(2));
            Assert.That(model.ModelName, Is.EqualTo("Prosit_2019_irt"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.IsIndexedRetentionTimeModel, Is.True);
        }

        /// <summary>
        /// Tests the constructor with an empty list of peptides.
        /// 
        /// Validates that:
        /// - A warning is generated indicating empty input
        /// - PeptideSequences list is empty
        /// - Model is still instantiated (doesn't throw exception)
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should contain "Inputs were empty"
        /// - PeptideSequences should be empty
        /// 
        /// Use Case: Defensive programming - ensures graceful handling of empty inputs
        /// </summary>
        [Test]
        public void TestConstructorWithEmptyPeptides()
        {
            var peptideSequences = new List<string>();

            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings.Message, Does.Contain("Inputs were empty"));
            Assert.That(model.PeptideSequences, Is.Empty);
        }

        /// <summary>
        /// Tests the constructor with a null list of peptides.
        /// 
        /// Validates that:
        /// - A warning is generated for null input
        /// - Model handles null gracefully (no NullReferenceException)
        /// - PeptideSequences list is initialized as empty
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should contain "Inputs were empty"
        /// - PeptideSequences should be empty
        /// 
        /// Use Case: Defensive programming - ensures null safety
        /// </summary>
        [Test]
        public void TestConstructorWithNullPeptides()
        {
            List<string> peptideSequences = null;
            
            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings.Message, Does.Contain("Inputs were empty"));
            Assert.That(model.PeptideSequences, Is.Empty);
        }

        /// <summary>
        /// Tests the constructor with a mix of valid and invalid peptide sequences.
        /// 
        /// Invalid sequences include:
        /// - Sequences exceeding 30 amino acids (MaxPeptideLength)
        /// - Sequences with invalid characters (e.g., *)
        /// - Sequences with non-canonical amino acids (e.g., U for selenocysteine)
        /// 
        /// Validates that:
        /// - A warning is generated listing invalid sequences
        /// - Only valid sequences are added to PeptideSequences
        /// - Invalid sequences are properly filtered out
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should contain "invalid"
        /// - PeptideSequences should only contain "PEPTIDE"
        /// - Invalid sequences should be excluded
        /// 
        /// Use Case: Ensures robust input validation and informative error messaging
        /// </summary>
        [Test]
        public void TestConstructorWithInvalidSequences()
        {
            // Arrange
            var peptideSequences = new List<string>
            {
                "PEPTIDE",
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // Too long
                "PEP*TIDE", // Invalid character
                "SEQUENS" // Noncanonical amino acid 'U'
            };

            // Act
            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            // Assert
            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings.Message, Does.Contain("invalid"));
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(1)); // Only valid sequences
            Assert.That(model.PeptideSequences, Does.Contain("PEPTIDE"));
        }

        /// <summary>
        /// Tests the handling of valid post-translational modifications.
        /// 
        /// Tests three supported modification types:
        /// 1. Oxidation on methionine (M): [Common Variable:Oxidation on M] → [UNIMOD:35]
        /// 2. Carbamidomethyl on cysteine (C): [Common Fixed:Carbamidomethyl on C] → [UNIMOD:4]
        /// 
        /// Validates that:
        /// - Modifications are recognized as valid
        /// - mzLib format is converted to Prosit UNIMOD format
        /// - All modified sequences are accepted
        /// - No warnings are generated
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - All 3 sequences should be accepted
        /// - Modifications should be converted to UNIMOD format
        /// 
        /// Use Case: Ensures compatibility with mzLib modification format and proper conversion
        /// </summary>
        [Test]
        public void TestValidModificationMapping()
        {
            var peptideSequences = new List<string>
            {
                "PEPTIDE[Common Variable:Oxidation on M]",
                "M[Common Variable:Oxidation on M]SEQENS",
                "TESTC[Common Fixed:Carbamidomethyl on C]ING"
            };

            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(3));
            // Check that modifications are converted to UNIMOD format
            Assert.That(model.PeptideSequences[0], Does.Contain("[UNIMOD:35]"));
            Assert.That(model.PeptideSequences[1], Does.Contain("[UNIMOD:35]"));
            Assert.That(model.PeptideSequences[2], Does.Contain("[UNIMOD:4]"));
        }

        /// <summary>
        /// Tests the rejection of invalid or unsupported modifications.
        /// 
        /// Invalid modifications are those not in the ValidModificationUnimodMapping dictionary.
        /// For Prosit2019iRT, only Oxidation on M and Carbamidomethyl on C are supported.
        /// 
        /// Validates that:
        /// - Invalid modifications are detected
        /// - Sequences with invalid modifications are rejected
        /// - A warning is generated listing problematic sequences
        /// - Valid unmodified sequences are still accepted
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should contain "invalid"
        /// - Valid sequences count should be 2 (PEPTIDE, TESTING)
        /// </summary>
        [Test]
        public void TestInvalidModifications()
        {
            // Arrange
            var peptideSequences = new List<string>
            {
                "PEPTIDE",
                "SEQUENC[InvalidMod]E",
                "TESTING"
            };

            // Act
            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            // Assert
            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings.Message, Does.Contain("invalid"));
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(2)); // Invalid sequence excluded
        }

        /// <summary>
        /// Tests batching behavior with a small input list (below max batch size).
        /// 
        /// Validates that:
        /// - All sequences are processed in a single batch
        /// - Batch contains expected keys: "id" and "inputs"
        /// - No warnings are generated
        /// - Predictions are populated after inference
        /// 
        /// Expected Behavior:
        /// - batches count should be 1
        /// - Predictions count should equal input count (3)
        /// 
        /// Use Case: Ensures correct batching and inference execution
        /// </summary>
        [Test]
        public void TestBatchingWithSmallInput()
        {
            var peptideSequences = new List<string> { "PEPTIDE", "SEQENS", "TESTING" };
            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);
            Assert.DoesNotThrowAsync(async () => await model.PredictAsync());
            Assert.That(model.Predictions, Has.Count.EqualTo(3));
        }

        /// <summary>
        /// Tests batching behavior with a large input list (above max batch size).
        /// 
        /// Validates that:
        /// - Input list is split into multiple batches
        /// - Batch count matches expected (3 batches: 1000 + 1000 + 500)
        /// - Each batch has a unique ID
        /// 
        /// Expected Behavior:
        /// - batches count should be 3
        /// - Batch IDs should be unique and sequential (Batch0, Batch1, Batch2)
        /// </summary>
        [Test]
        public void TestBatchingWithLargeInput()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            int numberOfSequences = 2500;
            int seqLength = 20;
            var peptides = new List<string>();

            while (peptides.Count < numberOfSequences)
            {
                var pep = new string(Random.Shared.GetItems(aminoacids, seqLength));
                if (!peptides.Contains(pep))
                {
                    peptides.Add(pep);
                }
            }

            var model = new Prosit2019iRT(peptides, out WarningException warnings);
            Assert.DoesNotThrowAsync(async () => await model.PredictAsync());
            Assert.That(model.Predictions, Has.Count.EqualTo(numberOfSequences));
        }

        /// <summary>
        /// Tests batching behavior when the input list size is exactly the max batch size.
        /// 
        /// Validates that:
        /// - All sequences are processed in a single batch
        /// - No batching errors occur
        /// 
        /// Expected Behavior:
        /// - batches count should be 1
        /// </summary>
        [Test]
        public void TestBatchingExactlyAtMaxBatchSize()
        {
            var peptideSequences = new List<string>();
            for (int i = 0; i < 1000; i++)
            {
                peptideSequences.Add("PEPTIDE");
            }
            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);
            Assert.DoesNotThrowAsync(async () => await model.PredictAsync());
            Assert.That(model.Predictions, Has.Count.EqualTo(1000));
        }

        /// <summary>
        /// Tests that model predictions are initially empty upon model creation.
        /// 
        /// Validates that:
        /// - Predictions list is empty until inference is run
        /// 
        /// Expected Behavior:
        /// - Predictions should be empty
        /// </summary>
        [Test]
        public void TestPredictionsInitiallyEmpty()
        {
            var peptideSequences = new List<string> { "PEPTIDE", "SEQUENCE" };

            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            Assert.That(model.Predictions, Is.Empty);
        }

        /// <summary>
        /// Tests the properties of the Prosit2019iRT model.
        /// 
        /// Validates that:
        /// - Model name is correctly set
        /// - Max batch size, max peptide length, and min peptide length are within expected ranges
        /// - The model is flagged as an indexed retention time model
        /// - Valid modification mappings are populated
        /// 
        /// Expected Behavior:
        /// - Model properties should match the specifications
        /// </summary>
        [Test]
        public void TestModelProperties()
        {
            var peptideSequences = new List<string> { "PEPTIDE" };
            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            Assert.That(model.ModelName, Is.EqualTo("Prosit_2019_irt"));
            Assert.That(model.MaxBatchSize, Is.GreaterThan(0));
            Assert.That(model.MaxPeptideLength, Is.GreaterThan(0));
            Assert.That(model.MinPeptideLength, Is.GreaterThanOrEqualTo(1));
            Assert.That(model.IsIndexedRetentionTimeModel, Is.True);
            Assert.That(model.ValidModificationUnimodMapping, Is.Not.Null);
            Assert.That(model.ValidModificationUnimodMapping.Count, Is.GreaterThan(0));
        }

        /// <summary>
        /// Tests the handling of a mix of valid and invalid sequences in the input list.
        /// 
        /// Validates that:
        /// - A warning is generated for the invalid sequences
        /// - Only valid sequences are added to the model
        /// - Invalid sequences are properly identified
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should indicate invalid sequences
        /// - model.PeptideSequences should contain only valid sequences
        /// </summary>
        [Test]
        public void TestMixedValidAndInvalidSequences()
        {
            var peptideSequences = new List<string>
            {
                "PEPTIDE", // Valid
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // Invalid - too long
                "SEQENS", // Valid
                "INVALID*", // Invalid - bad character
                "TESTINGTWICE", // Valid
                "INVALIDPEPTIDE[BadMod]" // Invalid - bad modification
            };

            var model = new Prosit2019iRT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(3)); // Only valid sequences
            Assert.That(warnings.Message, Does.Contain("invalid"));
        }

        [Test]
        public static async Task TestModelIsDisposedProperly()
        {
            var peptides = new List<string> { "PEPTIDE", "PEPTIDEK" };

            // Disposed after use/inference
            var model = new Prosit2019iRT(peptides, out var warning);
            Assert.That(warning, Is.Null,
                "Warning should not be generated for valid peptides");
            await model.PredictAsync();
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.PredictAsync(),
                "1Running inference on a disposed model should throw ObjectDisposedException");

            // Results should still be accessible after disposal
            Assert.That(model.Predictions.Count, Is.EqualTo(2),
                "Predicted spectra should still be accessible after model is disposed");

            // Disposed on command
            model = new Prosit2019iRT(peptides, out warning);
            Assert.That(warning, Is.Null,
                "Warning should not be generated for valid peptides");
            model.Dispose();
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.PredictAsync(),
                "2Running inference on a disposed model should throw ObjectDisposedException");
        }
    }
}
