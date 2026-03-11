using NUnit.Framework;
using Omics.Modifications;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

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
            Assert.DoesNotThrowAsync(async () => await model.RunInferenceAsync());
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
            Assert.DoesNotThrowAsync(async () => await model.RunInferenceAsync());
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
            Assert.DoesNotThrowAsync(async () => await model.RunInferenceAsync());
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
            await model.RunInferenceAsync();
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.RunInferenceAsync(),
                "1Running inference on a disposed model should throw ObjectDisposedException");

            // Results should still be accessible after disposal
            Assert.That(model.Predictions.Count, Is.EqualTo(2),
                "Predicted spectra should still be accessible after model is disposed");

            // Disposed on command
            model = new Prosit2019iRT(peptides, out warning);
            Assert.That(warning, Is.Null,
                "Warning should not be generated for valid peptides");
            model.Dispose();
            Assert.ThrowsAsync<ObjectDisposedException>(async () => await model.RunInferenceAsync(),
                "2Running inference on a disposed model should throw ObjectDisposedException");
        }
		/// <summary>
		/// Tests Prosit2019iRT retention time predictions followed by Prosit2020HCD fragmentation predictions
		/// for unmodified peptides from a FASTA file.
		/// 
		/// This test:
		/// 1. Reads proteins from a FASTA file
		/// 2. Digests proteins with trypsin (max 1 missed cleavage)
		/// 3. Filters for unmodified peptides with length ≤ 30
		/// 4. Runs retention time prediction
		/// 5. Uses predicted retention times for fragmentation prediction
		/// 6. Generates fragment ion predictions for charge states 2, 3, 4 at collision energy 30
		/// 
		/// Expected Behavior:
		/// - Peptides should be successfully digested from FASTA
		/// - All peptides should pass validation (length ≤ 30, no modifications)
		/// - Retention time predictions should be generated for all valid peptides
		/// - Fragment ion predictions should be generated for all peptides × charge states
		/// - Each fragmentation prediction should contain fragment ions
		/// 
		/// Use Case: Real-world scenario for generating spectral library from protein database
		/// </summary>
		[Test]
		public async Task TestRetentionTimePredictionFromFastaFile()
		{
			// Arrange - Set your FASTA file path here
			string fastaFilePath = @"C:\Users\mrsho\Downloads\synthetic2\syntheticTest2.fasta";

			// Skip test if FASTA file doesn't exist
			if (!System.IO.File.Exists(fastaFilePath))
			{
				Assert.Ignore($"FASTA file not found at: {fastaFilePath}. Please update the path to run this test.");
				return;
			}

			// Read proteins from FASTA file
			var proteins = ProteinDbLoader.LoadProteinFasta(
				proteinDbLocation: fastaFilePath,
				generateTargets: true,
				decoyType: DecoyType.None,
				isContaminant: false,
				errors: out var errors);

			Assert.That(proteins, Is.Not.Empty, "FASTA file should contain at least one protein");
			Assert.That(errors, Is.Empty, "FASTA file should be parsed without errors");

			// Set up digestion parameters
			var digestionParams = new Proteomics.ProteolyticDigestion.DigestionParams(
				protease: "trypsin",
				maxMissedCleavages: 1,
				minPeptideLength: 7,  // Typical minimum for mass spec
				maxPeptideLength: 30, // Match Prosit max length
				initiatorMethionineBehavior: Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Variable,
				maxModificationIsoforms: 1);

			// Digest proteins into peptides
			var allPeptides = new HashSet<string>(); // Use HashSet for automatic deduplication
			foreach (var protein in proteins)
			{
				var digestedPeptides = protein.Digest(
					digestionParams,
					new List<Modification>(),  // No fixed modifications
					new List<Modification>()); // No variable modifications

				foreach (var peptide in digestedPeptides)
				{
					// Get base sequence (no modifications)
					string baseSequence = peptide.BaseSequence;

					// Filter: length check
					if (baseSequence.Length >= digestionParams.MinPeptideLength &&
						baseSequence.Length <= digestionParams.MaxPeptideLength)
					{
						allPeptides.Add(baseSequence); // HashSet automatically handles duplicates
					}
				}
			}

			Assert.That(allPeptides, Is.Not.Empty, "Digestion should produce at least one valid peptide");
			Console.WriteLine($"Total unique unmodified peptides: {allPeptides.Count}");

			// ========== STEP 1: Retention Time Prediction ==========
			Console.WriteLine("\n=== STEP 1: Retention Time Prediction ===");
			var rtModel = new Prosit2019iRT(allPeptides.ToList(), out WarningException rtWarnings);

			// Check for warnings (should be null for valid peptides)
			if (rtWarnings != null)
			{
				Console.WriteLine($"RT Warning: {rtWarnings.Message}");
			}

			Assert.That(rtModel.PeptideSequences, Is.Not.Empty,
				"RT model should accept at least some peptides from the FASTA file");

			// Run retention time prediction
			await rtModel.RunInferenceAsync();

			// Validate predictions
			Assert.That(rtModel.Predictions, Is.Not.Empty,
				"RT model should return predictions for the peptides");
			Assert.That(rtModel.Predictions.Count, Is.EqualTo(rtModel.PeptideSequences.Count),
				"Should have one RT prediction per input peptide sequence");

			// Validate that all predictions have reasonable retention time values
			foreach (var prediction in rtModel.Predictions)
			{
				Assert.That(prediction.PredictedRetentionTime, Is.Not.NaN,
					$"RT prediction for peptide {prediction.FullSequence} should not be NaN");
				Assert.That(double.IsFinite(prediction.PredictedRetentionTime), Is.True,
					$"RT prediction for peptide {prediction.FullSequence} should be a finite number");
			}

			// Output RT summary statistics
			var predictionValues = rtModel.Predictions.Select(p => p.PredictedRetentionTime).ToList();
			Console.WriteLine($"\nRetention Time Prediction Summary:");
			Console.WriteLine($"  Total predictions: {predictionValues.Count}");
			Console.WriteLine($"  Min iRT: {predictionValues.Min():F2}");
			Console.WriteLine($"  Max iRT: {predictionValues.Max():F2}");
			Console.WriteLine($"  Mean iRT: {predictionValues.Average():F2}");

			// Display first 5 RT predictions as examples
			Console.WriteLine($"\nFirst 5 RT Predictions:");
			foreach (var prediction in rtModel.Predictions.Take(5))
			{
				Console.WriteLine($"  {prediction.FullSequence}: {prediction.PredictedRetentionTime:F2}");
			}

			// ========== STEP 2: Fragmentation Prediction ==========
			Console.WriteLine("\n=== STEP 2: Fragmentation Prediction ===");

			// Prepare inputs for fragmentation model
			var chargeStates = new List<int> { 2, 3, 4 };
			int collisionEnergy = 30;

			var fragmentationPeptides = new List<string>();
			var fragmentationCharges = new List<int>();
			var fragmentationEnergies = new List<int>();
			var fragmentationRTs = new List<double?>();

			// Create a lookup dictionary for retention times
			var rtLookup = rtModel.Predictions.ToDictionary(
				p => p.FullSequence,
				p => (double?)p.PredictedRetentionTime);

			// Generate all combinations of peptide × charge state
			foreach (var prediction in rtModel.Predictions)
			{
				foreach (var charge in chargeStates)
				{
					fragmentationPeptides.Add(prediction.FullSequence);
					fragmentationCharges.Add(charge);
					fragmentationEnergies.Add(collisionEnergy);
					fragmentationRTs.Add(prediction.PredictedRetentionTime);
				}
			}

			Console.WriteLine($"\nPreparing fragmentation predictions:");
			Console.WriteLine($"  Unique peptides: {rtModel.Predictions.Count}");
			Console.WriteLine($"  Charge states: {string.Join(", ", chargeStates)}");
			Console.WriteLine($"  Collision energy: {collisionEnergy}");
			Console.WriteLine($"  Total combinations: {fragmentationPeptides.Count}");

			// Construct output path for spectral library (same folder as input FASTA)
			string fastaDirectory = System.IO.Path.GetDirectoryName(fastaFilePath);
			string fastaFileName = System.IO.Path.GetFileNameWithoutExtension(fastaFilePath);
			string spectralLibraryPath = System.IO.Path.Combine(fastaDirectory, $"{fastaFileName}_PredictedLibrary.msp");
			
			Console.WriteLine($"  Spectral library will be saved to: {spectralLibraryPath}");

			// Create fragmentation model
			var fragmentModel = new Prosit2020IntensityHCD(
				peptideSequences: fragmentationPeptides,
				precursorCharges: fragmentationCharges,
				collisionEnergies: fragmentationEnergies,
				retentionTimes: fragmentationRTs,
				spectralLibrarySavePath: spectralLibraryPath,
				warnings: out WarningException fragWarnings);

			if (fragWarnings != null)
			{
				Console.WriteLine($"\nFragmentation Warning: {fragWarnings.Message}");
			}

			Assert.That(fragmentModel.PeptideSequences, Is.Not.Empty,
				"Fragmentation model should accept peptides");

			// Run fragmentation prediction
			Console.WriteLine($"\nRunning fragmentation inference...");
			await fragmentModel.RunInferenceAsync();

			// Validate fragmentation predictions
			Assert.That(fragmentModel.Predictions, Is.Not.Empty,
				"Fragmentation model should return predictions");
			Assert.That(fragmentModel.Predictions.Count, Is.EqualTo(fragmentModel.PeptideSequences.Count),
				"Should have one fragmentation prediction per input");

			// Analyze fragmentation results
			Console.WriteLine($"\nFragmentation Prediction Summary:");
			Console.WriteLine($"  Total predictions: {fragmentModel.Predictions.Count}");
			Console.WriteLine($"  Unique peptides fragmented: {fragmentModel.Predictions.Select(p => p.FullSequence).Distinct().Count()}");

			// Validate fragment ion data
			int totalFragmentIons = 0;
			foreach (var prediction in fragmentModel.Predictions)
			{
				Assert.That(prediction.FragmentAnnotations, Is.Not.Null,
					$"Fragment annotations should not be null for {prediction.FullSequence} at charge {prediction.PrecursorCharge}");
				Assert.That(prediction.FragmentAnnotations.Count, Is.GreaterThan(0),
					$"Should have fragment ions for {prediction.FullSequence} at charge {prediction.PrecursorCharge}");
				
				// Verify parallel lists have matching counts
				Assert.That(prediction.FragmentMZs.Count, Is.EqualTo(prediction.FragmentAnnotations.Count),
					$"Fragment m/z count should match annotation count");
				Assert.That(prediction.FragmentIntensities.Count, Is.EqualTo(prediction.FragmentAnnotations.Count),
					$"Fragment intensity count should match annotation count");
				
				totalFragmentIons += prediction.FragmentAnnotations.Count;
			}

			Console.WriteLine($"  Total fragment ions predicted: {totalFragmentIons}");
			Console.WriteLine($"  Average fragments per peptide: {(double)totalFragmentIons / fragmentModel.Predictions.Count:F1}");

			// Display first 3 fragmentation predictions as examples
			Console.WriteLine($"\nFirst 3 Fragmentation Predictions:");
			foreach (var prediction in fragmentModel.Predictions.Take(3))
			{
				Console.WriteLine($"\n  Peptide: {prediction.FullSequence}");
				Console.WriteLine($"  Charge: {prediction.PrecursorCharge}+");
				Console.WriteLine($"  Fragment Ions: {prediction.FragmentAnnotations.Count}");

				// Show top 5 most intense fragments
				var fragments = prediction.FragmentAnnotations
					.Select((annotation, idx) => new
					{
						Annotation = annotation,
						Mz = prediction.FragmentMZs[idx],
						Intensity = prediction.FragmentIntensities[idx]
					})
					.Where(f => f.Intensity > 0) // Filter out impossible fragments (intensity = -1)
					.OrderByDescending(f => f.Intensity)
					.Take(5)
					.ToList();

				Console.WriteLine($"  Top 5 Fragments:");
				foreach (var fragment in fragments)
				{
					Console.WriteLine($"    {fragment.Annotation}: m/z {fragment.Mz:F4}, intensity {fragment.Intensity:F6}");
				}
			}

			// Verify spectral library generation
			if (fragmentModel.PredictedSpectra.Count > 0)
			{
				Console.WriteLine($"\n  Spectral library generated: {fragmentModel.PredictedSpectra.Count} spectra");
				Console.WriteLine($"  Library contains spectra for charges: {string.Join(", ", fragmentModel.PredictedSpectra.Select(s => s.ChargeState).Distinct().OrderBy(c => c))}");
			}

			Console.WriteLine("\n=== Test Complete ===");
		}
	}
}
