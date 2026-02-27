using Easy.Common.Extensions;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.Koina.Util;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class Prosit2020iRTTMTTests
    {
        [Test]
        public void TestProsit2020iRTTMTPrediction()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]SEQENS"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTING")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3));
            Assert.That(model.ValidInputsMask, Is.All.True, "All inputs should be valid with N-terminal modifications");
            Assert.That(predictions.All(p => p.IsIndexed == true));
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null && p.PredictedRetentionTime.Value.IsFinite()));
            Assert.That(modelInputs.All(p => p.SequenceWarning == null), "Valid inputs should not have warnings");
            Assert.That(predictions.All(p => p.Warning == null), "Predictions from perfectly valid inputs should not have warnings");
            for (int i = 0; i < modelInputs.Count; i++)
            {
                Assert.That(predictions[i].FullSequence, Is.EqualTo(modelInputs[i].FullSequence), $"Full sequence in prediction should match input for index {i}");
                Assert.That(predictions[i].FullSequence, Is.EqualTo(modelInputs[i].ValidatedFullSequence), $"Validated full sequence in prediction should match input for index {i}");
            }
        }

        /// <summary>
        /// Tests that the model correctly processes valid peptide sequences with required N-terminal modifications.
        /// All sequences MUST have an N-terminal TMT or iTRAQ modification for this model.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelAcceptsValidPeptides()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]TESTING")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(2), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All valid peptides should have predictions");
            Assert.That(predictions.All(p => p.Warning == null), Is.True, "Valid peptides should not produce warnings");

            // Verify model properties
            Assert.That(model.ModelName, Is.EqualTo("Prosit_2020_irt_TMT"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.IsIndexedRetentionTimeModel, Is.True);

            // Make sure ThrowException mode also accepts valid peptides without throwing
            model = new Prosit2020iRTTMT(modHandlingMode: IncompatibleModHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.Predict(modelInputs), "Model should not throw exception for valid peptides with required N-terminal modifications when ModHandlingMode is set to ThrowException");
        }

        /// <summary>
        /// Tests that sequences WITHOUT N-terminal modifications are REJECTED.
        /// ModHandlingMode.ReturnNull rejects sequences with invalid/missing modifications.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelRejectsSequencesWithoutNTerminalModification()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"), // Valid - has N-term
                new RetentionTimePredictionInput("PEPTIDE"), // Invalid - no N-term mod
                new RetentionTimePredictionInput("PEPTK[Common Fixed:TMT6plex on K]IDE"), // Invalid - only K-TMT, no N-term
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]TESTING"), // Valid - has N-term
                new RetentionTimePredictionInput("M[Common Variable:Oxidation on M]SEQUENCE") // Invalid - only oxidation, no N-term
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(5), "Should return entries for all inputs");

            // First peptide: valid with N-term
            Assert.That(predictions[0].PredictedRetentionTime, Is.Not.Null, "Valid peptide with N-term should have predictions");
            Assert.That(predictions[0].Warning, Is.Null, "Valid peptide should not have warnings");

            // Second peptide: no N-term mod
            Assert.That(predictions[1].PredictedRetentionTime, Is.Null, "Peptide without N-term should be rejected");
            Assert.That(predictions[1].Warning, Is.Not.Null, "Peptide without N-term should have warning");

            // Third peptide: K-TMT but no N-term
            Assert.That(predictions[2].PredictedRetentionTime, Is.Null, "Peptide with only K-TMT should be rejected");
            Assert.That(predictions[2].Warning, Is.Not.Null, "Peptide without N-term should have warning");

            // Fourth peptide: valid with N-term
            Assert.That(predictions[3].PredictedRetentionTime, Is.Not.Null, "Valid peptide with N-term should have predictions");
            Assert.That(predictions[3].Warning, Is.Null, "Valid peptide should not have warnings");

            // Fifth peptide: oxidation but no N-term
            Assert.That(predictions[4].PredictedRetentionTime, Is.Null, "Peptide without N-term should be rejected");
            Assert.That(predictions[4].Warning, Is.Not.Null, "Peptide without N-term should have warning");

            model = new Prosit2020iRTTMT(modHandlingMode: IncompatibleModHandlingMode.ThrowException);
            Assert.Throws<ArgumentException>(() => model.Predict(modelInputs), "Model should throw exception for peptides without required N-terminal modifications when ModHandlingMode is set to ThrowException");
        }

        /// <summary>
        /// Tests handling of empty input lists.
        /// Empty lists should not throw exceptions and should return empty predictions.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelEmptyInputHandling()
        {
            var emptyInputs = new List<RetentionTimePredictionInput>();

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(emptyInputs);

            Assert.That(predictions.Count, Is.EqualTo(0), "Empty input should result in no predictions");
            Assert.DoesNotThrow(() => model.Predict(emptyInputs), "Empty input should not throw exception");

            model = new Prosit2020iRTTMT(modHandlingMode: IncompatibleModHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.Predict(emptyInputs), "Model should throw exception for peptides without required N-terminal modifications when ModHandlingMode is set to ThrowException");
        }

        /// <summary>
        /// Tests handling of null input list.
        /// Null list should throw ArgumentNullException.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelNullInput()
        {
            List<RetentionTimePredictionInput> nullList = null;

            var model = new Prosit2020iRTTMT();

            Assert.DoesNotThrow(() => model.Predict(nullList));
            Assert.That(model.Predictions.Count, Is.EqualTo(0), "Empty input should result in no predictions");
            Assert.That(model.ValidInputsMask.Count, Is.EqualTo(0), "Empty input should result in empty valid inputs mask");

            model = new Prosit2020iRTTMT(modHandlingMode: IncompatibleModHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.Predict(nullList), "Model should not throw exception for null input list");
        }

        /// <summary>
        /// Tests handling of mixed valid and invalid peptide sequences.
        /// Invalid sequences include too long, invalid characters, and sequences without N-terminal mods.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelInvalidSequences()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"), // Valid
                new RetentionTimePredictionInput("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), // Invalid - too long (and no N-term)
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]PEP*TIDE"), // Invalid - has N-term but invalid character
                new RetentionTimePredictionInput("SEQUENS"), // Invalid - no N-term mod (and noncanonical U)
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTING") // Valid
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(5), "Should return entries for all inputs");
            Assert.That(predictions[0].PredictedRetentionTime, Is.Not.Null, "First valid peptide should have predictions");
            Assert.That(predictions[1].PredictedRetentionTime, Is.Null, "Too long peptide should be rejected");
            Assert.That(predictions[2].PredictedRetentionTime, Is.Null, "Invalid character peptide should be rejected");
            Assert.That(predictions[3].PredictedRetentionTime, Is.Null, "Peptide without N-term should be rejected");
            Assert.That(predictions[4].PredictedRetentionTime, Is.Not.Null, "Second valid peptide should have predictions");

            model = new Prosit2020iRTTMT(modHandlingMode: IncompatibleModHandlingMode.ThrowException);
            Assert.Throws<ArgumentException>(() => model.Predict(modelInputs), "Model should throw exception for invalid peptides when ModHandlingMode is set to ThrowException");
        }

        /// <summary>
        /// Tests handling of all supported N-terminal modifications.
        /// The model supports TMT6plex, TMTpro, iTRAQ4plex, and iTRAQ8plex.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelAllNTerminalModificationTypes()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]SEQENCE"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTING"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ8plex on N-terminus]ANTHER")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(4), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All should have predictions");

            // Verify each N-terminal modification type is correctly converted
            Assert.That(predictions[0].FullSequence, Does.StartWith("[Common Fixed:TMT6plex on N-terminus]")); // TMT6plex
            Assert.That(predictions[1].FullSequence, Does.StartWith("[Common Fixed:TMTpro on N-terminus]")); // TMTpro
            Assert.That(predictions[2].FullSequence, Does.StartWith("[Common Fixed:iTRAQ4plex on N-terminus]")); // iTRAQ4plex
            Assert.That(predictions[3].FullSequence, Does.StartWith("[Common Fixed:iTRAQ8plex on N-terminus]")); // iTRAQ8plex
            Assert.That(model.ValidInputsMask, Is.All.True);
        }

        /// <summary>
        /// Tests handling of valid TMT/iTRAQ/SILAC modifications with required N-terminal mods.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelValidModificationMapping()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]M[Common Variable:Oxidation on M]SEQENS"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTC[Common Fixed:Carbamidomethyl on C]ING"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ8plex on N-terminus]PEPTK[Common Fixed:TMT6plex on K]IDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTK[Common Fixed:TMTpro on K]IDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]PEPTK[Common Fixed:iTRAQ4plex on K]IDE"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]PEPTK[Common Fixed:iTRAQ8plex on K]IDE"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ8plex on N-terminus]PEPTK[Common Variable:Label:13C(6)15N(2) on K]IDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTR[Common Variable:Label:13C(6)15N(4) on R]IDE")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(9), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All should have predictions");

            // Check that all sequences have N-terminal modifications
            foreach (var pred in predictions)
            {
                var modEnding = pred.FullSequence.IndexOf(']') + 1;
                Assert.That(model.ValidModificationUnimodMapping.Keys.Contains(pred.FullSequence.Substring(0, modEnding)), $"All sequences should have valid N-terminal modification: {pred.FullSequence}");
            }

            // Check that side-chain modifications are kept in the validated sequence and correctly identified as valid
            Assert.That(modelInputs[1].ValidatedFullSequence, Does.Contain("[Common Variable:Oxidation on M]")); // Oxidation
            Assert.That(modelInputs[2].ValidatedFullSequence, Does.Contain("[Common Fixed:Carbamidomethyl on C]"));  // Carbamidomethyl
            Assert.That(modelInputs[3].ValidatedFullSequence, Does.Contain("[Common Fixed:TMT6plex on K]")); // TMT6plex on K
            Assert.That(modelInputs[4].ValidatedFullSequence, Does.Contain("[Common Fixed:TMTpro on K]")); // TMTpro on K
            Assert.That(modelInputs[5].ValidatedFullSequence, Does.Contain("[Common Fixed:iTRAQ4plex on K]")); // iTRAQ4plex on K
            Assert.That(modelInputs[6].ValidatedFullSequence, Does.Contain("[Common Fixed:iTRAQ8plex on K]")); // iTRAQ8plex on K
            Assert.That(modelInputs[7].ValidatedFullSequence, Does.Contain("[Common Variable:Label:13C(6)15N(2) on K]")); // SILAC K
            Assert.That(modelInputs[8].ValidatedFullSequence, Does.Contain("[Common Variable:Label:13C(6)15N(4) on R]")); // SILAC R

            Assert.That(predictions[1].FullSequence, Does.Contain("[Common Variable:Oxidation on M]")); // Oxidation
            Assert.That(predictions[2].FullSequence, Does.Contain("[Common Fixed:Carbamidomethyl on C]"));  // Carbamidomethyl
            Assert.That(predictions[3].FullSequence, Does.Contain("[Common Fixed:TMT6plex on K]")); // TMT6plex on K
            Assert.That(predictions[4].FullSequence, Does.Contain("[Common Fixed:TMTpro on K]")); // TMTpro on K
            Assert.That(predictions[5].FullSequence, Does.Contain("[Common Fixed:iTRAQ4plex on K]")); // iTRAQ4plex on K
            Assert.That(predictions[6].FullSequence, Does.Contain("[Common Fixed:iTRAQ8plex on K]")); // iTRAQ8plex on K
            Assert.That(predictions[7].FullSequence, Does.Contain("[Common Variable:Label:13C(6)15N(2) on K]")); // SILAC K
            Assert.That(predictions[8].FullSequence, Does.Contain("[Common Variable:Label:13C(6)15N(4) on R]")); // SILAC R


        }

        /// <summary>
        /// Tests rejection of invalid or unsupported modifications.
        /// ModHandlingMode.ReturnNull rejects sequences with unsupported modifications.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelInvalidModifications()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"), // Valid
                new RetentionTimePredictionInput("SEQUENC[InvalidMod]E"), // Invalid - no N-term AND invalid mod
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]TESTING"), // Valid
                new RetentionTimePredictionInput("PEPTK[Common Fixed:Acetyl on K]IDE"), // Invalid - no N-term (acetylation also not supported)
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]PEPTK[Common Fixed:Acetyl on K]IDE") // Invalid - has N-term but acetylation not supported
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(5), "Should return entries for all inputs");
            Assert.That(predictions[0].PredictedRetentionTime, Is.Not.Null, "First valid peptide should have predictions");
            Assert.That(predictions[1].PredictedRetentionTime, Is.Null, "Invalid mod peptide should be rejected");
            Assert.That(predictions[2].PredictedRetentionTime, Is.Not.Null, "Second valid peptide should have predictions");
            Assert.That(predictions[3].PredictedRetentionTime, Is.Null, "Peptide without N-term should be rejected");
            Assert.That(predictions[4].PredictedRetentionTime, Is.Null, "Peptide with unsupported mod should be rejected");
        }

        /// <summary>
        /// Tests N-terminal TMT/iTRAQ modifications are correctly identified and converted.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelNTerminalTMTModifications()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ8plex on N-terminus]PEPTIDE")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(4), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All should have predictions");
            Assert.That(model.ValidInputsMask, Is.All.True, "All inputs should be valid with N-terminal modifications");

            // Check that N-terminal modifications are converted correctly
            for (int i = 0; i < modelInputs.Count; i++)
            {
                var modEnding = modelInputs[i].FullSequence.IndexOf(']') + 1;
                Assert.That(model.ValidModificationUnimodMapping.Keys.Contains(modelInputs[i].FullSequence.Substring(0, modEnding)), $"Input sequence should have valid N-terminal modification: {modelInputs[i].FullSequence}");
            }
        }

        /// <summary>
        /// Tests SILAC-labeled peptides with required N-terminal modifications.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelSILACModifications()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTK[Common Variable:Label:13C(6)15N(2) on K]IDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]PEPTR[Common Variable:Label:13C(6)15N(4) on R]IDE"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTK[Common Variable:Label:13C(6)15N(2) on K]R[Common Variable:Label:13C(6)15N(4) on R]ING")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), "Should return predictions for all inputs");
            // Only if silac mods are valid would retention time not be null.
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All should have predictions");
            Assert.That(model.ValidInputsMask, Is.All.True, "All inputs should be valid with N-terminal modifications");
        }

        /// <summary>
        /// Tests complex peptides with multiple modification types.
        /// All must have N-terminal modifications.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelComplexModificationCombinations()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTK[Common Fixed:TMT6plex on K]IDEC[Common Fixed:Carbamidomethyl on C]"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]M[Common Variable:Oxidation on M]EPTK[Common Fixed:TMTpro on K]IDEC[Common Fixed:Carbamidomethyl on C]"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTK[Common Variable:Label:13C(6)15N(2) on K]R[Common Variable:Label:13C(6)15N(4) on R]C[Common Fixed:Carbamidomethyl on C]ING")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All should have predictions");
            foreach (var input in modelInputs)
            {
                // since all inputs are valid and have N-terminal modifications, the validated sequence should match the original full sequence
                Assert.That(input.FullSequence == input.ValidatedFullSequence, $"Input sequence should be validated and match original for: {input.FullSequence}");
                Assert.That(input.SequenceWarning, Is.Null, $"Valid input should not have warnings for: {input.FullSequence}");
            }
        }

        /// <summary>
        /// Tests batching behavior with a small input list (below max batch size).
        /// All sequences must have N-terminal modifications.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelSmallBatchProcessing()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"),
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]SEQENS"),
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTING")
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(3), "Should return predictions for all inputs");
            Assert.That(predictions.All(p => p.PredictedRetentionTime != null), Is.True, "All should have predictions");
        }

        /// <summary>
        /// Tests batching behavior with a large input list (above max batch size).
        /// Generates 2500 random sequences with N-terminal TMT modifications.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelRequestBatching()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            int numberOfSequences = 2500;
            int seqLength = 20;
            var peptides = new HashSet<string>();

            while (peptides.Count < numberOfSequences)
            {
                var pep = "[Common Fixed:TMT6plex on N-terminus]" + new string(Random.Shared.GetItems(aminoacids, seqLength));
                peptides.Add(pep);
            }

            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2020iRTTMT();

            Assert.DoesNotThrow(() => model.Predict(modelInputs));

            var predictions = model.Predict(modelInputs);
            Assert.That(predictions.Count, Is.EqualTo(numberOfSequences), "Should return predictions for all inputs");
        }

        /// <summary>
        /// Tests batching when input size equals max batch size.
        /// All sequences must have N-terminal modifications.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelExactlyAtMaxBatchSize()
        {
            var modelInputs = new List<RetentionTimePredictionInput>();
            for (int i = 0; i < 1000; i++)
            {
                modelInputs.Add(new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"));
            }

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(1000), "Should return predictions for all inputs");
        }

        /// <summary>
        /// Tests the properties of the Prosit2020iRTTMT model.
        /// Validates model metadata and configuration.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelProperties()
        {
            var model = new Prosit2020iRTTMT();

            Assert.That(model.ModelName, Is.EqualTo("Prosit_2020_irt_TMT"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.IsIndexedRetentionTimeModel, Is.True);
            Assert.That(model.ModHandlingMode, Is.EqualTo(IncompatibleModHandlingMode.ReturnNull));
            Assert.That(model.ValidModificationUnimodMapping, Is.Not.Null);
            Assert.That(model.ValidModificationUnimodMapping.Count, Is.EqualTo(12), "Should support 12 modification types");
        }

        /// <summary>
        /// Tests handling of mixed valid and invalid sequences.
        /// Invalid sequences should have null predictions with warnings.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelMixedValidAndInvalidSequences()
        {
            var modelInputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]PEPTIDE"), // Valid
                new RetentionTimePredictionInput("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), // Invalid - too long AND no N-term
                new RetentionTimePredictionInput("[Common Fixed:TMTpro on N-terminus]SEQENS"), // Valid
                new RetentionTimePredictionInput("INVALID*"), // Invalid - bad character AND no N-term
                new RetentionTimePredictionInput("[Common Fixed:iTRAQ4plex on N-terminus]TESTINGTWICE"), // Valid
                new RetentionTimePredictionInput("INVALIDPEPTIDE[BadMod]"), // Invalid - bad mod AND no N-term
                new RetentionTimePredictionInput("VALIDSEQENCE") // Invalid - no N-term mod
            };

            var model = new Prosit2020iRTTMT();
            var predictions = model.Predict(modelInputs);

            Assert.That(predictions.Count, Is.EqualTo(7), "Should return entries for all inputs");
            Assert.That(model.ValidInputsMask.Count(v => v), Is.EqualTo(3), "Should have 3 valid inputs with N-term modifications");
            Assert.That(predictions.Count(p => p.PredictedRetentionTime != null), Is.EqualTo(3), "Only valid sequences with N-term should have predictions");
            Assert.That(predictions.Count(p => p.Warning != null), Is.EqualTo(4), "Invalid sequences should have warnings");

            // Verify all accepted sequences have N-terminal modifications
            foreach (var pred in predictions.Where(p => p.PredictedRetentionTime != null))
            {
                int modEnding = pred.FullSequence.IndexOf(']') + 1;
                Assert.That(model.ValidModificationUnimodMapping.Keys.Contains(pred.FullSequence.Substring(0, modEnding)), $"Accepted sequence should have N-terminal modification: {pred.FullSequence}");
            }
        }

        /// <summary>
        /// Tests custom constructor parameters for batching and throttling.
        /// </summary>
        [Test]
        public void TestProsit2020iRTTMTModelCustomConstructorParameters()
        {
            var model = new Prosit2020iRTTMT(
                modHandlingMode: IncompatibleModHandlingMode.ReturnNull,
                maxNumberOfBatchesPerRequest: 100,
                throttlingDelayInMilliseconds: 50
            );

            Assert.That(model.MaxNumberOfBatchesPerRequest, Is.EqualTo(100));
            Assert.That(model.ThrottlingDelayInMilliseconds, Is.EqualTo(50));
            Assert.That(model.ModHandlingMode, Is.EqualTo(IncompatibleModHandlingMode.ReturnNull));
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark test for the Prosit2020iRTTMT model with a large number of unique peptide sequences.
        /// This test is meant to evaluate the model's ability to handle large batch sizes and the efficiency of request batching logic.
        /// 
        /// Notes: 
        ///  - 1 MaxBatchSize is 1000 peptides and PFly2024FineTuned.Predict() takes under 1.5 seconds to run.
        ///  - With 4 million unique peptides of length 30, the test takes approximately 2 minutes to run. 
        ///  - Timing is linearly proportional to peptide number due to the throttled approach.
        ///  - The default MaxNumberOfBatchesPerRequest is currently set to 500 which was tested up to 4 million peptides without hitting client issues. 
        ///  - All peptides must include a required N-terminal modification (TMT6plex used here) for this model.
        /// </summary>
        public static void BenchmarkModelInputCountPerformance()
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            var modelInputs = new List<RetentionTimePredictionInput>();
            var seqLength = 30; // max length increases combinatorial space to ensure we always get unique sequences and benchmark response length handling
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
                modelInputs.Add(new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]" + peptide));
            }
            var model = new Prosit2020iRTTMT();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(numberOfSequences));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict {numberOfSequences:N0} peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }
    }
}

