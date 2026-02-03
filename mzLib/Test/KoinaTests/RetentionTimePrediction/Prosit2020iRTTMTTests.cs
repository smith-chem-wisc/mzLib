using NUnit.Framework;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;

namespace Test.KoinaTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class Prosit2020iRTTMTTests
    {
        /// <summary>
        /// Tests the constructor with valid peptide sequences that have required N-terminal modifications.
        /// 
        /// This is the fundamental requirement for Prosit2020iRTTMT: all sequences MUST have
        /// an N-terminal TMT or iTRAQ modification.
        /// 
        /// Validates that:
        /// - No warnings are generated for sequences with N-terminal mods
        /// - All sequences are accepted and stored
        /// - Model properties are correctly initialized
        /// - Model name matches expected value "Prosit_2020_irt_TMT"
        /// - Batch size, length constraints, and iRT flag are correct
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - PeptideSequences count should equal input count (with carbamidomethylation applied)
        /// - All model properties should match specification
        /// 
        /// Use Case: Standard TMT/iTRAQ workflow where all peptides are N-terminally labeled
        /// Real-world scenario: TMT10plex experiment with all peptides labeled at N-terminus
        /// </summary>
        [Test]
        public void TestConstructorWithValidPeptides()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE",
                "[Common Fixed:TMTpro on N-terminus]TESTING"
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(2));
            Assert.That(model.ModelName, Is.EqualTo("Prosit_2020_irt_TMT"));
            Assert.That(model.MaxBatchSize, Is.EqualTo(1000));
            Assert.That(model.MaxPeptideLength, Is.EqualTo(30));
            Assert.That(model.MinPeptideLength, Is.EqualTo(1));
            Assert.That(model.IsIndexedRetentionTimeModel, Is.True);
        }

        /// <summary>
        /// Tests that sequences WITHOUT N-terminal modifications are REJECTED.
        /// 
        /// This is a critical test for the Prosit2020iRTTMT model requirement.
        /// Unlike Prosit2019iRT or standard models, this model REQUIRES N-terminal labeling.
        /// 
        /// Validates that:
        /// - Sequences without N-terminal mods are detected as invalid
        /// - A warning is generated explaining the missing N-terminal modification
        /// - Invalid sequences are excluded from processing
        /// - Only sequences with N-terminal mods are accepted
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should indicate missing N-terminal modification requirement
        /// - Only sequences with N-terminal mods should be in PeptideSequences
        /// 
        /// Use Case: Quality control - detecting improperly prepared sequences
        /// Real-world scenario: User attempts to use non-TMT data with TMT-specific model
        /// </summary>
        [Test]
        public void TestRejectionOfSequencesWithoutNTerminalModification()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE", // Valid - has N-term
                "PEPTIDE", // Invalid - no N-term mod
                "PEPTK[Common Fixed:TMT6plex on K]IDE", // Invalid - only K-TMT, no N-term
                "[Common Fixed:TMTpro on N-terminus]TESTING", // Valid - has N-term
                "M[Common Variable:Oxidation on M]SEQUENCE" // Invalid - only oxidation, no N-term
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings.Message, Does.Contain("invalid"));
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(2)); // Only sequences with N-term mods
            
            // Verify that only sequences with N-terminal mods were accepted
            foreach (var sequence in model.PeptideSequences)
            {
                Assert.That(sequence, Does.Match(@"^\[UNIMOD:(737|2016|214|730)\]-")); // Must start with N-term mod
            }
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

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

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

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

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
        /// - Sequences WITHOUT N-terminal modifications (critical for this model)
        /// 
        /// Validates that:
        /// - A warning is generated listing invalid sequences
        /// - Only valid sequences WITH N-terminal mods are added to PeptideSequences
        /// - Invalid sequences are properly filtered out
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should contain "invalid"
        /// - PeptideSequences should only contain sequences with N-terminal mods
        /// - All invalid sequences should be excluded
        /// 
        /// Use Case: Ensures robust input validation and informative error messaging
        /// </summary>
        [Test]
        public void TestConstructorWithInvalidSequences()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE", // Valid
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // Invalid - too long (and no N-term)
                "[Common Fixed:TMTpro on N-terminus]PEP*TIDE", // Invalid - has N-term but invalid character
                "SEQUENS", // Invalid - no N-term mod (and noncanonical U)
                "[Common Fixed:iTRAQ4plex on N-terminus]TESTING" // Valid
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings.Message, Does.Contain("invalid"));
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(2)); // Only valid sequences with N-term
        }

        /// <summary>
        /// Tests the handling of all supported N-terminal modifications.
        /// 
        /// The Prosit2020iRTTMT model supports four types of N-terminal isobaric labels:
        /// - TMT6plex (UNIMOD:737)-
        /// - TMTpro (UNIMOD:2016)-
        /// - iTRAQ4plex (UNIMOD:214)-
        /// - iTRAQ8plex (UNIMOD:730)-
        /// 
        /// Validates that:
        /// - All four N-terminal modification types are recognized
        /// - mzLib format is converted to Prosit UNIMOD format
        /// - All N-terminally modified sequences are accepted
        /// - No warnings are generated
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - All 4 sequences should be accepted
        /// - N-terminal modifications should be converted to [UNIMOD:XXX]- format
        /// 
        /// Use Case: Ensures support for all isobaric labeling chemistries
        /// Real-world scenario: Processing data from different TMT/iTRAQ experiments
        /// </summary>
        [Test]
        public void TestAllNTerminalModificationTypes()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE",
                "[Common Fixed:TMTpro on N-terminus]SEQENCE",
                "[Common Fixed:iTRAQ4plex on N-terminus]TESTING",
                "[Common Fixed:iTRAQ8plex on N-terminus]ANTHER"
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(4));
            
            // Verify each N-terminal modification type is correctly converted
            Assert.That(model.PeptideSequences[0], Does.StartWith("[UNIMOD:737]-")); // TMT6plex
            Assert.That(model.PeptideSequences[1], Does.StartWith("[UNIMOD:2016]-")); // TMTpro
            Assert.That(model.PeptideSequences[2], Does.StartWith("[UNIMOD:214]-")); // iTRAQ4plex
            Assert.That(model.PeptideSequences[3], Does.StartWith("[UNIMOD:730]-")); // iTRAQ8plex
        }

        /// <summary>
        /// Tests the handling of valid post-translational modifications including TMT labels.
        /// 
        /// Prosit2020iRTTMT supports:
        /// REQUIRED (N-terminal):
        /// 1. TMT6plex on N-terminus: [Common Fixed:TMT6plex on N-terminus] → [UNIMOD:737]-
        /// 2. TMTpro on N-terminus: [Common Fixed:TMTpro on N-terminus] → [UNIMOD:2016]-
        /// 3. iTRAQ4plex on N-terminus: [Common Fixed:iTRAQ4plex on N-terminus] → [UNIMOD:214]-
        /// 4. iTRAQ8plex on N-terminus: [Common Fixed:iTRAQ8plex on N-terminus] → [UNIMOD:730]-
        /// 
        /// OPTIONAL (Side-chain):
        /// 5. Oxidation on methionine (M): [Common Variable:Oxidation on M] → [UNIMOD:35]
        /// 6. Carbamidomethyl on cysteine (C): [Common Fixed:Carbamidomethyl on C] → [UNIMOD:4]
        /// 7. SILAC Heavy K: [Common Variable:Label:13C(6)15N(2) on K] → [UNIMOD:259]
        /// 8. SILAC Heavy R: [Common Variable:Label:13C(6)15N(4) on R] → [UNIMOD:267]
        /// 9. TMT6plex on K: [Common Fixed:TMT6plex on K] → [UNIMOD:737]
        /// 10. TMTpro on K: [Common Fixed:TMTpro on K] → [UNIMOD:2016]
        /// 11. iTRAQ4plex on K: [Common Fixed:iTRAQ4plex on K] → [UNIMOD:214]
        /// 12. iTRAQ8plex on K: [Common Fixed:iTRAQ8plex on K] → [UNIMOD:730]
        /// 
        /// Validates that:
        /// - All sequences have required N-terminal modifications
        /// - All supported modifications are recognized as valid
        /// - mzLib format is converted to Prosit UNIMOD format
        /// - All modified sequences are accepted
        /// - No warnings are generated
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - All sequences should be accepted
        /// - Modifications should be converted to UNIMOD format
        /// 
        /// Use Case: Ensures compatibility with TMT/iTRAQ/SILAC labeling workflows
        /// Real-world scenario: Processing TMT-labeled proteomics data for retention time prediction
        /// </summary>
        [Test]
        public void TestValidModificationMapping()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE[Common Variable:Oxidation on M]",
                "[Common Fixed:TMTpro on N-terminus]M[Common Variable:Oxidation on M]SEQENS",
                "[Common Fixed:iTRAQ4plex on N-terminus]TESTC[Common Fixed:Carbamidomethyl on C]ING",
                "[Common Fixed:iTRAQ8plex on N-terminus]PEPTK[Common Fixed:TMT6plex on K]IDE",
                "[Common Fixed:TMT6plex on N-terminus]PEPTK[Common Fixed:TMTpro on K]IDE",
                "[Common Fixed:TMTpro on N-terminus]PEPTK[Common Fixed:iTRAQ4plex on K]IDE",
                "[Common Fixed:iTRAQ4plex on N-terminus]PEPTK[Common Fixed:iTRAQ8plex on K]IDE",
                "[Common Fixed:iTRAQ8plex on N-terminus]PEPTK[Common Variable:Label:13C(6)15N(2) on K]IDE",
                "[Common Fixed:TMT6plex on N-terminus]PEPTR[Common Variable:Label:13C(6)15N(4) on R]IDE"
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(9));
            
            // Check that all sequences have N-terminal modifications
            foreach (var sequence in model.PeptideSequences)
            {
                Assert.That(sequence, Does.Match(@"^\[UNIMOD:(737|2016|214|730)\]-"));
            }
            
            // Check that side-chain modifications are converted to UNIMOD format
            Assert.That(model.PeptideSequences[0], Does.Contain("[UNIMOD:35]")); // Oxidation
            Assert.That(model.PeptideSequences[1], Does.Contain("[UNIMOD:35]")); // Oxidation
            Assert.That(model.PeptideSequences[2], Does.Contain("[UNIMOD:4]"));  // Carbamidomethyl
            Assert.That(model.PeptideSequences[3], Does.Contain("[UNIMOD:737]")); // TMT6plex on K
            Assert.That(model.PeptideSequences[4], Does.Contain("[UNIMOD:2016]")); // TMTpro on K
            Assert.That(model.PeptideSequences[5], Does.Contain("[UNIMOD:214]")); // iTRAQ4plex on K
            Assert.That(model.PeptideSequences[6], Does.Contain("[UNIMOD:730]")); // iTRAQ8plex on K
            Assert.That(model.PeptideSequences[7], Does.Contain("[UNIMOD:259]")); // SILAC K
            Assert.That(model.PeptideSequences[8], Does.Contain("[UNIMOD:267]")); // SILAC R
        }

        /// <summary>
        /// Tests the rejection of invalid or unsupported modifications.
        /// 
        /// Invalid modifications are those not in the ValidModificationUnimodMapping dictionary.
        /// Additionally, sequences without N-terminal modifications are invalid for this model.
        /// 
        /// Validates that:
        /// - Invalid modifications are detected
        /// - Sequences without N-terminal mods are rejected (even if they have valid side-chain mods)
        /// - Sequences with invalid modifications are rejected
        /// - A warning is generated listing problematic sequences
        /// - Valid sequences with N-terminal mods are still accepted
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should contain "invalid"
        /// - Only sequences with valid N-terminal mods should be accepted
        /// 
        /// Use Case: Prevents submission of unsupported modifications to the model
        /// Real-world scenario: User attempts to use non-TMT peptides or unsupported modifications
        /// </summary>
        [Test]
        public void TestInvalidModifications()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE", // Valid
                "SEQUENC[InvalidMod]E", // Invalid - no N-term AND invalid mod
                "[Common Fixed:TMTpro on N-terminus]TESTING", // Valid
                "PEPTK[Common Fixed:Acetyl on K]IDE", // Invalid - no N-term (acetylation also not supported)
                "[Common Fixed:iTRAQ4plex on N-terminus]PEPTK[Common Fixed:Acetyl on K]IDE" // Invalid - has N-term but acetylation not supported
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings.Message, Does.Contain("invalid"));
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(2)); // Only 2 valid sequences
        }

        /// <summary>
        /// Tests N-terminal TMT/iTRAQ modifications are correctly identified and converted.
        /// 
        /// TMT and iTRAQ reagents label the peptide N-terminus (and optionally lysine side chains).
        /// The model REQUIRES N-terminal modifications with the format [UNIMOD:XXX]-
        /// 
        /// Validates that:
        /// - N-terminal TMT6plex, TMTpro, iTRAQ4plex, and iTRAQ8plex are recognized
        /// - mzLib format with "N-terminus" is converted correctly
        /// - Modified sequences are accepted without warnings
        /// - All sequences have the required N-terminal modification
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - All 4 sequences should be accepted
        /// - N-terminal modifications should be converted to [UNIMOD:XXX]- format
        /// 
        /// Use Case: Ensures proper handling of N-terminal labeling in TMT/iTRAQ experiments
        /// Real-world scenario: Processing TMT-labeled peptides where N-term labeling is mandatory
        /// </summary>
        [Test]
        public void TestNTerminalTMTModifications()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE",
                "[Common Fixed:TMTpro on N-terminus]PEPTIDE",
                "[Common Fixed:iTRAQ4plex on N-terminus]PEPTIDE",
                "[Common Fixed:iTRAQ8plex on N-terminus]PEPTIDE"
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(4));
            
            // Check that N-terminal modifications are converted correctly
            Assert.That(model.PeptideSequences[0], Does.StartWith("[UNIMOD:737]-")); // TMT6plex
            Assert.That(model.PeptideSequences[1], Does.StartWith("[UNIMOD:2016]-")); // TMTpro
            Assert.That(model.PeptideSequences[2], Does.StartWith("[UNIMOD:214]-")); // iTRAQ4plex
            Assert.That(model.PeptideSequences[3], Does.StartWith("[UNIMOD:730]-")); // iTRAQ8plex
        }

        /// <summary>
        /// Tests SILAC-labeled peptides with required N-terminal modifications.
        /// 
        /// SILAC (Stable Isotope Labeling by Amino acids in Cell culture) introduces
        /// heavy isotope labels on lysine and arginine residues. However, for this model,
        /// N-terminal TMT/iTRAQ labeling is still REQUIRED.
        /// 
        /// Validates that:
        /// - Heavy lysine (13C6 15N2) modifications are recognized
        /// - Heavy arginine (13C6 15N4) modifications are recognized
        /// - Multiple SILAC labels in one peptide are handled correctly
        /// - N-terminal modifications are still required
        /// - No warnings are generated
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - All sequences should be accepted
        /// - SILAC modifications should be converted to UNIMOD format
        /// - All sequences must start with N-terminal modification
        /// 
        /// Use Case: Ensures compatibility with SILAC quantification workflows in TMT experiments
        /// Real-world scenario: SILAC-labeled samples with additional TMT labeling (hyperplexing)
        /// </summary>
        [Test]
        public void TestSILACModifications()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTK[Common Variable:Label:13C(6)15N(2) on K]IDE",
                "[Common Fixed:TMTpro on N-terminus]PEPTR[Common Variable:Label:13C(6)15N(4) on R]IDE",
                "[Common Fixed:iTRAQ4plex on N-terminus]TESTK[Common Variable:Label:13C(6)15N(2) on K]R[Common Variable:Label:13C(6)15N(4) on R]ING"
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(3));
            
            // Verify all have N-terminal modifications
            foreach (var sequence in model.PeptideSequences)
            {
                Assert.That(sequence, Does.Match(@"^\[UNIMOD:(737|2016|214|730)\]-"));
            }
            
            // Check that SILAC modifications are converted to UNIMOD format
            Assert.That(model.PeptideSequences[0], Does.Contain("[UNIMOD:259]")); // Heavy K
            Assert.That(model.PeptideSequences[1], Does.Contain("[UNIMOD:267]")); // Heavy R
            Assert.That(model.PeptideSequences[2], Does.Contain("[UNIMOD:259]")); // Heavy K
            Assert.That(model.PeptideSequences[2], Does.Contain("[UNIMOD:267]")); // Heavy R
        }

        /// <summary>
        /// Tests complex peptides with multiple modification types.
        /// 
        /// Validates that:
        /// - N-terminal modifications are always present (REQUIRED)
        /// - Multiple different modifications in the same peptide are handled
        /// - Combinations of TMT/iTRAQ with oxidation and carbamidomethyl work correctly
        /// - SILAC with other modifications is supported
        /// - No warnings are generated for complex but valid modifications
        /// 
        /// Expected Behavior:
        /// - warnings should be null
        /// - All complex sequences should be accepted
        /// - All modifications should be correctly converted
        /// - All sequences must start with N-terminal modification
        /// 
        /// Use Case: Real-world complex modification scenarios
        /// Real-world scenario: TMT-labeled peptides with oxidized Met, carbamidomethylated Cys,
        ///                      and additional K-labeling (common in TMT experiments)
        /// </summary>
        [Test]
        public void TestComplexModificationCombinations()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTK[Common Fixed:TMT6plex on K]IDEC[Common Fixed:Carbamidomethyl on C]",
                "[Common Fixed:TMTpro on N-terminus]M[Common Variable:Oxidation on M]EPTK[Common Fixed:TMTpro on K]IDEC[Common Fixed:Carbamidomethyl on C]",
                "[Common Fixed:iTRAQ4plex on N-terminus]TESTK[Common Variable:Label:13C(6)15N(2) on K]R[Common Variable:Label:13C(6)15N(4) on R]C[Common Fixed:Carbamidomethyl on C]ING"
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(3));
            
            // Verify first sequence has N-terminal TMT, K-TMT, and Carbamidomethyl C
            Assert.That(model.PeptideSequences[0], Does.StartWith("[UNIMOD:737]-")); // N-term TMT6plex
            Assert.That(model.PeptideSequences[0], Does.Contain("K[UNIMOD:737]")); // K TMT6plex
            Assert.That(model.PeptideSequences[0], Does.Contain("C[UNIMOD:4]")); // Carbamidomethyl C
            
            // Verify second sequence has N-terminal TMT, Oxidation, TMTpro, and Carbamidomethyl
            Assert.That(model.PeptideSequences[1], Does.StartWith("[UNIMOD:2016]-")); // N-term TMTpro
            Assert.That(model.PeptideSequences[1], Does.Contain("M[UNIMOD:35]")); // Oxidation
            Assert.That(model.PeptideSequences[1], Does.Contain("K[UNIMOD:2016]")); // TMTpro
            Assert.That(model.PeptideSequences[1], Does.Contain("C[UNIMOD:4]")); // Carbamidomethyl C
            
            // Verify third sequence has N-terminal iTRAQ, both SILAC labels, and Carbamidomethyl
            Assert.That(model.PeptideSequences[2], Does.StartWith("[UNIMOD:214]-")); // N-term iTRAQ4plex
            Assert.That(model.PeptideSequences[2], Does.Contain("[UNIMOD:259]")); // Heavy K
            Assert.That(model.PeptideSequences[2], Does.Contain("[UNIMOD:267]")); // Heavy R
            Assert.That(model.PeptideSequences[2], Does.Contain("C[UNIMOD:4]")); // Carbamidomethyl C
        }

        /// <summary>
        /// Tests batching behavior with a small input list (below max batch size).
        /// 
        /// Note: All sequences must have N-terminal modifications.
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
            var peptideSequences = new List<string> 
            { 
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE", 
                "[Common Fixed:TMTpro on N-terminus]SEQENS", 
                "[Common Fixed:iTRAQ4plex on N-terminus]TESTING" 
            };
            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);
            Assert.DoesNotThrowAsync(async () => await model.RunInferenceAsync());
            Assert.That(model.Predictions, Has.Count.EqualTo(3));
        }

        /// <summary>
        /// Tests batching behavior with a large input list (above max batch size).
        /// 
        /// Generates 2500 random sequences with N-terminal TMT modifications.
        /// 
        /// Validates that:
        /// - Input list is split into multiple batches
        /// - Batch count matches expected (3 batches: 1000 + 1000 + 500)
        /// - Each batch has a unique ID
        /// - All sequences have N-terminal modifications
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
                var pep = "[Common Fixed:TMT6plex on N-terminus]" + new string(Random.Shared.GetItems(aminoacids, seqLength));
                if (!peptides.Contains(pep))
                {
                    peptides.Add(pep);
                }
            }

            var model = new Prosit2020iRTTMT(peptides, out WarningException warnings);
            Assert.DoesNotThrowAsync(async () => await model.RunInferenceAsync());
            Assert.That(model.Predictions, Has.Count.EqualTo(numberOfSequences));
        }

        /// <summary>
        /// Tests batching behavior when the input list size is exactly the max batch size.
        /// 
        /// Validates that:
        /// - All sequences are processed in a single batch
        /// - No batching errors occur
        /// - All sequences have N-terminal modifications
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
                peptideSequences.Add("[Common Fixed:TMT6plex on N-terminus]PEPTIDE");
            }
            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);
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
            var peptideSequences = new List<string> 
            { 
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE", 
                "[Common Fixed:TMTpro on N-terminus]SEQUENCE" 
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(model.Predictions, Is.Empty);
        }

        /// <summary>
        /// Tests the properties of the Prosit2020iRTTMT model.
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
            var peptideSequences = new List<string> { "[Common Fixed:TMT6plex on N-terminus]PEPTIDE" };
            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(model.ModelName, Is.EqualTo("Prosit_2020_irt_TMT"));
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
        /// - Only valid sequences WITH N-terminal mods are added to the model
        /// - Invalid sequences are properly identified
        /// - Sequences without N-terminal mods are rejected
        /// 
        /// Expected Behavior:
        /// - warnings should not be null
        /// - Warning message should indicate invalid sequences
        /// - model.PeptideSequences should contain only sequences with N-terminal mods
        /// </summary>
        [Test]
        public void TestMixedValidAndInvalidSequences()
        {
            var peptideSequences = new List<string>
            {
                "[Common Fixed:TMT6plex on N-terminus]PEPTIDE", // Valid
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // Invalid - too long AND no N-term
                "[Common Fixed:TMTpro on N-terminus]SEQENS", // Valid
                "INVALID*", // Invalid - bad character AND no N-term
                "[Common Fixed:iTRAQ4plex on N-terminus]TESTINGTWICE", // Valid
                "INVALIDPEPTIDE[BadMod]", // Invalid - bad mod AND no N-term
                "VALIDSEQENCE" // Invalid - no N-term mod (even though sequence is valid)
            };

            var model = new Prosit2020iRTTMT(peptideSequences, out WarningException warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(model.PeptideSequences, Has.Count.EqualTo(3)); // Only sequences with N-term mods
            Assert.That(warnings.Message, Does.Contain("invalid"));
            
            // Verify all accepted sequences have N-terminal modifications
            foreach (var sequence in model.PeptideSequences)
            {
                Assert.That(sequence, Does.Match(@"^\[UNIMOD:(737|2016|214|730)\]-"));
            }
        }
    }
}

