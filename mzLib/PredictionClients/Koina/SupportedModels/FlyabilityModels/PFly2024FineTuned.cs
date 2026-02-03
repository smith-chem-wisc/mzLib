using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using Easy.Common.Extensions;
using MzLibUtil;
using System.Text.RegularExpressions;
using System.ComponentModel;


namespace PredictionClients.Koina.SupportedModels.FlyabilityModels
{
    /// <summary>
    /// Represents a peptide detectability prediction result with probability scores for each detectability class.
    /// Contains the original sequence and probability scores across four detectability categories.
    /// </summary>
    /// <param name="PeptideSequence">The peptide sequence used for detectability prediction</param>
    /// <param name="NotDetectable">Probability score (0-1) that the peptide is not detectable by MS</param>
    /// <param name="LowDetectability">Probability score (0-1) that the peptide has low detectability</param>
    /// <param name="IntermediateDetectability">Probability score (0-1) that the peptide has intermediate detectability</param>
    /// <param name="HighDetectability">Probability score (0-1) that the peptide has high detectability</param>
    /// <remarks>
    /// The four probability scores sum to 1.0, representing a complete probability distribution
    /// across detectability classes. Higher scores indicate greater likelihood for that class.
    /// </remarks>
    public record PeptideDetectabilityPrediction(
        string PeptideSequence,
        double NotDetectable,
        double LowDetectability,
        double IntermediateDetectability,
        double HighDetectability
    );

    /// <summary>
    /// Implementation of the pFly 2024 fine-tuned peptide detectability prediction model.
    /// Predicts the likelihood of peptide detection in mass spectrometry experiments using machine learning.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length up to 40 amino acids
    /// - Processes up to 128 peptides per batch
    /// - Predicts detectability across 4 classes: Not Detectable, Low, Intermediate, High
    /// - Returns probability scores for each detectability class
    /// - Trained on experimental MS data to predict peptide flyability (detectability)
    /// 
    /// Detectability classes:
    /// - Not Detectable: Peptides unlikely to be observed in MS experiments
    /// - Low Detectability: Peptides with poor ionization/detection characteristics
    /// - Intermediate Detectability: Peptides with moderate detection likelihood
    /// - High Detectability: Peptides with excellent ionization/detection characteristics
    /// 
    /// Use cases:
    /// - Peptide selection for targeted proteomics experiments
    /// - Proteome coverage optimization
    /// - Method development and peptide filtering
    /// - In silico peptide screening prior to synthesis
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/pfly_2024_fine_tuned/infer
    /// </remarks>
    public class PFly2024FineTuned : IKoinaModelIO
    {
        /// <summary>The Koina API model name identifier for peptide detectability prediction</summary>
        public string ModelName => "pfly_2024_fine_tuned";
        /// <summary>Maximum number of peptides that can be processed in a single API request</summary>
        public int MaxBatchSize => 128;
        /// <summary>
        /// Number of detectability classes predicted by the model.
        /// Fixed at 4 classes: Not Detectable, Low, Intermediate, High.
        /// </summary>
        public int NumberOfDetectabilityClasses => 4;
        /// <summary>
        /// Human-readable names for the detectability classes in prediction order.
        /// Corresponds to the probability scores returned by the model.
        /// </summary>
        public List<string> DetectabilityClasses => new() { "Not Detectable", "Low Detectability", "Intermediate Detectability", "High Detectability" };
        /// <summary>Maximum allowed peptide sequence length in amino acids</summary>
        public int MaxPeptideLength => 40;
        /// <summary>
        /// Regex pattern matching valid canonical amino acid sequences.
        /// Includes the 20 standard proteinogenic amino acids.
        /// </summary>
        public string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        /// <summary>
        /// Regex pattern for detecting modifications in peptide sequences.
        /// Matches content within square brackets: [modification_name]
        /// </summary>
        public string ModificationPattern => @"\[[^\]]+\]";
        /// <summary>
        /// List of validated peptide sequences for detectability prediction.
        /// Contains sequences that passed validation for length and amino acid composition.
        /// </summary>
        public List<string> PeptideSequences { get; } = new();
        /// <summary>
        /// Collection of detectability prediction results after inference completion.
        /// Each prediction contains probability scores for all four detectability classes.
        /// </summary>
        public List<PeptideDetectabilityPrediction> Predictions = new();

        /// <summary>
        /// Initializes a new instance of the PFly2024FineTuned model with input validation and filtering.
        /// Validates sequences against model requirements and prepares valid sequences for prediction.
        /// </summary>
        /// <param name="peptideSequences">
        /// List of peptide sequences for detectability prediction.
        /// Sequences can contain modifications in square bracket notation, but only sequence length 
        /// and amino acid composition are validated (modifications are ignored for this model).
        /// </param>
        /// <param name="warnings">
        /// Output parameter containing details about any invalid sequences that were filtered out.
        /// Will be null if all input sequences are valid.
        /// </param>
        /// <exception cref="WarningException">
        /// Returned via warnings parameter when invalid sequences are encountered or input is empty.
        /// </exception>
        /// <remarks>
        /// Validation criteria:
        /// - Sequence length: ≤ 40 amino acids (after removing modification annotations)
        /// - Only canonical amino acids (20 standard proteinogenic amino acids)
        /// - Empty or null input results in warning
        /// - Modifications are stripped for validation but preserved in output
        /// 
        /// Processing steps:
        /// 1. Validates input is not empty
        /// 2. Strips modification annotations for length/composition validation
        /// 3. Filters sequences that meet model requirements
        /// 4. Preserves original sequences (with modifications) for prediction
        /// 5. Collects invalid sequences for warning message
        /// </remarks>
        public PFly2024FineTuned(List<string> peptideSequences, out WarningException? warnings)
        {
            // Handle empty input case early
            if (peptideSequences.IsNullOrEmpty())
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            var invalidSequences = new List<string>();
            // Validate each sequence against model requirements
            foreach (var seq in peptideSequences)
            {
                if (IsValidPeptideSequence(seq))
                {
                    PeptideSequences.Add(seq); // Keep original sequence with modifications
                }
                else
                {
                    invalidSequences.Add(seq);
                }
            }

            // Generate warning message for invalid sequences if any were found
            warnings = invalidSequences.Count > 0
                ? new WarningException($"The following peptide sequences were invalid and will be skipped: {string.Join(", ", invalidSequences)}")
                : null;
        }

        /// <summary>
        /// Creates batched API requests formatted for the pFly 2024 fine-tuned detectability model.
        /// Each batch contains up to MaxBatchSize peptide sequences for optimal API performance.
        /// </summary>
        /// <returns>
        /// List of request dictionaries compatible with Koina API format.
        /// Each request contains peptide sequences for detectability classification.
        /// </returns>
        /// <remarks>
        /// Request structure follows Koina API specification for pfly_2024_fine_tuned:
        /// - peptide_sequences: BYTES array containing peptide sequences (modifications preserved)
        /// - Shape: [batch_size, 1] for tensor compatibility
        /// - Datatype: BYTES for string sequence data
        /// 
        /// Batching strategy:
        /// - Splits input sequences into chunks of MaxBatchSize (128)
        /// - Each batch gets a unique identifier for tracking
        /// - Smaller batch size optimized for detectability model's computational requirements
        /// - Enables concurrent processing of large peptide sets
        /// </remarks>
        protected List<Dictionary<string, object>> ToBatchedRequests()
        {
            // Split sequences into smaller batches optimized for detectability prediction
            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
            var batchedRequests = new List<Dictionary<string, object>>();

            for (int i = 0; i < batchedPeptides.Count; i++)
            {
                // Create API request structure following Koina specification
                var request = new Dictionary<string, object>
                {
                    { "id", $"Batch{i}_" + Guid.NewGuid()}, // Unique identifier for batch tracking
                    { "inputs", new List<object>
                        {
                            new {
                                name = "peptide_sequences",           // Model input parameter name
                                shape = new[]{ batchedPeptides[i].Length, 1 }, // Tensor shape [batch_size, 1]
                                datatype = "BYTES",                   // String data type for sequences
                                data = batchedPeptides[i]            // Actual sequence data
                            },
                        }
                    }
                };
                batchedRequests.Add(request);
            }
            return batchedRequests;
        }
        public async Task RunInferenceAsync()
        {
            // Dynamic timeout: ~2 minutes per batch + 2 minute buffer for network/processing overhead. Typically a 
            // batch takes less than a minute. 
            int numBatches = (int)Math.Ceiling((double)PeptideSequences.Count / MaxBatchSize);
            using var _http = new HTTP(timeoutInMinutes: numBatches * 2 + 2);
            
            var responses = await Task.WhenAll(ToBatchedRequests().Select(request => _http.InferenceRequest(ModelName, request)));
            ResponseToPredictions(responses);
            
        }

        /// <summary>
        /// Converts Koina API responses into structured PeptideDetectabilityPrediction objects.
        /// Expects responses with detectability probability scores for each peptide across 4 classes.
        /// </summary>
        /// <param name="responses">Array of JSON response strings from Koina API</param>
        /// <exception cref="Exception">
        /// Thrown when:
        /// - Response deserialization fails
        /// - Response format is unexpected
        /// - Number of predictions doesn't match input sequences
        /// </exception>
        /// <remarks>
        /// Response processing steps:
        /// 1. Deserializes JSON responses from all batches
        /// 2. Extracts detectability probability arrays from the single output
        /// 3. Chunks probability data into groups of 4 (one per detectability class)
        /// 4. Creates PeptideDetectabilityPrediction objects with sequence mapping
        /// 5. Validates that probability arrays have correct dimensions
        /// 
        /// Expected response format:
        /// - Single output containing flattened probability scores
        /// - 4 consecutive values per peptide (Not Detectable, Low, Intermediate, High)
        /// - Probability scores should sum to 1.0 for each peptide
        /// </remarks>
        protected void ResponseToPredictions(string[] responses)
        {
            if (PeptideSequences.Count == 0)
            {
                return; // No input sequences to process
            }

            // Deserialize all batch responses
            var deserializedResponses = responses.Select(r => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(r)).ToList();

            // Validate successful deserialization
            if (deserializedResponses.IsNullOrEmpty() || deserializedResponses.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }

            // Extract and chunk detectability predictions (4 probabilities per peptide)
            var detectabilityPredictions = deserializedResponses
                .SelectMany(batch => batch!.Outputs[0].Data.Chunk(NumberOfDetectabilityClasses))
                .ToList();

            // Create prediction objects with probability scores for each detectability class
            for (int i = 0; i < PeptideSequences.Count; i++)
            {
                var probs = detectabilityPredictions[i];
                Predictions.Add(new PeptideDetectabilityPrediction(
                    PeptideSequences[i],               // Original sequence
                    Convert.ToDouble(probs[0]),        // Not Detectable probability
                    Convert.ToDouble(probs[1]),        // Low Detectability probability
                    Convert.ToDouble(probs[2]),        // Intermediate Detectability probability
                    Convert.ToDouble(probs[3])         // High Detectability probability
                ));
            }
        }

        /// <summary>
        /// Validates that a peptide sequence meets model constraints for length and amino acid composition.
        /// Removes modification annotations before validation to check only the base amino acid sequence.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential modifications</param>
        /// <returns>True if sequence meets model requirements; otherwise false</returns>
        /// <remarks>
        /// Validation criteria:
        /// 1. Strips modification annotations using ModificationPattern regex
        /// 2. Checks unmodified sequence length against MaxPeptideLength (40 amino acids)
        /// 3. Validates amino acid composition against CanonicalAminoAcidPattern (20 standard amino acids)
        /// 
        /// Note: This model focuses on sequence-based detectability prediction, so modifications
        /// are preserved in the sequence but not used for validation constraints.
        /// The base amino acid sequence determines detectability characteristics.
        /// </remarks>
        protected bool IsValidPeptideSequence(string sequence)
        {
            // Remove modification annotations to get base amino acid sequence
            var unmodifiedSequence = Regex.Replace(sequence, ModificationPattern, string.Empty);

            // Check length constraint
            if (unmodifiedSequence.Length > MaxPeptideLength)
            {
                return false;
            }

            // Check amino acid composition (canonical amino acids only)
            if (!Regex.IsMatch(unmodifiedSequence, CanonicalAminoAcidPattern))
            {
                return false;
            }

            return true;
        }
    }
}
