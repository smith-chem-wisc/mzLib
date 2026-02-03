using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace PredictionClients.Koina.AbstractClasses
{
    /// <summary>
    /// Represents a retention time prediction result for a single peptide sequence.
    /// Contains the original sequence, predicted retention time, and indexing information.
    /// </summary>
    /// <param name="PeptideSequence">The peptide sequence used for prediction (in UNIMOD format)</param>
    /// <param name="PredictedRetentionTime">Predicted retention time value (units depend on model - typically minutes or indexed RT)</param>
    /// <param name="IsIndexed">True if the model predicts indexed retention time (iRT); false for absolute retention time</param>
    public record PeptideRTPrediction(
        string PeptideSequence,
        double PredictedRetentionTime,
        bool IsIndexed
    );

    /// <summary>
    /// Abstract base class for retention time prediction models using the Koina API.
    /// Provides common functionality for sequence validation, request batching, and response parsing.
    /// Derived classes implement model-specific details including name, constraints, and request formatting.
    /// </summary>
    /// <remarks>
    /// Retention time models can predict either:
    /// - Absolute retention times (in minutes)
    /// - Indexed retention times (iRT values for relative comparison)
    /// 
    /// The model automatically handles:
    /// - Input validation for sequence length and amino acid composition
    /// - Modification format conversion from mzLib to UNIMOD format
    /// - Request batching for optimal API performance
    /// - Response parsing and error handling
    /// 
    /// Derived classes need to implement:
    /// - Model metadata (name, batch size, constraints)
    /// - Constructor that handles input data properties (sequences)
    /// - Request formatting method (ToBatchedRequests)
    /// - Model-specific modification handling if needed
    /// </remarks>
    public abstract class RetentionTimeModel : IKoinaModelIO
    {
        #region Model Metadata from Koina
        /// <summary>The Koina API model name identifier for retention time prediction</summary>
        public abstract string ModelName { get; }
        /// <summary>Maximum number of peptides that can be processed in a single API request</summary>
        public abstract int MaxBatchSize { get; }
        /// <summary>
        /// Indicates whether this model predicts indexed retention time (iRT) or absolute retention time.
        /// True = indexed retention time (relative scale), False = absolute retention time (minutes)
        /// </summary>
        public abstract bool IsIndexedRetentionTimeModel { get; }
        #endregion

        #region Input Constraints
        /// <summary>Maximum allowed peptide sequence length in amino acids</summary>
        public abstract int MaxPeptideLength { get; }
        /// <summary>Minimum allowed peptide sequence length in amino acids</summary>
        public abstract int MinPeptideLength { get; }
        /// <summary>
        /// Maps mzLib modification format to UNIMOD format for model compatibility.
        /// Example: "[Common Variable:Oxidation on M]" -> "[UNIMOD:35]"
        /// </summary>
        public virtual Dictionary<string, string> ValidModificationUnimodMapping => new();
        #endregion

        #region Validation Patterns and Filters
        /// <summary>
        /// Regex pattern for detecting modifications in peptide sequences.
        /// Default matches content within square brackets: [modification_name]
        /// </summary>
        public virtual string ModificationPattern => @"\[[^\]]+\]";
        /// <summary>
        /// Regex pattern matching valid canonical amino acid sequences.
        /// Default includes the 20 standard proteinogenic amino acids.
        /// </summary>
        public virtual string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        #endregion

        #region Input Data
        /// <summary>
        /// List of peptide sequences for retention time prediction.
        /// Sequences should be in a format compatible with the specific model requirements.
        /// </summary>
        public abstract List<string> PeptideSequences { get; }
        #endregion

        #region Output Data
        /// <summary>
        /// Collection of retention time prediction results after inference completion.
        /// Each prediction contains the sequence, predicted RT, and indexing information.
        /// </summary>
        public abstract List<PeptideRTPrediction> Predictions { get; protected set; }
        #endregion

        #region Querying Methods for the Koina API
        /// <summary>
        /// Converts input peptide sequences into batched requests formatted for the specific model's API endpoint.
        /// Each batch should not exceed MaxBatchSize peptides and must include all required model inputs.
        /// </summary>
        /// <returns>List of request dictionaries, each containing a batch of model inputs in Koina API format</returns>
        /// <remarks>
        /// Implementation should handle:
        /// - Model-specific input formatting and data types
        /// - Proper batching based on MaxBatchSize constraints
        /// - Required tensor shapes and metadata for the API
        /// </remarks>
        protected abstract List<Dictionary<string, object>> ToBatchedRequests();

        /// <summary>
        /// Executes retention time prediction by sending batched requests to the Koina API.
        /// Handles HTTP client lifecycle, request batching, and response parsing automatically.
        /// </summary>
        /// <returns>Task representing the asynchronous inference operation</returns>
        /// <remarks>
        /// Process flow:
        /// 1. Creates HTTP client with dynamic timeout based on input size
        /// 2. Sends all batched requests concurrently using Task.WhenAll
        /// 3. Processes responses and populates Predictions collection
        /// 4. Ensures proper resource cleanup regardless of success/failure
        /// </remarks>
        /// <exception cref="Exception">Thrown when API responses cannot be deserialized or processed</exception>
        public virtual async Task RunInferenceAsync()
        {
            // Dynamic timeout: ~1 minute per batch + 2 minute buffer for network/processing overhead
            var _http = new HTTP(timeoutInMinutes: PeptideSequences.Count / MaxBatchSize * 2 + 2);

            try
            {
                var responses = await Task.WhenAll(ToBatchedRequests().Select(request => _http.InferenceRequest(ModelName, request)));
                ResponseToPredictions(responses);
            }
            finally
            {
                _http.Dispose(); // Ensure HTTP client is properly disposed
            }
        }

        /// <summary>
        /// Converts Koina API responses into structured PeptideRTPrediction objects.
        /// Expects responses with a single output containing retention time predictions.
        /// </summary>
        /// <param name="responses">Array of JSON response strings from Koina API</param>
        /// <exception cref="Exception">
        /// Thrown when:
        /// - Response deserialization fails
        /// - Number of predictions doesn't match input peptides
        /// - Response format is unexpected
        /// </exception>
        /// <remarks>
        /// Processing steps:
        /// 1. Deserializes JSON responses from all batches
        /// 2. Extracts retention time values from the first (and only) output
        /// 3. Validates prediction count matches input sequence count
        /// 4. Creates PeptideRTPrediction objects with sequence mapping
        /// </remarks>
        protected virtual void ResponseToPredictions(string[] responses)
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

            // Extract retention time predictions from all batches (flattened)
            var rtOutputs = deserializedResponses.SelectMany(r => r.Outputs[0].Data).ToList();

            // Ensure prediction count matches input count
            if (rtOutputs.Count != PeptideSequences.Count)
            {
                throw new Exception("The number of predictions does not match the number of input peptides.");
            }

            // Create prediction objects with sequence-to-prediction mapping
            Predictions = PeptideSequences
                .Select((seq, index) => new PeptideRTPrediction(
                    PeptideSequence: seq,
                    PredictedRetentionTime: Convert.ToDouble(rtOutputs[index]),
                    IsIndexed: IsIndexedRetentionTimeModel
                ))
                .ToList();
        }
        #endregion

        #region Validation and Modification Handling
        /// <summary>
        /// Validates that all modifications in a peptide sequence are supported by the model.
        /// Returns true for sequences without modifications or when all modifications are recognized.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential modifications in mzLib format</param>
        /// <returns>True if all modifications are valid or no modifications present; otherwise false</returns>
        /// <remarks>
        /// Modification validation process:
        /// 1. Uses ModificationPattern regex to find all modification annotations
        /// 2. Checks each modification against ValidModificationUnimodMapping
        /// 3. Empty modification list (unmodified peptides) is considered valid
        /// </remarks>
        protected virtual bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            if (matches.Count == 0)
            {
                return true; // No modifications found - valid for all models
            }

            // Check if any modifications are not in the valid mapping
            return matches.Where(m => !ValidModificationUnimodMapping.ContainsKey(m.Value)).Count() == 0;
        }

        /// <summary>
        /// Validates that a peptide sequence meets model constraints for length and amino acid composition.
        /// Removes modifications before validation to check only the base amino acid sequence.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential modifications</param>
        /// <returns>True if sequence meets model length and composition requirements; otherwise false</returns>
        /// <remarks>
        /// Validation criteria:
        /// - Strips modification annotations using ModificationPattern
        /// - Checks against CanonicalAminoAcidPattern for valid amino acids
        /// - Validates length is within [MinPeptideLength, MaxPeptideLength] range
        /// </remarks>
        protected virtual bool IsValidSequence(string sequence)
        {
            // Remove modification annotations to get base sequence
            var baseSequence = Regex.Replace(sequence, ModificationPattern, "");

            return Regex.IsMatch(baseSequence, CanonicalAminoAcidPattern) // Valid amino acids only
                && baseSequence.Length <= MaxPeptideLength                 // Within max length
                && baseSequence.Length >= MinPeptideLength;                // Above min length (implicit from abstract property)
        }
        #endregion

        #region Full Sequence Modification Conversion Methods
        /// <summary>
        /// Converts peptide sequence from mzLib modification format to UNIMOD format required by the model.
        /// Base implementation performs standard mzLib to UNIMOD conversion using the ValidModificationUnimodMapping.
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib modification format</param>
        /// <returns>Sequence converted to UNIMOD format</returns>
        /// <remarks>
        /// Conversion process:
        /// 1. Replaces mzLib modification names with UNIMOD identifiers using ValidModificationUnimodMapping
        /// 
        /// Example transformations:
        /// - "PEPT[Common Variable:Oxidation on M]IDE" -> "PEPT[UNIMOD:35]IDE" (oxidation converted)
        /// - "C[Common Fixed:Carbamidomethyl on C]PEPTIDE" -> "C[UNIMOD:4]PEPTIDE" (carbamidomethyl converted)
        /// 
        /// Derived classes may override this method to implement model-specific modification handling,
        /// such as automatic carbamidomethylation of cysteines or other required modifications.
        /// </remarks>
        protected virtual string ConvertMzLibModificationsToUnimod(string sequence)
        {
            // Apply custom modification mappings first
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Key, mod.Value);
            }
            return sequence;
        }
        #endregion
    }
}
