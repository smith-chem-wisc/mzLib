using PredictionClients.Koina.Client;
using MzLibUtil;

namespace PredictionClients.Koina.AbstractClasses
{
    /// <summary>
    /// Represents a retention time prediction result for a single peptide sequence.
    /// Contains the original sequence, predicted retention time, and indexing information.
    /// </summary>
    /// <param name="FullSequence">The peptide sequence used for prediction (in UNIMOD format)</param>
    /// <param name="PredictedRetentionTime">Predicted retention time value (units depend on model - typically minutes or indexed RT)</param>
    /// <param name="IsIndexed">True if the model predicts indexed retention time (iRT); false for absolute retention time</param>
    public record PeptideRTPrediction(
        string FullSequence,
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
    public abstract class RetentionTimeModel : KoinaModelBase
    {
        #region Model Metadata from Koina
        /// <summary>
        /// Indicates whether this model predicts indexed retention time (iRT) or absolute retention time.
        /// True = indexed retention time (relative scale), False = absolute retention time (minutes)
        /// </summary>
        public abstract bool IsIndexedRetentionTimeModel { get; }
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
        public override async Task RunInferenceAsync()
        {
            if (_disposed)
            {
                throw new ObjectDisposedException(nameof(RetentionTimeModel), "Cannot run inference on a disposed model instance. The model is meant to be used only on initialized peptides. The results are still accessible.");
            }
            // Dynamic timeout: ~2 minutes per batch + 2 minute buffer for network/processing overhead. Typically a 
            // batch takes less than a minute. 
            int numBatches = (int)Math.Ceiling((double)PeptideSequences.Count / MaxBatchSize);
            using var _http = new HTTP(timeoutInMinutes: numBatches * 2 + 2);
            
            var responses = await Task.WhenAll(ToBatchedRequests().Select(request => _http.InferenceRequest(ModelName, request)));
            ResponseToPredictions(responses);
            Dispose();
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
        protected override void ResponseToPredictions(string[] responses)
        {
            if (PeptideSequences.Count == 0)
            {
                Predictions = new List<PeptideRTPrediction>();
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
                    FullSequence: seq,
                    PredictedRetentionTime: Convert.ToDouble(rtOutputs[index]),
                    IsIndexed: IsIndexedRetentionTimeModel
                ))
                .ToList();
        }
        #endregion

    }
}
