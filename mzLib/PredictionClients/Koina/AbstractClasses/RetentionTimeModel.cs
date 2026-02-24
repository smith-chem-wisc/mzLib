using PredictionClients.Koina.Client;
using MzLibUtil;
using System.ComponentModel;
using PredictionClients.Koina.Interfaces;
using Easy.Common.Extensions;

namespace PredictionClients.Koina.AbstractClasses
{
    /// <summary>
    /// Represents a retention time prediction result for a single peptide sequence.
    /// Contains the original sequence, predicted retention time, and indexing information.
    /// </summary>
    /// <param name="FullSequence">The peptide sequence used for prediction (in UNIMOD format)</param>
    /// <param name="PredictedRetentionTime">Predicted retention time value (units depend on model - typically minutes or indexed RT)</param>
    /// <param name="IsIndexed">True if the model predicts indexed retention time (iRT); false for absolute retention time</param>
    /// <param name="Warning">Warning message if any issues occurred during prediction</param>
    public record PeptideRTPrediction(
        string FullSequence,
        double? PredictedRetentionTime,
        bool? IsIndexed,
        WarningException? Warning = null
    );

    /// <summary>
    /// Represents the input parameters for retention time prediction models from the Koina API.
    /// This record captures the input information required for peptide retention time prediction.
    /// </summary>
    /// <param name="FullSequence">Peptide sequence with modifications in mzLib format</param>
    public record RetentionTimePredictionInput(
        string FullSequence
    )
    {
        public string? ValidatedFullSequence { get; set; }
        public WarningException? Warning { get; set; }
    }

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
    public abstract class RetentionTimeModel : KoinaModelBase<RetentionTimePredictionInput, PeptideRTPrediction>, IPredictor<RetentionTimePredictionInput, PeptideRTPrediction>
    {
        #region Model-Specific Properties
        /// <summary>
        /// Indicates whether this model predicts indexed retention time (iRT) or absolute retention time.
        /// True = indexed retention time (relative scale), False = absolute retention time (minutes)
        /// </summary>
        public abstract bool IsIndexedRetentionTimeModel { get; }
        #endregion

        #region Inputs and Outputs for internal processing
        /// <summary>
        /// Inputs provided to the model for the LATEST prediction session. This list is populated during the prediction workflow 
        /// and is used to keep track of the original input parameters for each prediction, especially when batching is involved.
        /// </summary>
        public List<RetentionTimePredictionInput> ModelInputs { get; protected set; } = new();

        /// <summary>
        /// Boolean mask indicating which inputs from the original list were valid for prediction after applying model-specific validation criteria.
        /// This is used for realigning predictions back to the original input list and for filtering out invalid inputs from the prediction results.
        /// </summary>
        public bool[] ValidInputsMask { get; protected set; } = Array.Empty<bool>();

        /// <summary>
        /// Collection of retention time prediction results after inference completion for the LATEST prediction session.
        /// Each prediction contains the sequence, predicted RT, and indexing information.
        /// </summary>
        public List<PeptideRTPrediction> Predictions { get; protected set; } = new();
        #endregion

        /// <summary>
        /// Executes retention time prediction by sending batched requests to the Koina API.
        /// The method performs the following steps:
        /// 1. Validates and cleans input sequences according to model constraints, populating the ValidInputsMask to keep track of which inputs are valid for prediction.
        /// 2. Converts valid inputs into batched request payloads formatted for the specific model using the ToBatchedRequests method.
        /// 3. Sends batched requests to the Koina API with throttling between batches to avoid overwhelming the server, and processes responses to extract predictions.
        /// 4. Realigns predictions back to the original input list using the ValidInputsMask, ensuring that the output list corresponds to the original input order and includes placeholders for invalid inputs with appropriate warnings.
        /// </summary>
        /// <returns>Task representing the asynchronous inference operation</returns>
        protected virtual async Task<List<PeptideRTPrediction>> AsyncThrottledPredictor(List<RetentionTimePredictionInput> modelInputs)
        {
            #region Input Validation and Cleaning
            ModelInputs = modelInputs;
            ValidInputsMask = new bool[ModelInputs.Count];
            var validInputs = new List<RetentionTimePredictionInput>();

            for (int i = 0; i < ModelInputs.Count; i++)
            {
                WarningException? warning = null;
                var cleanedSequence = TryCleanSequence(ModelInputs[i].FullSequence, out var modHandlingWarning);
                warning = modHandlingWarning;

                if (cleanedSequence != null &&
                    HasValidModifications(cleanedSequence) &&
                    IsValidBaseSequence(cleanedSequence))
                {
                    ModelInputs[i] = ModelInputs[i] with { ValidatedFullSequence = cleanedSequence, Warning = warning };
                    ValidInputsMask[i] = true;
                    validInputs.Add(ModelInputs[i]);
                }
                else
                {
                    ModelInputs[i] = ModelInputs[i] with { ValidatedFullSequence = null, Warning = warning };
                    ValidInputsMask[i] = false;
                }
            }
            #endregion

            #region Request Batching, Throttling Setup
            var batchedRequests = ToBatchedRequests(validInputs);
            var batchChunks = batchedRequests.Chunk(MaxNumberOfBatchesPerRequest).ToList();
            int sessionTimeoutInMinutes = batchedRequests.Count * 2 + (int)(ThrottlingDelayInMilliseconds / 6000 * batchChunks.Count) + 2; // Dynamic timeout: ~2 minutes per batch + throttle time between batches + 2 minute buffer for network/processing overhead.
            #endregion

            #region Throttled API Requests and Response Processing
            var predictions = new List<PeptideRTPrediction>();
            var responses = new List<string>();
            using var _http = new HTTP(timeoutInMinutes: sessionTimeoutInMinutes);

            for (int i = 0; i < batchChunks.Count; i++)
            {
                var batchChunk = batchChunks[i];
                var responseChunk = await Task.WhenAll(batchChunk.Select(request => _http.InferenceRequest(ModelName, request)));
                responses.AddRange(responseChunk);

                if (i < batchChunks.Count - 1) // No need to throttle after the last batch
                {
                    await Task.Delay(ThrottlingDelayInMilliseconds);
                }
            }

            predictions = ResponseToPredictions(responses.ToArray(), validInputs);
            #endregion

            #region Realign Predictions to Original Input List
            // Realign predictions back to the original input list using the ValidInputsMask
            var realignedPredictions = new List<PeptideRTPrediction>();
            int predictionIndex = 0;

            for (int i = 0; i < ValidInputsMask.Length; i++)
            {
                if (ValidInputsMask[i])
                {
                    realignedPredictions.Add(predictions[predictionIndex]);
                    predictionIndex++;
                }
                else
                {
                    // For invalid inputs, add a placeholder prediction with a warning
                    realignedPredictions.Add(new PeptideRTPrediction(
                        FullSequence: ModelInputs[i].FullSequence,
                        PredictedRetentionTime: null,
                        IsIndexed: null,
                        Warning: ModelInputs[i].Warning ?? new WarningException("Input was invalid and skipped during prediction.")
                    ));
                }
            }
            #endregion

            Predictions = realignedPredictions;
            return Predictions;
        }

        public List<PeptideRTPrediction> Predict(List<RetentionTimePredictionInput> modelInputs)
        {
            return AsyncThrottledPredictor(modelInputs).GetAwaiter().GetResult();
        }

        /// <summary>
        /// Converts Koina API responses into structured PeptideRTPrediction objects.
        /// Expects responses with a single output containing retention time predictions.
        /// </summary>
        /// <param name="responses">Array of JSON response strings from Koina API</param>
        /// <param name="requestInputs">List of input parameters that were sent to the API</param>
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
        protected virtual List<PeptideRTPrediction> ResponseToPredictions(string[] responses, List<RetentionTimePredictionInput> requestInputs)
        {
            var predictions = new List<PeptideRTPrediction>();
            if (requestInputs.IsNullOrEmpty())
            {
                return predictions;
            }

            // Deserialize all batch responses
            var deserializedResponses = responses.Select(r => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(r)).ToList();

            // Validate successful deserialization
            if (deserializedResponses.IsNullOrEmpty() || deserializedResponses.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }

            // Extract retention time predictions from all batches (flattened)
            var rtOutputs = deserializedResponses.SelectMany(r => r!.Outputs[0].Data).ToList();

            // Ensure prediction count matches input count
            if (rtOutputs.Count != requestInputs.Count)
            {
                throw new Exception("The number of predictions does not match the number of input peptides.");
            }

            // Create prediction objects with sequence-to-prediction mapping
            for (int i = 0; i < requestInputs.Count; i++)
            {
                predictions.Add(new PeptideRTPrediction(
                    FullSequence: requestInputs[i].ValidatedFullSequence!,
                    PredictedRetentionTime: Convert.ToDouble(rtOutputs[i]),
                    IsIndexed: IsIndexedRetentionTimeModel,
                    Warning: requestInputs[i].Warning
                ));
            }

            return predictions;
        }
    }
}

