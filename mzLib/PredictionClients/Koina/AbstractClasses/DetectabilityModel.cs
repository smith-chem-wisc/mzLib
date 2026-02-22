using PredictionClients.Koina.Client;
using Easy.Common.Extensions;
using System.ComponentModel;
using MzLibUtil;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Interfaces;

namespace PredictionClients.Koina.SupportedModels.FlyabilityModels
{
    /// <summary>
    /// Represents the prediction results for a single peptide, containing detectability probability scores
    /// for each detectability class from a detectability prediction model.
    /// </summary>
    /// <param name="FullSequence">The peptide sequence (with modifications in UNIMOD format)</param>
    /// <param name="DetectabilityProbabilities">Probability scores for each detectability class (Not Detectable, Low, Intermediate, High)</param>
    /// <param name="Warning">Warning message if any issues occurred during prediction</param>
    public record PeptideDetectabilityPrediction(
        string FullSequence,
        (double NotDetectable,
         double LowDetectability,
         double IntermediateDetectability,
         double HighDetectability) DetectabilityProbabilities,
        WarningException? Warning = null
    );

    /// <summary>
    /// Represents the input parameters for detectability prediction models from the Koina API.
    /// This record captures the input information required for peptide detectability prediction.
    /// </summary>
    /// <param name="FullSequence">Peptide sequence with modifications in mzLib format</param>
    public record DetectabilityPredictionInput(
        string FullSequence
    )
    {
        public string? ValidatedFullSequence { get; set; }
        public WarningException? Warning { get; set; }
    }

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
    public abstract class DetectabilityModel : KoinaModelBase<DetectabilityPredictionInput, PeptideDetectabilityPrediction>, IPredictor<DetectabilityPredictionInput, PeptideDetectabilityPrediction>
    {
        #region Model-Specific Properties
        /// <summary>
        /// Number of detectability classes predicted by the model.
        /// Fixed at 4 classes: Not Detectable, Low, Intermediate, High.
        /// </summary>
        public abstract int NumberOfDetectabilityClasses { get; }

        /// <summary>
        /// Human-readable names for the detectability classes in prediction order.
        /// Corresponds to the probability scores returned by the model.
        /// </summary>
        public abstract List<string> DetectabilityClasses { get; }
        #endregion

        #region Inputs and Outputs for internal processing
        /// <summary>
        /// Inputs provided to the model for the LATEST prediction session. This list is populated during the prediction workflow 
        /// and is used to keep track of the original input parameters for each prediction, especially when batching is involved.
        /// </summary>
        public List<DetectabilityPredictionInput> ModelInputs { get; protected set; } = new();

        /// <summary>
        /// Boolean mask indicating which inputs from the original list were valid for prediction after applying model-specific validation criteria.
        /// This is used for realigning predictions back to the original input list and for filtering out invalid inputs from the prediction results.
        /// </summary>
        public bool[] ValidInputsMask { get; protected set; } = Array.Empty<bool>();

        /// <summary>
        /// Collection of detectability prediction results after inference completion for the LATEST prediction session. 
        /// Each PeptideDetectabilityPrediction contains probability scores for all four detectability classes.
        /// </summary>
        public List<PeptideDetectabilityPrediction> Predictions { get; protected set; } = new();
        #endregion

        /// <summary>
        /// Executes peptide detectability prediction by sending batched requests to the Koina API.
        /// The method performs the following steps:
        /// 1. Validates and cleans input sequences according to model constraints, populating the ValidInputsMask to keep track of which inputs are valid for prediction.
        /// 2. Converts valid inputs into batched request payloads formatted for the specific model using the ToBatchedRequests method.
        /// 3. Sends batched requests to the Koina API with throttling between batches to avoid overwhelming the server, and processes responses to extract predictions.
        /// 4. Realigns predictions back to the original input list using the ValidInputsMask, ensuring that the output list corresponds to the original input order and includes placeholders for invalid inputs with appropriate warnings.
        /// </summary>
        /// <returns>Task representing the asynchronous inference operation</returns>
        protected virtual async Task<List<PeptideDetectabilityPrediction>> AsyncThrottledPredictor(List<DetectabilityPredictionInput> modelInputs)
        {
            #region Input Validation and Cleaning
            ModelInputs = modelInputs;
            ValidInputsMask = new bool[ModelInputs.Count];
            var validInputs = new List<DetectabilityPredictionInput>();

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
            var predictions = new List<PeptideDetectabilityPrediction>();
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
            var realignedPredictions = new List<PeptideDetectabilityPrediction>();
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
                    realignedPredictions.Add(new PeptideDetectabilityPrediction(
                        FullSequence: ModelInputs[i].FullSequence,
                        DetectabilityProbabilities: (0.0, 0.0, 0.0, 0.0),
                        Warning: ModelInputs[i].Warning ?? new WarningException("Input was invalid and skipped during prediction.")
                    ));
                }
            }
            #endregion

            Predictions = realignedPredictions;
            return Predictions;
        }

        public List<PeptideDetectabilityPrediction> Predict(List<DetectabilityPredictionInput> modelInputs)
        {
            return AsyncThrottledPredictor(modelInputs).GetAwaiter().GetResult();
        }

        /// <summary>
        /// Converts Koina API responses into structured PeptideDetectabilityPrediction objects.
        /// Expects responses with detectability probability scores for each peptide across 4 classes.
        /// </summary>
        /// <param name="responses">Array of JSON response strings from Koina API</param>
        /// <param name="requestInputs">List of input parameters that were sent to the API</param>
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
        protected virtual List<PeptideDetectabilityPrediction> ResponseToPredictions(string[] responses, List<DetectabilityPredictionInput> requestInputs)
        {
            var predictions = new List<PeptideDetectabilityPrediction>();
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

            // Extract and chunk detectability predictions (4 probabilities per peptide)
            var detectabilityPredictions = deserializedResponses
                .SelectMany(batch => batch!.Outputs[0].Data.Chunk(NumberOfDetectabilityClasses))
                .ToList();

            // Create prediction objects with probability scores for each detectability class
            for (int i = 0; i < requestInputs.Count; i++)
            {
                var peptideFlyabilityClassProbs = detectabilityPredictions[i].Select(p => (double)p).ToList();
                predictions.Add(new PeptideDetectabilityPrediction(
                    requestInputs[i].ValidatedFullSequence!,
                    (
                        NotDetectable: peptideFlyabilityClassProbs[0],
                        LowDetectability: peptideFlyabilityClassProbs[1],
                        IntermediateDetectability: peptideFlyabilityClassProbs[2],
                        HighDetectability: peptideFlyabilityClassProbs[3]
                    ),
                    Warning: requestInputs[i].Warning
                ));
            }

            return predictions;
        }
    }
}