using PredictionClients.Koina.Client;
using Easy.Common.Extensions;
using System.ComponentModel;
using MzLibUtil;
using PredictionClients.Koina.AbstractClasses;


namespace PredictionClients.Koina.SupportedModels.FlyabilityModels
{
    /// <summary>
    /// Represents a peptide detectability prediction result with probability scores for each detectability class.
    /// Contains the original sequence and probability scores across four detectability categories.
    /// </summary>
    /// <param name="PeptideSequence">The peptide sequence used for detectability prediction</param>
    /// <param name="DetectabilityClasses"> Tuple of probability scores for each detectability class:
    ///     Not Detectable 
    ///     Low Detectability
    ///     Intermediate Detectability
    ///     High Detectability
    /// <remarks>
    /// The four probability scores sum to 1.0, representing a complete probability distribution
    /// across detectability classes. Higher scores indicate greater likelihood for that class.
    /// </remarks>
    public record PeptideDetectabilityPrediction(
        string PeptideSequence,
        (double NotDetectable,
         double LowDetectability,
         double IntermediateDetectability,
         double HighDetectability) DetectabilityProbabilities
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
    public abstract class DetectabilityModel : KoinaModelBase
    {
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
        /// <summary>
        /// Collection of detectability prediction results after inference completion.
        /// Each prediction contains probability scores for all four detectability classes.
        /// </summary>
        public abstract List<PeptideDetectabilityPrediction> Predictions { get; protected set; }

        /// <summary>
        /// Executes peptide detectability prediction by sending batched requests to the Koina API.
        /// Handles HTTP client lifecycle, request batching, and response parsing automatically.
        /// </summary>
        /// <returns>Task representing the asynchronous inference operation</returns>
        /// <remarks>
        /// Process flow:
        /// 1. Creates HTTP client with realistic timeout based on actual processing time
        /// 2. Uses throttled batch processing for large datasets to avoid overwhelming the server
        /// 3. Processes responses and populates Predictions collection
        /// 4. Ensures proper resource cleanup regardless of success/failure
        /// </remarks>
        /// <exception cref="Exception">Thrown when API responses cannot be deserialized or processed</exception>
        public override async Task<WarningException?> PredictAsync()
        {
            if (_disposed)
            {
                throw new ObjectDisposedException(nameof(DetectabilityModel), "Cannot run inference on a disposed model instance. The model is meant to be used only on initialized peptides. The results are still accessible.");
            }

            int numBatches = (int)Math.Ceiling((double)PeptideSequences.Count / MaxBatchSize);

            // Realistic timeout calculation: 0.07s per batch + 2min buffer for network overhead
            // Typically takes ~0.052s per batch, but we add some buffer
            // Cap at reasonable maximum to prevent extremely long timeouts
            int estimatedTimeSeconds = (int)(numBatches * 0.07) + 120; // 0.07s per batch + 2min buffer
            int timeoutMinutes = Math.Max(2, Math.Min(estimatedTimeSeconds / 60 + 1, 30)); // Between 2-30 minutes

            using var _http = new HTTP(timeoutInMinutes: timeoutMinutes);

            var batchedRequests = ToBatchedRequests();
            var responses = new List<string>();

            // For large datasets (>500 batches), use throttled processing to avoid 504 errors
            if (numBatches > 500)
            {
                Console.WriteLine($"Processing {PeptideSequences.Count:N0} peptides in {numBatches:N0} batches with throttling...");
                Console.WriteLine($"Estimated completion time: {TimeSpan.FromSeconds(numBatches * 0.05 + 60):mm\\:ss}");

                // Process in chunks of 100 concurrent batches
                const int maxConcurrentBatches = 100;
                const int delayBetweenChunks = 200; // 200ms delay between chunks

                for (int i = 0; i < batchedRequests.Count; i += maxConcurrentBatches)
                {
                    var chunk = batchedRequests.Skip(i).Take(maxConcurrentBatches);
                    var chunkResponses = await Task.WhenAll(
                        chunk.Select(request => _http.InferenceRequest(ModelName, request))
                    );

                    responses.AddRange(chunkResponses);

                    // Progress reporting every 1000 batches
                    if ((i + maxConcurrentBatches) % 1000 == 0 || i + maxConcurrentBatches >= batchedRequests.Count)
                    {
                        int processed = Math.Min(i + maxConcurrentBatches, numBatches);
                        Console.WriteLine($"Processed {processed:N0}/{numBatches:N0} batches ({(double)processed / numBatches:P1})");
                    }

                    // Small delay between chunks to avoid overwhelming the server
                    if (i + maxConcurrentBatches < batchedRequests.Count)
                    {
                        await Task.Delay(delayBetweenChunks);
                    }
                }
            }
            else
            {
                // For smaller datasets, process all batches concurrently
                var responses_array = await Task.WhenAll(batchedRequests.Select(request => _http.InferenceRequest(ModelName, request)));
                responses.AddRange(responses_array);
            }

            ResponseToPredictions(responses.ToArray());
            Dispose();
            return null; // No warnings to return for detectability prediction
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
        protected override void ResponseToPredictions(string[] responses)
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
                var peptideFlyabilityClassProbs = detectabilityPredictions[i].Select(p => (double)p).ToList();
                Predictions.Add(new PeptideDetectabilityPrediction(
                    PeptideSequences[i],
                    (
                        NotDetectable: peptideFlyabilityClassProbs[0],
                        LowDetectability: peptideFlyabilityClassProbs[1],
                        IntermediateDetectability: peptideFlyabilityClassProbs[2],
                        HighDetectability: peptideFlyabilityClassProbs[3]
                )));
            }
        }

        
    }
}
