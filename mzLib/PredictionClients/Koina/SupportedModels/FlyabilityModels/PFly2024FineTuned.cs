using Easy.Common.Extensions;
using MzLibUtil;
using System.ComponentModel;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.FlyabilityModels
{
    /// <summary>
    /// Implementation of the pFly 2024 fine-tuned peptide detectability prediction model.
    /// Predicts the likelihood of peptide detection in mass spectrometry experiments using machine learning.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length up to 40 amino acids
    /// - Does not support modified peptides for detectability prediction
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
    public class PFly2024FineTuned : DetectabilityModel
    {
        /// <summary>The Koina API model name identifier for peptide detectability prediction</summary>
        public override string ModelName => "pfly_2024_fine_tuned";

        /// <summary>Maximum number of peptides that can be processed in a single API request</summary>
        public override int MaxBatchSize => 128;

        /// <summary>
        /// Maximum number of batches that should be processed in a single API request. This is necessary 
        /// to prevent overwhelming the server with too many concurrent requests, which can lead to timeouts 
        /// or rate limiting. Adjust this value based on the expected number of peptides and server capacity.
        /// </summary>
        public override int MaxNumberOfBatchesPerRequest { get; init; }

        /// <summary>
        /// Throttle time between batches to avoid overwhelming the server. 
        /// Adjust as needed based on model performance and server capacity.
        /// </summary> 
        public override int ThrottlingDelayInMilliseconds { get; init; }

        /// <summary>
        /// Number of detectability classes predicted by the model.
        /// Fixed at 4 classes: Not Detectable, Low, Intermediate, High.
        /// </summary>
        public override int NumberOfDetectabilityClasses => 4;

        /// <summary>
        /// Human-readable names for the detectability classes in prediction order.
        /// Corresponds to the probability scores returned by the model.
        /// </summary>
        public override List<string> DetectabilityClasses => new() { "Not Detectable", "Low Detectability", "Intermediate Detectability", "High Detectability" };

        /// <summary>Maximum allowed peptide sequence length in amino acids</summary>
        public override int MaxPeptideLength => 40;

        /// <summary>Minimum allowed peptide sequence length in amino acids</summary>
        public override int MinPeptideLength => 1;

        public PFly2024FineTuned(int maxNumberOfBatchesPerRequest = 500, int throttlingDelayInMilliseconds = 100)
        {
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
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
        protected override List<Dictionary<string, object>> ToBatchedRequests(List<DetectabilityPredictionInput> validInputs)
        {
            // Split sequences into smaller batches optimized for detectability prediction
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence).Chunk(MaxBatchSize).ToList();
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
    }
}