using PredictionClients.Koina.Util;
using Easy.Common.Extensions;
using MzLibUtil;
using PredictionClients.Koina.AbstractClasses;
using System.ComponentModel;

namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    /// <summary>
    /// Implementation of the Prosit 2019 indexed retention time (iRT) prediction model.
    /// Predicts indexed retention times for peptide sequences using the Koina API.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Processes up to 1000 peptides per batch
    /// - Predicts indexed retention time (iRT) values for relative comparison
    /// - Supports carbamidomethylation on cysteine and oxidation on methionine
    /// 
    /// iRT values provide relative retention time measurements that are independent of 
    /// chromatographic conditions, enabling cross-laboratory comparison of retention times.
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Prosit_2019_irt/infer
    /// </remarks>
    public class Prosit2019iRT : RetentionTimeModel
    {
        /// <summary>The Koina API model name identifier</summary>
        public override string ModelName => "Prosit_2019_irt";

        /// <summary>Maximum number of peptides that can be processed in a single API request</summary>
        public override int MaxBatchSize => 1000;

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
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;

        /// <summary>Maximum allowed peptide sequence length in amino acids</summary>
        public override int MaxPeptideLength => 30;

        /// <summary>Minimum allowed peptide sequence length in amino acids</summary>
        public override int MinPeptideLength => 1;

        /// <summary>
        /// Indicates this model predicts indexed retention time (iRT) values.
        /// iRT values are relative measurements independent of chromatographic conditions.
        /// </summary>
        public override bool IsIndexedRetentionTimeModel => true;

        /// <summary>
        /// Supported modifications mapping from mzLib format to UNIMOD format.
        /// Only carbamidomethylation on cysteine and oxidation on methionine are supported.
        /// </summary>
        public override Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            {"[Common Variable:Oxidation on M]", "[UNIMOD:35]"},
            {"[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]"}
        };
        public override IncompatibleModHandlingMode ModHandlingMode { get; init; }

        public Prosit2019iRT(IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.RemoveIncompatibleMods, int maxNumberOfBatchesPerRequest = 500, int throttlingDelayInMilliseconds = 100)
        {
            ModHandlingMode = modHandlingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
        }

        /// <summary>
        /// Creates batched API requests formatted for the Prosit 2019 iRT model.
        /// Each batch contains up to MaxBatchSize peptide sequences for optimal API performance.
        /// </summary>
        /// <returns>
        /// List of request dictionaries compatible with Koina API format.
        /// Each request contains peptide sequences in UNIMOD format with carbamidomethylated cysteines.
        /// </returns>
        /// <remarks>
        /// Request structure follows Koina API specification for Prosit_2019_irt:
        /// - peptide_sequences: BYTES array containing UNIMOD-formatted sequences
        /// - Shape: [batch_size, 1] for tensor compatibility
        /// - Datatype: BYTES for string sequence data
        /// 
        /// Batching strategy:
        /// - Splits input sequences into chunks of MaxBatchSize (1000)
        /// - Each batch gets a unique identifier for tracking
        /// - Enables concurrent processing of large sequence sets
        /// </remarks>
        protected override List<Dictionary<string, object>> ToBatchedRequests(List<RetentionTimePredictionInput> validInputs)
        {
            // Split sequences into batches for optimal API performance
            // ValidatedFullSequence should not be null at this point due to prior validation steps
            var batchedPeptides = validInputs.Select(p => ConvertMzLibModificationsToUnimod(p.ValidatedFullSequence!)).Chunk(MaxBatchSize).ToList();
            var batchedRequests = new List<Dictionary<string, object>>();

            for (int i = 0; i < batchedPeptides.Count; i++)
            {
                // Create API request structure following Koina specification
                var request = new Dictionary<string, object>
                {
                    { "id", $"Batch{i}_" + Guid.NewGuid()}, // Unique identifier for request tracking
                    { "inputs", new List<object>
                        {
                            new {
                                name = "peptide_sequences",           // Model input parameter name
                                shape = new[]{ batchedPeptides[i].Length, 1 }, // Tensor shape [batch_size, 1]
                                datatype = "BYTES",                   // String data type
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

