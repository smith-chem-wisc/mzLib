using System.ComponentModel;
using MzLibUtil;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    /// <summary>
    /// Implementation of the Prosit 2020 indexed retention time (iRT) prediction model with TMT labeling support.
    /// Predicts indexed retention times for TMT-labeled peptide sequences using the Koina API.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Processes up to 1000 peptides per batch
    /// - Predicts indexed retention time (iRT) values for relative comparison
    /// - Specialized for TMT (Tandem Mass Tag) and iTRAQ labeled peptides
    /// - Supports various isobaric labeling strategies including TMT6plex, TMTpro, iTRAQ 4-plex and 8-plex
    /// - Also supports SILAC labeling and standard modifications
    /// 
    /// Supported labeling methods:
    /// - TMT6plex and TMTpro labeling on lysine and N-terminus
    /// - iTRAQ 4-plex and 8-plex labeling on lysine and N-terminus
    /// - SILAC heavy labeling (13C6 15N2 on K, 13C6 15N4 on R)
    /// - Standard modifications (oxidation on M, carbamidomethyl on C)
    /// 
    /// iRT values provide relative retention time measurements that are independent of 
    /// chromatographic conditions, enabling cross-laboratory comparison of retention times
    /// for labeled peptides.
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Prosit_2020_irt_TMT/infer
    /// </remarks>
    public class Prosit2020iRTTMT : RetentionTimeModel
    {
        /// <summary>The Koina API model name identifier for TMT-capable iRT prediction</summary>
        public override string ModelName => "Prosit_2020_irt_TMT";

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
        /// iRT values are relative measurements independent of chromatographic conditions,
        /// specifically calibrated for TMT-labeled peptides.
        /// </summary>
        public override bool IsIndexedRetentionTimeModel => true;

        private static readonly IReadOnlySet<int> AllowedMods = new HashSet<int>
        {
            35, 4, 259, 267, 737, 2016, 214, 730, 2016
        };

        private static readonly UnimodSequenceFormatSchema TmtSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");

        public override IReadOnlySet<int> AllowedUnimodIds => AllowedMods;
        protected override UnimodSequenceFormatSchema UnimodSchema => TmtSchema;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; } 

        // Labeling a sequence as invalid when it contains modifications that are not supported by the model seems better than removing unsupported
        // mods and sending a sequence without the required TMT/iTRAQ labels, which would likely lead to crashes.
        public Prosit2020iRTTMT(SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull, int maxNumberOfBatchesPerRequest = 500, int throttlingDelayInMilliseconds = 100)
        {
            ModHandlingMode = modHandlingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
        }

        /// <summary>
        /// Creates batched API requests formatted for the Prosit 2020 iRT TMT model.
        /// Each batch contains up to MaxBatchSize TMT-labeled peptide sequences for optimal API performance.
        /// </summary>
        /// <returns>
        /// List of request dictionaries compatible with Koina API format.
        /// Each request contains peptide sequences in UNIMOD format with preserved isobaric labeling.
        /// </returns>
        /// <remarks>
        /// Request structure follows Koina API specification for Prosit_2020_irt_TMT:
        /// - peptide_sequences: BYTES array containing UNIMOD-formatted sequences with TMT labels
        /// - Shape: [batch_size, 1] for tensor compatibility
        /// - Datatype: BYTES for string sequence data with modification annotations
        /// 
        /// TMT-specific considerations:
        /// - Preserves N-terminal and lysine labeling information
        /// - Maintains isobaric tag consistency across batch
        /// - Handles various TMT/iTRAQ labeling formats uniformly
        /// - Enables accurate retention time prediction for labeled peptides
        /// 
        /// Batching strategy:
        /// - Splits input sequences into chunks of MaxBatchSize (1000)
        /// - Each batch gets a unique identifier for tracking
        /// - Optimized for concurrent processing of large TMT datasets
        /// </remarks>
        protected override List<Dictionary<string, object>> ToBatchedRequests(List<RetentionTimePredictionInput> validInputs)
        {
            // Split TMT-labeled sequences into batches for optimal API performance
            // ValidatedFullSequence should not be null at this point due to prior validation steps
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToList();
            var batchedRequests = new List<Dictionary<string, object>>();

            for (int i = 0; i < batchedPeptides.Count; i++)
            {
                // Create API request structure following Koina specification for TMT model
                var request = new Dictionary<string, object>
                {
                    { "id", $"Batch{i}_" + Guid.NewGuid()}, // Unique identifier for batch tracking
                    { "inputs", new List<object>
                        {
                            new {
                                name = "peptide_sequences",           // TMT model input parameter name
                                shape = new[]{ batchedPeptides[i].Length, 1 }, // Tensor shape [batch_size, 1]
                                datatype = "BYTES",                   // String data type for labeled sequences
                                data = batchedPeptides[i]            // TMT-labeled sequence data in UNIMOD format
                            },
                        }
                    }
                };
                batchedRequests.Add(request);
            }
            return batchedRequests;
        }

        protected override string? TryCleanSequence(string sequence, out string? apiSequence, out WarningException? warning)
        {
            var sanitized = base.TryCleanSequence(sequence, out apiSequence, out warning);
            if (sanitized == null || apiSequence == null)
            {
                return sanitized;
            }

            if (!HasAllowedNTerminalLabel(apiSequence))
            {
                var message = "Sequence must contain a supported N-terminal TMT/iTRAQ label.";
                switch (ModHandlingMode)
                {
                    case SequenceConversionHandlingMode.ThrowException:
                        throw new ArgumentException(message);
                    case SequenceConversionHandlingMode.ReturnNull:
                        warning = new WarningException(message);
                        return null;
                    case SequenceConversionHandlingMode.RemoveIncompatibleElements:
                    case SequenceConversionHandlingMode.UsePrimarySequence:
                        warning = new WarningException(message);
                        return null;
                }
            }

            return apiSequence;
        }

        private static bool HasAllowedNTerminalLabel(string apiSequence)
        {
            return apiSequence.StartsWith("[UNIMOD:737]-")
                   || apiSequence.StartsWith("[UNIMOD:2016]-")
                   || apiSequence.StartsWith("[UNIMOD:214]-")
                   || apiSequence.StartsWith("[UNIMOD:730]-");
        }
    }
}


