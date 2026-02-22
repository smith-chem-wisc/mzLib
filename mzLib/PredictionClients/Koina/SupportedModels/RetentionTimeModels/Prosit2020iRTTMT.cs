using MzLibUtil;
using PredictionClients.Koina.AbstractClasses;
using System.ComponentModel;
using System.Text;
using System.Text.RegularExpressions;

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

        /// <summary>
        /// Comprehensive modification mapping supporting TMT, iTRAQ, SILAC, and standard modifications.
        /// Maps mzLib format to UNIMOD format for various isobaric labeling strategies.
        /// </summary>
        /// <remarks>
        /// Supported modification categories:
        /// - Standard: Oxidation on methionine, carbamidomethylation on cysteine
        /// - SILAC: Heavy lysine (13C6 15N2) and arginine (13C6 15N4) labeling
        /// - TMT: TMT6plex and TMTpro labeling on lysine and N-terminus
        /// - iTRAQ: 4-plex and 8-plex labeling on lysine and N-terminus
        /// 
        /// N-terminal modifications are denoted with "-" suffix in UNIMOD format.
        /// </remarks>
        public override Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            // Standard modifications
            {"[Common Variable:Oxidation on M]", "[UNIMOD:35]"},
            {"[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]"},
            
            // SILAC modifications for quantitative proteomics
            {"[Common Variable:Label:13C(6)15N(2) on K]", "[UNIMOD:259]"},
            {"[Common Variable:Label:13C(6)15N(4) on R]", "[UNIMOD:267]"},
            
            // TMT6plex modifications (6-channel multiplexing)
            {"[Common Fixed:TMT6plex on K]", "[UNIMOD:737]"},
            {"[Common Fixed:TMT6plex on N-terminus]", "[UNIMOD:737]-"},
            
            // TMTpro modifications (16+ channel multiplexing)
            {"[Common Fixed:TMTpro on K]", "[UNIMOD:2016]"},
            {"[Common Fixed:TMTpro on N-terminus]", "[UNIMOD:2016]-"},
            
            // iTRAQ 4-plex modifications (4-channel multiplexing)
            {"[Common Fixed:iTRAQ4plex on K]", "[UNIMOD:214]"},
            {"[Common Fixed:iTRAQ4plex on N-terminus]", "[UNIMOD:214]-"},
            
            // iTRAQ 8-plex modifications (8-channel multiplexing)
            {"[Common Fixed:iTRAQ8plex on K]", "[UNIMOD:730]"},
            {"[Common Fixed:iTRAQ8plex on N-terminus]", "[UNIMOD:730]-"}
        };

        public Prosit2020iRTTMT(int maxNumberOfBatchesPerRequest = 500, int throttlingDelayInMilliseconds = 100)
        {
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
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence).Chunk(MaxBatchSize).ToList();
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

        /// <summary>
        /// Validates modifications in TMT-labeled peptide sequences with enhanced N-terminal labeling validation.
        /// Overrides base method to provide TMT-specific validation logic for isobaric labeling patterns.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential TMT/iTRAQ modifications in mzLib format</param>
        /// <returns>True if all modifications are valid and N-terminal labeling is properly applied; otherwise false</returns>
        /// <remarks>
        /// TMT-specific validation enhancements:
        /// 1. Validates all modifications against the comprehensive TMT modification mapping
        /// 2. Specifically checks N-terminal modification validity for isobaric labeling
        /// 3. Ensures proper labeling patterns for TMT/iTRAQ experiments
        /// 
        /// Validation process:
        /// - Uses base ModificationPattern regex to find all modifications
        /// - Extracts and validates N-terminal modification (if present at index 0)
        /// - Checks all modifications against ValidModificationUnimodMapping
        /// - Ensures N-terminal labeling consistency with overall labeling strategy
        /// 
        /// Common TMT labeling patterns validated:
        /// - N-terminal TMT6plex/TMTpro labeling
        /// - Lysine TMT6plex/TMTpro labeling  
        /// - N-terminal iTRAQ 4-plex/8-plex labeling
        /// - Lysine iTRAQ 4-plex/8-plex labeling
        /// </remarks>
        protected override bool HasValidModifications(string sequence)
        {
            // Find all modification annotations in the sequence
            var matches = Regex.Matches(sequence, ModificationPattern);

            // Extract N-terminal modification for special validation (TMT labeling often requires N-terminal tags)
            // nTermMod will be empty string if no N-terminal mod found in the sequence
            // firstModIsValid checks if the N-terminal modification is in the valid mapping
            var nTermMod = matches.FirstOrDefault(m => m.Index == 0)?.Value ?? string.Empty;
            var firstModIsValid = ValidModificationUnimodMapping.TryGetValue(nTermMod, out var _);

            // Validate all modifications are in the supported mapping AND N-terminal mod exists and is valid
            return matches.All(m => ValidModificationUnimodMapping.ContainsKey(m.Value))
                && firstModIsValid;
        }
    }
}


