using Easy.Common.Extensions;
using MzLibUtil;
using PredictionClients.Koina.AbstractClasses;
using System.ComponentModel;
using System.Text;


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
        /// <summary>
        /// List of validated peptide sequences formatted for model prediction.
        /// Sequences are converted to UNIMOD format with automatic cysteine carbamidomethylation.
        /// </summary>
        public override List<string> PeptideSequences { get; } = new();
        /// <summary>
        /// Collection of iRT prediction results after inference completion.
        /// Each prediction contains the sequence and predicted iRT value.
        /// </summary>
        public override List<PeptideRTPrediction> Predictions { get; protected set; } = new();

        /// <summary>
        /// Initializes a new instance of the Prosit2019iRT model with input validation and filtering.
        /// Validates sequences against model requirements and converts valid sequences to UNIMOD format.
        /// </summary>
        /// <param name="peptideSequences">
        /// List of peptide sequences in mzLib format for iRT prediction.
        /// Valid modifications: oxidation on methionine, carbamidomethylation on cysteine.
        /// Unmodified cysteines will be automatically carbamidomethylated.
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
        /// - Sequence length: 1-30 amino acids
        /// - Only canonical amino acids (20 standard)
        /// - Only supported modifications (oxidation on M, carbamidomethyl on C)
        /// - Empty or null input results in warning
        /// 
        /// Processing steps:
        /// 1. Validates input sequences against model constraints
        /// 2. Converts valid sequences from mzLib to UNIMOD format
        /// 3. Collects invalid sequences for warning message
        /// </remarks>
        /// <example>
        /// <code>
        /// var sequences = new List&lt;string&gt; 
        /// { 
        ///     "PEPTIDER", 
        ///     "PEPTC[Common Fixed:Carbamidomethyl on C]IDE",
        ///     "PEPTM[Common Variable:Oxidation on M]IDE" 
        /// };
        /// 
        /// var model = new Prosit2019iRT(sequences, out var warnings);
        /// if (warnings != null)
        /// {
        ///     Console.WriteLine(warnings.Message);
        /// }
        /// 
        /// await model.RunInferenceAsync();
        /// foreach (var prediction in model.Predictions)
        /// {
        ///     Console.WriteLine($"Sequence: {prediction.PeptideSequence}, iRT: {prediction.PredictedRetentionTime}");
        /// }
        /// </code>
        /// </example>
        public Prosit2019iRT(List<string> peptideSequences, out WarningException? warnings)
        {
            // Handle empty input case early
            if (peptideSequences.IsNullOrEmpty())
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            // Validate each sequence and collect invalid ones for reporting
            var invalidSequences = new List<string>();
            foreach (var seq in peptideSequences)
            {
                // Apply comprehensive validation: sequence structure, length, and modifications
                if (!IsValidSequence(seq) || !HasValidModifications(seq))
                {
                    invalidSequences.Add(seq);
                }
                else
                {
                    // Convert to UNIMOD format
                    PeptideSequences.Add(ConvertMzLibModificationsToUnimod(seq));
                }
            }

            // Generate warning message for invalid sequences if any were found
            warnings = null;
            if (invalidSequences.Count > 0)
            {
                var sb = new StringBuilder();
                sb.AppendLine("The following peptide sequences were invalid and will be skipped:");
                sb.AppendLine($"Model requirements: Length 1-{MaxPeptideLength}, canonical amino acids only,");
                sb.AppendLine("supported modifications: oxidation on M, carbamidomethyl on C");
                sb.AppendLine();
                foreach (var invalid in invalidSequences)
                {
                    sb.AppendLine($"  - {invalid}");
                }
                warnings = new WarningException(sb.ToString());
            }
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
        protected override List<Dictionary<string, object>> ToBatchedRequests()
        {
            // Split sequences into batches for optimal API performance
            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
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
