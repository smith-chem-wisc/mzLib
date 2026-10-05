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
    /// - Specialized for TMT (Tandem Mass Tag) and iTRAQ labeled peptides; requires an N-terminal label,
    ///   so UsePrimarySequence is rejected at construction
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
        private static readonly UnimodSequenceFormatSchema TmtSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        // Koina: ALPHABET_MOD in models/Prosit/Prosit_Preprocess_peptide_2020_TMT/1/sequence_conversion.py
        private static readonly IReadOnlySet<string> SupportedModificationTokens = new HashSet<string>
        {
            "M[UNIMOD:35]", "C[UNIMOD:4]", "K[UNIMOD:259]", "R[UNIMOD:267]",
            "K[UNIMOD:737]", "K[UNIMOD:2016]", "K[UNIMOD:214]", "K[UNIMOD:730]",
            "[UNIMOD:737]-", "[UNIMOD:2016]-", "[UNIMOD:214]-", "[UNIMOD:730]-"
        };
        private static readonly IReadOnlySet<int> SupportedUnimodIds = UnimodIdsOf(SupportedModificationTokens);
        private static readonly IReadOnlySet<int> NTerminalLabelIds = UnimodIdsOf(SupportedModificationTokens.Where(t => t.EndsWith('-')));
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(TmtSchema, SupportedUnimodIds);

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

        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override IReadOnlySet<string>? AllowedModificationTokens => SupportedModificationTokens;
        public override IReadOnlySet<int>? RequiredNTerminalUnimodIds => NTerminalLabelIds;
        private readonly SequenceConversionHandlingMode _modHandlingMode;
        public override SequenceConversionHandlingMode ModHandlingMode
        {
            get => _modHandlingMode;
            init => _modHandlingMode = value == SequenceConversionHandlingMode.UsePrimarySequence
                ? throw new ArgumentException($"{ModelName} requires an N-terminal TMT/iTRAQ label, which UsePrimarySequence would strip from every sequence. Use ReturnNull, ThrowException or RemoveIncompatibleElements.", nameof(ModHandlingMode))
                : value;
        }

        // Labeling a sequence as invalid when it contains modifications that are not supported by the model seems better than removing unsupported
        // mods and sending a sequence without the required TMT/iTRAQ labels, which would likely lead to crashes.
        public Prosit2020iRTTMT(SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull, int maxNumberOfBatchesPerRequest = 500, int throttlingDelayInMilliseconds = 100)
            : base(Converter)
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
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToArray();
            var batchedRequests = new List<Dictionary<string, object>>(batchedPeptides.Length);
            for (int i = 0; i < batchedPeptides.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences", "BYTES", batchedPeptides[i])));
            }
            return batchedRequests;
        }
    }
}


