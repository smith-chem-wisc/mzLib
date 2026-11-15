using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    /// <summary>
    /// Implementation of the Chronologer retention time prediction model.
    /// Predicts absolute retention times (not iRT) using a TorchSharp-based model.
    /// </summary>
    public class ChronologerRT : RetentionTimeModel
    {
        // Koina: amino_acid_unimod_map and n_terminal_mod_map in models/Chronologer/Chronologer_Preprocess/1/sequence_conversion.py
        // Left out: K[UNIMOD:125], which mzLib's Unimod lacks, and [UNIMOD:385]-, which Koina's table lists only before its own C[UNIMOD:385] token.
        private static readonly IReadOnlySet<string> SupportedModificationTokens = new HashSet<string>
        {
            "C[UNIMOD:4]", "M[UNIMOD:35]", "S[UNIMOD:21]", "T[UNIMOD:21]", "Y[UNIMOD:21]",
            "K[UNIMOD:1]", "K[UNIMOD:121]", "K[UNIMOD:34]", "K[UNIMOD:36]", "K[UNIMOD:37]", "R[UNIMOD:34]", "R[UNIMOD:36]",
            "K[UNIMOD:739]", "K[UNIMOD:737]", "C[UNIMOD:385]", "E[UNIMOD:27]", "Q[UNIMOD:28]",
            "[UNIMOD:1]-", "[UNIMOD:27]-", "[UNIMOD:28]-", "[UNIMOD:737]-", "[UNIMOD:739]-"
        };
        private static readonly IReadOnlySet<int> SupportedUnimodIds = UnimodIdsOf(SupportedModificationTokens);
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(
            UnimodSequenceFormatSchema.DashedNTerminalInstance, SupportedUnimodIds);

        public override string ModelName => "Chronologer_RT";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => false;
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override IReadOnlySet<string>? AllowedModificationTokens => SupportedModificationTokens;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }

        public ChronologerRT(
            SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull,
            int maxNumberOfBatchesPerRequest = 500,
            int throttlingDelayInMilliseconds = 100)
            : base(Converter)
        {
            ModHandlingMode = modHandlingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<RetentionTimePredictionInput> validInputs)
        {
            var batchedPeptides = validInputs.Select(p => GetKoinaSequence(p)).Chunk(MaxBatchSize).ToArray();
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
