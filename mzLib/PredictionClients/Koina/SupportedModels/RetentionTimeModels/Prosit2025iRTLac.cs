using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    /// <summary>
    /// Implementation of the Prosit 2025 lactylation iRT prediction model.
    /// Specialized for lactylated peptides.
    /// </summary>
    public class Prosit2025iRTLac : RetentionTimeModel
    {
        private static readonly UnimodSequenceFormatSchema LacSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        // Koina: ALPHABET_MOD in models/Prosit/Prosit_Preprocess_peptide_lac/1/sequence_conversion.py
        private static readonly IReadOnlySet<string> SupportedModificationTokens = new HashSet<string> { "M[UNIMOD:35]", "C[UNIMOD:4]", "K[UNIMOD:2114]", "[UNIMOD:2114]-" };
        private static readonly IReadOnlySet<int> SupportedUnimodIds = UnimodIdsOf(SupportedModificationTokens);
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(LacSchema, SupportedUnimodIds);

        public override string ModelName => "Prosit_2025_irt_lac";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => true;
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override IReadOnlySet<string>? AllowedModificationTokens => SupportedModificationTokens;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }

        public Prosit2025iRTLac(
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
