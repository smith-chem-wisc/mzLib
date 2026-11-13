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
        // Koina's Chronologer reads an N-terminal modification only in the dashed form, [UNIMOD:1]-PEPTIDE.
        private static readonly UnimodSequenceFormatSchema NTermDashSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        private static readonly ISequenceConverter Converter = CreateUnimodConverterAcceptAll(NTermDashSchema);

        public override string ModelName => "Chronologer_RT";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => false;
        public override IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>();
        public override bool AcceptsAllUnimodModifications => true;
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
