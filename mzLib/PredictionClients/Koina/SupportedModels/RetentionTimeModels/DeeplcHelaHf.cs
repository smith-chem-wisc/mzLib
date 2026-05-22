using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    /// <summary>
    /// Implementation of the DeepLC HeLa HF retention time prediction model.
    /// Uses deep learning with extensive feature engineering for accurate RT prediction.
    /// </summary>
    public class DeeplcHelaHf : RetentionTimeModel
    {
        private static readonly ISequenceConverter Converter = CreateUnimodConverterAcceptAll(
            UnimodSequenceFormatSchema.Instance);

        public override string ModelName => "Deeplc_hela_hf";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 2000;
        public override int MaxPeptideLength => 60;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => true;
        public override IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>(); // Accepts all UNIMOD modifications
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }

        public DeeplcHelaHf(
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
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToList();
            var batchedRequests = new List<Dictionary<string, object>>();

            for (int i = 0; i < batchedPeptides.Count; i++)
            {
                var request = new Dictionary<string, object>
                {
                    { "id", $"Batch{i}_" + Guid.NewGuid()},
                    { "inputs", new List<object>
                        {
                            new {
                                name = "peptide_sequences",
                                shape = new[]{ batchedPeptides[i].Length, 1 },
                                datatype = "BYTES",
                                data = batchedPeptides[i]
                            }
                        }
                    }
                };
                batchedRequests.Add(request);
            }
            return batchedRequests;
        }
    }
}
