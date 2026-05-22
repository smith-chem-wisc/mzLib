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
        private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int> { 35, 4 };
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(
            UnimodSequenceFormatSchema.Instance, SupportedUnimodIds);

        public override string ModelName => "Chronologer_RT";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => false;
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
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
