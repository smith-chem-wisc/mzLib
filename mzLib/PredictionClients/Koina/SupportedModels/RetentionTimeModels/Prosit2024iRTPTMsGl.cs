using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    /// <summary>
    /// Implementation of the Prosit 2024 generalized PTM iRT prediction model.
    /// Supports a wide range of post-translational modifications.
    /// </summary>
    public class Prosit2024iRTPTMsGl : RetentionTimeModel
    {
        private static readonly UnimodSequenceFormatSchema PTMSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        // Source of truth: koina repo, models/Prosit/Prosit_Preprocess_ac_gain/1/model.py
        // and Prosit_Preprocess_ac_loss/1/model.py — union of UNIMOD ids referenced in both dicts.
        private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int>
        {
            1, 4, 7, 21, 27, 28, 34, 35, 36, 37, 43, 56, 58, 59, 121, 214, 267, 312, 411,
            535, 730, 737, 1263, 1289, 1293, 1848, 1990, 2016, 2062, 5634, 12118, 19903, 129317
        };
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(PTMSchema, SupportedUnimodIds);

        public override string ModelName => "Prosit_2024_irt_PTMs_gl";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => true;
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }

        public Prosit2024iRTPTMsGl(
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
