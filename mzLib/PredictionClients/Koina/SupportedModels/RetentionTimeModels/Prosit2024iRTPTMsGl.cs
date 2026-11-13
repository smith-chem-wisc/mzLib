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
        // Koina runs the 2025 PTM models' preprocessing for these, which places any UNIMOD id it can encode.
        private static readonly ISequenceConverter Converter = CreateUnimodConverterAcceptAll(PTMSchema);

        public override string ModelName => "Prosit_2024_irt_PTMs_gl";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => true;
        public override bool AcceptsAllUnimodModifications => true;
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
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToArray();
            var batchedRequests = new List<Dictionary<string, object>>(batchedPeptides.Length);
            for (int i = 0; i < batchedPeptides.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences", "BYTES", batchedPeptides[i])));
            }
            return batchedRequests;
        }

        /// <summary>
        /// Refuses the pairs Koina cannot place; see <see cref="PrositAtomCountModificationRules"/>.
        /// </summary>
        protected override string? DecideModification(CanonicalSequence sequence, CanonicalModification modification)
            => PrositAtomCountModificationRules.Decide(sequence, modification);
    }
}
