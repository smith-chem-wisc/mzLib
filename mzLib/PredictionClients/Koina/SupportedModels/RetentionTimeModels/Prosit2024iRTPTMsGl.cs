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
        // Koina: the keys of the two atom-count dictionaries in models/Prosit/Prosit_Preprocess_ac_gain/1/model.py and
        // Prosit_Preprocess_ac_loss/1/model.py, where "K_737" is K[UNIMOD:737] and "_737" is the N-terminal [UNIMOD:737]-.
        // Four of those keys are left out because Koina answers 400 to them: R_267 has an empty atom count the
        // preprocessing cannot read, and K_12118, K_19903 and K_129317 are ids its UNIMOD table does not have.
        private static readonly IReadOnlySet<string> SupportedModificationTokens = new HashSet<string>
        {
            "K[UNIMOD:1]", "C[UNIMOD:4]", "K[UNIMOD:4]", "N[UNIMOD:7]", "Q[UNIMOD:7]", "R[UNIMOD:7]",
            "H[UNIMOD:21]", "P[UNIMOD:21]", "S[UNIMOD:21]", "T[UNIMOD:21]", "Y[UNIMOD:21]", "E[UNIMOD:27]",
            "Q[UNIMOD:28]", "C[UNIMOD:34]", "D[UNIMOD:34]", "E[UNIMOD:34]", "H[UNIMOD:34]", "I[UNIMOD:34]",
            "K[UNIMOD:34]", "L[UNIMOD:34]", "N[UNIMOD:34]", "Q[UNIMOD:34]", "R[UNIMOD:34]", "C[UNIMOD:35]",
            "H[UNIMOD:35]", "K[UNIMOD:35]", "M[UNIMOD:35]", "P[UNIMOD:35]", "W[UNIMOD:35]", "K[UNIMOD:36]",
            "R[UNIMOD:36]", "K[UNIMOD:37]", "S[UNIMOD:43]", "T[UNIMOD:43]", "K[UNIMOD:56]", "K[UNIMOD:58]",
            "K[UNIMOD:59]", "K[UNIMOD:121]", "K[UNIMOD:214]", "C[UNIMOD:312]", "K[UNIMOD:535]",
            "K[UNIMOD:730]", "K[UNIMOD:737]", "C[UNIMOD:1263]", "K[UNIMOD:1263]", "K[UNIMOD:1289]",
            "K[UNIMOD:1293]", "K[UNIMOD:1848]", "K[UNIMOD:1990]", "K[UNIMOD:2016]", "C[UNIMOD:2062]",
            "K[UNIMOD:5634]",
            "[UNIMOD:1]-", "[UNIMOD:58]-", "[UNIMOD:59]-", "[UNIMOD:214]-", "[UNIMOD:411]-", "[UNIMOD:730]-",
            "[UNIMOD:737]-", "[UNIMOD:2016]-"
        };
        private static readonly IReadOnlySet<int> SupportedUnimodIds = UnimodIdsOf(SupportedModificationTokens);
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
        public override IReadOnlySet<string>? AllowedModificationTokens => SupportedModificationTokens;
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
