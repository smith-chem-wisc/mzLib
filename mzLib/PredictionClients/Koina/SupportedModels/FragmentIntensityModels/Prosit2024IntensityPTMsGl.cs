using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// Implementation of the Prosit 2024 generalized PTM intensity prediction model.
    /// Supports a wide range of post-translational modifications for fragment intensity prediction.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Handles precursor charges 1-6
    /// - Predicts up to 174 fragment ions per peptide
    /// - Supports many PTM types via generalized representation
    /// - Requires fragmentation type input
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Prosit_2024_intensity_PTMs_gl/infer
    /// </remarks>
    public class Prosit2024IntensityPTMsGl : FragmentIntensityModel
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

        public override string ModelName => "Prosit_2024_intensity_PTMs_gl";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override HashSet<int>? AllowedCollisionEnergies => new HashSet<int>(); // Koina accepts any FP32 collision energy
        public override HashSet<string>? AllowedFragmentationTypes => new() { "HCD", "CID" };
        public override int NumberOfPredictedFragmentIons => 174;
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override IReadOnlySet<string>? AllowedModificationTokens => SupportedModificationTokens;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public Prosit2024IntensityPTMsGl(
            SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull,
            IncompatibleParameterHandlingMode parameterHandlingMode = IncompatibleParameterHandlingMode.ReturnNull,
            FragmentIonMappingMode fragmentIonMappingMode = FragmentIonMappingMode.MapToValidatedFullSequence,
            int maxNumberOfBatchesPerRequest = 250,
            int throttlingDelayInMilliseconds = 100)
            : base(Converter)
        {
            ModHandlingMode = modHandlingMode;
            ParameterHandlingMode = parameterHandlingMode;
            FragmentIonMappingMode = fragmentIonMappingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
        {
            var batchedPeptides = validInputs.Select(p => GetKoinaSequence(p)).Chunk(MaxBatchSize).ToArray();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToArray();
            var batchedEnergies = validInputs.Select(p => (float)p.CollisionEnergy!).Chunk(MaxBatchSize).ToArray();
            var batchedFragTypes = validInputs.Select(p => p.FragmentationType ?? "HCD").Chunk(MaxBatchSize).ToArray();

            var batchedRequests = new List<Dictionary<string, object>>(batchedPeptides.Length);
            for (int i = 0; i < batchedPeptides.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences", "BYTES", batchedPeptides[i]),
                    new InputField("precursor_charges", "INT32", batchedCharges[i]),
                    new InputField("collision_energies", "FP32", batchedEnergies[i]),
                    new InputField("fragmentation_types", "BYTES", batchedFragTypes[i])));
            }
            return batchedRequests;
        }
    }
}
