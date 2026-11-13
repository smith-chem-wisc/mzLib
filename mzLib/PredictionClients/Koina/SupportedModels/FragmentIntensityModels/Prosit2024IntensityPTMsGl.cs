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
        // Koina runs the 2025 PTM models' preprocessing for these, which places any UNIMOD id it can encode.
        private static readonly ISequenceConverter Converter = CreateUnimodConverterAcceptAll(PTMSchema);

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
        public override bool AcceptsAllUnimodModifications => true;
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
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToArray();
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

        /// <summary>
        /// Refuses the pairs Koina cannot place; see <see cref="PrositAtomCountModificationRules"/>.
        /// </summary>
        protected override string? DecideModification(CanonicalSequence sequence, CanonicalModification modification)
            => PrositAtomCountModificationRules.Decide(sequence, modification);
    }
}
