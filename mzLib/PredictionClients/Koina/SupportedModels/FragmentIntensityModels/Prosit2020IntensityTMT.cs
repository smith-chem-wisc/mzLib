using MzLibUtil;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// Implementation of the Prosit 2020 TMT intensity prediction model.
    /// Predicts fragment ion intensities for TMT-labeled peptides using HCD fragmentation.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Handles precursor charges 1-6
    /// - Predicts up to 174 fragment ions per peptide
    /// - Requires N-terminal TMT/iTRAQ labeling, so UsePrimarySequence is rejected at construction
    /// - Supports TMT6plex, TMTpro, iTRAQ4/8plex, SILAC, oxidation, carbamidomethyl
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_TMT/infer
    /// </remarks>
    public class Prosit2020IntensityTMT : FragmentIntensityModel
    {
        private static readonly UnimodSequenceFormatSchema TmtSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        // Koina: ALPHABET_MOD in models/Prosit/Prosit_Preprocess_peptide_2020_TMT/1/sequence_conversion.py
        private static readonly IReadOnlySet<string> SupportedModificationTokens = new HashSet<string>
        {
            "M[UNIMOD:35]", "C[UNIMOD:4]", "K[UNIMOD:259]", "R[UNIMOD:267]",
            "K[UNIMOD:737]", "K[UNIMOD:2016]", "K[UNIMOD:214]", "K[UNIMOD:730]",
            "[UNIMOD:737]-", "[UNIMOD:2016]-", "[UNIMOD:214]-", "[UNIMOD:730]-"
        };
        private static readonly IReadOnlySet<int> SupportedUnimodIds = UnimodIdsOf(SupportedModificationTokens);
        private static readonly IReadOnlySet<int> NTerminalLabelIds = UnimodIdsOf(SupportedModificationTokens.Where(t => t.EndsWith('-')));
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(TmtSchema, SupportedUnimodIds);

        public override string ModelName => "Prosit_2020_intensity_TMT";
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
        public override IReadOnlySet<int>? RequiredNTerminalUnimodIds => NTerminalLabelIds;
        private readonly SequenceConversionHandlingMode _modHandlingMode;
        public override SequenceConversionHandlingMode ModHandlingMode
        {
            get => _modHandlingMode;
            init => _modHandlingMode = value == SequenceConversionHandlingMode.UsePrimarySequence
                ? throw new ArgumentException($"{ModelName} requires an N-terminal TMT/iTRAQ label, which UsePrimarySequence would strip from every sequence. Use ReturnNull, ThrowException or RemoveIncompatibleElements.", nameof(ModHandlingMode))
                : value;
        }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public Prosit2020IntensityTMT(
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
    }
}
