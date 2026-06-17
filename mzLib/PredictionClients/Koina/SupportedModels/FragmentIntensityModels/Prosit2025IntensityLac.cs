using System.ComponentModel;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// Implementation of the Prosit 2025 lactylation intensity prediction model.
    /// Specialized for predicting fragment ion intensities of lactylated peptides.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Handles precursor charges 1-6
    /// - Predicts up to 174 fragment ions per peptide
    /// - Specialized for lactylation
    /// - Requires fragmentation type and instrument type inputs
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Prosit_2025_intensity_lac/infer
    /// </remarks>
    public class Prosit2025IntensityLac : FragmentIntensityModel
    {
        private static readonly UnimodSequenceFormatSchema LacSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int> { 35, 4, 2114 };
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(LacSchema, SupportedUnimodIds);

        public override string ModelName => "Prosit_2025_intensity_lac";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override HashSet<int>? AllowedCollisionEnergies => new HashSet<int>(); // Koina accepts any FP32 collision energy
        public override HashSet<string>? AllowedFragmentationTypes => new() { "HCD", "CID" };
        // Koina's Prosit_Preprocess_instrument_types matches these byte-for-byte (uppercase) and
        // silently defaults anything else to LUMOS, so the client must send them uppercased.
        public override HashSet<string>? AllowedInstrumentTypes => new() { "ECLIPSE", "ASTRAL", "LUMOS" };
        public override int NumberOfPredictedFragmentIons => 174;
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public Prosit2025IntensityLac(
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

        protected override bool ValidateModelSpecificInputs(FragmentIntensityPredictionInput input, out WarningException? warning)
        {
            // Accept instrument types case-insensitively; the request is sent uppercased.
            input = input with { InstrumentType = input.InstrumentType?.ToUpperInvariant() };
            return base.ValidateModelSpecificInputs(input, out warning);
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
        {
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToArray();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToArray();
            var batchedEnergies = validInputs.Select(p => (float)p.CollisionEnergy!).Chunk(MaxBatchSize).ToArray();
            var batchedFragTypes = validInputs.Select(p => p.FragmentationType ?? "HCD").Chunk(MaxBatchSize).ToArray();
            var batchedInstTypes = validInputs.Select(p => (p.InstrumentType ?? "LUMOS").ToUpperInvariant()).Chunk(MaxBatchSize).ToArray();

            var batchedRequests = new List<Dictionary<string, object>>(batchedPeptides.Length);
            for (int i = 0; i < batchedPeptides.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences", "BYTES", batchedPeptides[i]),
                    new InputField("precursor_charges", "INT32", batchedCharges[i]),
                    new InputField("collision_energies", "FP32", batchedEnergies[i]),
                    new InputField("fragmentation_types", "BYTES", batchedFragTypes[i]),
                    new InputField("instrument_types", "BYTES", batchedInstTypes[i])));
            }
            return batchedRequests;
        }
    }
}
