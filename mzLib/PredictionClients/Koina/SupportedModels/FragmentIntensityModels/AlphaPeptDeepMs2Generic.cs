using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// Implementation of the AlphaPeptDeep MS2 generic fragment intensity prediction model.
    /// Predicts fragment ion intensities with dynamic ion count based on peptide properties.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with variable length
    /// - Handles precursor charges 1-6
    /// - Predicts dynamic number of fragment ions
    /// - Requires instrument type input
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/AlphaPeptDeep_ms2_generic/infer
    /// </remarks>
    public class AlphaPeptDeepMs2Generic : FragmentIntensityModel
    {
        private static readonly ISequenceConverter Converter = CreateUnimodConverterAcceptAll(
            UnimodSequenceFormatSchema.Instance);

        public override string ModelName => "AlphaPeptDeep_ms2_generic";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
        public override int MaxPeptideLength => 500; // Model has no upper limit; using 500 as practical bound to simplify validation
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override HashSet<int>? AllowedCollisionEnergies => new HashSet<int>(); // Koina accepts any FP32 collision energy
        public override HashSet<string>? AllowedInstrumentTypes => new() { "QE", "LUMOS", "TIMSTOF", "SCIEXTOF" };
        public override IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>(); // Accepts all UNIMOD modifications
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public AlphaPeptDeepMs2Generic(
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
            var batchedEnergies = validInputs.Select(p => p.CollisionEnergy).Chunk(MaxBatchSize).ToArray();
            var batchedInstTypes = validInputs.Select(p => p.InstrumentType ?? "QE").Chunk(MaxBatchSize).ToArray();

            var batchedRequests = new List<Dictionary<string, object>>(batchedPeptides.Length);
            for (int i = 0; i < batchedPeptides.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences", "BYTES", batchedPeptides[i]),
                    new InputField("precursor_charges", "INT32", batchedCharges[i]),
                    new InputField("collision_energies", "FP32", batchedEnergies[i]),
                    new InputField("instrument_types", "BYTES", batchedInstTypes[i])));
            }
            return batchedRequests;
        }
    }
}
