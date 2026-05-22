using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels
{
    /// <summary>
    /// Implementation of the Prosit 2024 crosslink intensity prediction model for NMS2 crosslinks.
    /// Takes two peptide sequences (alpha and beta chains) and predicts fragment intensities.
    /// </summary>
    public class Prosit2024IntensityXLNMS2 : CrosslinkFragmentIntensityModel
    {
        private static readonly IReadOnlySet<int> Nms2UnimodIds = new HashSet<int> { 4, 35, 1898 };
        private static readonly ISequenceConverter Nms2Converter = CreateUnimodConverter(CrosslinkSchema, Nms2UnimodIds);

        public override string ModelName => "Prosit_2024_intensity_XL_NMS2";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override int NumberOfPredictedFragmentIons => 174;
        public override IReadOnlySet<int> AllowedUnimodIds => Nms2UnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public Prosit2024IntensityXLNMS2(
            SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull,
            IncompatibleParameterHandlingMode parameterHandlingMode = IncompatibleParameterHandlingMode.ReturnNull,
            FragmentIonMappingMode fragmentIonMappingMode = FragmentIonMappingMode.MapToValidatedFullSequence,
            int maxNumberOfBatchesPerRequest = 250,
            int throttlingDelayInMilliseconds = 100)
            : base(Nms2Converter)
        {
            ModHandlingMode = modHandlingMode;
            ParameterHandlingMode = parameterHandlingMode;
            FragmentIonMappingMode = fragmentIonMappingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<CrosslinkIntensityPredictionInput> validInputs)
        {
            var batchedAlpha = validInputs.Select(p => p.ValidatedAlphaSequence!).Chunk(MaxBatchSize).ToArray();
            var batchedBeta = validInputs.Select(p => p.ValidatedBetaSequence!).Chunk(MaxBatchSize).ToArray();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToArray();
            var batchedEnergies = validInputs.Select(p => p.CollisionEnergy).Chunk(MaxBatchSize).ToArray();

            var batchedRequests = new List<Dictionary<string, object>>(batchedAlpha.Length);
            for (int i = 0; i < batchedAlpha.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences_1", "BYTES", batchedAlpha[i]),
                    new InputField("peptide_sequences_2", "BYTES", batchedBeta[i]),
                    new InputField("precursor_charges", "INT32", batchedCharges[i]),
                    new InputField("collision_energies", "FP32", batchedEnergies[i])));
            }
            return batchedRequests;
        }
    }
}
