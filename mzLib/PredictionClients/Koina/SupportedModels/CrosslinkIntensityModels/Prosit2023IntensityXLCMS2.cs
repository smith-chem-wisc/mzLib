using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels
{
    /// <summary>
    /// Implementation of the Prosit 2023 crosslink intensity prediction model for CMS2 crosslinks.
    /// Takes two peptide sequences (alpha and beta chains) and predicts fragment intensities.
    /// </summary>
    public class Prosit2023IntensityXLCMS2 : CrosslinkFragmentIntensityModel
    {
        private static readonly IReadOnlySet<int> Cms2UnimodIds = new HashSet<int> { 4, 35, 1896 };
        private static readonly ISequenceConverter Cms2Converter = CreateUnimodConverter(CrosslinkSchema, Cms2UnimodIds);

        public override string ModelName => "Prosit_2023_intensity_XL_CMS2";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override int NumberOfPredictedFragmentIons => 348;
        public override IReadOnlySet<int> AllowedUnimodIds => Cms2UnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public Prosit2023IntensityXLCMS2(
            SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull,
            IncompatibleParameterHandlingMode parameterHandlingMode = IncompatibleParameterHandlingMode.ReturnNull,
            FragmentIonMappingMode fragmentIonMappingMode = FragmentIonMappingMode.MapToValidatedFullSequence,
            int maxNumberOfBatchesPerRequest = 250,
            int throttlingDelayInMilliseconds = 100)
            : base(Cms2Converter)
        {
            ModHandlingMode = modHandlingMode;
            ParameterHandlingMode = parameterHandlingMode;
            FragmentIonMappingMode = fragmentIonMappingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<CrosslinkIntensityPredictionInput> validInputs)
        {
            var batchedAlpha = validInputs.Select(p => p.ValidatedAlphaSequence!).Chunk(MaxBatchSize).ToList();
            var batchedBeta = validInputs.Select(p => p.ValidatedBetaSequence!).Chunk(MaxBatchSize).ToList();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToList();
            var batchedEnergies = validInputs.Select(p => p.CollisionEnergy).Chunk(MaxBatchSize).ToList();

            var batchedRequests = new List<Dictionary<string, object>>();

            for (int i = 0; i < batchedAlpha.Count; i++)
            {
                var request = new Dictionary<string, object>
                {
                    { "id", $"Batch{i}_" + Guid.NewGuid()},
                    { "inputs", new List<object>
                        {
                            new {
                                name = "peptide_sequences_1",
                                shape = new[]{ batchedAlpha[i].Length, 1 },
                                datatype = "BYTES",
                                data = batchedAlpha[i]
                            },
                            new {
                                name = "peptide_sequences_2",
                                shape = new[]{ batchedBeta[i].Length, 1 },
                                datatype = "BYTES",
                                data = batchedBeta[i]
                            },
                            new {
                                name = "precursor_charges",
                                shape = new[]{ batchedCharges[i].Length, 1 },
                                datatype = "INT32",
                                data = batchedCharges[i]
                            },
                            new {
                                name = "collision_energies",
                                shape = new[]{ batchedEnergies[i].Length, 1 },
                                datatype = "FP32",
                                data = batchedEnergies[i]
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
