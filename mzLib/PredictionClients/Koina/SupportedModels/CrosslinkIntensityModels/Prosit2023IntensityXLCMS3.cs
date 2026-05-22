using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels
{
    /// <summary>
    /// Implementation of the Prosit 2023 crosslink intensity prediction model for CMS3 crosslinks.
    /// Takes a single combined peptide sequence (preprocessed internally) and predicts fragment intensities.
    /// Uses fixed NCE of 35 (no collision energy input required).
    /// </summary>
    public class Prosit2023IntensityXLCMS3 : CrosslinkFragmentIntensityModel
    {
        private static readonly IReadOnlySet<int> Cms3UnimodIds = new HashSet<int> { 4, 35, 1881 };
        private static readonly ISequenceConverter Cms3Converter = CreateUnimodConverter(CrosslinkSchema, Cms3UnimodIds);

        public override string ModelName => "Prosit_2023_intensity_XL_CMS3";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override HashSet<int> AllowedCollisionEnergies => new HashSet<int>(); // Fixed NCE=35, no CE input
        // CMS3 takes a single combined sequence in the alpha slot; beta is unused.
        public override bool RequiresBetaSequence => false;
        public override int NumberOfPredictedFragmentIons => 174;
        public override IReadOnlySet<int> AllowedUnimodIds => Cms3UnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public Prosit2023IntensityXLCMS3(
            SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull,
            IncompatibleParameterHandlingMode parameterHandlingMode = IncompatibleParameterHandlingMode.ReturnNull,
            FragmentIonMappingMode fragmentIonMappingMode = FragmentIonMappingMode.MapToValidatedFullSequence,
            int maxNumberOfBatchesPerRequest = 250,
            int throttlingDelayInMilliseconds = 100)
            : base(Cms3Converter)
        {
            ModHandlingMode = modHandlingMode;
            ParameterHandlingMode = parameterHandlingMode;
            FragmentIonMappingMode = fragmentIonMappingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<CrosslinkIntensityPredictionInput> validInputs)
        {
            var batchedPeptides = validInputs.Select(p => p.ValidatedAlphaSequence!).Chunk(MaxBatchSize).ToList();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToList();

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
                            },
                            new {
                                name = "precursor_charges",
                                shape = new[]{ batchedCharges[i].Length, 1 },
                                datatype = "INT32",
                                data = batchedCharges[i]
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
