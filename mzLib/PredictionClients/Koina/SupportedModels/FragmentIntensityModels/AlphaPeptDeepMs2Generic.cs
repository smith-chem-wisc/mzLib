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
        public override HashSet<int> AllowedCollisionEnergies => new HashSet<int>(); // Koina accepts any FP32 collision energy
        public override HashSet<string> AllowedInstrumentTypes => new() { "QE", "LUMOS", "TIMSTOF", "SCIEXTOF" };
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
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToList();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToList();
            var batchedEnergies = validInputs.Select(p => p.CollisionEnergy).Chunk(MaxBatchSize).ToList();
            var batchedInstTypes = validInputs.Select(p => p.InstrumentType ?? "QE").Chunk(MaxBatchSize).ToList();

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
                            },
                            new {
                                name = "collision_energies",
                                shape = new[]{ batchedEnergies[i].Length, 1 },
                                datatype = "FP32",
                                data = batchedEnergies[i]
                            },
                            new {
                                name = "instrument_types",
                                shape = new[]{ batchedInstTypes[i].Length, 1 },
                                datatype = "BYTES",
                                data = batchedInstTypes[i]
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
