using System.ComponentModel;
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
    /// - Requires N-terminal TMT/iTRAQ labeling
    /// - Supports TMT6plex, TMTpro, iTRAQ4/8plex, SILAC, oxidation, carbamidomethyl
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_TMT/infer
    /// </remarks>
    public class Prosit2020IntensityTMT : FragmentIntensityModel
    {
        private static readonly UnimodSequenceFormatSchema TmtSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int>
        {
            35, 4, 259, 267, 737, 2016, 214, 730
        };
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(TmtSchema, SupportedUnimodIds);

        public override string ModelName => "Prosit_2020_intensity_TMT";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override HashSet<int> AllowedCollisionEnergies => new HashSet<int>(); // Koina accepts any FP32 collision energy
        public override HashSet<string> AllowedFragmentationTypes => new() { "HCD", "CID" };
        public int NumberOfPredictedFragmentIons => 174;
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
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
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToList();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToList();
            var batchedEnergies = validInputs.Select(p => p.CollisionEnergy).Chunk(MaxBatchSize).ToList();
            var batchedFragTypes = validInputs.Select(p => p.FragmentationType ?? "HCD").Chunk(MaxBatchSize).ToList();

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
                                name = "fragmentation_types",
                                shape = new[]{ batchedFragTypes[i].Length, 1 },
                                datatype = "BYTES",
                                data = batchedFragTypes[i]
                            }
                        }
                    }
                };
                batchedRequests.Add(request);
            }

            return batchedRequests;
        }

        protected override string? TryCleanSequence(string sequence, out string? apiSequence, out WarningException? warning)
        {
            var sanitized = base.TryCleanSequence(sequence, out apiSequence, out warning);
            if (sanitized == null || apiSequence == null)
            {
                return sanitized;
            }

            if (!HasAllowedNTerminalLabel(apiSequence))
            {
                var message = "Sequence must contain a supported N-terminal TMT/iTRAQ label.";
                switch (ModHandlingMode)
                {
                    case SequenceConversionHandlingMode.ThrowException:
                        throw new ArgumentException(message);
                    case SequenceConversionHandlingMode.ReturnNull:
                        warning = new WarningException(message);
                        return null;
                    case SequenceConversionHandlingMode.RemoveIncompatibleElements:
                    case SequenceConversionHandlingMode.UsePrimarySequence:
                        warning = new WarningException(message);
                        return null;
                }
            }

            return apiSequence;
        }

        private static bool HasAllowedNTerminalLabel(string apiSequence)
        {
            return apiSequence.StartsWith("[UNIMOD:737]-")
                   || apiSequence.StartsWith("[UNIMOD:2016]-")
                   || apiSequence.StartsWith("[UNIMOD:214]-")
                   || apiSequence.StartsWith("[UNIMOD:730]-");
        }
    }
}
