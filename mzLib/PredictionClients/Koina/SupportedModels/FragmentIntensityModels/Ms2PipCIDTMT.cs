using System.ComponentModel;
using MzLibUtil;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// MS2PIP CID TMT intensity prediction model.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Handles precursor charges 1-6
    /// - Predicts up to 58 fragment ions per peptide
    /// - Requires N-terminal TMT/iTRAQ labeling
    /// - Supports TMT6plex, TMTpro, iTRAQ4/8plex
    /// </remarks>
    public class Ms2PipCIDTMT : FragmentIntensityModel
    {
        private static readonly UnimodSequenceFormatSchema TmtSchema = new(UnimodLabelStyle.UpperCase, '[', ']', "-", "-");
        private static readonly ISequenceConverter Converter = CreateUnimodConverterAcceptAll(TmtSchema);

        public override string ModelName => "ms2pip_CID_TMT";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 750;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override HashSet<int>? AllowedCollisionEnergies => null; // Fixed CID, no CE input
        public override IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>(); // Accepts all UNIMOD (mods only affect m/z)
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }
        public override int NumberOfPredictedFragmentIons => 58;

        public Ms2PipCIDTMT(
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

            var batchedRequests = new List<Dictionary<string, object>>(batchedPeptides.Length);
            for (int i = 0; i < batchedPeptides.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences", "BYTES", batchedPeptides[i]),
                    new InputField("precursor_charges", "INT32", batchedCharges[i])));
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

            return sanitized;
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
