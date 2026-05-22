using System.ComponentModel;
using System.Text.RegularExpressions;
using Chemistry;
using MzLibUtil;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// Implementation of the UniSpec fragment intensity prediction model.
    /// Predicts fragment ion intensities with dynamic ion count based on peptide properties.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides up to length 40
    /// - Handles precursor charges 1-8 (with instrument-dependent narrowing: QE/QEHFX/ELITE 2-5, VELOS 2-4)
    /// - Predicts dynamic number of fragment ions
    /// - Requires instrument type input (case-insensitive on the client; normalised to upper before send)
    ///
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/UniSpec/infer
    /// </remarks>
    public class UniSpec : FragmentIntensityModel
    {
        private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int> { 1, 4, 21, 26, 27, 28, 35 };
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(
            UnimodSequenceFormatSchema.Instance, SupportedUnimodIds);

        public override string ModelName => "UniSpec";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
        public override int MaxPeptideLength => 40;
        public override int MinPeptideLength => 1;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6, 7, 8};
        public override HashSet<int> AllowedCollisionEnergies => new HashSet<int>(); // Koina accepts any FP32 collision energy
        public override HashSet<string> AllowedInstrumentTypes => new() { "QE", "QEHFX", "LUMOS", "ELITE", "VELOS", "NONE"};
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        public UniSpec(
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
            // Base validation guarantees PrecursorCharge is in AllowedPrecursorCharges,
            // InstrumentType is non-null, and InstrumentType is in AllowedInstrumentTypes (exact case).
            if (!base.ValidateModelSpecificInputs(input, out warning))
                return false;

            var instrument = input.InstrumentType!;
            var charge = input.PrecursorCharge;

            bool valid = instrument switch
            {
                "QE" or "QEHFX" or "ELITE" => charge >= 2 && charge <= 5,
                "VELOS" => charge >= 2 && charge <= 4,
                "NONE" or "LUMOS" => true,
                _ => true
            };

            if (!valid)
            {
                string message = $"Charge {charge} is not supported for instrument '{instrument}'. " +
                    "QE/QEHFX/ELITE support charges 2-5; VELOS supports charges 2-4. " +
                    "Use instrument 'NONE' or 'LUMOS' to bypass instrument-specific charge constraints.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(message);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(message);
                        return false;
                }
            }

            return true;
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
        {
            var batchedPeptides = validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize).ToList();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToList();
            var batchedEnergies = validInputs.Select(p => p.CollisionEnergy).Chunk(MaxBatchSize).ToList();
            var batchedInstTypes = validInputs.Select(p => p.InstrumentType!).Chunk(MaxBatchSize).ToList();

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

        protected override ParsedFragmentAnnotation ParseFragmentAnnotation(string annotation)
        {
            if (string.IsNullOrEmpty(annotation))
            {
                throw new Exception($"UniSpec model returned empty fragment annotation.");
            }

            var plusIndex = annotation.LastIndexOf('+');
            if (plusIndex > 0)
            {
                var fragmentIdentifier = annotation.Substring(0, plusIndex);
                var chargeStr = annotation.Substring(plusIndex + 1);
                if (!int.TryParse(chargeStr, out int charge))
                {
                    throw new Exception($"Cannot parse charge from UniSpec fragment annotation '{annotation}'. Charge value: '{chargeStr}'.");
                }

                if (TryExtractNeutralLoss(fragmentIdentifier, out var nlFormula, out var nlMass))
                {
                    return new ParsedFragmentAnnotation(fragmentIdentifier, charge, nlFormula, nlMass);
                }

                return new ParsedFragmentAnnotation(fragmentIdentifier, charge);
            }

            if (Regex.IsMatch(annotation, @"^[a-zA-Z]+\d+$"))
            {
                // Check for neutral loss in charge-1 annotations (though rare, handle for consistency)
                if (TryExtractNeutralLoss(annotation, out var nlFormula, out var nlMass))
                {
                    return new ParsedFragmentAnnotation(annotation, 1, nlFormula, nlMass);
                }
                return new ParsedFragmentAnnotation(annotation, 1);
            }

            throw new Exception($"UniSpec fragment annotation '{annotation}' is not a standard fragment ion (e.g., internal fragment or neutral loss) and cannot be mapped to theoretical products. Expected format: 'type+charge' (e.g., 'y8+1') or 'type' (e.g., 'a2').");
        }

        private static bool TryExtractNeutralLoss(string fragmentIdentifier, out string ? nlFormula, out double ? nlMass)
        {
            nlFormula = null;
            nlMass = null;

            var dashIndex = fragmentIdentifier.IndexOf('-');
            if (dashIndex <= 0)
            {
                return false;
            }

            var formula = fragmentIdentifier.Substring(dashIndex + 1);
            try
            {
                nlMass = ChemicalFormula.ParseFormula(formula).MonoisotopicMass;
                nlFormula = formula;
                return true;
            }
            catch
            {
                return false;
            }
        }
    }
}
