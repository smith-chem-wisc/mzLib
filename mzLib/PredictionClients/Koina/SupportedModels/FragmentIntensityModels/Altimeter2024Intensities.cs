using Chemistry;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// Altimeter 2024 intensity prediction model.
    /// Predicts fragment ion intensities using spline-based representation for B, Y, P, and immonium ions.
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Handles precursor charges 1-6
    /// - Predicts up to 380 fragment ions per peptide
    /// - Returns B, Y, P, and immonium ions
    /// - Configurable via parameters: return_b_ions, return_y_ions, return_p_ions, return_imm_ions
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Altimeter_2024_intensities/infer
    /// </remarks>
    public class Altimeter2024Intensities : FragmentIntensityModel
    {
        private static readonly IReadOnlySet<int> SupportedUnimodIds = new HashSet<int> { 35, 4 };
        private static readonly ISequenceConverter Converter = CreateUnimodConverter(
            UnimodSequenceFormatSchema.Instance, SupportedUnimodIds);

        public override string ModelName => "Altimeter_2024_intensities";
        public override int MaxBatchSize => 1000;
        public override int MaxNumberOfBatchesPerRequest { get; init; }
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
        public override int MaxPeptideLength => 40; // Koina docs: 6-40 (6-30 recommended)
        public override int MinPeptideLength => 6;
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6, 7 };
        public override HashSet<int>? AllowedCollisionEnergies => Enumerable.Range(20, 21).ToHashSet(); // NCE 20-40
        public override IReadOnlySet<int> AllowedUnimodIds => SupportedUnimodIds;
        public override SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public override FragmentIonMappingMode FragmentIonMappingMode { get; init; }
        public override int NumberOfPredictedFragmentIons => 380;

        public Altimeter2024Intensities(
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

            var batchedRequests = new List<Dictionary<string, object>>(batchedPeptides.Length);
            for (int i = 0; i < batchedPeptides.Length; i++)
            {
                batchedRequests.Add(BuildBatchedRequest(i,
                    new InputField("peptide_sequences", "BYTES", batchedPeptides[i]),
                    new InputField("precursor_charges", "INT32", batchedCharges[i]),
                    new InputField("collision_energies", "FP32", batchedEnergies[i])));
            }
            return batchedRequests;
        }

        protected override ParsedFragmentAnnotation ParseFragmentAnnotation(string annotation)
        {
            if (string.IsNullOrEmpty(annotation))
            {
                throw new Exception($"Altimeter model returned empty fragment annotation.");
            }

            var caretIndex = annotation.IndexOf('^');
            if (caretIndex > 0)
            {
                var fragmentIdentifier = annotation.Substring(0, caretIndex);
                var chargeStr = annotation.Substring(caretIndex + 1);
                if (!int.TryParse(chargeStr, out int charge))
                {
                    throw new Exception($"Cannot parse charge from Altimeter fragment annotation '{annotation}'. Charge value: '{chargeStr}'.");
                }

                if (TryExtractNeutralLoss(fragmentIdentifier, out var baseId, out var nlFormula, out var nlMass))
                {
                    return new ParsedFragmentAnnotation(fragmentIdentifier, charge, nlFormula, nlMass);
                }

                return new ParsedFragmentAnnotation(fragmentIdentifier, charge);
            }

            // Altimeter uses '+' for internal fragment and neutral gain notation (e.g., "IN+CO"), not for charge.
            // These cannot be mapped to theoretical products.
            if (annotation.Contains('+'))
            {
                throw new Exception($"Altimeter fragment annotation '{annotation}' uses internal fragment or neutral gain notation which is not supported for m/z recalculation. Expected format: 'type^charge' (e.g., 'y8^2') or 'type' (e.g., 'y8').");
            }

            // No charge delimiter — default to charge 1, but check for neutral loss (e.g., "IY-NH3")
            if (TryExtractNeutralLoss(annotation, out var baseIdNoCharge, out var nlFormulaNoCharge, out var nlMassNoCharge))
            {
                return new ParsedFragmentAnnotation(annotation, 1, nlFormulaNoCharge, nlMassNoCharge);
            }

            return new ParsedFragmentAnnotation(annotation, 1);
        }

        private static bool TryExtractNeutralLoss(string fragmentIdentifier, out string baseId, out string? nlFormula, out double? nlMass)
        {
            baseId = fragmentIdentifier;
            nlFormula = null;
            nlMass = null;

            var dashIndex = fragmentIdentifier.IndexOf('-');
            if (dashIndex <= 0)
            {
                return false;
            }

            baseId = fragmentIdentifier.Substring(0, dashIndex);
            var formula = fragmentIdentifier.Substring(dashIndex + 1);
            try
            {
                nlMass = ChemicalFormula.ParseFormula(formula).MonoisotopicMass;
                nlFormula = formula;
                return true;
            }
            catch
            {
                // Formula not parseable — reset to defaults
                baseId = fragmentIdentifier;
                return false;
            }
        }
    }
}
