using Omics.SpectrumMatch;
using System.ComponentModel;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    public class Prosit2020IntensityHCD : FragmentIntensityModel
    {
        public override string ModelName => "Prosit_2020_intensity_HCD";
        public override int MaxBatchSize => 1000;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public int NumberOfPredictedFragmentIons => 174;
        public override Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            {"[Common Variable:Oxidation on M]", "[UNIMOD:35]"},
            {"[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]"}
        };
        public override Dictionary<string, double> ValidModificationsMonoisotopicMasses => new()
        {
            {"[Common Variable:Oxidation on M]", 15.994915 },
            {"[Common Fixed:Carbamidomethyl on C]", 57.021464 }
        };
        public override List<string> PeptideSequences { get; } = new();
        public override List<int> PrecursorCharges { get; } = new();
        public List<int> CollisionEnergies { get; } = new();
        public override List<double?> RetentionTimes { get; } = new();
        public override List<PeptideFragmentIntensityPrediction> Predictions { get; protected set; } = new();
        public override List<LibrarySpectrum> PredictedSpectra { get; protected set; } = new();
        public override double MinIntensityFilter { get; protected set; } = 1e-4;

        /// <summary>
        /// Client for the Prosit 2020 HCD Intensity Model 
        /// Model details: https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_HCD/infer
        /// It is expected that all of the input lists are of the same length.
        /// Designed to perform optimally with large input lists, rather than requesting predictions per peptide.
        /// Large inputs are automatically batched into smaller requests to the Koina server, and a single
        /// HTTP client is used for all requests to improve performance. If the requests are made for each individual 
        /// peptide for large amounts of peptides, socket exhaustion may occur.
        /// 
        /// If any of the input entries are invalid, they will be skipped and a WarningException will be returned.
        /// If the inputs are all invalid or empty, no predictions will be made.
        /// 
        /// All cysteines will be carbamidomethylated as per model input requirements.
        /// The valid modifications are carbamidomethylation on C and oxidation on M.
        /// 
        /// </summary>
        /// <param name="peptideSequences">Peptide sequences with valid modifications for fragment intensity predictions.</param>
        /// <param name="precursorCharges">Precursor charge states. Valid charge states are 1-6.</param>
        /// <param name="collisionEnergies">HCD collision energies. The model performs best for collision energies 20, 23, 25, 28, 30, and 35.</param>
        /// <param name="retentionTimes">Retention time for the peptide. This is used for LibrarySpectrum creation.</param>
        /// <param name="warnings">String message with the peptides omitted from predictions. These peptides are filtered out based on model input requirements. </param>
        /// <param name="minIntensityFilter">This is the intensity filter for fragment ions kept and included in the library spectrum.</param>
        /// <exception cref="ArgumentException"></exception>
        public Prosit2020IntensityHCD(List<string> peptideSequences, List<int> precursorCharges, List<int> collisionEnergies, List<double?> retentionTimes, out WarningException? warnings, double minIntensityFilter=1e-4)
        {
            // Verify input lists are of the same length
            if (peptideSequences.Count != precursorCharges.Count
                || precursorCharges.Count != collisionEnergies.Count
                || collisionEnergies.Count != retentionTimes.Count)
            {
                throw new ArgumentException("Input lists must have the same length.");
            }

            if (peptideSequences.Count == 0)
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            // Ensure that the inputs meet the model requirements.
            List<string> invalidArguments = new();
            for (int i = 0; i < peptideSequences.Count; i++)
            {
                var peptide = peptideSequences[i];
                var charge = precursorCharges[i];
                var energy = collisionEnergies[i];
                var retentionTime = retentionTimes[i];

                // Skip invalid peptides
                if (!HasValidModifications(peptide) ||
                    !IsValidPeptideSequence(peptide) ||
                    !AllowedPrecursorCharges.Contains(charge) ||
                    energy <= 0
                    )
                {
                    invalidArguments.Add($"Index {i}: Peptide '{peptide}' (Length: {peptide.Length}), Charge: {charge}, Collision Energy: {energy}");
                }
                else
                {
                    PeptideSequences.Add(ConvertToPrositModificationFormat(peptide));
                    PrecursorCharges.Add(charge);
                    CollisionEnergies.Add(energy);
                    RetentionTimes.Add(retentionTime);
                }
            }

            warnings = null;
            if (invalidArguments.Count > 0)
            {
                string warningMessage = "The following input entries are invalid and will be skipped:\n"
                    + string.Join("\n", invalidArguments)
                    + "\nModel Requirements:\n"
                    + $"- Peptide length <= {MaxPeptideLength}\n"
                    + "- Peptide sequence is not empty\n"
                    + "- Peptide has valid modifications\n"
                    + $"- Precursor charge in [{string.Join(", ", AllowedPrecursorCharges)}]\n"
                    + "- Collision energy > 0\n"
                    + "- Retention time is null or >= 0 (for Spectrum Library only)";
                warnings = new WarningException(warningMessage);
            }

            MinIntensityFilter = minIntensityFilter;
        }


        // Construct from a list of LibrarySpectrum objects
        public Prosit2020IntensityHCD(List<LibrarySpectrum> spectralLibrary)
        {
            throw new NotImplementedException();
        }

        // Construct from a a spectral library file path
        public Prosit2020IntensityHCD(string filePath)
        {
            throw new NotImplementedException();
        }

        protected override List< Dictionary<string, object> > ToBatchedRequests()
        {

            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
            var batchedCharges = PrecursorCharges.Chunk(MaxBatchSize).ToList();
            var batchedEnergies = CollisionEnergies.Chunk(MaxBatchSize).ToList();

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
