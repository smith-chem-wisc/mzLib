using Omics.SpectrumMatch;
using System.ComponentModel;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.SupportedModels.FragmentIntensityModels
{
    /// <summary>
    /// Implementation of the Prosit 2020 HCD intensity prediction model for fragment ion intensity prediction.
    /// This model predicts fragment ion intensities for peptides using Higher-energy Collisional Dissociation (HCD).
    /// </summary>
    /// <remarks>
    /// Model specifications:
    /// - Supports peptides with length 1-30 amino acids
    /// - Handles precursor charges 1-6
    /// - Predicts up to 174 fragment ions per peptide
    /// - Supports carbamidomethylation on cysteine (required) and oxidation on methionine (optional)
    /// - Optimal collision energies: 20, 23, 25, 28, 30, 35
    /// 
    /// API Documentation: https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_HCD/infer
    /// </remarks>
    public class Prosit2020IntensityHCD : FragmentIntensityModel
    {
        /// <summary>The Koina API model name identifier</summary>
        public override string ModelName => "Prosit_2020_intensity_HCD";
        /// <summary>Maximum number of peptides that can be processed in a single API request</summary>
        public override int MaxBatchSize => 1000;
        /// <summary>Maximum allowed peptide sequence length in amino acids</summary>
        public override int MaxPeptideLength => 30;
        /// <summary>Minimum allowed peptide sequence length in amino acids</summary>
        public override int MinPeptideLength => 1;
        /// <summary>Set of supported precursor charge states for this model</summary>
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        /// <summary>Total number of fragment ions predicted by this model per peptide</summary>
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
        /// <summary>
        /// HCD collision energies for fragmentation prediction.
        /// Model performs optimally with energies: 20, 23, 25, 28, 30, 35.
        /// </summary>
        public List<int> CollisionEnergies { get; } = new();
        public override List<double?> RetentionTimes { get; } = new();
        public override List<PeptideFragmentIntensityPrediction> Predictions { get; protected set; } = new();
        public override List<LibrarySpectrum> PredictedSpectra { get; protected set; } = new();
        public override string? SpectralLibrarySavePath { get; } = null;
        /// <summary>
        /// Minimum intensity threshold for fragment ions to be included in spectral library generation.
        /// </summary>
        public override double MinIntensityFilter { get; protected set; } = 1e-4;

        /// <summary>
        /// Initializes a new instance of the Prosit2020IntensityHCD model with input validation and filtering.
        /// </summary>
        /// <param name="peptideSequences">
        /// Peptide sequences with modifications in mzLib format. 
        /// Valid modifications: "[Common Variable:Oxidation on M]" and "[Common Fixed:Carbamidomethyl on C]"
        /// </param>
        /// <param name="precursorCharges">
        /// Precursor charge states (1-6). Must match the length of peptideSequences.
        /// </param>
        /// <param name="collisionEnergies">
        /// HCD collision energies in eV. Must be positive values. 
        /// Optimal performance at: 20, 23, 25, 28, 30, 35 eV.
        /// Must match the length of peptideSequences.
        /// </param>
        /// <param name="retentionTimes">
        /// Optional retention times used for spectral library generation.
        /// Can contain null values. Must match the length of peptideSequences.
        /// </param>
        /// <param name="warnings">
        /// Output parameter containing details about any invalid entries that were filtered out.
        /// Will be null if all inputs are valid.
        /// </param>
        /// <param name="minIntensityFilter">
        /// Minimum intensity threshold for including fragment ions in spectral library.
        /// Default: 1e-4. Lower values include more low-intensity fragments.
        /// </param>
        /// <exception cref="ArgumentException">
        /// Thrown when input lists have different lengths or when all inputs are empty.
        /// </exception>
        /// <example>
        /// <code>
        /// var sequences = new List&lt;string&gt; { "PEPTIDE", "SEQUENCE[Common Variable:Oxidation on M]" };
        /// var charges = new List&lt;int&gt; { 2, 3 };
        /// var energies = new List&lt;int&gt; { 25, 30 };
        /// var times = new List&lt;double?&gt; { 12.5, null };
        /// 
        /// var model = new Prosit2020IntensityHCD(sequences, charges, energies, times, out var warnings);
        /// if (warnings != null) Console.WriteLine(warnings.Message);
        /// </code>
        /// </example>
        public Prosit2020IntensityHCD(List<string> peptideSequences, List<int> precursorCharges, List<int> collisionEnergies, List<double?> retentionTimes, out WarningException? warnings, 
            string? spectralLibrarySavePath = null, double minIntensityFilter = 1e-4)
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
                    !IsValidSequence(peptide) ||
                    !AllowedPrecursorCharges.Contains(charge) ||
                    energy <= 0
                    )
                {
                    invalidArguments.Add($"Index {i}: Peptide '{peptide}' (Length: {peptide.Length}), Charge: {charge}, Collision Energy: {energy}");
                }
                else
                {
                    PeptideSequences.Add(ConvertMzLibModificationsToUnimod(peptide));
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
            SpectralLibrarySavePath = spectralLibrarySavePath;
        }



        // TODO: Construct from a list of LibrarySpectrum objects
        public Prosit2020IntensityHCD(List<LibrarySpectrum> spectralLibrary)
        {
            throw new NotImplementedException();
        }

        // TODO:Construct from a a spectral library file path
        public Prosit2020IntensityHCD(string filePath)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Creates batched API requests formatted for the Prosit 2020 HCD intensity model.
        /// Each batch contains up to MaxBatchSize peptides with their associated parameters.
        /// </summary>
        /// <returns>
        /// List of request dictionaries compatible with Koina API format.
        /// Each request contains peptide sequences, precursor charges, and collision energies.
        /// </returns>
        /// <remarks>
        /// Request structure follows Koina API specification:
        /// - peptide_sequences: BYTES array with UNIMOD-formatted sequences
        /// - precursor_charges: INT32 array with charge states
        /// - collision_energies: FP32 array with HCD collision energies
        /// 
        /// Batching improves performance by reducing HTTP overhead and enables
        /// efficient processing of large peptide datasets.
        /// </remarks>
        protected override List< Dictionary<string, object> > ToBatchedRequests()
        {
            // Split inputs into batches
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
                                name = "peptide_sequences", // Model input name
                                shape = new[]{ batchedPeptides[i].Length, 1 }, // Tensor shape
                                datatype = "BYTES", // Data type for string sequences
                                data = batchedPeptides[i] // Actual peptide sequences
                            },
                            new {
                                name = "precursor_charges",
                                shape = new[]{ batchedCharges[i].Length, 1 },
                                datatype = "INT32", // Data type for integer charges
                                data = batchedCharges[i]
                            },
                            new {
                                name = "collision_energies",
                                shape = new[]{ batchedEnergies[i].Length, 1 },
                                datatype = "FP32", // Data type for floating-point energies
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
