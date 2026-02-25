using Omics.SpectrumMatch;
using System.ComponentModel;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;

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

        /// <summary>
        /// Maximum number of batches that should be processed in a single API request. This is necessary 
        /// to prevent overwhelming the server with too many concurrent requests, which can lead to timeouts 
        /// or rate limiting. Adjust this value based on the expected number of peptides and server capacity.
        /// </summary>
        public override int MaxNumberOfBatchesPerRequest { get; init; }

        /// <summary>
        /// Throttle time between batches to avoid overwhelming the server. 
        /// Adjust as needed based on model performance and server capacity.
        /// </summary> 
        public override int ThrottlingDelayInMilliseconds { get; init; }
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;

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
        public override IncompatibleModHandlingMode ModHandlingMode { get; init; }

        public Prosit2020IntensityHCD(IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.RemoveIncompatibleMods, int maxNumberOfBatchesPerRequest = 250, int throttlingDelayInMilliseconds = 100)
        {
            ModHandlingMode = modHandlingMode;
            MaxNumberOfBatchesPerRequest = maxNumberOfBatchesPerRequest;
            ThrottlingDelayInMilliseconds = throttlingDelayInMilliseconds;
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
        protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
        {
            // Split inputs into batches
            // ValidatedFullSequence should not be null at this point due to prior validation steps
            var batchedPeptides = validInputs.Select(p => ConvertMzLibModificationsToUnimod(p.ValidatedFullSequence!)).Chunk(MaxBatchSize).ToList();
            var batchedCharges = validInputs.Select(p => p.PrecursorCharge).Chunk(MaxBatchSize).ToList();
            var batchedEnergies = validInputs.Select(p => p.CollisionEnergy).Chunk(MaxBatchSize).ToList();

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
