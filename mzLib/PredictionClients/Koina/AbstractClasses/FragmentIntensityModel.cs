using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using Readers.SpectralLibrary;
using System.ComponentModel;
using System.Text.RegularExpressions;
using Chemistry;
using Proteomics.AminoAcidPolymer;

namespace PredictionClients.Koina.AbstractClasses
{
    /// <summary>
    /// Represents the prediction results for a single peptide, containing fragment annotations,
    /// m/z values, and predicted intensities from a fragment intensity model.
    /// </summary>
    /// <param name="FullSequence">The peptide sequence (with modifications in UNIMOD format)</param>
    /// <param name="FragmentAnnotations">Fragment ion annotations (e.g., "b5+1", "y3+2")</param>
    /// <param name="FragmentMZs">Theoretical m/z values for each fragment ion</param>
    /// <param name="FragmentIntensities">Predicted relative intensities (0-1 scale, -1 indicates impossible ions)</param>
    public record PeptideFragmentIntensityPrediction(
        string FullSequence,
        int PrecursorCharge,
        List<string> FragmentAnnotations,
        List<double> FragmentMZs,
        List<double> FragmentIntensities
    );

    /// <summary>
    /// Abstract base class for fragment intensity prediction models using the Koina API.
    /// Derived classes must implement model-specific details (name, batch size, and input constraints).
    /// Default behaviors are provided for sequence/mod validation, request batching, and response parsing.
    /// Spectral library generation from predictions is also included.
    /// 
    /// Designed to perform optimally with large input lists, rather than requesting predictions per peptide.
    /// Large inputs are automatically batched into smaller requests to the Koina server, and a single
    /// HTTP client is used for all requests to improve performance. If the requests are made for each individual 
    /// peptide for large amounts of peptides, socket exhaustion may occur.
    /// 
    /// In most cases, users will only need to implement a constructor to properly set up the model parameters
    /// and a ToBatchedRequests method to batch requests according to the specific model's input format.
    /// </summary>
    public abstract class FragmentIntensityModel : IKoinaModelIO
    {
        #region Model Metadata from Koina
        public abstract string ModelName { get; }
        public abstract int MaxBatchSize { get; }
        #endregion

        #region Input Constraints
        /// <summary>Maximum allowed peptide length (amino acids) for the model</summary>
        public abstract int MaxPeptideLength { get; }
        /// <summary>Minimum allowed peptide length (amino acids) for the model</summary>
        public abstract int MinPeptideLength { get; }
        /// <summary>Set of precursor charge states supported by the model (e.g., {2, 3, 4})</summary>
        public abstract HashSet<int> AllowedPrecursorCharges { get; }
        /// <summary>Maps mzLib modification format to UNIMOD format (e.g., "[Common Fixed:Carbamidomethyl on C]" -> "[UNIMOD:4]")</summary>
        public virtual Dictionary<string, string> ValidModificationUnimodMapping { get; protected set; } = new();
        /// <summary>
        /// Maps mzLib modification format to monoisotopic mass differences for mass-only conversions. 
        /// The monoisotopic masses should be obtained from the UNIMOD database for accuracy (https://www.unimod.org/).
        /// While UNIMOD monoisotopic masses already account for net mass effects (such as loss of H+ for Na+ mod),
        /// this is something that should be kept in mind when adding custom modifications.
        /// </summary>
        public virtual Dictionary<string, double> ValidModificationsMonoisotopicMasses { get; protected set; } = new();
        #endregion

        #region Validation Patterns and Filters
        /// <summary>Regex pattern matching valid canonical amino acid sequences (default: 20 standard amino acids)</summary>
        public virtual string CanonicalAminoAcidPattern { get; protected set; } = @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        /// <summary>
        /// Regex pattern for detecting modifications in peptide sequences (default: matches content within square brackets)
        /// ex. "PEC[Common Fixed:Carbamidomethyl on C]TIDE"
        /// </summary>
        public virtual string ModificationPattern { get; protected set; } = @"\[[^\]]+\]";
        /// <summary>Minimum intensity threshold for including predicted fragments in spectral library generation</summary>
        public virtual double MinIntensityFilter { get; protected set; } = 1e-6;
        #endregion

        #region Input Data
        // These two inputs are shared across all fragment intensity models.
        // Other models may require additional inputs, which should be implemented in the derived classes.
        public abstract List<string> PeptideSequences { get; }
        public abstract List<int> PrecursorCharges { get; }

        // Retention times are not used for fragment intensity prediction. They are useful for library spectrum generation.
        public abstract List<double?> RetentionTimes { get; }
        #endregion

        #region Output Data
        public abstract List<PeptideFragmentIntensityPrediction> Predictions { get; protected set; }
        public virtual List<LibrarySpectrum> PredictedSpectra { get; protected set; } = new();
        #endregion

        #region Querying Methods for the Koina API
        /// <summary>
        /// Converts input data into batched requests formatted for the specific model's API endpoint.
        /// Each batch should not exceed MaxBatchSize peptides and must include all required model inputs.
        /// </summary>
        /// <returns>List of request dictionaries, each containing a batch of model inputs</returns>
        /// <remarks>
        /// Implementation should handle model-specific input formatting and ensure proper batching
        /// based on the model's constraints and API requirements.
        /// </remarks>
        protected abstract List<Dictionary<string, object>> ToBatchedRequests();

        /// <summary>
        /// Executes fragment intensity prediction by sending batched requests to the Koina API.
        /// Handles HTTP client lifecycle, request batching, and response parsing automatically.
        /// </summary>
        /// <returns>Task representing the asynchronous inference operation</returns>
        /// <exception cref="Exception">Thrown when API responses cannot be deserialized or have unexpected format</exception>
        public virtual async Task RunInferenceAsync()
        {
            // Dynamic timeout: ~2 minutes per batch + 2 minute buffer for network/processing overhead. Typically a 
            // batch takes less than a minute. 
            int numBatches = (int)Math.Ceiling((double)PeptideSequences.Count / MaxBatchSize);
            using var _http = new HTTP(timeoutInMinutes: numBatches * 2 + 2);

            var responses = await Task.WhenAll(ToBatchedRequests().Select(request => _http.InferenceRequest(ModelName, request)));
            ResponseToPredictions(responses);
        }

        /// <summary>
        /// Converts Koina API responses into structured prediction objects.
        /// Expects responses with exactly 3 outputs: annotations, m/z values, and intensities.
        /// </summary>
        /// <param name="responses">Array of JSON response strings from Koina API</param>
        /// <exception cref="Exception">Thrown when responses are malformed or contain unexpected number of outputs</exception>
        protected virtual void ResponseToPredictions(string[] responses)
        {
            if (PeptideSequences.Count == 0)
            {
                Predictions = new List<PeptideFragmentIntensityPrediction>();
                return;
            }
            var deserializedResponses = responses.Select(response => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(response)).ToList();
            var numBatches = deserializedResponses.Count;
            if (deserializedResponses.IsNullOrEmpty() || deserializedResponses.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }
            for (int batchIndex = 0; batchIndex < numBatches; batchIndex++)
            {
                var response = deserializedResponses[batchIndex];

                /// Each response is excpected to have an outputs list with
                /// 1) Fragment Annotations
                /// 2) Fragment MZs
                /// 3) Fragment Intensities
                if (response == null || response.Outputs.Count != 3)
                {
                    throw new Exception($"API response is not in the expected format. Expected 3 outputs, got {response?.Outputs.Count}.");
                }

                var outputAnnotations = response.Outputs[0].Data;
                var outputMZs = response.Outputs[1].Data;
                var outputIntensities = response.Outputs[2].Data;
                var batchPeptides = PeptideSequences.Skip(batchIndex * MaxBatchSize).Take(MaxBatchSize).ToList();
                // Assuming outputData is structured such that each peptide's data is sequential
                var fragmentCount = outputAnnotations.Count / batchPeptides.Count;
                for (int i = 0; i < batchPeptides.Count; i++)
                {
                    var peptideSequence = batchPeptides[i];
                    var fragmentIons = new List<string>();
                    var fragmentMZs = new List<double>();
                    var predictedIntensities = new List<double>();
                    for (int j = 0; j < fragmentCount; j++)
                    {
                        fragmentIons.Add(outputAnnotations[i * fragmentCount + j].ToString()!);
                        fragmentMZs.Add(Convert.ToDouble(outputMZs[i * fragmentCount + j]));
                        predictedIntensities.Add(Convert.ToDouble(outputIntensities[i * fragmentCount + j]));
                    }
                    Predictions.Add(new PeptideFragmentIntensityPrediction(
                        peptideSequence,
                        PrecursorCharges[batchIndex * MaxBatchSize + i],
                        fragmentIons,
                        fragmentMZs,
                        predictedIntensities
                    ));
                }
            }
        }

        /// <summary>
        /// Validates that a peptide sequence meets model constraints for length and amino acid composition.
        /// Strips modifications before validation to check only the base sequence.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential modifications</param>
        /// <returns>True if sequence is valid for the model; otherwise false</returns>
        protected virtual bool IsValidPeptideSequence(string sequence)
        {
            var baseSequence = Regex.Replace(sequence, ModificationPattern, "");
            return Regex.IsMatch(baseSequence, CanonicalAminoAcidPattern)
                && baseSequence.Length <= MaxPeptideLength
                && baseSequence.Length >= MinPeptideLength;
        }

        /// <summary>
        /// Checks if all modifications in a sequence are supported by the model.
        /// Base implementation only allows sequences without modifications.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential modifications</param>
        /// <returns>True if all modifications are valid; otherwise false</returns>
        protected virtual bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            return matches.All(m => ValidModificationUnimodMapping.ContainsKey(m.Value));
        }
        #endregion

        #region Full Sequence Modification Conversion Methods
        /// <summary>
        /// Converts a peptide sequence from the mzLib modification format to the UNIMOD format.
        /// Example: "[Common Fixed:Carbamidomethyl on C]" becomes "[UNIMOD:4]"
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib modification format.</param>
        protected virtual string ConvertMzLibModificationsToUnimod(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Key, mod.Value);
            }
            return sequence;
        }

        /// <summary>
        /// Converts peptide sequence from UNIMOD format back to mzLib modification format.
        /// Inverse operation of ConvertMzLibModificationsToUnimod.
        /// </summary>
        protected virtual string ConvertUnimodToMzLibModifications(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Value, mod.Key);
            }
            return sequence;
        }

        /// <summary>
        /// Converts mzLib modifications to mass-only format with 6 decimal precision.
        /// Example: "[Common Fixed:Carbamidomethyl on C]" becomes "[57.021464]"
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib format</param>
        /// <returns>Peptide sequence with modifications as mass values</returns>
        protected virtual string ConvertMzLibModificationsToMassesOnly(string sequence)
        {
            foreach (var mod in ValidModificationsMonoisotopicMasses)
            {
                sequence = sequence.Replace(mod.Key, $"[{mod.Value.ToString("F6")}]");
            }
            return sequence;
        }
        #endregion

        #region Spectral Library Generation
        /// <summary>
        /// Transforms prediction results into LibrarySpectrum objects suitable for spectral library creation.
        /// Filters out impossible ions (intensity = -1) and low-intensity fragments below MinIntensityFilter.
        /// Parses fragment annotations (e.g., "b5+1") to extract ion type, fragment number, and charge state.
        /// </summary>
        /// <remarks>
        /// The conversion process:
        /// 1. Converts sequence format: UNIMOD -> mzLib -> mass-only for Peptide object creation
        /// 2. Parses each fragment annotation to determine ion properties
        /// 3. Creates MatchedFragmentIon objects with experimental m/z and predicted intensities
        /// 4. Builds LibrarySpectrum with precursor information and fragment data
        /// 5. Validates uniqueness of generated spectra by name
        /// </remarks>
        /// <exception cref="WarningException">Thrown when duplicate spectra are detected in predictions</exception>
        public void GenerateLibrarySpectraFromPredictions()
        {
            if (Predictions.Count == 0)
            {
                return; // No predictions to process
            }

            for (int predictionIndex = 0; predictionIndex < Predictions.Count; predictionIndex++)
            {
                var prediction = Predictions[predictionIndex];

                var peptide = new Peptide(
                    ConvertMzLibModificationsToMassesOnly(
                        ConvertUnimodToMzLibModifications(prediction.FullSequence)
                        )
                    );

                List<MatchedFragmentIon> fragmentIons = new();
                for (int fragmentIndex = 0; fragmentIndex < prediction.FragmentAnnotations.Count; fragmentIndex++)
                {
                    if ((int)prediction.FragmentIntensities[fragmentIndex] == -1 || prediction.FragmentIntensities[fragmentIndex] < MinIntensityFilter)
                    {
                        // Skip impossible ions and peaks with near zero intensity. The model uses -1 to indicate impossible ions.
                        continue;
                    }

                    var annotation = prediction.FragmentAnnotations[fragmentIndex];
                    // Parse the annotation to get ion type, number and charge from something like 'b5+1'
                    var ionType = annotation.First().ToString(); // 'b' or 'y'
                    var plusIndex = annotation.IndexOf('+');
                    var fragmentNumber = int.Parse(annotation.Substring(1, plusIndex - 1));
                    var fragmentIonCharge = int.Parse(annotation.Substring(plusIndex + 1));

                    // Create a new MatchedFragmentIon for each output
                    var fragmentIon = new MatchedFragmentIon
                    (
                        neutralTheoreticalProduct: new Product
                        (
                            productType: Enum.Parse<ProductType>(ionType),
                            terminus: ionType == "b" ? FragmentationTerminus.N : FragmentationTerminus.C,
                            // neutralMass is not directly provided by models, and it is not necessary here. If needed, 
                            // compute it from the peptide sequence and fragment information as shown below.
                            // neutralMass: peptide.Fragment(Enum.Parse<FragmentTypes>(ionType), fragmentNumber).First().MonoisotopicMass,
                            neutralMass: 0.0, // Placeholder, not used in this context
                            fragmentNumber: fragmentNumber,
                            residuePosition: fragmentNumber, // For b / y ions, the fragment number corresponds to the residue count from the respective terminus.
                            neutralLoss: 0 // Annotations like "b5+1" do not encode neutral losses, so we explicitly assume no loss as placeholder.
                        ),
                        experMz: prediction.FragmentMZs[fragmentIndex],
                        experIntensity: prediction.FragmentIntensities[fragmentIndex],
                        charge: fragmentIonCharge
                    );

                    fragmentIons.Add(fragmentIon);
                }

                var spectrum = new LibrarySpectrum
                (
                    sequence: prediction.FullSequence,
                    precursorMz: peptide.ToMz(prediction.PrecursorCharge),
                    chargeState: prediction.PrecursorCharge,
                    peaks: fragmentIons,
                    rt: RetentionTimes[predictionIndex]
                );

                PredictedSpectra.Add(spectrum);
            }
            var unique = PredictedSpectra.DistinctBy(p => p.Name).ToList();
            if (unique.Count != PredictedSpectra.Count)
            {
                throw new WarningException($"Duplicate spectra found in predictions. Reduced from {PredictedSpectra.Count} predicted spectra to {unique.Count} unique spectra.");
            }
        }

        /// <summary>
        /// Saves predicted spectra to a spectral library file.
        /// Automatically calls GenerateLibrarySpectraFromPredictions() if not already done.
        /// </summary>
        /// <param name="filePath">Output file path for the spectral library</param>
        /// <example>
        /// Usage:
        /// model.SavePredictedSpectralLibrary("predicted_library.msp");
        /// </example>
        public void SavePredictedSpectralLibrary(string filePath)
        {
            var spectralLibrary = new SpectralLibrary();
            if (PredictedSpectra.Count == 0)
            {
                GenerateLibrarySpectraFromPredictions();
            }
            spectralLibrary.Results = PredictedSpectra;
            spectralLibrary.WriteResults(filePath);
        }
        #endregion
    }
}