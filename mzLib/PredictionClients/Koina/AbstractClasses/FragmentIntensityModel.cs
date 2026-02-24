using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.Client;
using Readers.SpectralLibrary;
using System.ComponentModel;
using Chemistry;
using Proteomics.AminoAcidPolymer;
using PredictionClients.Koina.Interfaces;
using System.Data.SQLite;
using Easy.Common.Extensions;

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
        List<string>? FragmentAnnotations,
        List<double>? FragmentMZs,
        List<double>? FragmentIntensities,
        WarningException? Warning = null
    );

    /// <summary>
    /// Represents all input parameters for the fragment intensity prediction models from the Koina API.
    /// For each model, the required and optional parameters may differ, but this record captures the common 
    /// set of inputs that are typically used across different fragment intensity models.
    /// Each model will look for specific parameters within this record and may ignore others, but this provides 
    /// a standardized way to pass all relevant information to the models.
    /// </summary>
    /// <param name="FullSequence">Peptide sequence with modifications in UNIMOD format (used in every model)</param>
    /// <param name="PrecursorCharge">Charge state of the precursor ion (used in every model)</param>
    /// <param name="CollisionEnergy">Collision energy used for fragmentation (not used by some models)</param>
    /// <param name="InstrumentType">Type of mass spectrometer instrument (not used by some models)</param>
    /// <param name="FragmentationType">Fragmentation method used (not used by some models)</param>
    public record FragmentIntensityPredictionInput(
        string FullSequence,
        int PrecursorCharge,
        int? CollisionEnergy,
        string? InstrumentType,
        string? FragmentationType
    )
    {
        public string? ValidatedFullSequence { get; set; }
        public WarningException? SequenceWarning { get; set; }
        public WarningException? ParameterWarning { get; set; }
    }

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
    public abstract class FragmentIntensityModel : KoinaModelBase<FragmentIntensityPredictionInput, PeptideFragmentIntensityPrediction>, IPredictor<FragmentIntensityPredictionInput, PeptideFragmentIntensityPrediction>
    {

        #region Additional Model-Type Constraints
        /// <summary>Set of precursor charge states supported by the model (e.g., {2, 3, 4})</summary>
        public abstract HashSet<int> AllowedPrecursorCharges { get; }
        public virtual HashSet<int> AllowedCollisionEnergies => new HashSet<int>(); 
        public virtual HashSet<string> AllowedInstrumentTypes => new HashSet<string>();
        public virtual HashSet<string> AllowedFragmentationTypes => new HashSet<string>();
        /// <summary>
        /// Maps mzLib modification format to monoisotopic mass differences for mass-only conversions. 
        /// The monoisotopic masses should be obtained from the UNIMOD database for accuracy (https://www.unimod.org/).
        /// While UNIMOD monoisotopic masses already account for net mass effects (such as loss of H+ for Na+ mod),
        /// this is something that should be kept in mind when adding custom modifications.
        /// </summary>
        public virtual Dictionary<string, double> ValidModificationsMonoisotopicMasses { get; protected set; } = new();
        #endregion

        #region Validation Patterns and Filters
        /// <summary>Minimum intensity threshold for including predicted fragments in spectral library generation</summary>
        public virtual double MinIntensityFilter { get; protected set; } = 1e-6;
        #endregion

        #region Inputs and Outputs for internal processing 
        // NOTE: Since the model is expected to be reusable for multiple predictions, these properties are not static and are intended to be set during the prediction workflow.
        // Additionally, these properties are recorded for each prediction session to allow for easy access to the original inputs and outputs when batching, response parsing,
        // and spectral library generation. 

        /// <summary>
        /// Inputs provided to the model for the LATEST prediction session. This list is populated during the PredictAsync workflow and is used to keep track of the original 
        /// input parameters for each prediction, especially when batching is involved. Each FragmentIntensityPredictionInput contains all relevant information (sequence, charge, 
        /// collision energy, etc.) that was used for that specific prediction. This allows for better traceability and debugging, as well as providing necessary context when 
        /// generating spectral libraries from the predictions.
        /// </summary>
        public List<FragmentIntensityPredictionInput> ModelInputs { get; protected set; } = new();
        /// <summary>
        /// Boolean mask indicating which inputs from the original list were valid for prediction after applying model-specific validation criteria.
        /// This is used for realigning predictions back to the original input list and for filtering out invalid inputs from the prediction results.
        /// The mask is populated during the PredictAsync workflow after validating each input against the model's constraints (e.g., allowed precursor charges, collision energies, etc.). 
        /// A value of 'true' at index i indicates that the input at index i in ModelInputs was valid and included in the prediction process, 
        /// while 'false' indicates that it was filtered out due to incompatibility with the model. 
        /// This allows for better handling of mixed input lists where some entries may not meet the model's requirements, without losing track of their original positions in the input list.
        /// </summary>
        public bool[] ValidInputsMask { get; protected set; } = Array.Empty<bool>();
        /// <summary>
        /// Collection of fragment intensity prediction results after inference completion for the LATEST prediction session. Each PeptideFragmentIntensityPrediction contains the original input sequence,
        /// precursor charge, fragment annotations, m/z values, and predicted intensities for that specific peptide. This list is populated during the PredictAsync workflow after parsing the API responses 
        /// and realigning the predictions back to the original input list using the ValidInputsMask.
        /// </summary>
        public List<PeptideFragmentIntensityPrediction> Predictions { get; protected set; } = new();
        #endregion
        // TODO: Implement Caches to optimize performance by avoiding redundant computations during sequence validation and modification conversions.


        #region Querying Methods for the Koina API

        /// <summary>
        /// Executes fragment intensity prediction by sending batched requests to the Koina API.
        /// The method performs the following steps:
        /// 1. Validates and cleans input sequences according to model constraints, populating the ValidInputsMask to keep track of which inputs are valid for prediction.
        /// 2. Converts valid inputs into batched request payloads formatted for the specific model using the ToBatchedRequests method.
        /// 3. Sends batched requests to the Koina API with throttling between batches to avoid overwhelming the server, and processes responses to extract predictions.
        /// 4. Realigns predictions back to the original input list using the ValidInputsMask, ensuring that the output list corresponds to the original input order and includes placeholders for invalid inputs with appropriate warnings.
        /// </summary>
        /// <returns>Task representing the asynchronous inference operation</returns>
        protected virtual async Task<List<PeptideFragmentIntensityPrediction>> AsyncThrottledPredictor(List<FragmentIntensityPredictionInput> modelInputs)
        {
            #region Input Validation and Cleaning
            ModelInputs = modelInputs;
            ValidInputsMask = new bool[ModelInputs.Count];
            var validInputs = new List<FragmentIntensityPredictionInput>();
            for (int i = 0; i < ModelInputs.Count; i++)
            {
                var cleanedSequence = TryCleanSequence(ModelInputs[i].FullSequence, out var modHandlingWarning); // mod handling happens here
                var validModelParams = ValidateModelSpecificInputs(ModelInputs[i], out var modelParametersWarning);
                if (cleanedSequence != null && validModelParams)
                {
                    ModelInputs[i] = ModelInputs[i] with { ValidatedFullSequence = cleanedSequence, SequenceWarning = modHandlingWarning, ParameterWarning = modelParametersWarning };
                    ValidInputsMask[i] = true;
                    validInputs.Add(ModelInputs[i]);
                }
                else
                {
                    ModelInputs[i] = ModelInputs[i] with { ValidatedFullSequence = null, SequenceWarning = modHandlingWarning, ParameterWarning = modelParametersWarning };
                    ValidInputsMask[i] = false;
                }
            }
            #endregion

            #region Request Batching, Throttling Setup
            var batchedRequests = ToBatchedRequests(validInputs);
            var batchChunks = batchedRequests.Chunk(MaxNumberOfBatchesPerRequest).ToList();
            int sessionTimeoutInMinutes = batchedRequests.Count * 2 + (int)(ThrottlingDelayInMilliseconds/6000*batchChunks.Count) + 2; // Dynamic timeout: ~2 minutes per batch + throttle time between batches + 2 minute buffer for network/processing overhead.
            #endregion

            #region Throttled API Requests and Response Processing
            var predictions = new List<PeptideFragmentIntensityPrediction>();
            var responses = new string[validInputs.Count];
            using var _http = new HTTP(timeoutInMinutes: sessionTimeoutInMinutes); // Set a reasonable timeout for each batch chunk
            for (int i = 0; i < batchChunks.Count; i++)
            {
                var batchChunk = batchChunks[i];
                var responseChunk = await Task.WhenAll(batchChunk.Select(request => _http.InferenceRequest(ModelName, request)));
                responses.Add(responseChunk);
                if (i < batchChunks.Count - 1) // No need to throttle after the last batch
                {
                    await Task.Delay((int)ThrottlingDelayInMilliseconds);
                }
            }

            predictions = ResponseToPredictions(responses, validInputs);
            #endregion

            #region Realign Predictions to Original Input List
            // Realign predictions back to the original input list using the ValidInputsMask
            var realignedPredictions = new List<PeptideFragmentIntensityPrediction>();
            int predictionIndex = 0;
            for (int i = 0; i < ValidInputsMask.Length; i++)
            {
                if (ValidInputsMask[i])
                {
                    realignedPredictions.Add(predictions[predictionIndex]);
                    predictionIndex++;
                }
                else
                {
                    // For invalid inputs, we can choose to add a placeholder prediction with a warning, or simply skip them. Here we add a placeholder with a warning for traceability.
                    realignedPredictions.Add(new PeptideFragmentIntensityPrediction(
                        FullSequence: ModelInputs[i].FullSequence,
                        PrecursorCharge: ModelInputs[i].PrecursorCharge,
                        FragmentAnnotations: null,
                        FragmentMZs: null,
                        FragmentIntensities: null,
                        Warning: ModelInputs[i].ParameterWarning ?? ModelInputs[i].SequenceWarning ?? new WarningException("Input was invalid and skipped during prediction.")
                    ));
                }
            }
            #endregion

            Predictions = realignedPredictions;
            return Predictions;
        }

        public List<PeptideFragmentIntensityPrediction> Predict(List<FragmentIntensityPredictionInput> modelInputs)
        {
            return AsyncThrottledPredictor(modelInputs).GetAwaiter().GetResult();
        }

        /// <summary>
        /// Converts Koina API responses into structured prediction objects.
        /// Expects responses with exactly 3 outputs: annotations, m/z values, and intensities.
        /// The only filtering done at this stage is removing impossible ions (intensity = -1). 
        /// MinIntensityFilter is applied later during spectral library generation.
        /// </summary>
        /// <param name="responses">Array of JSON response strings from Koina API</param>
        /// <exception cref="Exception">Thrown when responses are malformed or contain unexpected number of outputs</exception>
        protected virtual List<PeptideFragmentIntensityPrediction> ResponseToPredictions(string[] responses, List<FragmentIntensityPredictionInput> requestInputs)
        {
            var predictions = new List<PeptideFragmentIntensityPrediction>();
            if (requestInputs.IsNullOrEmpty())
            {
                return predictions;
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
                var batchPeptides = requestInputs.Skip(batchIndex * MaxBatchSize).Take(MaxBatchSize).ToList();
                // Assuming outputData is structured such that each peptide's data is sequential
                var fragmentCount = outputAnnotations.Count / batchPeptides.Count;
                for (int i = 0; i < batchPeptides.Count; i++)
                {
                    var peptideSequence = batchPeptides[i].ValidatedFullSequence;
                    var fragmentIons = new List<string>();
                    var fragmentMZs = new List<double>();
                    var predictedIntensities = new List<double>();
                    for (int j = 0; j < fragmentCount; j++)
                    {
                        fragmentIons.Add(outputAnnotations[i * fragmentCount + j].ToString()!);
                        fragmentMZs.Add(Convert.ToDouble(outputMZs[i * fragmentCount + j]));
                        predictedIntensities.Add(Convert.ToDouble(outputIntensities[i * fragmentCount + j]));
                    }
                    predictions.Add(new PeptideFragmentIntensityPrediction(
                        peptideSequence!,
                        ModelInputs[batchIndex * MaxBatchSize + i].PrecursorCharge,
                        fragmentIons,
                        fragmentMZs,
                        predictedIntensities
                    ) with
                    { Warning = batchPeptides[i].SequenceWarning });
                }
            }
            return predictions;
        }
        #endregion

        #region Full Sequence Modification Conversion Methods

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

        /// <summary>
        /// Validates the input parameters against the model's specific constraints (e.g., allowed precursor charges, collision energies, instrument types, fragmentation types).
        /// Returns true if the input is valid for this model, false otherwise. If invalid, a WarningException is provided with details about the incompatibility.
        /// </summary>
        /// <param name="input"> FragmentIntensityPredictionInput object containing all relevant input parameters for the model</param>
        /// <param name="warning"> Output parameter that will contain a WarningException with details if the input is invalid, or null if the input is valid. This allows the calling code to log or handle warnings without throwing exceptions for expected incompatibilities. </param>
        /// <returns></returns>
        protected virtual bool ValidateModelSpecificInputs(FragmentIntensityPredictionInput input, out WarningException? warning)
        {
            warning = null;
            if (!AllowedPrecursorCharges.IsNullOrEmpty() && !AllowedPrecursorCharges.Contains(input.PrecursorCharge))
            {
                warning = new WarningException($"Precursor charge {input.PrecursorCharge} is not supported by this model.");
                return false;
            }
            if (!AllowedCollisionEnergies.IsNullOrEmpty() && input.CollisionEnergy.HasValue && !AllowedCollisionEnergies.Contains(input.CollisionEnergy.Value))
            {
                warning = new WarningException($"Collision energy {input.CollisionEnergy.Value} is not supported by this model.");
                return false;
            }
            if (!AllowedInstrumentTypes.IsNullOrEmpty() && input.InstrumentType != null && !AllowedInstrumentTypes.Contains(input.InstrumentType))
            {
                warning = new WarningException($"Instrument type '{input.InstrumentType}' is not supported by this model.");
                return false;
            }
            if (!AllowedFragmentationTypes.IsNullOrEmpty() && input.FragmentationType != null && !AllowedFragmentationTypes.Contains(input.FragmentationType))
            {
                warning = new WarningException($"Fragmentation type '{input.FragmentationType}' is not supported by this model.");
                return false;
            }
            return true;
        }

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
        /// <exception cref="WarningException">Recorded in the out parameter when duplicate spectra are detected in predictions</exception>
        public List<LibrarySpectrum> GenerateLibrarySpectraFromPredictions(double[] alignedRetentionTimes, out WarningException? warning, string? filepath=null, double minIntensityFilter=1e-4)
        {
            warning = null;
            if (Predictions.Count == 0)
            {
                return new List<LibrarySpectrum>(); // No predictions to process
            }

            if (Predictions.Count != alignedRetentionTimes.Length)
            {
                throw new ArgumentException("The number of predictions must match the number of aligned retention times.");
            }

            var predictedSpectra = new List<LibrarySpectrum>();
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
                    if ((int)prediction.FragmentIntensities[fragmentIndex] == -1 || prediction.FragmentIntensities[fragmentIndex] < minIntensityFilter)
                    {
                        // Skip impossible ions and peaks with near zero intensity. The model uses -1 to indicate impossible ions.
                        continue;
                    }

                    var annotation = prediction.FragmentAnnotations[fragmentIndex];

                    // This check is to ensure that the annotation contains the expected format (e.g., "b5+1") before attempting to parse it.
                    // The API WILL always contain it in the expected format, but this is a safeguard against any malformed annotations that could cause exceptions during parsing.
                    if (annotation == null || !annotation.Contains('+'))
                    {
                        // Skip malformed annotations that do not contain expected ion type and charge information.
                        continue;
                    }

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
                    rt: alignedRetentionTimes[predictionIndex]
                );

                predictedSpectra.Add(spectrum);
            }
            var warningString = $"Generated {predictedSpectra.Count} spectra from predictions.\n";
            var unique = predictedSpectra.DistinctBy(p => p.Name).ToList();
            if (unique.Count != predictedSpectra.Count)
            {
                warning = new WarningException($"Duplicate spectra found in predictions. Reduced from {predictedSpectra.Count} predicted spectra to {unique.Count} unique spectra.");
                predictedSpectra = unique;
            }

            if (filepath == null)
            {
                warningString += "No file path provided for spectral library output. Generated spectra will not be saved to disk.\n";
            }
            else
            {
                var spectralLibrary = new SpectralLibrary();
                spectralLibrary.Results = predictedSpectra;
                spectralLibrary.WriteResults(filepath);
            }

            return predictedSpectra;
        }
        #endregion
    }
}