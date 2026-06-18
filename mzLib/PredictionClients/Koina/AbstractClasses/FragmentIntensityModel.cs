using MzLibUtil;
using Omics.SequenceConversion;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.Client;
using Readers.SpectralLibrary;
using System.ComponentModel;
using Chemistry;
using PredictionClients.Koina.Interfaces;
using Proteomics.ProteolyticDigestion;
using Easy.Common.Extensions;
using PredictionClients.Koina.Util;
namespace PredictionClients.Koina.AbstractClasses
{
    /// <summary>
    /// Represents parsed fragment annotation data extracted from model output strings.
    /// </summary>
    /// <param name="FragmentIdentifier">Full fragment identifier including neutral loss suffix (e.g., "b5", "y3-H2O")</param>
    /// <param name="Charge">Fragment charge state</param>
    /// <param name="NeutralLossFormula">Neutral loss chemical formula if present (e.g., "H2O", "NH3"); null if no neutral loss</param>
    /// <param name="NeutralLossMass">Monoisotopic mass of the neutral loss; null if no neutral loss</param>
    public record ParsedFragmentAnnotation(string FragmentIdentifier, int Charge, string? NeutralLossFormula = null, double? NeutralLossMass = null)
    {
        /// <summary>
        /// The base fragment identifier with the neutral loss suffix stripped (e.g., "y3" for "y3-H2O").
        /// If there is no neutral loss, this is identical to <see cref="FragmentIdentifier"/>.
        /// </summary>
        public string BaseFragmentIdentifier
        {
            get
            {
                if (NeutralLossFormula != null)
                {
                    var dashIndex = FragmentIdentifier.IndexOf('-');
                    if (dashIndex > 0)
                        return FragmentIdentifier.Substring(0, dashIndex);
                }
                return FragmentIdentifier;
            }
        }
    }

    /// <summary>
    /// Represents the prediction results for a single peptide, containing fragment annotations,
    /// m/z values, and predicted intensities from a fragment intensity model.
    /// </summary>
    /// <param name="FullSequence">Original peptide sequence provided by the user (mzLib format)</param>
    /// <param name="ValidatedFullSequence">Validated and cleaned peptide sequence that was actually used for prediction (Unimod format). This may also differ from the original FullSequence if modifications were removed or if the sequence was deemed invalid for the model. This is the sequence that reflects the actual input to the model.</param>
    /// <param name="PrecursorCharge">Charge state of the precursor ion used for prediction</param>
    /// <param name="FragmentAnnotations">Fragment ion annotations (e.g., "b5+1", "y3+2")</param>
    /// <param name="FragmentMZs">Theoretical m/z values for each fragment ion</param>
    /// <param name="FragmentIntensities">Predicted relative intensities (0-1 scale, -1 indicates impossible ions)</param>
    public record PeptideFragmentIntensityPrediction(
        string FullSequence,
        string ValidatedFullSequence,
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
    /// <param name="FullSequence">Peptide sequence with modifications in mzLib format (original input)</param>
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
    ///
    /// Thread safety: instances are NOT thread-safe. Predict and related methods
    /// mutate instance state (ModelInputs, ValidInputsMask, Predictions); callers must not invoke
    /// these methods concurrently on the same instance, nor read Predictions while a call is in flight.
    /// Use one instance per concurrent caller (or serialize externally) when sharing across pipelines.
    /// </summary>
    public abstract class FragmentIntensityModel : KoinaModelBase<FragmentIntensityPredictionInput, PeptideFragmentIntensityPrediction>, IPredictor<FragmentIntensityPredictionInput, PeptideFragmentIntensityPrediction>
    {
        protected FragmentIntensityModel(ISequenceConverter sequenceConverter)
            : base(sequenceConverter)
        {
        }


        #region Additional Model-Type Constraints
        /// <summary>
        /// Set of precursor charge states supported by the model (e.g., {2, 3, 4}).
        /// null = charge is not applicable to this model (skip validation).
        /// empty = charge IS required but any value is accepted.
        /// populated = only listed charge values are accepted.
        /// </summary>
        public abstract HashSet<int> AllowedPrecursorCharges { get; }
        /// <summary>
        /// Set of collision energies (NCE) accepted by the model.
        /// null = collision energy is not applicable (skip validation).
        /// empty = collision energy IS required but any FP32 value is accepted.
        /// populated = only listed values are accepted (e.g., 20-40 for Altimeter).
        /// </summary>
        public virtual HashSet<int>? AllowedCollisionEnergies => null;
        /// <summary>
        /// Set of instrument types accepted by the model.
        /// null = instrument type is not applicable (skip validation).
        /// empty = instrument type IS required but any value is accepted.
        /// populated = only listed values are accepted.
        /// </summary>
        public virtual HashSet<string>? AllowedInstrumentTypes => null;
        /// <summary>
        /// Set of fragmentation types accepted by the model.
        /// null = fragmentation type is not applicable (skip validation).
        /// empty = fragmentation type IS required but any value is accepted.
        /// populated = only listed values are accepted (e.g., {"HCD", "CID"} for Prosit).
        /// </summary>
        public virtual HashSet<string>? AllowedFragmentationTypes => null;
        /// <summary>
        /// Maps mzLib modification format to monoisotopic mass differences for mass-only conversions. 
        /// The monoisotopic masses should be obtained from the UNIMOD database for accuracy (https://www.unimod.org/).
        /// While UNIMOD monoisotopic masses already account for net mass effects (such as loss of H+ for Na+ mod),
        /// this is something that should be kept in mind when adding custom modifications.
        /// </summary>
        public virtual Dictionary<string, double> ValidModificationsMonoisotopicMasses { get; protected set; } = new();
        #endregion

        #region Validation Patterns and Filters
        /// <summary>Total number of fragment ions predicted per peptide by this model. -1 if not known or dynamic.</summary>
        public virtual int NumberOfPredictedFragmentIons => -1;

        /// <summary>Minimum intensity threshold for including predicted fragments in spectral library generation</summary>
        public virtual double MinIntensityFilter { get; protected set; } = 1e-6;

        public abstract IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public abstract FragmentIonMappingMode FragmentIonMappingMode { get; init; }
        /// <summary>
        /// Reusable static instance — no M ions, no M-ion losses. This is a safeguard when creating PeptideWithSetModifications objects,
        /// such as in ResponseToPredictions() with MapToInputFullSequence mapping mode selected, since Koina does not return M-ions and 
        /// any potential issues with generated M-ions would be inherited here.
        /// </summary>
        private static readonly FragmentationParams _noMIonParams = new() { GenerateMIon = false };
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
        /// Extracts annotation, m/z, and intensity outputs from a Koina API response by matching output names.
        /// Models may return outputs in different orders; this method resolves them by name.
        /// </summary>
        protected virtual (List<object> annotations, List<object> mz, List<object> intensities) ExtractOutputs(ResponseJSONStruct response)
        {
            List<object>? annotations = null;
            List<object>? mz = null;
            List<object>? intensities = null;

            foreach (var output in response.Outputs)
            {
                if (output.Name == "annotation" || output.Name == "annotations")
                    annotations = output.Data ?? throw new Exception($"Output '{output.Name}' has null data.");
                else if (output.Name == "mz")
                    mz = output.Data ?? throw new Exception($"Output '{output.Name}' has null data.");
                else if (output.Name == "intensities")
                    intensities = output.Data ?? throw new Exception($"Output '{output.Name}' has null data.");
            }

            if (annotations == null || mz == null || intensities == null)
            {
                throw new Exception($"API response is missing expected outputs. Found: {string.Join(", ", response.Outputs.Select(o => o.Name))}. Expected: annotations, mz, intensities.");
            }

            return (annotations, mz, intensities);
        }

        /// <summary>
        /// Parses a fragment annotation string from model output into its identifier and charge components.
        /// Default implementation expects the format "type+charge" (e.g., "b5+1", "y10+2").
        /// Override in derived classes for models that use different annotation formats.
        /// </summary>
        /// <param name="annotation">Raw annotation string from model output</param>
        /// <returns>Parsed fragment identifier and charge</returns>
        /// <exception cref="Exception">Thrown when the annotation cannot be parsed</exception>
        protected virtual ParsedFragmentAnnotation ParseFragmentAnnotation(string annotation)
        {
            var plusIndex = annotation.LastIndexOf('+');
            if (plusIndex <= 0)
            {
                throw new Exception($"Cannot parse fragment annotation '{annotation}'. Expected format: 'type+charge' (e.g., 'b5+1', 'y10+2').");
            }
            var fragmentIdentifier = annotation.Substring(0, plusIndex);
            var chargeStr = annotation.Substring(plusIndex + 1);
            if (!int.TryParse(chargeStr, out int charge))
            {
                throw new Exception($"Cannot parse charge from fragment annotation '{annotation}'. Charge value: '{chargeStr}'.");
            }

            // Check for neutral loss suffix (e.g., "y3-H2O+1" -> base "y3", NL "H2O")
            var dashIndex = fragmentIdentifier.IndexOf('-');
            if (dashIndex > 0)
            {
                var baseId = fragmentIdentifier.Substring(0, dashIndex);
                var nlFormula = fragmentIdentifier.Substring(dashIndex + 1);
                try
                {
                    var nlMass = ChemicalFormula.ParseFormula(nlFormula).MonoisotopicMass;
                    return new ParsedFragmentAnnotation(fragmentIdentifier, charge, nlFormula, nlMass);
                }
                catch
                {
                    // If the formula cannot be parsed (e.g., numeric mass like "18.01", or unknown notation),
                    // fall through to return the annotation without neutral loss information.
                }
            }

            return new ParsedFragmentAnnotation(fragmentIdentifier, charge);
        }

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
            if (modelInputs.IsNullOrEmpty())
            {
                Predictions = new List<PeptideFragmentIntensityPrediction>();
                return Predictions;
            }

            ModelInputs = modelInputs;
            ValidInputsMask = new bool[ModelInputs.Count];
            var validInputs = new List<FragmentIntensityPredictionInput>();
            for (int i = 0; i < ModelInputs.Count; i++)
            {
                var cleanedSequence = TryCleanSequence(ModelInputs[i].FullSequence, out var apiSequence, out var modHandlingWarning); // mod handling happens here
                var validModelParams = ValidateModelSpecificInputs(ModelInputs[i], out var modelParametersWarning);
                if (cleanedSequence != null && apiSequence != null && validModelParams)
                {
                    ModelInputs[i] = ModelInputs[i] with { ValidatedFullSequence = apiSequence, SequenceWarning = modHandlingWarning, ParameterWarning = modelParametersWarning };
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

            var predictions = new List<PeptideFragmentIntensityPrediction>();
            if (validInputs.Count > 0)
            {
                #region Request Batching, Throttling Setup
                var batchedRequests = ToBatchedRequests(validInputs);
                var batchChunks = batchedRequests.Chunk(MaxNumberOfBatchesPerRequest).ToList();
                // We calculate a dynamic timeout based on the number of batches at (BenchmarkedTimeForOneMaxBatchSizeInMilliseconds x 2)ms/batch
                // for buffer to ensure we don't hit timeouts during processing plus throttling time.
                // Note: the time per batch is benchmarked for the entire Predict() method, so it includes some overhead beyond just the API call. Large peptide
                // requests will not necessarily scale linearly, so this is a rough estimate to provide a reasonable timeout and is an aggressive 
                // upper bound to avoid timeouts.
                int sessionTimeoutInMinutes = (int)Math.Ceiling((batchedRequests.Count * 2 * BenchmarkedTimeForOneMaxBatchSizeInMilliseconds + ThrottlingDelayInMilliseconds * batchChunks.Count) / 6e4); // 60000ms/min
                #endregion

                #region Throttled API Requests and Response Processing
                var responses = new List<string>();
                using var cts = new CancellationTokenSource(TimeSpan.FromMinutes(sessionTimeoutInMinutes)); // Bound the whole session
                for (int i = 0; i < batchChunks.Count; i++)
                {
                    var batchChunk = batchChunks[i];
                    var responseChunk = await Task.WhenAll(batchChunk.Select(request => SendInferenceRequestAsync(ModelName, request, cts.Token)));
                    responses.AddRange(responseChunk);
                    if (i < batchChunks.Count - 1) // No need to throttle after the last batch
                    {
                        await Task.Delay((int)ThrottlingDelayInMilliseconds);
                    }
                }

                predictions = ResponseToPredictions(responses, validInputs);
                #endregion
            }

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
                        ValidatedFullSequence: ModelInputs[i].ValidatedFullSequence ?? null,
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
        protected virtual List<PeptideFragmentIntensityPrediction> ResponseToPredictions(
            IReadOnlyList<string> responses, 
            List<FragmentIntensityPredictionInput> requestInputs)
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

                var (outputAnnotations, outputMZs, outputIntensities) = ExtractOutputs(response);
                var batchPeptides = requestInputs.Skip(batchIndex * MaxBatchSize).Take(MaxBatchSize).ToList();
                // Assuming outputData is structured such that each peptide's data is sequential
                if (outputAnnotations.Count % batchPeptides.Count != 0)
                {
                    throw new Exception($"Fragment annotation count ({outputAnnotations.Count}) is not evenly divisible by peptide count ({batchPeptides.Count}).");
                }
                var fragmentCount = outputAnnotations.Count / batchPeptides.Count;
                if (FragmentIonMappingMode == FragmentIonMappingMode.MapToValidatedFullSequence)
                {
                    for (int i = 0; i < batchPeptides.Count; i++)
                    {
                        var peptide = batchPeptides[i];
                        var fragmentIons = new List<string>();
                        var fragmentMZs = new List<double>();
                        var predictedIntensities = new List<double>();
                        for (int j = 0; j < fragmentCount; j++)
                        {
                            double intensity = Convert.ToDouble(outputIntensities[i * fragmentCount + j]);
                            if (intensity == -1)
                            {
                                // Skip impossible ions as indicated by the model with an intensity of -1. This is a convention used by the models to indicate that a particular fragment ion cannot be formed from the given peptide sequence and fragmentation conditions.
                                continue;
                            }

                            fragmentIons.Add(outputAnnotations[i * fragmentCount + j].ToString()!);
                            fragmentMZs.Add(Convert.ToDouble(outputMZs[i * fragmentCount + j]));
                            predictedIntensities.Add(intensity);
                        }
                        predictions.Add(new PeptideFragmentIntensityPrediction(
                            peptide.FullSequence,
                            peptide.ValidatedFullSequence!,
                            peptide.PrecursorCharge,
                            fragmentIons,
                            fragmentMZs,
                            predictedIntensities
                        ) with
                        { Warning = peptide.SequenceWarning });
                    }
                }
                else // FragmentIonMappingMode == FragmentIonMappingMode.MapToInputFullSequence
                {
                    for (int i = 0; i < batchPeptides.Count; i++)
                    {
                        var peptide = batchPeptides[i];
                        var pwsm = new PeptideWithSetModifications(peptide.FullSequence);
                        List<Product> theoreticalProducts = new();
                        pwsm.Fragment(MassSpectrometry.DissociationType.HCD, FragmentationTerminus.Both, theoreticalProducts, fragmentationParams: _noMIonParams);
                        Dictionary<string, Product> tpLookup = theoreticalProducts.DistinctBy(tp => tp.Annotation).ToDictionary(tp => tp.Annotation);

                        var fragmentIons = new List<string>();
                        var fragmentMZs = new List<double>();
                        var predictedIntensities = new List<double>();
                        for (int j = 0; j < fragmentCount; j++)
                        {
                            double intensity = Convert.ToDouble(outputIntensities[i * fragmentCount + j]);
                            if (intensity == -1)
                            {
                                // Skip impossible ions as indicated by the model with an intensity of -1. This is a convention used by the models to indicate that a particular fragment ion cannot be formed from the given peptide sequence and fragmentation conditions.
                                continue;
                            }
                            var fragmentIon = outputAnnotations[i * fragmentCount + j].ToString()!;
                            // Use the model's annotation parser (handles non-'+' charge delimiters and
                            // neutral losses) and recompute m/z like GenerateLibrarySpectraFromPredictions.
                            ParsedFragmentAnnotation parsed;
                            try
                            {
                                parsed = ParseFragmentAnnotation(fragmentIon);
                            }
                            catch
                            {
                                continue;
                            }
                            if (!tpLookup.TryGetValue(parsed.BaseFragmentIdentifier, out var tp))
                            {
                                continue;
                            }
                            double mz = parsed.NeutralLossMass.HasValue
                                ? (tp.NeutralMass - parsed.NeutralLossMass.Value).ToMz(parsed.Charge)
                                : tp.ToMz(parsed.Charge);
                            fragmentIons.Add(fragmentIon);
                            fragmentMZs.Add(mz);
                            predictedIntensities.Add(intensity);
                        }
                        predictions.Add(new PeptideFragmentIntensityPrediction(
                            peptide.FullSequence,
                            peptide.ValidatedFullSequence!,
                            peptide.PrecursorCharge,
                            fragmentIons,
                            fragmentMZs,
                            predictedIntensities
                        ) with
                        { Warning = peptide.SequenceWarning });
                    }
                }
            }
            return predictions;
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

            // Check precursor charge constraints.
            if (AllowedPrecursorCharges != null)
            {
                if (!AllowedPrecursorCharges.IsNullOrEmpty() && !AllowedPrecursorCharges.Contains(input.PrecursorCharge))
                {
                    string exceptionMessage = $"Precursor charge {input.PrecursorCharge} is not supported by this model. Allowed precursor charges: {string.Join(", ", AllowedPrecursorCharges)}.";
                    switch (ParameterHandlingMode)
                    {
                        case IncompatibleParameterHandlingMode.ThrowException:
                            throw new ArgumentException(exceptionMessage);

                        case IncompatibleParameterHandlingMode.ReturnNull:
                            warning = new WarningException(exceptionMessage);
                            return false;
                        default:
                            throw new ArgumentException($"Unhandled ParameterHandlingMode: {ParameterHandlingMode}");
                    }
                }
            }

            // Check that input has defined all required additional parameters for the model.
            if (AllowedCollisionEnergies != null && input.CollisionEnergy == null)
            {
                string exceptionMessage = "Input is missing required parameter CollisionEnergy for this model.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                    default:
                        throw new ArgumentException($"Unhandled ParameterHandlingMode: {ParameterHandlingMode}");
                }
            }

            if (AllowedInstrumentTypes != null && input.InstrumentType == null)
            {
                string exceptionMessage = "Input is missing required parameter InstrumentType for this model.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                    default:
                        throw new ArgumentException($"Unhandled ParameterHandlingMode: {ParameterHandlingMode}");
                }
            }

            if (AllowedFragmentationTypes != null && input.FragmentationType == null)
            {
                string exceptionMessage = "Input is missing required parameter FragmentationType for this model.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                    default:
                        throw new ArgumentException($"Unhandled ParameterHandlingMode: {ParameterHandlingMode}");
                }
            }

            // Check collision energy constraints
            if (AllowedCollisionEnergies != null && input.CollisionEnergy != null && !AllowedCollisionEnergies.IsNullOrEmpty() && !AllowedCollisionEnergies.Contains(input.CollisionEnergy.Value))
            {
                string exceptionMessage = $"Collision energy {input.CollisionEnergy} is not supported by this model. Allowed collision energies: {string.Join(", ", AllowedCollisionEnergies)}.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                    default:
                        throw new ArgumentException($"Unhandled ParameterHandlingMode: {ParameterHandlingMode}");
                }
            }

            // Check instrument type constraints
            if (AllowedInstrumentTypes != null && input.InstrumentType != null && !AllowedInstrumentTypes.IsNullOrEmpty() && !AllowedInstrumentTypes.Contains(input.InstrumentType))
            {
                string exceptionMessage = $"Instrument type '{input.InstrumentType}' is not supported by this model. Allowed instrument types: {string.Join(", ", AllowedInstrumentTypes)}.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                    default:
                        throw new ArgumentException($"Unhandled ParameterHandlingMode: {ParameterHandlingMode}");
                }
            }

            // Check fragmentation type constraints
            if (AllowedFragmentationTypes != null && input.FragmentationType != null && !AllowedFragmentationTypes.IsNullOrEmpty() && !AllowedFragmentationTypes.Contains(input.FragmentationType))
            {
                string exceptionMessage = $"Fragmentation type '{input.FragmentationType}' is not supported by this model. Allowed fragmentation types: {string.Join(", ", AllowedFragmentationTypes)}.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                    default:
                        throw new ArgumentException($"Unhandled ParameterHandlingMode: {ParameterHandlingMode}");
                }
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
        public List<LibrarySpectrum> GenerateLibrarySpectraFromPredictions(double?[] alignedRetentionTimes, out WarningException? warning, string? filepath=null, double minIntensityFilter=1e-4)
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
            List<PeptideFragmentIntensityPrediction> predictions = new();
            List<double?> rts = new();

            for (int i = 0; i < Predictions.Count; i++)
            {
                if (ValidInputsMask[i])
                {
                    predictions.Add(Predictions[i]);
                    rts.Add(alignedRetentionTimes[i]);
                }
            }

            var predictedSpectra = new List<LibrarySpectrum>();
            for (int predictionIndex = 0; predictionIndex < predictions.Count; predictionIndex++)
            {
                var prediction = predictions[predictionIndex];

                PeptideWithSetModifications peptide = FragmentIonMappingMode == FragmentIonMappingMode.MapToValidatedFullSequence 
                    ? new PeptideWithSetModifications(prediction.ValidatedFullSequence) 
                    : new PeptideWithSetModifications(prediction.FullSequence);
                List<MatchedFragmentIon> fragmentIons = new();

                List<Product> theoreticalProducts = new();
                peptide.Fragment(MassSpectrometry.DissociationType.HCD, FragmentationTerminus.Both, theoreticalProducts, fragmentationParams: _noMIonParams); // Generate theoretical fragments to get the m/z values for the input sequence
                Dictionary<string, double> predictionAnnotationIntensityLookup = new();
                Dictionary<string, Product> tpLookup = theoreticalProducts.DistinctBy(tp => tp.Annotation).ToDictionary(tp => tp.Annotation);

                for (int i = 0; i < prediction.FragmentAnnotations.Count; i++)
                {
                    if (prediction.FragmentIntensities[i] == -1 || 
                        prediction.FragmentIntensities[i] < minIntensityFilter ||
                        prediction.FragmentAnnotations[i] == null)
                    {
                        // Skip impossible ions (intensity == -1) and peaks with near zero intensity.
                        // The model uses -1 to indicate impossible ions.
                        continue;
                    }
                    predictionAnnotationIntensityLookup[prediction.FragmentAnnotations[i]] = prediction.FragmentIntensities[i];
                }

                foreach (var pa in predictionAnnotationIntensityLookup.Keys)
                {
                    ParsedFragmentAnnotation parsed;
                    try
                    {
                        parsed = ParseFragmentAnnotation(pa);
                    }
                    catch
                    {
                        continue;
                    }
                    if (!tpLookup.TryGetValue(parsed.BaseFragmentIdentifier, out var tp))
                    {
                        continue;
                    }
                    var charge = parsed.Charge;
                    double experMz = parsed.NeutralLossMass.HasValue
                        ? (tp.NeutralMass - parsed.NeutralLossMass.Value).ToMz(charge)
                        : tp.ToMz(charge);

                    var fragmentIon = new MatchedFragmentIon
                    (
                        neutralTheoreticalProduct: tp,
                        experMz: experMz,
                        experIntensity: predictionAnnotationIntensityLookup[pa],
                        charge: charge
                    );

                    fragmentIons.Add(fragmentIon);
                }

                var spectrum = new LibrarySpectrum
                (
                    sequence: prediction.FullSequence,
                    precursorMz: peptide.ToMz(prediction.PrecursorCharge),
                    chargeState: prediction.PrecursorCharge,
                    peaks: fragmentIons,
                    rt: rts[predictionIndex]
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
