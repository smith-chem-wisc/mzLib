using System.Text.RegularExpressions;
using Omics.SequenceConversion;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using PredictionClients.Koina.Util;
using System.ComponentModel;
using MzLibUtil;
using Easy.Common.Extensions;

namespace PredictionClients.Koina.AbstractClasses
{
    /// <summary>
    /// Represents crosslink fragment intensity prediction results.
    /// <see cref="BetaSequence"/> and <see cref="ValidatedBetaSequence"/> are null for
    /// single-sequence models (e.g. Prosit_2023_intensity_XL_CMS3) where the alpha sequence
    /// already encodes the full crosslinked pair.
    /// </summary>
    public record CrosslinkFragmentIntensityPrediction(
        string AlphaSequence,
        string? BetaSequence,
        string ValidatedAlphaSequence,
        string? ValidatedBetaSequence,
        int PrecursorCharge,
        List<string>? FragmentAnnotations,
        List<double>? FragmentMZs,
        List<double>? FragmentIntensities,
        WarningException? Warning = null
    );

    /// <summary>
    /// Input parameters for crosslink fragment intensity prediction.
    /// <see cref="BetaSequence"/> may be null for single-sequence models that take a
    /// pre-combined crosslinked sequence in <see cref="AlphaSequence"/>.
    /// </summary>
    public record CrosslinkIntensityPredictionInput(
        string AlphaSequence,
        string? BetaSequence,
        int PrecursorCharge,
        int? CollisionEnergy
    )
    {
        public string? ValidatedAlphaSequence { get; set; }
        public string? ValidatedBetaSequence { get; set; }
        public WarningException? AlphaSequenceWarning { get; set; }
        public WarningException? BetaSequenceWarning { get; set; }
        public WarningException? ParameterWarning { get; set; }
    }

    /// <summary>
    /// Abstract base class for crosslink fragment intensity prediction models using the Koina API.
    /// </summary>
    public abstract class CrosslinkFragmentIntensityModel : KoinaModelBase<CrosslinkIntensityPredictionInput, CrosslinkFragmentIntensityPrediction>, IPredictor<CrosslinkIntensityPredictionInput, CrosslinkFragmentIntensityPrediction>
    {
        protected CrosslinkFragmentIntensityModel(ISequenceConverter sequenceConverter)
            : base(sequenceConverter)
        {
        }

        public virtual HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public virtual HashSet<int> AllowedCollisionEnergies => new HashSet<int>(); // Koina accepts any FP32 collision energy
        /// <summary>
        /// Whether this Koina model requires a separate beta sequence input.
        /// Paired-sequence models (CMS2, NMS2) leave this true; single-sequence helpers
        /// (CMS3) override to false and consume only <see cref="CrosslinkIntensityPredictionInput.AlphaSequence"/>.
        /// </summary>
        public virtual bool RequiresBetaSequence => true;
        public override abstract IReadOnlySet<int> AllowedUnimodIds { get; }
        public override abstract SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public abstract IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }
        public abstract FragmentIonMappingMode FragmentIonMappingMode { get; init; }

        protected static readonly UnimodSequenceFormatSchema CrosslinkSchema = UnimodSequenceFormatSchema.Instance;

        public List<CrosslinkIntensityPredictionInput> ModelInputs { get; protected set; } = new();
        public bool[] ValidInputsMask { get; protected set; } = Array.Empty<bool>();
        public List<CrosslinkFragmentIntensityPrediction> Predictions { get; protected set; } = new();

        protected virtual async Task<List<CrosslinkFragmentIntensityPrediction>> AsyncThrottledPredictor(List<CrosslinkIntensityPredictionInput> modelInputs)
        {
            if (modelInputs.IsNullOrEmpty())
            {
                Predictions = new List<CrosslinkFragmentIntensityPrediction>();
                return Predictions;
            }

            ModelInputs = modelInputs;
            ValidInputsMask = new bool[ModelInputs.Count];
            var validInputs = new List<CrosslinkIntensityPredictionInput>();

            for (int i = 0; i < ModelInputs.Count; i++)
            {
                var cleanedAlpha = TryCleanSequence(ModelInputs[i].AlphaSequence, out var apiAlpha, out var alphaWarning);

                string? cleanedBeta = null;
                string? apiBeta = null;
                WarningException? betaWarning = null;
                bool betaOk;
                if (RequiresBetaSequence)
                {
                    if (ModelInputs[i].BetaSequence == null)
                    {
                        betaOk = false;
                        betaWarning = new WarningException("Beta sequence is required by this model but was null.");
                    }
                    else
                    {
                        cleanedBeta = TryCleanSequence(ModelInputs[i].BetaSequence!, out apiBeta, out betaWarning);
                        betaOk = cleanedBeta != null && apiBeta != null;
                    }
                }
                else
                {
                    // Single-sequence model: beta is not consumed; accept anything (including null).
                    betaOk = true;
                }

                var validModelParams = ValidateModelSpecificInputs(ModelInputs[i], out var parameterWarning);

                if (cleanedAlpha != null && apiAlpha != null && betaOk && validModelParams)
                {
                    ModelInputs[i] = ModelInputs[i] with
                    {
                        ValidatedAlphaSequence = apiAlpha,
                        ValidatedBetaSequence = apiBeta,
                        AlphaSequenceWarning = alphaWarning,
                        BetaSequenceWarning = betaWarning,
                        ParameterWarning = parameterWarning
                    };
                    ValidInputsMask[i] = true;
                    validInputs.Add(ModelInputs[i]);
                }
                else
                {
                    ModelInputs[i] = ModelInputs[i] with
                    {
                        ValidatedAlphaSequence = null,
                        ValidatedBetaSequence = null,
                        AlphaSequenceWarning = alphaWarning,
                        BetaSequenceWarning = betaWarning,
                        ParameterWarning = parameterWarning
                    };
                    ValidInputsMask[i] = false;
                }
            }

            var predictions = new List<CrosslinkFragmentIntensityPrediction>();
            if (validInputs.Count > 0)
            {
                var batchedRequests = ToBatchedRequests(validInputs);
                var batchChunks = batchedRequests.Chunk(MaxNumberOfBatchesPerRequest).ToList();
                int sessionTimeoutInMinutes = (int)Math.Ceiling((batchedRequests.Count * 2 * BenchmarkedTimeForOneMaxBatchSizeInMilliseconds + ThrottlingDelayInMilliseconds * batchChunks.Count) / 6e4);
                sessionTimeoutInMinutes = Math.Max(sessionTimeoutInMinutes, 1);

                var responses = new List<string>();
                using var _http = new HTTP(timeoutInMinutes: sessionTimeoutInMinutes);

                for (int i = 0; i < batchChunks.Count; i++)
                {
                    var batchChunk = batchChunks[i];
                    var responseChunk = await Task.WhenAll(batchChunk.Select(request => _http.InferenceRequest(ModelName, request)));
                    responses.AddRange(responseChunk);

                    if (i < batchChunks.Count - 1)
                    {
                        await Task.Delay(ThrottlingDelayInMilliseconds);
                    }
                }

                predictions = ResponseToPredictions(responses, validInputs);
            }

            var realignedPredictions = new List<CrosslinkFragmentIntensityPrediction>();
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
                    realignedPredictions.Add(new CrosslinkFragmentIntensityPrediction(
                        AlphaSequence: ModelInputs[i].AlphaSequence,
                        BetaSequence: ModelInputs[i].BetaSequence,
                        ValidatedAlphaSequence: null,
                        ValidatedBetaSequence: null,
                        PrecursorCharge: ModelInputs[i].PrecursorCharge,
                        FragmentAnnotations: null,
                        FragmentMZs: null,
                        FragmentIntensities: null,
                        Warning: ModelInputs[i].AlphaSequenceWarning ?? ModelInputs[i].BetaSequenceWarning ?? new WarningException("Input was invalid and skipped during prediction.")
                    ));
                }
            }

            Predictions = realignedPredictions;
            return Predictions;
        }

        public List<CrosslinkFragmentIntensityPrediction> Predict(List<CrosslinkIntensityPredictionInput> modelInputs)
        {
            return AsyncThrottledPredictor(modelInputs).GetAwaiter().GetResult();
        }

        /// <summary>
        /// Override to preserve crosslink modification annotations in the API sequence.
        /// The standard converter would strip the crosslink UNIMOD markers (e.g. UNIMOD:1881,
        /// UNIMOD:1896, UNIMOD:1898) because they aren't in the mzLib local mod database — but
        /// Koina's helper models require those markers in place to determine the crosslink
        /// position. Rather than bypass validation entirely, this override:
        ///   1. Validates the bare amino-acid sequence and length bounds.
        ///   2. Requires every bracketed annotation to be in UNIMOD:N notation.
        ///   3. Rejects any UNIMOD id not in <see cref="AllowedUnimodIds"/>.
        /// </summary>
        protected override string? TryCleanSequence(
            string sequence,
            out string? apiSequence,
            out WarningException? warning)
        {
            apiSequence = null;
            warning = null;

            var rawBase = BaseStripper.Replace(sequence, string.Empty);
            if (!Regex.IsMatch(rawBase, AllowedAminoAcidPattern))
            {
                HandleFailure(ModHandlingMode, "Invalid base sequence.");
                return null;
            }

            if (!IsValidBaseSequence(rawBase, AllowedAminoAcidPattern, MinPeptideLength, MaxPeptideLength))
            {
                HandleFailure(ModHandlingMode, "Invalid base sequence.");
                return null;
            }

            // Validate every bracketed annotation. Crosslink models expect Koina-format
            // UNIMOD:N notation (e.g. K[UNIMOD:1896]); mzLib-format notation
            // (e.g. K[Common Variable:Acetyl on K]) would otherwise be forwarded to Koina
            // verbatim and rejected with an opaque server error.
            var disallowed = new List<string>();
            foreach (Match m in BaseStripper.Matches(sequence))
            {
                var inner = m.Value.Substring(1, m.Value.Length - 2); // strip [ and ]
                if (!inner.StartsWith("UNIMOD:", StringComparison.OrdinalIgnoreCase)
                    || !int.TryParse(inner.AsSpan(7), out var id)
                    || !AllowedUnimodIds.Contains(id))
                {
                    disallowed.Add(inner);
                }
            }

            if (disallowed.Count > 0)
            {
                string message = $"Crosslink sequence contains unsupported modification(s): {string.Join(", ", disallowed)}. " +
                                 $"This model accepts only UNIMOD:N notation with ids in {{{string.Join(", ", AllowedUnimodIds)}}}.";
                HandleFailure(ModHandlingMode, message);
                warning = new WarningException(message);
                return null;
            }

            apiSequence = sequence;
            return apiSequence;
        }

        protected virtual bool ValidateModelSpecificInputs(CrosslinkIntensityPredictionInput input, out WarningException? warning)
        {
            warning = null;

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
                }
            }

            if (!AllowedCollisionEnergies.IsNullOrEmpty() && input.CollisionEnergy == null)
            {
                string exceptionMessage = "Input is missing required parameter CollisionEnergy for this model.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                }
            }

            if (!AllowedCollisionEnergies.IsNullOrEmpty() && input.CollisionEnergy != null && !AllowedCollisionEnergies.Contains(input.CollisionEnergy.Value))
            {
                string exceptionMessage = $"Collision energy {input.CollisionEnergy} is not supported by this model. Allowed collision energies: {string.Join(", ", AllowedCollisionEnergies)}.";
                switch (ParameterHandlingMode)
                {
                    case IncompatibleParameterHandlingMode.ThrowException:
                        throw new ArgumentException(exceptionMessage);
                    case IncompatibleParameterHandlingMode.ReturnNull:
                        warning = new WarningException(exceptionMessage);
                        return false;
                }
            }

            return true;
        }

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

        protected virtual List<CrosslinkFragmentIntensityPrediction> ResponseToPredictions(
            IReadOnlyList<string> responses,
            List<CrosslinkIntensityPredictionInput> requestInputs)
        {
            var predictions = new List<CrosslinkFragmentIntensityPrediction>();
            if (requestInputs.IsNullOrEmpty())
            {
                return predictions;
            }

            var deserializedResponses = responses.Select(r => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(r)).ToList();
            if (deserializedResponses.IsNullOrEmpty() || deserializedResponses.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }

            for (int batchIndex = 0; batchIndex < deserializedResponses.Count; batchIndex++)
            {
                var response = deserializedResponses[batchIndex];
                if (response == null || response.Outputs.Count != 3)
                {
                    throw new Exception($"API response is not in the expected format. Expected 3 outputs, got {response?.Outputs.Count}.");
                }

                var (outputAnnotations, outputMZs, outputIntensities) = ExtractOutputs(response);

                var batchInputs = requestInputs.Skip(batchIndex * MaxBatchSize).Take(MaxBatchSize).ToList();
                if (outputAnnotations.Count % batchInputs.Count != 0)
                {
                    throw new Exception($"Fragment annotation count ({outputAnnotations.Count}) is not evenly divisible by peptide count ({batchInputs.Count}).");
                }
                var fragmentCount = outputAnnotations.Count / batchInputs.Count;

                for (int i = 0; i < batchInputs.Count; i++)
                {
                    var input = batchInputs[i];
                    var fragmentIons = new List<string>();
                    var fragmentMZs = new List<double>();
                    var predictedIntensities = new List<double>();

                    for (int j = 0; j < fragmentCount; j++)
                    {
                        double intensity = Convert.ToDouble(outputIntensities[i * fragmentCount + j]);
                        if (intensity == -1)
                        {
                            continue;
                        }

                        fragmentIons.Add(outputAnnotations[i * fragmentCount + j].ToString()!);
                        fragmentMZs.Add(Convert.ToDouble(outputMZs[i * fragmentCount + j]));
                        predictedIntensities.Add(intensity);
                    }

                    predictions.Add(new CrosslinkFragmentIntensityPrediction(
                        AlphaSequence: input.AlphaSequence,
                        BetaSequence: input.BetaSequence,
                        ValidatedAlphaSequence: input.ValidatedAlphaSequence!,
                        ValidatedBetaSequence: input.ValidatedBetaSequence,
                        PrecursorCharge: input.PrecursorCharge,
                        FragmentAnnotations: fragmentIons,
                        FragmentMZs: fragmentMZs,
                        FragmentIntensities: predictedIntensities,
                        Warning: input.AlphaSequenceWarning ?? input.BetaSequenceWarning
                    ));
                }
            }

            return predictions;
        }

        public abstract int NumberOfPredictedFragmentIons { get; }
    }
}
