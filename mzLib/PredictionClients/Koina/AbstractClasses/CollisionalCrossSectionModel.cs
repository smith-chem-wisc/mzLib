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
    /// Represents a collisional cross section prediction result for a single peptide.
    /// </summary>
    public record PeptideCCSPrediction(
        string FullSequence,
        string ValidatedFullSequence,
        int PrecursorCharge,
        double? PredictedCCS,
        WarningException? Warning = null
    );

    /// <summary>
    /// Input parameters for CCS prediction models.
    /// </summary>
    public record CCSPredictionInput(
        string FullSequence,
        int PrecursorCharge
    )
    {
        public string? ValidatedFullSequence { get; set; }
        public WarningException? SequenceWarning { get; set; }
        public WarningException? ParameterWarning { get; set; }
    }

    /// <summary>
    /// Abstract base class for collisional cross section (CCS) prediction models using the Koina API.
    /// </summary>
    public abstract class CollisionalCrossSectionModel : KoinaModelBase<CCSPredictionInput, PeptideCCSPrediction>, IPredictor<CCSPredictionInput, PeptideCCSPrediction>
    {
        protected CollisionalCrossSectionModel(ISequenceConverter sequenceConverter)
            : base(sequenceConverter)
        {
        }

        public virtual HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };
        public override IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>();
        public override abstract SequenceConversionHandlingMode ModHandlingMode { get; init; }
        public virtual IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }

        public List<CCSPredictionInput> ModelInputs { get; protected set; } = new();
        public bool[] ValidInputsMask { get; protected set; } = Array.Empty<bool>();
        public List<PeptideCCSPrediction> Predictions { get; protected set; } = new();

        protected virtual async Task<List<PeptideCCSPrediction>> AsyncThrottledPredictor(List<CCSPredictionInput> modelInputs)
        {
            if (modelInputs.IsNullOrEmpty())
            {
                Predictions = new List<PeptideCCSPrediction>();
                return Predictions;
            }

            ModelInputs = modelInputs;
            ValidInputsMask = new bool[ModelInputs.Count];
            var validInputs = new List<CCSPredictionInput>();

            for (int i = 0; i < ModelInputs.Count; i++)
            {
                var cleanedSequence = TryCleanSequence(ModelInputs[i].FullSequence, out var apiSequence, out var modHandlingWarning);
                var validModelParams = ValidateModelSpecificInputs(ModelInputs[i], out var parameterWarning);
                if (cleanedSequence != null && apiSequence != null && validModelParams)
                {
                    ModelInputs[i] = ModelInputs[i] with { ValidatedFullSequence = apiSequence, SequenceWarning = modHandlingWarning, ParameterWarning = parameterWarning };
                    ValidInputsMask[i] = true;
                    validInputs.Add(ModelInputs[i]);
                }
                else
                {
                    ModelInputs[i] = ModelInputs[i] with { ValidatedFullSequence = null, SequenceWarning = modHandlingWarning, ParameterWarning = parameterWarning };
                    ValidInputsMask[i] = false;
                }
            }

            var predictions = new List<PeptideCCSPrediction>();
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

            var realignedPredictions = new List<PeptideCCSPrediction>();
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
                    realignedPredictions.Add(new PeptideCCSPrediction(
                        FullSequence: ModelInputs[i].FullSequence,
                        ValidatedFullSequence: null,
                        PrecursorCharge: ModelInputs[i].PrecursorCharge,
                        PredictedCCS: null,
                        Warning: ModelInputs[i].SequenceWarning ?? new WarningException("Input was invalid and skipped during prediction.")
                    ));
                }
            }

            Predictions = realignedPredictions;
            return Predictions;
        }

        public List<PeptideCCSPrediction> Predict(List<CCSPredictionInput> modelInputs)
        {
            return AsyncThrottledPredictor(modelInputs).GetAwaiter().GetResult();
        }

        protected virtual bool ValidateModelSpecificInputs(CCSPredictionInput input, out WarningException? warning)
        {
            warning = null;

            // Reject non-physical charges even when the model declares no charge constraint.
            // Koina rejects charge <= 0 server-side; catching it here gives a clearer client-side warning.
            if (input.PrecursorCharge <= 0)
            {
                string exceptionMessage = $"Precursor charge {input.PrecursorCharge} is not physical; precursor charge must be a positive integer.";
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

            return true;
        }

        protected virtual List<PeptideCCSPrediction> ResponseToPredictions(
            IReadOnlyList<string> responses,
            List<CCSPredictionInput> requestInputs)
        {
            var predictions = new List<PeptideCCSPrediction>();
            if (requestInputs.IsNullOrEmpty())
            {
                return predictions;
            }

            var deserializedResponses = responses.Select(r => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(r)).ToList();
            if (deserializedResponses.IsNullOrEmpty() || deserializedResponses.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }

            var ccsOutputs = deserializedResponses.SelectMany(r =>
            {
                var output = r!.Outputs.FirstOrDefault(o => o.Name == "ccs")
                    ?? throw new Exception($"API response is missing expected output 'ccs'. Found: {string.Join(", ", r.Outputs.Select(o => o.Name))}.");
                return output.Data.Select(d => Convert.ToDouble(d));
            }).ToList();

            if (ccsOutputs.Count != requestInputs.Count)
            {
                throw new Exception("The number of CCS predictions does not match the number of input peptides.");
            }

            for (int i = 0; i < requestInputs.Count; i++)
            {
                predictions.Add(new PeptideCCSPrediction(
                    FullSequence: requestInputs[i].FullSequence,
                    ValidatedFullSequence: requestInputs[i].ValidatedFullSequence,
                    PrecursorCharge: requestInputs[i].PrecursorCharge,
                    PredictedCCS: ccsOutputs[i],
                    Warning: requestInputs[i].SequenceWarning
                ));
            }

            return predictions;
        }
    }
}
