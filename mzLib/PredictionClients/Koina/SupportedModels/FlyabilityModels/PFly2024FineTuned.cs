using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using Easy.Common.Extensions;
using MzLibUtil;
using System.Text.RegularExpressions;
using System.ComponentModel;


namespace PredictionClients.Koina.SupportedModels.FlyabilityModels
{
    public record PeptideDetectabilityPrediction(
        string PeptideSequence, 
        double NotDetectable,
        double LowDetectability,
        double IntermediateDetectability,
        double HighDetectability
    );

    public class PFly2024FineTuned : IKoinaModelIO
    {
        public string ModelName => "pfly_2024_fine_tuned";
        public int MaxBatchSize => 128;
        public int NumberOfDetectabilityClasses => 4;
        public List<string> DetectabilityClasses => new() { "Not Detectable", "Low Detectability", "Intermediate Detectability", "High Detectability" };
        public int MaxPeptideLength => 40;
        public string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        public string ModificationPattern => @"\[[^\]]+\]";

        public readonly List<string> PeptideSequences = new();
        public List<PeptideDetectabilityPrediction> Predictions = new();
        public PFly2024FineTuned(List<string> peptideSequences, out WarningException warnings)
        {
            if (peptideSequences.IsNullOrEmpty())
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            var invalidSequences = new List<string>();
            //Ensure all peptides are valid
            foreach (var seq in peptideSequences)
            {
                if (IsValidPeptideSequence(seq))
                {
                    PeptideSequences.Add(seq);
                }
                else
                {
                    invalidSequences.Add(seq);
                }
            }

            if (invalidSequences.Count > 0)
            {
                warnings = new WarningException($"The following peptide sequences were invalid and will be skipped: {string.Join(", ", invalidSequences)}");
            }
            else
            {
                warnings = null;
            }
        }

        protected List<Dictionary<string, object>> ToBatchedRequests()
        {
            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
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
                        }
                    }
                };
                batchedRequests.Add(request);
            }
            return batchedRequests;
        }
        public async Task RunInferenceAsync()
        {
            var _http = new HTTP(timeoutInMinutes: PeptideSequences.Count / MaxBatchSize * 2 + 2); // Typically a full batch takes about a minute. Setting it to double that for safety.

            try
            {
                var responses = await Task.WhenAll(ToBatchedRequests().Select(request => _http.InferenceRequest(ModelName, request)));
                ResponseToPredictions(responses);
            }
            finally
            {
                _http.Dispose();
            }
        }
        protected void ResponseToPredictions(string[] responses)
        {
            if (PeptideSequences.Count == 0)
            {
                return;
            }
            
            var deserializedResponses = responses.Select(r => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(r)).ToList();

            if (deserializedResponses.IsNullOrEmpty() || deserializedResponses.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }

            var detectabilityPredictions = deserializedResponses.SelectMany(batch => batch!.Outputs[0].Data.Chunk(NumberOfDetectabilityClasses)).ToList();
            for (int i = 0; i < PeptideSequences.Count; i++)
            {
                var probs = detectabilityPredictions[i];
                Predictions.Add(new PeptideDetectabilityPrediction(
                    PeptideSequences[i],
                    Convert.ToDouble(probs[0]),
                    Convert.ToDouble(probs[1]),
                    Convert.ToDouble(probs[2]),
                    Convert.ToDouble(probs[3])
                ));
            }
        }
        protected bool IsValidPeptideSequence(string sequence)
        {
            var unmodifiedSequence = Regex.Replace(sequence, ModificationPattern, string.Empty);
            if (unmodifiedSequence.Length > MaxPeptideLength)
            {
                return false;
            }
            if (!Regex.IsMatch(unmodifiedSequence, CanonicalAminoAcidPattern))
            {
                return false;
            }
            return true;
        }
    }
}
