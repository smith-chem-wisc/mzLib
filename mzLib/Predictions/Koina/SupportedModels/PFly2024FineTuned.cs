using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Predictions.Koina.Interfaces;
using Predictions.Koina.Client;
using Easy.Common.Extensions;
using MzLibUtil;
using System.Text.RegularExpressions;
using System.ComponentModel;

namespace Predictions.Koina.SupportedModels
{
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
        public List<List<double>> DetectabilityProbabilityTable = new();
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

        public List<Dictionary<string, object>> ToBatchedRequests()
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
                ResponseToDetectabilityProbabilityTable(responses);
            }
            finally
            {
                _http.Dispose();
            }
        }
        public void ResponseToDetectabilityProbabilityTable(string[] responses)
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

            DetectabilityProbabilityTable = deserializedResponses.Select(batch => batch!.Outputs[0].Data)
                .Chunk(NumberOfDetectabilityClasses)
                .Select(classProbabilities => classProbabilities.Select(d=>Convert.ToDouble(d)).ToList())
                .ToList();
        }
        public bool IsValidPeptideSequence(string sequence)
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
