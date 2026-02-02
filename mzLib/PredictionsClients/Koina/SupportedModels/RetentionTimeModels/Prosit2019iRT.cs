using Easy.Common.Extensions;
using MzLibUtil;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;


namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    public class Prosit2019iRT : IKoinaModelIO
    {
        public string ModelName => "Prosit_2019_irt";
        public int MaxBatchSize => 1000;
        public int MaxPeptideLength => 30;
        public Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            {"[Common Variable:Oxidation on M]", "[UNIMOD:35]"},
            {"[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]"}
        };
        public string ModificationPattern => @"\[[^\]]+\]";
        public string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        public readonly List<string> PeptideSequences = new();
        public List<double> PredictedIndexedRetentionTimes = new();
        public Prosit2019iRT(List<string> peptideSequences, out WarningException? warnings)
        {
            if (peptideSequences.IsNullOrEmpty())
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            //Ensure all peptides are valid
            var invalidSequences = new List<string>();
            foreach (var seq in peptideSequences)
            {
                if (!IsValidSequence(seq) || !HasValidModifications(seq))
                {
                    invalidSequences.Add(seq);
                }
                else
                {
                    PeptideSequences.Add(ConvertToPrositModificationFormat(seq));
                }
            }

            warnings = null;
            if (invalidSequences.Count > 0)
            {
                var sb = new StringBuilder();
                sb.AppendLine("The following peptide sequences were invalid and will be skipped:");
                foreach (var invalid in invalidSequences)
                {
                    sb.AppendLine(invalid);
                }
                warnings = new WarningException(sb.ToString());
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
                ResponseToiRTs(responses);
            }
            finally
            {
                _http.Dispose();
            }
        }

        public void ResponseToiRTs(string[] responses)
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

            PredictedIndexedRetentionTimes = deserializedResponses
                .SelectMany(batch => batch!.Outputs[0].Data)
                .Select(irt => Convert.ToDouble(irt))
                .ToList();
        }

        public bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            if (matches.Count == 0)
            {
                return true; // No modifications, valid
            }

            return matches.Where(m => !ValidModificationUnimodMapping.ContainsKey(m.Value)).Count() == 0;
        }

        public bool IsValidSequence(string sequence)
        {
            var baseSequence = Regex.Replace(sequence, ModificationPattern, "");
            return Regex.IsMatch(baseSequence, CanonicalAminoAcidPattern)
                && baseSequence.Length <= MaxPeptideLength;
        }

        /// <summary>
        /// Converts a peptide sequence from the mzLib modification format to the Prosit UNIMOD format.
        /// By default, all unmodified cysteines are carbamidomethylated (UNIMOD:4) to match the expectations
        /// of the Prosit 2020 HCD intensity model.
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib modification format.</param>
        /// <returns>The sequence converted to Prosit UNIMOD format with cysteines carbamidomethylated.</returns>
        internal string ConvertToPrositModificationFormat(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Key, mod.Value);
            }

            // Carbamidomethylate all Cysteines if not already modified
            return Regex.Replace(sequence, @"C(?!\[UNIMOD:4\])", "C[UNIMOD:4]");
        }
    }
}
