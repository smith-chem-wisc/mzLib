using MzLibUtil;
using Omics.SpectrumMatch;
using Predictions.Koina.Interfaces;
using System.ComponentModel;
using System.Text;
using System.Text.RegularExpressions;
using Predictions.Koina.Client;

namespace Predictions.Koina.SupportedModels
{
    public class Prosit2020iRTTMT : IKoinaModelIO
    {
        public string ModelName => "Prosit_2020_irt_TMT";
        public int MaxBatchSize => 1000;
        public int MaxPeptideLength => 30;
        public Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            {"[Common Variable:Oxidation on M]", "[UNIMOD:35]"},
            {"[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]"},
            // SILAC modifications
            {"[Common Variable:Label:13C(6)15N(2) on K]", "[UNIMOD:259]"},
            {"[Common Variable:Label:13C(6)15N(4) on R]", "[UNIMOD:267]"},
            // TMT modifications
            {"[Common Fixed:TMT6plex on K]", "[UNIMOD:737]"},
            {"[Common Fixed:TMT6plex on N-terminus]", "[UNIMOD:737]-"},
            // TMTpro modifications
            {"[Common Fixed:TMTpro on K]", "[UNIMOD:2016]"},
            {"[Common Fixed:TMTpro on N-terminus]", "[UNIMOD:2016]-"},
            // iTRAQ 4-plex modifications
            {"[Common Fixed:iTRAQ4plex on K]", "[UNIMOD:214]"},
            {"[Common Fixed:iTRAQ4plex on N-terminus]", "[UNIMOD:214]-"},
            // iTRAQ 8-plex modifications
            {"[Common Fixed:iTRAQ8plex on K]", "[UNIMOD:730]"},
            {"[Common Fixed:iTRAQ8plex on N-terminus]", "[UNIMOD:730]-"}
        };
        public string ModificationPattern => @"\[[^\]]+\]";
        public string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        public readonly List<string> PeptideSequences = new();
        public List<double> PredictedIndexedRetentionTimes = new();

        public Prosit2020iRTTMT(List<string> peptideSequences, out WarningException? warnings)
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
            var firstModIsValid = ValidModificationUnimodMapping.TryGetValue(matches.Where(m => m.Index == 0).First().Value, out var nTermMod);

            return matches.All(m => ValidModificationUnimodMapping.ContainsKey(m.Value))
                && firstModIsValid;
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
