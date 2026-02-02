using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace PredictionClients.Koina.AbstractClasses
{
    public record PeptideRTPrediction(
        string PeptideSequence,
        double PredictedRetentionTime,
        bool IsIndexed
    );

    public abstract class RetentionTimeModel : IKoinaModelIO
    {
        public abstract string ModelName { get; }
        public abstract int MaxBatchSize { get; }
        public abstract int MaxPeptideLength { get; }
        public abstract int MinPeptideLength { get; }
        public abstract bool IsIndexedRetentionTimeModel { get; }
        public virtual Dictionary<string, string> ValidModificationUnimodMapping => new();
        public virtual string ModificationPattern => @"\[[^\]]+\]";
        public virtual string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        public abstract List<string> PeptideSequences { get; }
        public abstract List<PeptideRTPrediction> Predictions { get; protected set; }

        public abstract List<Dictionary<string, object>> ToBatchedRequests();

        public virtual async Task RunInferenceAsync()
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

        public virtual void ResponseToPredictions(string[] responses)
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

            var rtOutputs = deserializedResponses.SelectMany(r => r.Outputs[0].Data).ToList();
            if (rtOutputs.Count != PeptideSequences.Count)
            {
                throw new Exception("The number of predictions does not match the number of input peptides.");
            }
            Predictions = PeptideSequences
                .Select((seq, index) => new PeptideRTPrediction(
                    PeptideSequence: seq,
                    PredictedRetentionTime: Convert.ToDouble(rtOutputs[index]),
                    IsIndexed: IsIndexedRetentionTimeModel
                ))
                .ToList();
        }

        public virtual bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            if (matches.Count == 0)
            {
                return true; // No modifications, valid
            }

            return matches.Where(m => !ValidModificationUnimodMapping.ContainsKey(m.Value)).Count() == 0;
        }

        public virtual bool IsValidSequence(string sequence)
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
        public virtual string ConvertToPrositModificationFormat(string sequence)
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
