using MzLibUtil;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using System.Text.RegularExpressions;

namespace PredictionClients.Koina.AbstractClasses
{
    public record PeptideFragmentIntensityPrediction(
        string PeptideSequence, 
        List<string> FragmentAnnotations, 
        List<double> FragmentMZs,
        List<double> FragmentIntensities
    );

    public abstract class FragmentIntensityModel : IKoinaModelIO
    {
        public abstract string ModelName { get; }
        public abstract int MaxBatchSize { get; }
        public abstract int MaxPeptideLength { get; }
        public abstract int MinPeptideLength { get; }
        public virtual string CanonicalAminoAcidPattern { get; protected set; } = @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        public virtual string ModificationPattern { get; protected set; } = @"\[[^\]]+\]";
        public virtual Dictionary<string, string> ValidModificationUnimodMapping { get; protected set; } = new();
        public virtual Dictionary<string, double> ValidModificationsMonoisotopicMasses { get; protected set; } = new();
        public abstract List<string> PeptideSequences { get; }
        public abstract List<int> PrecursorCharges { get; }
        public abstract List<PeptideFragmentIntensityPrediction> Predictions { get; }


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
            var deserializedResponses = responses.Select(response => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(response)).ToList();
            var numBatches = deserializedResponses.Count;

            if (deserializedResponses.IsNullOrEmpty() || deserializedResponses.Any(r => r == null))
            {
                throw new Exception("Something went wrong during deserialization of responses.");
            }
            for (int batchIndex = 0; batchIndex < numBatches; batchIndex++)
            {
                var response = deserializedResponses[batchIndex];
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
                        fragmentIons,
                        fragmentMZs,
                        predictedIntensities
                    ));
                }
            }
        }

        /// <summary>
        /// Defaults to checking for canonical amino acids and length constraints
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public virtual bool IsValidPeptideSequence(string sequence)
        {
            var baseSequence = Regex.Replace(sequence, ModificationPattern, "");
            return Regex.IsMatch(baseSequence, CanonicalAminoAcidPattern)
                && baseSequence.Length <= MaxPeptideLength
                && baseSequence.Length >= MinPeptideLength;
        }

        /// <summary>
        /// Defaults to only allowing non-modified sequences
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public virtual bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            foreach (Match match in matches)
            {
                if (!ValidModificationUnimodMapping.ContainsKey(match.Value))
                {
                    return false;
                }
            }
            return true;
        }
    }
}
