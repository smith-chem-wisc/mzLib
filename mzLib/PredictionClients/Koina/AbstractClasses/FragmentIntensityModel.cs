using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.Interfaces;
using Readers.SpectralLibrary;
using System.ComponentModel;
using System.Text.RegularExpressions;
using Chemistry;
using Proteomics.AminoAcidPolymer;

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
        public virtual double MinIntensityFilter { get; protected set; } = 1e-6;
        public virtual string CanonicalAminoAcidPattern { get; protected set; } = @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        public virtual string ModificationPattern { get; protected set; } = @"\[[^\]]+\]";
        public virtual Dictionary<string, string> ValidModificationUnimodMapping { get; protected set; } = new();
        public virtual Dictionary<string, double> ValidModificationsMonoisotopicMasses { get; protected set; } = new();
        public abstract List<string> PeptideSequences { get; }
        public abstract List<int> PrecursorCharges { get; }
        public abstract List<double?> RetentionTimes { get; }
        public abstract List<PeptideFragmentIntensityPrediction> Predictions { get; protected set; }
        public virtual List<LibrarySpectrum> PredictedSpectra { get; protected set; } = new();


        protected abstract List<Dictionary<string, object>> ToBatchedRequests();

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

        protected virtual void ResponseToPredictions(string[] responses)
        {
            if (PeptideSequences.Count == 0)
            {
                return;
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
        protected virtual bool IsValidPeptideSequence(string sequence)
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
        protected virtual bool HasValidModifications(string sequence)
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

        /// <summary>
        /// Converts a peptide sequence from the mzLib modification format to the Prosit UNIMOD format.
        /// By default, all unmodified cysteines are carbamidomethylated (UNIMOD:4) to match the expectations
        /// of the Prosit 2020 HCD intensity model.
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib modification format.</param>
        /// <returns>The sequence converted to Prosit UNIMOD format with cysteines carbamidomethylated.</returns>
        protected virtual string ConvertToPrositModificationFormat(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Key, mod.Value);
            }

            // Carbamidomethylate all Cysteines if not already modified
            return Regex.Replace(sequence, @"C(?!\[UNIMOD:4\])", "C[UNIMOD:4]");
        }

        protected virtual string ConvertToMzLibModificationFormat(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Value, mod.Key);
            }
            return sequence;
        }

        protected virtual string ConvertToMzLibModificationFormatWithMassesOnly(string sequence)
        {
            foreach (var mod in ValidModificationsMonoisotopicMasses)
            {
                sequence = sequence.Replace(mod.Key, $"[{mod.Value.ToString("F6")}]");
            }
            return sequence;
        }

        public void GenerateLibrarySpectraFromPredictions()
        {
            if (Predictions.Count == 0)
            {
                return; // No predictions to process
            }

            for (int predictionIndex = 0; predictionIndex < Predictions.Count; predictionIndex++)
            {
                var prediction = Predictions[predictionIndex];

                var peptide = new Peptide(
                    ConvertToMzLibModificationFormatWithMassesOnly(
                        ConvertToMzLibModificationFormat(prediction.PeptideSequence)
                        )
                    );

                List<MatchedFragmentIon> fragmentIons = new();
                for (int fragmentIndex = 0; fragmentIndex < prediction.FragmentAnnotations.Count; fragmentIndex++)
                {
                    if (prediction.FragmentIntensities[fragmentIndex] == -1 || prediction.FragmentIntensities[fragmentIndex] < MinIntensityFilter)
                    {
                        // Skip impossible ions and peaks with near zero intensity. The model uses -1 to indicate impossible ions.
                        continue;
                    }

                    var annotation = prediction.FragmentAnnotations[fragmentIndex];
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
                            // neutralMass is not directly provided by Prosit, and it is not necessary here. If needed, 
                            // compute it from the peptide sequence and fragment information as shown below.
                            // neutralMass: peptide.Fragment(Enum.Parse<FragmentTypes>(ionType), fragmentNumber).First().MonoisotopicMass,
                            neutralMass: 0.0, // Placeholder, not used in this context
                            fragmentNumber: fragmentNumber,
                            residuePosition: fragmentNumber, // For b / y ions, the fragment number corresponds to the residue count from the respective terminus.
                            neutralLoss: 0 // Prosit annotations like "b5+1" do not encode neutral losses, so we explicitly assume no loss as placeholder.
                        ),
                        experMz: prediction.FragmentMZs[fragmentIndex],
                        experIntensity: prediction.FragmentIntensities[fragmentIndex],
                        charge: fragmentIonCharge
                    );

                    fragmentIons.Add(fragmentIon);
                }

                var spectrum = new LibrarySpectrum
                (
                    sequence: prediction.PeptideSequence,
                    precursorMz: peptide.ToMz(PrecursorCharges[predictionIndex]),
                    chargeState: PrecursorCharges[predictionIndex],
                    peaks: fragmentIons,
                    rt: RetentionTimes[predictionIndex]
                );

                PredictedSpectra.Add(spectrum);
            }
            var unique = PredictedSpectra.DistinctBy(p => p.Name).ToList();
            if (unique.Count != PredictedSpectra.Count)
            {
                throw new WarningException($"Duplicate spectra found in predictions. Reduced from {PredictedSpectra.Count} predicted spectra to {unique.Count} unique spectra.");
            }
        }

        public void SavePredictedSpectralLibrary(string filePath)
        {
            var spectralLibrary = new SpectralLibrary();
            GenerateLibrarySpectraFromPredictions();
            spectralLibrary.Results = PredictedSpectra;
            spectralLibrary.WriteResults(filePath);
        }
    }
}
