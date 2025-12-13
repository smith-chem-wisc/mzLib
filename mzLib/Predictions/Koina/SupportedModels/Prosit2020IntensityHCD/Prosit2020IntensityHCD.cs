using Chemistry;
using Predictions.Koina.Client;
using Predictions.Koina.Interfaces;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Proteomics.AminoAcidPolymer;
using Readers.SpectralLibrary;
using System.ComponentModel;
using TopDownProteomics;
using System.Text.RegularExpressions;

namespace Predictions.Koina.SupportedModels.Prosit2020IntensityHCD
{
    public class Prosit2020IntensityHCD : IKoinaModelIO
    {
        public string ModelName => "Prosit_2020_intensity_HCD";
        public int MaxBatchSize => 1000;
        public readonly int MaxPeptideLength = 30;
        public readonly HashSet<int> AllowedPrecursorCharges = new() { 1, 2, 3, 4, 5, 6 };
        public readonly int NumberOfPredictedFragmentIons = 174;
        public readonly static Dictionary<string, string> ValidModificationUnimodMapping = new()
        {
            {"[Common Variable:Oxidation on M]", "[UNIMOD:35]"},
            {"[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]"}
        };
        public readonly static Dictionary<string, double> ValidModificationsMonoisotopicMasses = new()
        {
            {"[Common Variable:Oxidation on M]", 15.994915 },
            {"[Common Fixed:Carbamidomethyl on C]", 57.021464 }
        };
        public static string ModificationPattern = @"\[[^\]]+\]";
        public static string CanonicalAminoAcidPattern = @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        public readonly List<string> PeptideSequences = new();
        public readonly List<int> PrecursorCharges = new();
        public readonly List<int> CollisionEnergies = new();
        public readonly List<double?> RetentionTimes = new();
        public List<LibrarySpectrum> PredictedSpectra = new();
        public double MinIntensityFilter; // Tolerance for intensity filtering of predicted peaks


        public Prosit2020IntensityHCD(List<string> peptideSequences, List<int> precursorCharges, List<int> collisionEnergies, List<double?> retentionTimes, out WarningException? warnings, double minIntensityFilter=1e-4)
        {
            // Verify input lists are of the same length
            if (peptideSequences.Count != precursorCharges.Count
                || precursorCharges.Count != collisionEnergies.Count
                || collisionEnergies.Count != retentionTimes.Count)
            {
                throw new ArgumentException("Input lists must have the same length.");
            }

            // Ensure that the inputs meet the model requirements.
            List<string> invalidArguments = new();
            for (int i = 0; i < peptideSequences.Count; i++)
            {
                var peptide = peptideSequences[i];
                var charge = precursorCharges[i];
                var energy = collisionEnergies[i];
                var retentionTime = retentionTimes[i];

                // Skip invalid peptides
                if (!HasValidModifications(peptide) ||
                    !IsValidSequence(peptide) ||
                    !AllowedPrecursorCharges.Contains(charge) ||
                    energy <= 0
                    )
                {
                    invalidArguments.Add($"Index {i}: Peptide '{peptide}' (Length: {peptide.Length}), Charge: {charge}, Collision Energy: {energy}");
                }
                else
                {
                    PeptideSequences.Add(ConvertToPrositModificationFormat(peptide));
                    PrecursorCharges.Add(charge);
                    CollisionEnergies.Add(energy);
                    RetentionTimes.Add(retentionTime);
                }
            }

            warnings = null;
            if (invalidArguments.Count > 0)
            {
                string warningMessage = "The following input entries are invalid and will be skipped:\n"
                    + string.Join("\n", invalidArguments)
                    + "\nModel Requirements:\n"
                    + $"- Peptide length <= {MaxPeptideLength}\n"
                    + "- Peptide sequence is not empty\n"
                    + "- Peptide has valid modifications\n"
                    + $"- Precursor charge in [{string.Join(", ", AllowedPrecursorCharges)}]\n"
                    + "- Collision energy > 0\n"
                    + "- Retention time is null or >= 0";
                warnings = new WarningException(warningMessage);
            }

            MinIntensityFilter = minIntensityFilter;
        }


        // Construct from a list of LibrarySpectrum objects
        public Prosit2020IntensityHCD(List<LibrarySpectrum> spectralLibrary)
        {
            throw new NotImplementedException();
        }


        // Construct from a a spectral library file path
        public Prosit2020IntensityHCD(string filePath)
        {
            throw new NotImplementedException();
        }


        public List< Dictionary<string, object> > ToBatchedRequests()
        {

            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
            var batchedCharges = PrecursorCharges.Chunk(MaxBatchSize).ToList();
            var batchedEnergies = CollisionEnergies.Chunk(MaxBatchSize).ToList();

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
                            new {
                                name = "precursor_charges",
                                shape = new[]{ batchedCharges[i].Length, 1 },
                                datatype = "INT32",
                                data = batchedCharges[i]
                            },
                            new {
                                name = "collision_energies",
                                shape = new[]{ batchedEnergies[i].Length, 1 },
                                datatype = "FP32",
                                data = batchedEnergies[i]
                            }
                        }
                    }
                };
                batchedRequests.Add(request);
            }

            return batchedRequests;
        }


        public async Task RunInferenceAsync()
        {
            var _http = new HTTP(timeoutInMinutes: (int)(PeptideSequences.Count/MaxBatchSize) * 2 + 1); // Typically a full batch takes about a minute. Setting it to double that for safety.

            var responses = await Task.WhenAll(ToBatchedRequests().Select(request => _http.InferenceRequest(ModelName, request)));
            _ResponseToLibrarySpectra(responses);
            _http.Client.Dispose();
        }


        private void _ResponseToLibrarySpectra(string[] responses)
        {
            PredictedSpectra = new List<LibrarySpectrum>();
            var deserializedResponses = responses.Select(response => Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(response)).ToList();
            var numBatches = deserializedResponses.Count;

            if (deserializedResponses == null)
            {
                throw new Exception("Failed to deserialize response from Koina.");
            }

            for (int batchIndex = 0; batchIndex < numBatches; batchIndex++)
            {
                var responseBatch = deserializedResponses[batchIndex];
                if (responseBatch == null || responseBatch.Outputs.Count != 3)
                {
                    throw new Exception($"API response is not in the expected format. Expected 3 outputs, got {responseBatch?.Outputs.Count}.");
                }

                var currentBatchSize = responseBatch.Outputs[0].Shape[0];
                string[] outputIonAnnotations = responseBatch.Outputs[0].Data.Select(d => (string)d);
                double[] outputMZs = responseBatch.Outputs[1].Data.Select(d => double.Parse(String.Join("", d)));
                double[] outputIntensities = responseBatch.Outputs[2].Data.Select(d => double.Parse(String.Join("", d)));

                // Iterate through each batch/peptide to construct a LibrarySpectrum
                for (int precursorIndex = 0; precursorIndex < currentBatchSize; precursorIndex++)
                {
                    int linearIndexBatchOffset = batchIndex * MaxBatchSize;
                    var peptide = new Peptide(
                        ConvertToMzLibModificationFormatWithMassesOnly(
                            ConvertToMzLibModificationFormat(PeptideSequences[linearIndexBatchOffset + precursorIndex])
                            )
                        );
                    List<MatchedFragmentIon> fragmentIons = new();

                    for (int fragmentIndex = precursorIndex * NumberOfPredictedFragmentIons; fragmentIndex < (precursorIndex + 1) * NumberOfPredictedFragmentIons; fragmentIndex++)
                    {
                        if (outputMZs[fragmentIndex] == -1 || outputIntensities[fragmentIndex] < MinIntensityFilter)
                        {
                            // Skip impossible ions and peaks with near zero intensity. The model uses -1 to indicate impossible ions.
                            continue;
                        }

                        var annotation = outputIonAnnotations[fragmentIndex];
                        // Parse the annotation to get ion type, number and charge from something like 'b5+1'
                        var ionType = annotation.First().ToString(); // 'b' or 'y'
                        var fragmentNumber = int.Parse(String.Join("", annotation.SubSequence(1, annotation.IndexOf('+') - 1)));
                        var fragmentIonCharge = int.Parse(annotation.Last().ToString());

                        // Create a new MatchedFragmentIon for each output
                        var fragmentIon = new MatchedFragmentIon
                        (
                            neutralTheoreticalProduct: new Product
                            (
                                productType: Enum.Parse<ProductType>(ionType),
                                terminus: ionType == "b" ? FragmentationTerminus.N : FragmentationTerminus.C,
                                neutralMass: 0, // Placeholder, would need to calculate theoretical mass
                                fragmentNumber: fragmentNumber,
                                residuePosition: 0, // Placeholder, would need to calculate position
                                neutralLoss: 0 // Placeholder, would need to calculate neutral loss if any

                            ),
                            experMz: outputMZs[fragmentIndex],
                            experIntensity: outputIntensities[fragmentIndex],
                            charge: fragmentIonCharge
                        );

                        fragmentIons.Add(fragmentIon);
                    }

                    var spectrum = new LibrarySpectrum
                    (
                        sequence: PeptideSequences[linearIndexBatchOffset + precursorIndex],
                        precursorMz: peptide.ToMz(PrecursorCharges[linearIndexBatchOffset + precursorIndex]),
                        chargeState: PrecursorCharges[linearIndexBatchOffset + precursorIndex],
                        peaks: fragmentIons,
                        rt: RetentionTimes[linearIndexBatchOffset + precursorIndex]
                    );

                    PredictedSpectra.Add(spectrum);
                }
            }
            var unique = PredictedSpectra.DistinctBy(p => p.Name).ToList();
            if (unique.Count != PredictedSpectra.Count)
            {
                throw new Exception($"Duplicate spectra found in predictions. Reduced from {PredictedSpectra.Count} to {unique.Count} unique spectra.");
            }
        }


        public void SavePredictedSpectralLibrary(string filePath)
        {
            var spectralLibrary = new SpectralLibrary();
            spectralLibrary.Results = PredictedSpectra;
            spectralLibrary.WriteResults(filePath);
        }


        internal bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            if (matches.Count == 0)
            {
                return true; // No modifications, valid
            }
            else
            {
                foreach (Match match in matches)
                {
                    if (!ValidModificationUnimodMapping.ContainsKey(match.Value))
                    {
                        return false; // Invalid modification found
                    }
                }
                return true; // All modifications are valid
            }
        }

        internal bool IsValidSequence(string sequence)
        {
            var baseSequence = Regex.Replace(sequence, ModificationPattern, "");
            if (Regex.IsMatch(baseSequence, CanonicalAminoAcidPattern) &&
                baseSequence.Length <= MaxPeptideLength)
            {
                return true; // Valid sequence
            }
            else
            {
                return false; // Invalid characters found
            }
        }

        internal string ConvertToPrositModificationFormat(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Key, mod.Value);
            }

            // Carbamidomethylate all Cysteines if not already modified
            sequence = Regex.Replace(sequence, @"C(?!\[UNIMOD:4\])", "C[UNIMOD:4]");

            return sequence;
        }

        internal string ConvertToMzLibModificationFormat(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Value, mod.Key);
            }
            return sequence;
        }

        internal string ConvertToMzLibModificationFormatWithMassesOnly(string sequence)
        {
            foreach (var mod in ValidModificationsMonoisotopicMasses)
            {
                sequence = sequence.Replace(mod.Key, $"[{mod.Value.ToString("F6")}]");
            }
            return sequence;
        }
    }
}
