using Chemistry;
using Koina.Client;
using Koina.Interfaces;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using System.Threading.Tasks.Dataflow;
using TopDownProteomics;
using static System.Net.WebRequestMethods;

namespace Koina.SupportedModels.Prosit2020IntensityHCD
{
    public class Prosit2020IntensityHCD : IKoinaModelIO
    {
        public string ModelName { get; } = "Prosit_2020_intensity_HCD";
        public readonly int MaxPeptideLength = 30;
        public readonly HashSet<int> AllowedPrecursorCharges = new() { 1, 2, 3, 4, 5, 6 };
        public readonly int NumberOfPredictedIons = 174;
        public int BatchSize { get; private set; }
        public List<string> PeptideSequences;
        public List<int> PrecursorCharges;
        public List<int> CollisionEnergies;
        public List<LibrarySpectrum> PredictedSpectra;
        public List<double?> RetentionTimes;


        public SpectralLibrary PredictedSpectralLibrary { get; internal set; }
        public Prosit2020IntensityHCD(List<string> peptideSequences, List<int> precursorCharges, List<int> collisionEnergies, List<double?> retentionTimes)
        {
            // Ensure that the inputs meet the model requirements.
            for (int i = 0; i < peptideSequences.Count; i++)
            {
                var peptide = peptideSequences[i];
                var charge = precursorCharges[i];
                var energy = collisionEnergies[i];
                if (!(peptide.Length <= MaxPeptideLength)
                    && !AllowedPrecursorCharges.Contains(charge))
                {
                    throw new ArgumentException($"Input at index {i} does not meet model requirements: Peptide Length = {peptide.Length}, Precursor Charge = {charge}");
                }
            }

            PeptideSequences = peptideSequences;
            PrecursorCharges = precursorCharges;
            CollisionEnergies = collisionEnergies;
            BatchSize = PeptideSequences.Count;
            PredictedSpectralLibrary = new();
            RetentionTimes = retentionTimes;
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

        public Dictionary<string, object> ToRequest(string? requestID = null)
        {
            var request = new Dictionary<string, object>
            {
                { "id", requestID ?? Guid.NewGuid().ToString() },
                { "inputs", new List<object>
                    {
                        new {
                            name = "peptide_sequences",
                            shape = new[]{ BatchSize, 1 },
                            datatype = "BYTES",
                            data = PeptideSequences
                        },
                        new {
                            name = "precursor_charges",
                            shape = new[]{ BatchSize, 1 },
                            datatype = "INT32",
                            data = PrecursorCharges
                        },
                        new {
                            name = "collision_energies",
                            shape = new[]{ BatchSize, 1 },
                            datatype = "FP32",
                            data = CollisionEnergies
                        }
                    }
                }
            };

            return request;
        }

        public async Task RunInferenceAsync()
        {
            var _http = new HTTP();

            var response = await _http.InferenceRequest(ModelName, ToRequest());
            _ResponseToLibrarySpectra(response);
        }

        private void _ResponseToLibrarySpectra(string response)
        {
            PredictedSpectra = new List<LibrarySpectrum>();
            var deserializedResponse = Newtonsoft.Json.JsonConvert.DeserializeObject<ResponseJSONStruct>(response);

            if (deserializedResponse == null)
            {
                throw new Exception("Failed to deserialize response from Koina.");
            }

            var shape = deserializedResponse.Outputs[0].Shape; // Assuming all outputs have the same batch size. True for this model.
            string[] outputIonAnnotations = deserializedResponse.Outputs[0].Data.Select(d => (string)d);
            double[] outputMZs = deserializedResponse.Outputs[1].Data.Select(d => double.Parse(String.Join("", d)));
            double[] outputIntensities = deserializedResponse.Outputs[2].Data.Select(d => double.Parse(String.Join("", d)));

            // Create Peptide objects from sequences to facilitate mass calculations
            List<Peptide> peptides = new List<Peptide>();
            foreach (var peptideSeq in PeptideSequences)
            {
                peptides.Add(new Peptide(peptideSeq));
            }

            // Iterate through each batch/peptide to construct a LibrarySpectrum
            for (int i = 0; i < BatchSize; i++)
            {
                List<MatchedFragmentIon> fragmentIons = new();

                for (int j = i * NumberOfPredictedIons; j < (i + 1) * NumberOfPredictedIons; j++)
                {
                    if (outputMZs[j] == -1 || outputIntensities[j] == -1)
                    {
                        // Skip ions with invalid m/z. The model uses -1 to indicate impossible ions.
                        continue;
                    }

                    var annotation = outputIonAnnotations[j];
                    // Parse the annotation to get ion type, number and charge from something like 'b5+1'
                    var ionType = annotation.First().ToString(); // 'b' or 'y'
                    var fragmentNumber = int.Parse(String.Join("", annotation.SubSequence(1, annotation.IndexOf('+') - 1)));
                    var ionCharge = int.Parse(annotation.Last().ToString());

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
                        experMz: outputMZs[j],
                        experIntensity: outputIntensities[j],
                        charge: ionCharge
                    );

                    fragmentIons.Add(fragmentIon);
                }

                var spectrum = new LibrarySpectrum
                (
                    sequence: peptides[i].BaseSequence,
                    precursorMz: peptides[i].ToMz(PrecursorCharges[i]),
                    chargeState: PrecursorCharges[i],
                    peaks: fragmentIons,
                    rt: RetentionTimes[i]
                );

                PredictedSpectra.Add(spectrum);
            }
        }

        public void SavePredictedSpectralLibrary(string filePath)
        {
            var spectralLibrary = new SpectralLibrary();
            spectralLibrary.Results = PredictedSpectra;
            spectralLibrary.WriteResults(filePath);
        }
    }
}
