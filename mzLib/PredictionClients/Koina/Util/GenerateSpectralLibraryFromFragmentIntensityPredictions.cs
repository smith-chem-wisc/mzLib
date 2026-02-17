using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PredictionClients.Koina.AbstractClasses;

namespace PredictionClients.Koina.Util
{
    public static class GenerateSpectralLibraryFromFragmentIntensityPredictions
    {
        #region Spectral Library Generation
        /// <summary>
        /// Transforms prediction results into LibrarySpectrum objects suitable for spectral library creation.
        /// Filters out impossible ions (intensity = -1) and low-intensity fragments below MinIntensityFilter.
        /// Parses fragment annotations (e.g., "b5+1") to extract ion type, fragment number, and charge state.
        /// </summary>
        /// <remarks>
        /// The conversion process:
        /// 1. Converts sequence format: UNIMOD -> mzLib -> mass-only for Peptide object creation
        /// 2. Parses each fragment annotation to determine ion properties
        /// 3. Creates MatchedFragmentIon objects with experimental m/z and predicted intensities
        /// 4. Builds LibrarySpectrum with precursor information and fragment data
        /// 5. Validates uniqueness of generated spectra by name
        /// </remarks>
        /// <exception cref="WarningException">Recorded in the out parameter when duplicate spectra are detected in predictions</exception>
        public static List<LibrarySpectrum> GenerateLibrarySpectraFromPredictions(List<PeptideFragmentIntensityPrediction> predictions, double[] alignedRetentionTimes)
        {
            if (predictions.Count != alignedRetentionTimes.Length)
            {
                throw new ArgumentException("The number of predictions must match the number of aligned retention times.");
            }
            var predictedSpectra = new List<LibrarySpectrum>();

            for (int predictionIndex = 0; predictionIndex < predictions.Count; predictionIndex++)
            {
                var prediction = predictions[predictionIndex];

                var peptide = new Peptide(
                    ConvertMzLibModificationsToMassesOnly(
                        ConvertUnimodToMzLibModifications(prediction.FullSequence)
                        )
                    );

                List<MatchedFragmentIon> fragmentIons = new();
                for (int fragmentIndex = 0; fragmentIndex < prediction.FragmentAnnotations.Count; fragmentIndex++)
                {
                    if ((int)prediction.FragmentIntensities[fragmentIndex] == -1 || prediction.FragmentIntensities[fragmentIndex] < MinIntensityFilter)
                    {
                        // Skip impossible ions and peaks with near zero intensity. The model uses -1 to indicate impossible ions.
                        continue;
                    }

                    var annotation = prediction.FragmentAnnotations[fragmentIndex];

                    // This check is to ensure that the annotation contains the expected format (e.g., "b5+1") before attempting to parse it.
                    // The API WILL always contain it in the expected format, but this is a safeguard against any malformed annotations that could cause exceptions during parsing.
                    if (annotation == null || !annotation.Contains('+'))
                    {
                        // Skip malformed annotations that do not contain expected ion type and charge information.
                        continue;
                    }

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
                            // neutralMass is not directly provided by models, and it is not necessary here. If needed, 
                            // compute it from the peptide sequence and fragment information as shown below.
                            // neutralMass: peptide.Fragment(Enum.Parse<FragmentTypes>(ionType), fragmentNumber).First().MonoisotopicMass,
                            neutralMass: 0.0, // Placeholder, not used in this context
                            fragmentNumber: fragmentNumber,
                            residuePosition: fragmentNumber, // For b / y ions, the fragment number corresponds to the residue count from the respective terminus.
                            neutralLoss: 0 // Annotations like "b5+1" do not encode neutral losses, so we explicitly assume no loss as placeholder.
                        ),
                        experMz: prediction.FragmentMZs[fragmentIndex],
                        experIntensity: prediction.FragmentIntensities[fragmentIndex],
                        charge: fragmentIonCharge
                    );

                    fragmentIons.Add(fragmentIon);
                }

                var spectrum = new LibrarySpectrum
                (
                    sequence: prediction.FullSequence,
                    precursorMz: peptide.ToMz(prediction.PrecursorCharge),
                    chargeState: prediction.PrecursorCharge,
                    peaks: fragmentIons,
                    rt: alignedRetentionTimes[predictionIndex]
                );

                predictedSpectra.Add(spectrum);
            }
            var unique = predictedSpectra.DistinctBy(p => p.Name).ToList();

            //TODO: if filepath given then set predicted spectra to unique spectra and save to file, if not just return the list spectra to keep alignment.
        }

        /// <summary>
        /// Saves predicted spectra to a spectral library file.
        /// Automatically calls GenerateLibrarySpectraFromPredictions() if not already done.
        /// </summary>
        /// <param name="filePath">Output file path for the spectral library</param>
        /// <example>
        /// Usage:
        /// model.SavePredictedSpectralLibrary("predicted_library.msp");
        /// </example>
        protected void SavePredictedSpectralLibrary(string filePath, out WarningException? warning)
        {
            warning = null;
            var spectralLibrary = new SpectralLibrary();
            if (PredictedSpectra.Count == 0)
            {
                GenerateLibrarySpectraFromPredictions(out warning);
            }
            spectralLibrary.Results = PredictedSpectra;
            spectralLibrary.WriteResults(filePath);
        }
        #endregion
    }
}
