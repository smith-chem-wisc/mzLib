using Chromatography;
using Chromatography.RetentionTimePrediction;
using Omics;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using Readers.SpectralLibrary;

namespace PredictionClients.SpectralLibraryGeneration
{
    public enum SpectralLibraryFormat
    {
        Msp,
        Msl
    }

    /// <summary>
    /// Builds predicted library spectra for peptides from a fragment intensity model, an optional
    /// retention time predictor and an optional detectability filter, and writes them as msp or msl.
    /// The generator keeps no record of what it generated.
    /// </summary>
    public class SpectralLibraryGenerator
    {
        public FragmentIntensityModel FragmentModel { get; }
        public IRetentionTimePredictor? RetentionTimePredictor { get; }
        public DetectabilityModel? DetectabilityModel { get; }

        /// <summary>
        /// A peptide is kept when its predicted probability of being detected (one minus the
        /// not-detectable probability) is at least this value. Only used with a detectability model.
        /// </summary>
        public double DetectabilityThreshold { get; }

        public SpectralLibraryGenerator(
            FragmentIntensityModel fragmentModel,
            IRetentionTimePredictor? retentionTimePredictor = null,
            DetectabilityModel? detectabilityModel = null,
            double detectabilityThreshold = 0.5)
        {
            // CZE predicts from the peptide's mass, which the inputs do not carry.
            if (retentionTimePredictor?.SeparationType == SeparationType.CZE)
            {
                throw new NotSupportedException(
                    $"{retentionTimePredictor.PredictorName} predicts from the peptide's mass, which the generator does not have.");
            }

            FragmentModel = fragmentModel;
            RetentionTimePredictor = retentionTimePredictor;
            DetectabilityModel = detectabilityModel;
            DetectabilityThreshold = detectabilityThreshold;
        }

        /// <summary>
        /// Predicts one library spectrum per peptide and charge the models accept. Inputs the detectability model
        /// judges undetectable are left out; peptides it cannot judge are kept. A peptide whose
        /// retention time cannot be predicted gets a spectrum with no retention time.
        /// </summary>
        public List<LibrarySpectrum> GenerateLibrarySpectra(List<FragmentIntensityPredictionInput> inputs)
        {
            inputs = KeepDetectable(inputs);
            Dictionary<string, double?> retentionTimes = PredictRetentionTimes(inputs);

            FragmentModel.Predict(inputs);

            // The warning only says that duplicates were collapsed and that nothing was written to disk.
            List<LibrarySpectrum> spectra = FragmentModel.GenerateLibrarySpectraFromPredictions(
                inputs.Select(input => retentionTimes[input.FullSequence]).ToArray(),
                out _,
                filepath: null,
                minIntensityFilter: FragmentModel.MinIntensityFilter);

            foreach (LibrarySpectrum spectrum in spectra)
            {
                spectrum.FragmentModelName = FragmentModel.ModelName;
                if (spectrum.RetentionTime.HasValue)
                {
                    spectrum.RetentionTimeModelName = RetentionTimePredictor!.PredictorName;
                }
            }

            return spectra;
        }

        /// <summary>
        /// Generates the spectra and writes them. See <see cref="WriteSpectralLibrary"/> for the filters.
        /// </summary>
        public List<LibrarySpectrum> GenerateAndWrite(
            List<FragmentIntensityPredictionInput> inputs,
            string path,
            SpectralLibraryFormat format,
            double? minMz = null,
            double? maxMz = null,
            double? minRelativeIntensity = null,
            int? maxFragmentCount = null)
        {
            return WriteSpectralLibrary(GenerateLibrarySpectra(inputs), path, format,
                minMz, maxMz, minRelativeIntensity, maxFragmentCount);
        }

        /// <summary>
        /// Writes the spectra after filtering their fragment ions, and returns the filtered spectra it
        /// wrote. The spectra passed in are not changed. A spectrum left with no fragment ions is not
        /// written. Read back, a spectrum can differ from the one returned: msp scales each spectrum to
        /// its strongest ion, and msl stores a missing retention time as 0 and has no model names.
        /// </summary>
        /// <param name="minMz">Fragment ions below this m/z are left out.</param>
        /// <param name="maxMz">Fragment ions above this m/z are left out.</param>
        /// <param name="minRelativeIntensity">
        /// Fragment ions weaker than this fraction of the spectrum's strongest ion are left out.
        /// </param>
        /// <param name="maxFragmentCount">Only this many of the strongest remaining ions are kept.</param>
        public static List<LibrarySpectrum> WriteSpectralLibrary(
            List<LibrarySpectrum> spectra,
            string path,
            SpectralLibraryFormat format,
            double? minMz = null,
            double? maxMz = null,
            double? minRelativeIntensity = null,
            int? maxFragmentCount = null)
        {
            // The msp writer takes the strongest ion of each spectrum, so one with none would throw.
            List<LibrarySpectrum> toWrite = spectra
                .Select(spectrum => FilterFragments(spectrum, minMz, maxMz, minRelativeIntensity, maxFragmentCount))
                .Where(spectrum => spectrum.MatchedFragmentIons.Count > 0)
                .ToList();

            switch (format)
            {
                case SpectralLibraryFormat.Msp:
                    new SpectralLibrary { Results = toWrite }.WriteResults(path);
                    break;
                case SpectralLibraryFormat.Msl:
                    MslWriter.WriteFromLibrarySpectra(path, toWrite);
                    break;
            }

            return toWrite;
        }

        private static LibrarySpectrum FilterFragments(
            LibrarySpectrum spectrum,
            double? minMz,
            double? maxMz,
            double? minRelativeIntensity,
            int? maxFragmentCount)
        {
            List<MatchedFragmentIon> ions = spectrum.MatchedFragmentIons;
            double strongest = ions.Count > 0 ? ions.Max(ion => ion.Intensity) : 0;

            List<int> kept = Enumerable.Range(0, ions.Count)
                .Where(i => (minMz == null || ions[i].Mz >= minMz)
                    && (maxMz == null || ions[i].Mz <= maxMz)
                    && (minRelativeIntensity == null || ions[i].Intensity >= strongest * minRelativeIntensity))
                .ToList();

            if (maxFragmentCount is { } count)
            {
                kept = kept.OrderByDescending(i => ions[i].Intensity).Take(count).Order().ToList();
            }

            return new LibrarySpectrum(spectrum.Sequence, spectrum.PrecursorMz, spectrum.ChargeState,
                kept.Select(i => ions[i]).ToList(), spectrum.RetentionTime, spectrum.IsDecoy)
            {
                FragmentModelName = spectrum.FragmentModelName,
                RetentionTimeModelName = spectrum.RetentionTimeModelName
            };
        }

        private List<FragmentIntensityPredictionInput> KeepDetectable(List<FragmentIntensityPredictionInput> inputs)
        {
            if (DetectabilityModel == null)
            {
                return inputs;
            }

            List<string> sequences = inputs.Select(input => input.FullSequence).Distinct(StringComparer.Ordinal).ToList();
            List<PeptideDetectabilityPrediction> predictions =
                DetectabilityModel.Predict(sequences.Select(sequence => new DetectabilityPredictionInput(sequence)).ToList());

            var undetectable = new HashSet<string>(StringComparer.Ordinal);
            for (int i = 0; i < sequences.Count; i++)
            {
                if (predictions[i].DetectabilityProbabilities is { } probabilities
                    && 1 - probabilities.NotDetectable < DetectabilityThreshold)
                {
                    undetectable.Add(sequences[i]);
                }
            }

            return inputs.Where(input => !undetectable.Contains(input.FullSequence)).ToList();
        }

        private Dictionary<string, double?> PredictRetentionTimes(List<FragmentIntensityPredictionInput> inputs)
        {
            var retentionTimes = new Dictionary<string, double?>(StringComparer.Ordinal);
            foreach (FragmentIntensityPredictionInput input in inputs)
            {
                retentionTimes[input.FullSequence] = null;
            }

            if (RetentionTimePredictor == null)
            {
                return retentionTimes;
            }

            var peptides = retentionTimes.Keys.Select(sequence => new SequenceOnlyPeptide(sequence)).ToList();
            foreach (var (value, peptide, _) in RetentionTimePredictor.PredictRetentionTimeEquivalents(peptides))
            {
                retentionTimes[peptide.FullSequence] = value is { } rt && double.IsFinite(rt) ? rt : null;
            }

            return retentionTimes;
        }

        /// <summary>
        /// A peptide known only by its full sequence. Building a PeptideWithSetModifications instead would
        /// look each modification up in mzLib's dictionary, which throws on modifications it does not hold.
        /// </summary>
        private sealed record SequenceOnlyPeptide(string FullSequence) : IRetentionPredictable
        {
            public string BaseSequence { get; } = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(FullSequence);
            public double MonoisotopicMass => double.NaN;
            public string FullSequenceWithMassShifts => string.Empty;
        }
    }
}
