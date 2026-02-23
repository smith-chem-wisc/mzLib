using Omics.SpectrumMatch;
using System.ComponentModel;

namespace PredictionClients.MixedModels
{
    /// <summary>
    /// The typed output from a single IMixedModelComponent.
    ///
    /// A result carries a dictionary of LibrarySpectrum objects (one per peptide, keyed by
    /// "Sequence/Charge") plus the ContributionType that tells the LibrarySpectrumMerger
    /// how to fold this data into the combined spectrum.
    ///
    /// For contribution types that don't produce LibrarySpectra (e.g. RetentionTime,
    /// IonMobility), the Spectra dictionary will be empty and the scalar data is carried
    /// in the ScalarData bag. The merger knows how to handle both cases.
    ///
    /// DESIGN NOTE
    /// -----------
    /// Using LibrarySpectrum as the per-model output type (rather than, say, a raw float[])
    /// means each component can use the full existing infrastructure for fragment ion
    /// construction and MSP serialization. The merger's job is then purely additive:
    /// union fragment lists, resolve RT conflicts, etc. — not re-parse raw model outputs.
    /// </summary>
    public class MixedModelResult
    {
        /// <summary>The component that produced this result, for tracing and warning messages.</summary>
        public string ComponentName { get; init; }

        /// <summary>What kind of data this result contains.</summary>
        public ContributionType ContributionType { get; init; }

        /// <summary>
        /// Per-peptide spectra keyed by LibrarySpectrum.Name ("Sequence/Charge").
        /// Populated for PrimaryFragmentIntensities and InternalFragmentIntensities.
        /// Empty for scalar contribution types (RetentionTime, IonMobility).
        /// </summary>
        public IReadOnlyDictionary<string, LibrarySpectrum> Spectra { get; init; }
            = new Dictionary<string, LibrarySpectrum>();

        /// <summary>
        /// Scalar data keyed by "Sequence/Charge", for contribution types that don't
        /// produce fragment ions. For example, RetentionTime would store:
        ///   ScalarData["PEPTIDEK/2"] = 45.3
        ///
        /// Empty for fragment intensity contribution types.
        /// </summary>
        public IReadOnlyDictionary<string, double> ScalarData { get; init; }
            = new Dictionary<string, double>();

        /// <summary>
        /// Any non-fatal warning from the underlying model (e.g. filtered peptides,
        /// duplicate spectra). Null if the component ran cleanly.
        /// The CombinedLibraryModel surfaces these warnings to the caller.
        /// </summary>
        public WarningException? Warning { get; init; }

        /// <summary>
        /// True if the component completed successfully (even with warnings).
        /// False if the component threw and was caught — in that case, Error is set.
        /// </summary>
        public bool Succeeded { get; init; } = true;

        /// <summary>
        /// Set when Succeeded = false. Recorded but does not abort the merge —
        /// the CombinedLibraryModel will still merge results from other components
        /// and surface this as a warning to the caller.
        /// </summary>
        public Exception? Error { get; init; }

        // ── Convenience factories ────────────────────────────────────────────

        public static MixedModelResult FromSpectra(
            string componentName,
            ContributionType contributionType,
            IEnumerable<LibrarySpectrum> spectra,
            WarningException? warning = null)
        {
            return new MixedModelResult
            {
                ComponentName = componentName,
                ContributionType = contributionType,
                Spectra = spectra.ToDictionary(s => s.Name),
                Warning = warning,
                Succeeded = true,
            };
        }

        public static MixedModelResult FromScalars(
            string componentName,
            ContributionType contributionType,
            IReadOnlyDictionary<string, double> scalars,
            WarningException? warning = null)
        {
            return new MixedModelResult
            {
                ComponentName = componentName,
                ContributionType = contributionType,
                ScalarData = scalars,
                Warning = warning,
                Succeeded = true,
            };
        }

        public static MixedModelResult FromError(
            string componentName,
            ContributionType contributionType,
            Exception error)
        {
            return new MixedModelResult
            {
                ComponentName = componentName,
                ContributionType = contributionType,
                Error = error,
                Succeeded = false,
            };
        }
    }
}
