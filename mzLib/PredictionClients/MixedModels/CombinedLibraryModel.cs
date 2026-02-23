using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System.ComponentModel;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.LocalModels;
using PredictionClients.MixedModels.Components;

namespace PredictionClients.MixedModels
{
    /// <summary>
    /// Orchestrates multiple IMixedModelComponent instances and merges their outputs
    /// into a single list of LibrarySpectrum objects.
    ///
    /// From the caller's perspective this looks and feels like a single model:
    /// construct it, call RunAsync(), read PredictedSpectra. The fact that several
    /// models ran under the hood is transparent.
    ///
    /// USAGE — Prosit primary ions + local internal ions
    /// --------------------------------------------------
    ///
    ///   var primary = new Prosit2020IntensityHCD(
    ///       peptides, charges, energies, rts, out _, collisionEnergy: 35);
    ///
    ///   var internal = new InternalFragmentIntensityModel(
    ///       peptides, charges, rts, out _, onnxModelPath);
    ///
    ///   var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(
    ///       primary, internal);
    ///
    ///   var warning = await combined.RunAsync();
    ///   var spectra = combined.PredictedSpectra;  // primary + internal ions merged
    ///
    /// USAGE — custom component list
    /// ------------------------------
    ///   var combined = new CombinedLibraryModel(new List&lt;IMixedModelComponent&gt;
    ///   {
    ///       new PrimaryIntensityComponent(prositModel),
    ///       new InternalIntensityComponent(internalModel),
    ///       // new RetentionTimeComponent(rtModel),   // future
    ///   });
    ///
    /// PARALLELISM
    /// -----------
    /// All components are started simultaneously via Task.WhenAll. Koina components
    /// make async HTTP calls; local components run CPU work on the thread pool.
    /// The merge step runs after all components finish.
    ///
    /// FAILURE HANDLING
    /// ----------------
    /// A component that throws does not abort the whole run. Its error is captured in
    /// MixedModelResult.FromError and surfaced in the returned WarningException.
    /// The merge proceeds with whatever components succeeded. This means a Koina outage
    /// will produce internal-only spectra rather than crashing, and vice versa.
    ///
    /// SAVE PATH
    /// ---------
    /// If SpectralLibrarySavePath is set, the merged library is written to disk after
    /// the merge completes, in MetaMorpheus MSP format.
    /// </summary>
    public class CombinedLibraryModel
    {
        private readonly IReadOnlyList<IMixedModelComponent> _components;

        /// <summary>
        /// Combined LibrarySpectrum objects after RunAsync() completes.
        /// Each spectrum contains the unioned fragment ions from all components.
        /// Empty until RunAsync() is called.
        /// </summary>
        public List<LibrarySpectrum> PredictedSpectra { get; private set; } = new();

        /// <summary>
        /// Optional path to write the merged library as an MSP file after inference.
        /// If null, the library is only available in memory via PredictedSpectra.
        /// </summary>
        public string? SpectralLibrarySavePath { get; }

        // ── Constructors ──────────────────────────────────────────────────────────

        /// <summary>
        /// General-purpose constructor. Accepts any combination of components.
        /// </summary>
        public CombinedLibraryModel(
            IReadOnlyList<IMixedModelComponent> components,
            string? spectralLibrarySavePath = null)
        {
            if (components == null || components.Count == 0)
                throw new ArgumentException("At least one component is required.", nameof(components));

            _components = components;
            SpectralLibrarySavePath = spectralLibrarySavePath;
        }

        // ── Named factory methods ────────────────────────────────────────────────

        /// <summary>
        /// The primary use case today: Prosit primary ions merged with local internal ions.
        ///
        /// Handles constructing the two components and wiring them into a CombinedLibraryModel.
        /// Both models must already be constructed (validated, filtered) before calling this.
        /// </summary>
        /// <param name="primaryModel">
        /// A Koina fragment intensity model (e.g. Prosit2020IntensityHCD), not yet run.
        /// </param>
        /// <param name="internalModel">
        /// The local ONNX internal fragment model, not yet run.
        /// </param>
        /// <param name="spectralLibrarySavePath">
        /// If set, writes the merged MSP library to this path after inference.
        /// </param>
        public static CombinedLibraryModel WithPrimaryAndInternalFragments(
            FragmentIntensityModel primaryModel,
            InternalFragmentIntensityModel internalModel,
            string? spectralLibrarySavePath = null)
        {
            return new CombinedLibraryModel(
                new List<IMixedModelComponent>
                {
                    new PrimaryIntensityComponent(primaryModel),
                    new InternalIntensityComponent(internalModel),
                },
                spectralLibrarySavePath);
        }

        // ── Inference ─────────────────────────────────────────────────────────────

        /// <summary>
        /// Runs all components in parallel, then merges their outputs.
        /// </summary>
        /// <returns>
        /// A WarningException summarising any per-component warnings or merge-time
        /// conflicts. Null if everything completed without issue.
        /// </returns>
        public async Task<WarningException?> RunAsync()
        {
            // Fire all components simultaneously
            var results = await Task.WhenAll(
                _components.Select(c => c.RunAsync()));

            // Merge
            var mergedDict = LibrarySpectrumMerger.Merge(results, out var mergeWarning);
            PredictedSpectra = mergedDict.Values.ToList();

            // Optionally persist
            if (SpectralLibrarySavePath is not null)
                SaveLibrary(SpectralLibrarySavePath);

            return mergeWarning;
        }

        // ── Persistence ───────────────────────────────────────────────────────────

        /// <summary>
        /// Writes PredictedSpectra to an MSP file using the standard MetaMorpheus format.
        /// Can be called manually after RunAsync() if SpectralLibrarySavePath was not set.
        /// </summary>
        public void SaveLibrary(string filePath)
        {
            var library = new SpectralLibrary { Results = PredictedSpectra };
            library.WriteResults(filePath);
        }
    }
}
