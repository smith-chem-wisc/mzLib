using System.ComponentModel;

namespace PredictionClients.MixedModels
{
    /// <summary>
    /// The single contract that every contributor to a CombinedLibraryModel must satisfy.
    ///
    /// A component wraps exactly one underlying model (Koina, local ONNX, or anything else)
    /// and returns its output as a MixedModelResult — a typed container that the
    /// LibrarySpectrumMerger knows how to fold into a combined LibrarySpectrum.
    ///
    /// IMPLEMENTING A NEW COMPONENT
    /// -----------------------------
    /// 1. Create a class in PredictionClients/MixedModels/Components/
    /// 2. Implement IMixedModelComponent
    /// 3. RunAsync() should:
    ///    a) call the underlying model's RunInferenceAsync()
    ///    b) wrap its output as a MixedModelResult with the appropriate ContributionType
    ///    c) propagate any WarningException via MixedModelResult.Warning
    /// 4. Register it in a CombinedLibraryModel constructor
    ///
    /// Components are intentionally thin — they do not merge, deduplicate, or reformat data.
    /// All of that is the responsibility of LibrarySpectrumMerger.
    /// </summary>
    public interface IMixedModelComponent
    {
        /// <summary>
        /// Human-readable name used in warnings and log messages (e.g. "Prosit2020HCD", "InternalFragmentV3").
        /// </summary>
        string ComponentName { get; }

        /// <summary>
        /// The kind of data this component contributes.
        /// The merger uses this to decide how to fold the result into the combined spectrum.
        /// </summary>
        ContributionType ContributionType { get; }

        /// <summary>
        /// Runs the underlying model and returns its typed output.
        /// Must not throw for recoverable errors — record them in MixedModelResult.Warning instead.
        /// </summary>
        Task<MixedModelResult> RunAsync();
    }

    /// <summary>
    /// Describes what kind of data a component contributes to the merged LibrarySpectrum.
    /// The LibrarySpectrumMerger dispatches on this value.
    ///
    /// Adding a new contribution type (e.g. IonMobility) requires:
    ///   1. Adding a value here
    ///   2. Adding a merge case in LibrarySpectrumMerger.Merge()
    ///   3. Implementing a component that produces it
    /// No other code needs to change.
    /// </summary>
    public enum ContributionType
    {
        /// <summary>Primary b/y fragment ion intensities (e.g. from Prosit, ms2pip, UniSpec).</summary>
        PrimaryFragmentIntensities,

        /// <summary>Internal fragment ion intensities from our local ONNX model.</summary>
        InternalFragmentIntensities,

        /// <summary>
        /// Retention time prediction (e.g. from Prosit RT, DeepLC, AlphaFold).
        /// Not implemented today — reserved for future use.
        /// </summary>
        RetentionTime,

        /// <summary>
        /// Ion mobility / collisional cross-section prediction.
        /// Not implemented today — reserved for future use.
        /// </summary>
        IonMobility,
    }
}
