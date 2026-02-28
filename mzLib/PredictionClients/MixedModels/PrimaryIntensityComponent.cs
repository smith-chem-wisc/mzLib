using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.LocalModels;
using System.ComponentModel;

namespace PredictionClients.MixedModels.Components
{
    /// <summary>
    /// Wraps any Koina FragmentIntensityModel (Prosit2020IntensityHCD, etc.)
    /// and contributes its predictions as PrimaryFragmentIntensities.
    ///
    /// The component owns the model's lifetime — it disposes the model after RunAsync()
    /// completes (which mirrors how the model already disposes itself after RunInferenceAsync).
    /// </summary>
    public class PrimaryIntensityComponent : IMixedModelComponent
    {
        private readonly FragmentIntensityModel _model;

        public string ComponentName { get; }
        public ContributionType ContributionType => ContributionType.PrimaryFragmentIntensities;

        /// <param name="model">
        /// A fully constructed, not-yet-run FragmentIntensityModel (e.g. Prosit2020IntensityHCD).
        /// The component takes ownership and will call RunInferenceAsync() exactly once.
        /// </param>
        /// <param name="componentName">
        /// Human-readable name for warnings and logs. Defaults to model.ModelName.
        /// </param>
        public PrimaryIntensityComponent(
            FragmentIntensityModel model,
            string? componentName = null)
        {
            _model = model;
            ComponentName = componentName ?? model.ModelName;
        }

        /// <summary>
        /// Runs the underlying Koina model and wraps its PredictedSpectra as a MixedModelResult.
        /// </summary>
        public async Task<MixedModelResult> RunAsync()
        {
            try
            {
                var warning = await _model.RunInferenceAsync();
                return MixedModelResult.FromSpectra(
                    ComponentName,
                    ContributionType,
                    _model.PredictedSpectra,
                    warning);
            }
            catch (Exception ex)
            {
                return MixedModelResult.FromError(ComponentName, ContributionType, ex);
            }
        }
    }

    /// <summary>
    /// Wraps InternalFragmentIntensityModel and contributes its predictions
    /// as InternalFragmentIntensities.
    ///
    /// Mirrors PrimaryIntensityComponent exactly — the only difference is the
    /// ContributionType, which tells the merger to add these as internal ions
    /// rather than replacing primary ions.
    /// </summary>
    public class InternalIntensityComponent : IMixedModelComponent
    {
        private readonly InternalFragmentIntensityModel _model;

        public string ComponentName { get; }
        public ContributionType ContributionType => ContributionType.InternalFragmentIntensities;

        /// <param name="model">
        /// A fully constructed, not-yet-run InternalFragmentIntensityModel.
        /// </param>
        /// <param name="componentName">
        /// Human-readable name for warnings and logs. Defaults to model.ModelName.
        /// </param>
        public InternalIntensityComponent(
            InternalFragmentIntensityModel model,
            string? componentName = null)
        {
            _model = model;
            ComponentName = componentName ?? model.ModelName;
        }

        /// <summary>
        /// Runs the local ONNX model and wraps its PredictedSpectra as a MixedModelResult.
        /// </summary>
        public async Task<MixedModelResult> RunAsync()
        {
            try
            {
                var warning = await _model.RunInferenceAsync();
                return MixedModelResult.FromSpectra(
                    ComponentName,
                    ContributionType,
                    _model.PredictedSpectra,
                    warning);
            }
            catch (Exception ex)
            {
                return MixedModelResult.FromError(ComponentName, ContributionType, ex);
            }
        }
    }

    // ── Future component stubs ───────────────────────────────────────────────────
    //
    // RetentionTimeComponent
    // ----------------------
    // Wraps a Koina RT model (e.g. Prosit2019iRTHCD or AlphaPept RT) and contributes
    // scalar retention time predictions keyed by "Sequence/Charge".
    // ContributionType: RetentionTime
    // ScalarData["PEPTIDEK/2"] = predictedRT
    //
    // Usage would look like:
    //
    //   public class RetentionTimeComponent : IMixedModelComponent
    //   {
    //       private readonly RetentionTimeModel _model;
    //       public ContributionType ContributionType => ContributionType.RetentionTime;
    //       ...
    //       public async Task<MixedModelResult> RunAsync()
    //       {
    //           var warning = await _model.RunInferenceAsync();
    //           var scalars = _model.PredictedRetentionTimes
    //               .ToDictionary(p => p.Key, p => p.Value);
    //           return MixedModelResult.FromScalars(ComponentName, ContributionType, scalars, warning);
    //       }
    //   }
    //
    // IonMobilityComponent
    // --------------------
    // Similar pattern. ContributionType.IonMobility.
    // LibrarySpectrumMerger would need a new merge case to write CCS into the spectrum.
    // (LibrarySpectrum does not yet have a CCS field — that would be a separate mzLib PR.)
}
