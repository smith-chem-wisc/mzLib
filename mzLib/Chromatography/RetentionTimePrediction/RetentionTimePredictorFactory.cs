using Chromatography.RetentionTimePrediction.Chronologer;
using Chromatography.RetentionTimePrediction.CZE;
using Chromatography.RetentionTimePrediction.SSRCalc;

namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Discriminator identifying which concrete retention time predictor to build
/// via <see cref="RetentionTimePredictorFactory.Create"/>.
/// </summary>
/// <remarks>
/// This enum was previously nested inside <see cref="IRetentionTimePredictor"/>
/// (as <c>IRetentionTimePredictor.PredictorType</c>). It was promoted to a top-level
/// type alongside the factory so the interface no longer needs to reference any
/// concrete predictor. Callers that used the old nested name should update to the
/// top-level <see cref="PredictorType"/>.
/// </remarks>
public enum PredictorType
{
    /// <summary>
    /// Krokhin/SSRCalc v3 — analytical hydrophobicity predictor for reversed-phase HPLC.
    /// Lightweight, no external model files or native dependencies.
    /// </summary>
    SSRCalc3,

    /// <summary>
    /// Capillary zone electrophoresis (CZE) migration-time predictor.
    /// Lightweight, no external model files or native dependencies.
    /// </summary>
    CZE,

    /// <summary>
    /// Chronologer deep-learning retention time predictor (TorchSharp-backed).
    /// Pulls in TorchSharp and large native runtime binaries; use only when
    /// model-based prediction is actually required.
    /// </summary>
    Chronologer
}

/// <summary>
/// Factory for <see cref="IRetentionTimePredictor"/> instances. Given a
/// <see cref="PredictorType"/>, returns a fully initialized predictor ready for use.
/// </summary>
/// <remarks>
/// <para>
/// <b>Why this is a standalone class rather than a static method on the interface:</b>
/// </para>
/// <para>
/// The factory must reference every concrete predictor type it can build — including
/// <see cref="ChronologerRetentionTimePredictor"/>, which transitively depends on
/// TorchSharp and its large native binaries. When that factory lived as a
/// <c>static</c> member of <see cref="IRetentionTimePredictor"/>, any project that
/// merely referenced the interface was forced to compile against — and ship — the
/// full predictor dependency graph, even if it only described or consumed predictors
/// (e.g. dependency-injection hosts, tests with fakes, lightweight tooling).
/// </para>
/// <para>
/// Moving construction into this dedicated class restores a clean split:
/// </para>
/// <list type="bullet">
///   <item><description>
///     <b>Consumers of the abstraction</b> (code that accepts or produces
///     <see cref="IRetentionTimePredictor"/> instances) depend only on the interface
///     and pay no transitive cost for predictors they never instantiate.
///   </description></item>
///   <item><description>
///     <b>Composition roots</b> (the one place in an application that actually
///     chooses which predictor to build) reference this factory explicitly and
///     accept the resulting dependency footprint.
///   </description></item>
/// </list>
/// <para>
/// <b>Migration note:</b> the old entry point <c>IRetentionTimePredictor.Create(...)</c>
/// has been removed. Replace calls with <see cref="Create"/> on this class, and replace
/// the nested enum <c>IRetentionTimePredictor.PredictorType</c> with the top-level
/// <see cref="PredictorType"/>.
/// </para>
/// <para>
/// <b>Extending:</b> to add a new predictor, add a value to <see cref="PredictorType"/>
/// and a corresponding branch in <see cref="Create"/>. No change to
/// <see cref="IRetentionTimePredictor"/> is required.
/// </para>
/// </remarks>
public static class RetentionTimePredictorFactory
{
    /// <summary>
    /// Creates a new <see cref="IRetentionTimePredictor"/> of the specified
    /// <paramref name="type"/>, constructed with default parameters.
    /// </summary>
    /// <param name="type">The predictor variant to instantiate.</param>
    /// <returns>A fully initialized predictor ready for prediction calls.</returns>
    /// <remarks>
    /// <para>
    /// Each call returns a fresh instance — instances are not cached or pooled.
    /// </para>
    /// <para>
    /// The caller owns the lifetime of the returned predictor and is responsible
    /// for disposing predictors that implement <see cref="IDisposable"/>. In
    /// particular, <see cref="ChronologerRetentionTimePredictor"/> holds an
    /// unmanaged TorchSharp model and must be disposed; SSRCalc3 and CZE hold
    /// no unmanaged resources but should still be disposed for consistency.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="type"/> is not a defined <see cref="PredictorType"/>
    /// value (for example, when an out-of-range integer is cast to the enum).
    /// </exception>
    public static IRetentionTimePredictor Create(PredictorType type) => type switch
    {
        PredictorType.SSRCalc3 => new SSRCalc3RetentionTimePredictor(),
        PredictorType.CZE => new CZERetentionTimePredictor(),
        PredictorType.Chronologer => new ChronologerRetentionTimePredictor(),
        _ => throw new ArgumentOutOfRangeException(nameof(type), type, null)
    };
}
