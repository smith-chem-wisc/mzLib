namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Contract for a retention time predictor. Implementations produce a predicted
/// retention-time-equivalent value (time, iRT, or hydrophobicity, depending on the
/// predictor) for a given peptide, individually or in batches.
/// </summary>
/// <remarks>
/// <para>
/// <b>Disposal:</b> the interface extends <see cref="IDisposable"/> so callers holding
/// an <see cref="IRetentionTimePredictor"/>-typed variable can write
/// <c>using var p = RetentionTimePredictorFactory.Create(...);</c> directly. Implementors
/// that hold unmanaged resources (e.g. the TorchSharp model in
/// <see cref="Chronologer.ChronologerRetentionTimePredictor"/>) release them in
/// <see cref="IDisposable.Dispose"/>; lightweight predictors rely on the base-class no-op.
/// </para>
/// <para>
/// <b>Construction is deliberately not part of this interface.</b> Use
/// <see cref="RetentionTimePredictorFactory.Create"/> with a
/// <see cref="PredictorType"/> value to build concrete predictors. This keeps the
/// abstraction free of references to concrete types — most importantly the
/// TorchSharp-backed <see cref="Chronologer.ChronologerRetentionTimePredictor"/>,
/// which would otherwise force every consumer of this interface to transitively
/// depend on TorchSharp and ship its large native binaries.
/// </para>
/// <para>
/// <b>Migration:</b> the previously-nested <c>Create</c> method and
/// <c>PredictorType</c> enum on this interface have been moved to the top-level
/// <see cref="RetentionTimePredictorFactory"/> and <see cref="PredictorType"/>.
/// </para>
/// </remarks>
public interface IRetentionTimePredictor : IDisposable
{
    /// <summary>
    /// Name/identifier for this predictor (e.g., "SSRCalc3", "Chronologer")
    /// </summary>
    string PredictorName { get; }
    /// <summary>
    /// Gets the separation type this predictor is designed for
    /// </summary>
    SeparationType SeparationType { get; }
    /// <summary>
    /// Predicts a retention time equivalent for a given peptide.
    /// The value may represent a time, iRT, or hydrophobicity depending on the predictor.
    /// Returns null if prediction cannot be made.
    /// </summary>
    /// <returns>Predicted value in predictor-specific units, or null if prediction not possible</returns>
    double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);
    /// <summary>
    /// Predicts retention time equivalents for a batch of peptides.
    /// Results are materialized and safe to enumerate multiple times.
    /// </summary>
    IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1);
    /// <summary>
    /// Streams retention time equivalents as they are produced in parallel.
    /// Results may arrive out of order. Prefer this for large batches where
    /// results can be consumed incrementally.
    /// </summary>
    IEnumerable<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> StreamRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1);
    /// <inheritdoc cref="PredictRetentionTimeEquivalent"/>
    [Obsolete("Use PredictRetentionTimeEquivalent instead.")]
    public double? PredictRetentionTime(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        => PredictRetentionTimeEquivalent(peptide, out failureReason);
    public string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);
}
